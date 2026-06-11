library(SummarizedExperiment)
library(sva)

# Set your project directory:
# setwd("/path/to/Kbhb")
if (!interactive()) setwd(normalizePath("."))

# ============================================================
# PRE-PROCESAMIENTO PARA ARACNE-AP: Basal BRCA
#
# Salidas en data/:
#   tcga_basal_matrix.txt    — entrada ARACNe-AP (genes × muestras, tab-sep)
#   mtbrc_basal_matrix.txt   — entrada ARACNe-AP (genes × muestras, tab-sep)
#   TFs_Lambert2018.txt      — lista de TFs (Lambert et al. 2018)
#   tcga_basal_expr.rds      — matriz expresión Basal TCGA  (para mra_kbhb.R)
#   tcga_normal_expr.rds     — matriz expresión Normal TCGA (contraste MRA)
#   mtbrc_basal_expr.rds     — matriz expresión Basal METABRIC (ComBat-corregida)
#   mtbrc_normal_expr.rds    — matriz expresión Normal METABRIC (ComBat-corregida)
# ============================================================


# ============================================================
# SECCIÓN 1: TCGA-BRCA
# ============================================================

# ── 1.1  Cargar SE y metadatos ──────────────────────────────

se <- readRDS("data/tcga_brca_rnaseq_se.rds")
cd <- as.data.frame(colData(se))
rd <- as.data.frame(rowData(se))

# ── 1.2  Mapa Ensembl ID → símbolo HGNC (solo protein_coding) ─

pc_idx <- rd$gene_type == "protein_coding"
message("Genes protein_coding: ", sum(pc_idx))

# Ensembl IDs duplicados al mismo símbolo: conservar el de mayor expresión media
se_pc <- se[pc_idx, ]
rd_pc <- rd[pc_idx, ]

tpm_all <- assay(se_pc, "tpm_unstrand")
rownames(tpm_all) <- rd_pc$gene_name

dup_syms <- rownames(tpm_all)[duplicated(rownames(tpm_all))]
message("Símbolos duplicados (protein_coding): ", length(unique(dup_syms)))

df_tpm <- as.data.frame(tpm_all)
df_tpm$.sym <- rownames(tpm_all)
tpm_dedup <- do.call(rbind,
  lapply(
    split(df_tpm[, colnames(df_tpm) != ".sym"], df_tpm$.sym),
    function(chunk) apply(chunk, 2, median, na.rm = TRUE)
  )
)
message("Genes únicos tras colapso: ", nrow(tpm_dedup))

# ── 1.3  Asignar barcodes TCGA como colnames ────────────────

# Usar los primeros 16 caracteres del barcode (incluye tipo de muestra)
colnames(tpm_dedup) <- substr(cd$barcode, 1, 16)

# ── 1.4  Filtrar muestras por subtipo ───────────────────────

# %in% trata NA como FALSE; unique() elimina aliquots duplicados del mismo sample
basal_bc  <- unique(substr(cd$barcode[cd$paper_BRCA_Subtype_PAM50 %in% "Basal" &
                                       cd$shortLetterCode %in% "TP"], 1, 16))
normal_bc <- unique(substr(cd$barcode[cd$shortLetterCode %in% "NT"], 1, 16))

cat("Muestras Basal TP:", length(basal_bc), "\n")
cat("Muestras Normal NT:", length(normal_bc), "\n")

# ── 1.5  log2(TPM + 1) ──────────────────────────────────────

tcga_basal  <- log2(tpm_dedup[, basal_bc]  + 1)
tcga_normal <- log2(tpm_dedup[, normal_bc] + 1)

# ── 1.6  Filtro de expresión: ≥20 % de Basales con log2TPM > 1 ─

expressed <- rowSums(tcga_basal > 1) >= 0.2 * ncol(tcga_basal)
tcga_basal  <- tcga_basal[expressed, ]
tcga_normal <- tcga_normal[rownames(tcga_basal), ]

message("Genes tras filtro expresión TCGA: ", nrow(tcga_basal))

# ── 1.7  Guardar ────────────────────────────────────────────

saveRDS(tcga_basal,  "data/tcga_basal_expr.rds")
saveRDS(tcga_normal, "data/tcga_normal_expr.rds")
message("data/tcga_basal_expr.rds  — ", nrow(tcga_basal),  " genes × ", ncol(tcga_basal),  " muestras")
message("data/tcga_normal_expr.rds — ", nrow(tcga_normal), " genes × ", ncol(tcga_normal), " muestras")

# ── 1.8  Escribir matriz ARACNe-AP ──────────────────────────
# Formato: primera fila = "gene" \t sample1 \t sample2 ...

write_aracne_matrix <- function(mat, path) {
  df <- cbind(gene = rownames(mat), as.data.frame(mat))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(path, " escrito — ", nrow(mat), " genes × ", ncol(mat), " muestras")
}

write_aracne_matrix(tcga_basal, "data/tcga_basal_matrix.txt")


# ============================================================
# SECCIÓN 2: METABRIC
# ============================================================

# ── 2.1  Cargar datos ───────────────────────────────────────

cl   <- readRDS("data/mtbrc_clinical_mrna.rds")
expr <- readRDS("data/mtbrc_expr.rds")

basal_ids  <- cl$SAMPLE_ID[cl$CLAUDIN_SUBTYPE == "Basal"]
normal_ids <- cl$SAMPLE_ID[cl$CLAUDIN_SUBTYPE == "Normal"]

cat("Muestras Basal METABRIC:", length(basal_ids), "\n")
cat("Muestras Normal METABRIC:", length(normal_ids), "\n")

# ── 2.2  Subconjunto de muestras ────────────────────────────
# METABRIC ya está en escala log2 (microarray Illumina)

mtbrc_basal  <- expr[, colnames(expr) %in% basal_ids]
mtbrc_normal <- expr[, colnames(expr) %in% normal_ids]

# ── 2.3  Filtro de expresión: señal log2 > 6 en ≥20 % de Basales ─
# Equivalente a intensidad ~64 en escala linear (por encima del background)

expressed <- rowSums(mtbrc_basal > 6) >= 0.2 * ncol(mtbrc_basal)
mtbrc_basal  <- mtbrc_basal[expressed, ]
mtbrc_normal <- mtbrc_normal[rownames(mtbrc_basal), ]

message("Genes tras filtro expresión METABRIC: ", nrow(mtbrc_basal))

# ── 2.4  ComBat: corrección de efecto de cohorte ─────────────
# METABRIC tiene 5 cohortes (hospitales UK/CA). ComBat preserva
# la señal biológica Basal vs Normal mientras elimina el batch.

# Combinar Basal + Normal para corrección conjunta
combined      <- cbind(mtbrc_basal, mtbrc_normal)
group_labels  <- c(rep("Basal",  ncol(mtbrc_basal)),
                   rep("Normal", ncol(mtbrc_normal)))

# Recuperar COHORT para estas muestras
cohort_map <- setNames(cl$COHORT, cl$SAMPLE_ID)
batch_vec  <- cohort_map[colnames(combined)]

cat("Distribución de cohortes en Basal + Normal:\n")
print(table(batch_vec, group_labels))

# mod protege la señal grupo biológico durante la corrección
mod <- model.matrix(~ group_labels)

message("Aplicando ComBat (batch = COHORT, mod = ~grupo)...")
combined_corrected <- ComBat(
  dat   = combined,
  batch = batch_vec,
  mod   = mod
)

# Separar de nuevo en Basal y Normal
mtbrc_basal  <- combined_corrected[, group_labels == "Basal"]
mtbrc_normal <- combined_corrected[, group_labels == "Normal"]

message("ComBat aplicado")

# ── 2.5  Guardar ────────────────────────────────────────────

saveRDS(mtbrc_basal,  "data/mtbrc_basal_expr.rds")
saveRDS(mtbrc_normal, "data/mtbrc_normal_expr.rds")
message("data/mtbrc_basal_expr.rds  — ", nrow(mtbrc_basal),  " genes × ", ncol(mtbrc_basal),  " muestras")
message("data/mtbrc_normal_expr.rds — ", nrow(mtbrc_normal), " genes × ", ncol(mtbrc_normal), " muestras")

write_aracne_matrix(mtbrc_basal, "data/mtbrc_basal_matrix.txt")


# ============================================================
# SECCIÓN 3: TFs Lambert 2018
# ============================================================

lambert_csv <- "data/lambert2018_DatabaseExtract_v1.01.csv"
if (!file.exists(lambert_csv)) {
    download.file(
        url      = "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv",
        destfile = lambert_csv,
        mode     = "wb"
    )
}
lambert_db      <- read.csv(lambert_csv, check.names = FALSE)
lambert_symbols <- unique(na.omit(lambert_db[["HGNC symbol"]][lambert_db[["Is TF?"]] == "Yes"]))
cat("Lambert et al. TFs:", length(lambert_symbols), "\n")

writeLines(lambert_symbols, "data/TFs_Lambert2018.txt")
message("data/TFs_Lambert2018.txt escrito")


message("\n====  PRE-PROCESAMIENTO COMPLETO  ====")
message("Siguientes pasos:")
message("  1. bash run_aracne_ap.sh   — construir redes")
message("  2. Rscript mra_kbhb.R      — master regulator analysis")
