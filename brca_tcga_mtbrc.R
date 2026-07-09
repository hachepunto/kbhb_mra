library(TCGAbiolinks)
library(SummarizedExperiment)

# Set your project directory:
# setwd("/path/to/Kbhb")
if (!interactive()) setwd(normalizePath("."))

# ============================================================
# MÓDULO BASE: BRCA — TCGA + METABRIC
#
# Salidas en data/:
#   tcga_brca_subtypes.rds   — subtipos PAM50 + TNBC por IHC (TCGA)
#   tcga_brca_rnaseq_se.rds  — SummarizedExperiment RNA-seq STAR Counts
#   mtbrc_clinical.rds       — clínica completa METABRIC
#   mtbrc_clinical_mrna.rds  — clínica METABRIC, solo muestras con expresión
#   mtbrc_expr.rds           — matriz expresión METABRIC (genes x muestras)
#
# Datos descargados (fuera de data/):
#   GDCdata/                 — archivos crudos TCGA (TCGAbiolinks)
#   brca_metabric/           — archivos crudos METABRIC (cBioPortal)
# ============================================================

dir.create("data", showWarnings = FALSE)


# ============================================================
# SECCIÓN 1: TCGA-BRCA
# ============================================================

# ── 1.1  Descarga de datos clínicos ────────────────────────

query_clin <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)
GDCdownload(query_clin)
clinical_data <- GDCprepare(query_clin)

names(clinical_data)   # verificar clave antes de extraer
patient_table <- clinical_data[["clinical_patient_brca"]]


# ── 1.2  Subtipos PAM50 ─────────────────────────────────────

brca_subtypes <- TCGAquery_subtype(tumor = "BRCA")
table(brca_subtypes$BRCA_Subtype_PAM50, useNA = "always")


# ── 1.3  Cruzar PAM50 con estado TNBC por IHC ──────────────

patient_table$patient_id_short <- substr(patient_table$bcr_patient_barcode, 1, 12)

tcga_brca_subtypes <- merge(
  patient_table[, c("patient_id_short", "er_status_by_ihc",
                    "pr_status_by_ihc",  "her2_status_by_ihc")],
  brca_subtypes[, c("patient", "BRCA_Subtype_PAM50")],
  by.x = "patient_id_short",
  by.y = "patient"
)

tcga_brca_subtypes$is_TNBC <-
  tcga_brca_subtypes$er_status_by_ihc   %in% "Negative" &
  tcga_brca_subtypes$pr_status_by_ihc   %in% "Negative" &
  tcga_brca_subtypes$her2_status_by_ihc %in% "Negative"

table(tcga_brca_subtypes$BRCA_Subtype_PAM50, tcga_brca_subtypes$is_TNBC, useNA = "always")

saveRDS(tcga_brca_subtypes, "data/tcga_brca_subtypes.rds")
message("data/tcga_brca_subtypes.rds guardado")


# ── 1.4  Descarga RNA-seq STAR Counts ───────────────────────

message(format(Sys.time()), " — Iniciando descarga TCGA-BRCA RNA-seq")

query_rnaseq <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

message(format(Sys.time()), " — Query OK: ", nrow(getResults(query_rnaseq)), " archivos")

GDCdownload(query_rnaseq, files.per.chunk = 100)

message(format(Sys.time()), " — Descarga completa. Preparando SummarizedExperiment...")

tcga_brca_rnaseq_se <- GDCprepare(query_rnaseq)

message(format(Sys.time()), " — SE creado: ",
        nrow(tcga_brca_rnaseq_se), " genes × ", ncol(tcga_brca_rnaseq_se), " muestras")

saveRDS(tcga_brca_rnaseq_se, "data/tcga_brca_rnaseq_se.rds")
message("data/tcga_brca_rnaseq_se.rds guardado")


# ============================================================
# SECCIÓN 2: METABRIC
# ============================================================

# ── 2.1  Descarga desde cBioPortal ──────────────────────────

if (!file.exists("brca_metabric/data_clinical_patient.txt")) {
  download.file(
    url      = "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz",
    destfile = "brca_metabric.tar.gz",
    mode     = "wb"
  )
  untar("brca_metabric.tar.gz", exdir = ".")
  message("METABRIC extraído en brca_metabric/")
} else {
  message("brca_metabric/ already present — skipping download")
}


# ── 2.2  Clínica ─────────────────────────────────────────────

mtbrc_patient <- read.delim(
  "brca_metabric/data_clinical_patient.txt",
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  comment.char     = "#"
)

mtbrc_sample <- read.delim(
  "brca_metabric/data_clinical_sample.txt",
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  comment.char     = "#"
)

mtbrc_clinical <- merge(
  mtbrc_patient,
  mtbrc_sample[, c("PATIENT_ID", "SAMPLE_ID", "ER_STATUS", "HER2_STATUS", "PR_STATUS")],
  by    = "PATIENT_ID",
  all.x = TRUE
)

table(mtbrc_clinical$CLAUDIN_SUBTYPE, useNA = "always")
table(mtbrc_clinical$INTCLUST,        useNA = "always")

saveRDS(mtbrc_clinical, "data/mtbrc_clinical.rds")
message("data/mtbrc_clinical.rds guardado — ", nrow(mtbrc_clinical), " pacientes")


# ── 2.3  Clínica restringida a muestras con expresión ────────

mrna_ids_raw  <- readLines("brca_metabric/case_lists/cases_mRNA.txt")
mrna_ids_line <- grep("^case_list_ids:", mrna_ids_raw, value = TRUE)
mrna_ids      <- strsplit(sub("^case_list_ids:\\s*", "", mrna_ids_line), "\t")[[1]]

mtbrc_clinical_mrna <- mtbrc_clinical[mtbrc_clinical$SAMPLE_ID %in% mrna_ids, ]

table(mtbrc_clinical_mrna$CLAUDIN_SUBTYPE, useNA = "always")

saveRDS(mtbrc_clinical_mrna, "data/mtbrc_clinical_mrna.rds")
message("data/mtbrc_clinical_mrna.rds guardado — ", nrow(mtbrc_clinical_mrna), " muestras con expresión")


# ── 2.4  Matriz de expresión (microarray Illumina) ───────────

message("Leyendo data_mrna_illumina_microarray.txt (~658 MB, ~1 min)...")

mtbrc_expr_raw <- read.delim(
  "brca_metabric/data_mrna_illumina_microarray.txt",
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# Convertir a matriz numérica (columnas 3 en adelante son muestras)
expr_mat <- as.matrix(mtbrc_expr_raw[, -(1:2)])
rownames(expr_mat) <- mtbrc_expr_raw$Hugo_Symbol

# Resolver duplicados: colapsar por mediana de columnas
dup_genes <- unique(rownames(expr_mat)[duplicated(rownames(expr_mat))])
message(length(dup_genes), " símbolos duplicados — colapsando por mediana")

df_expr <- as.data.frame(expr_mat)
df_expr$.symbol <- rownames(expr_mat)
mtbrc_expr <- do.call(rbind,
  lapply(
    split(df_expr[, colnames(df_expr) != ".symbol"], df_expr$.symbol),
    function(chunk) apply(chunk, 2, median, na.rm = TRUE)
  )
)

message("Dimensiones finales: ", nrow(mtbrc_expr), " genes × ", ncol(mtbrc_expr), " muestras")

saveRDS(mtbrc_expr, "data/mtbrc_expr.rds")
message("data/mtbrc_expr.rds guardado")


message("\n====  MÓDULO BASE BRCA COMPLETO  ====")
message("Salidas en data/:")
message("  tcga_brca_subtypes.rds")
message("  tcga_brca_rnaseq_se.rds")
message("  mtbrc_clinical.rds")
message("  mtbrc_clinical_mrna.rds")
message("  mtbrc_expr.rds")
