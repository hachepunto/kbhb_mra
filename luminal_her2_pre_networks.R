library(SummarizedExperiment)
library(sva)

if (!interactive()) setwd(normalizePath("."))
dir.create("data", showWarnings = FALSE)

# ============================================================
# Extension of the MRA to Luminal A/B and HER2-enriched
#
# Same processing as basal_pre_networks.R (expression filters, ComBat in
# METABRIC, ARACNe-AP format), but for PAM50 LumA, LumB, Her2 instead of
# Basal. The SAME Normal group (TCGA-NT / METABRIC
# CLAUDIN_SUBTYPE=="Normal") is reused as the comparator for all 3
# signatures, exactly as was done for Basal.
#
# Outputs in data/ (per subtype x cohort):
#   tcga_{subtype}_matrix.txt / mtbrc_{subtype}_matrix.txt  — ARACNe-AP input
#   tcga_{subtype}_expr.rds   / mtbrc_{subtype}_expr.rds
#   tcga_{subtype}_normal_expr.rds / mtbrc_{subtype}_normal_expr.rds
# ============================================================

subtypes <- c(LumA = "LumA", LumB = "LumB", Her2 = "Her2")

write_aracne_matrix <- function(mat, path) {
  df <- cbind(gene = rownames(mat), as.data.frame(mat))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(path, " written — ", nrow(mat), " genes x ", ncol(mat), " samples")
}

# ============================================================
# SECTION 1: TCGA-BRCA — recompute the deduplicated protein-coding TPM
# matrix (same procedure as basal_pre_networks.R Section 1.1-1.5)
# ============================================================

tcga_outputs_exist <- all(file.exists(sprintf("data/tcga_%s_matrix.txt", names(subtypes))))

if (tcga_outputs_exist) {

  message("data/tcga_{", paste(names(subtypes), collapse=","), "}_matrix.txt already exist — skipping TCGA recompute")

} else {

message("Loading TCGA SummarizedExperiment...")
se <- readRDS("data/tcga_brca_rnaseq_se.rds")
cd <- as.data.frame(colData(se))
rd <- as.data.frame(rowData(se))

pc_idx <- rd$gene_type == "protein_coding"
se_pc  <- se[pc_idx, ]
rd_pc  <- rd[pc_idx, ]

tpm_all <- assay(se_pc, "tpm_unstrand")
rownames(tpm_all) <- rd_pc$gene_name

dup_syms <- rownames(tpm_all)[duplicated(rownames(tpm_all))]
df_tpm <- as.data.frame(tpm_all)
df_tpm$.sym <- rownames(tpm_all)
tpm_dedup <- do.call(rbind,
  lapply(
    split(df_tpm[, colnames(df_tpm) != ".sym"], df_tpm$.sym),
    function(chunk) apply(chunk, 2, median, na.rm = TRUE)
  )
)
colnames(tpm_dedup) <- substr(cd$barcode, 1, 16)
message("TCGA — deduplicated protein_coding genes: ", nrow(tpm_dedup))

normal_bc <- unique(substr(cd$barcode[cd$shortLetterCode %in% "NT"], 1, 16))
cat("TCGA — Normal NT samples (shared across subtypes):", length(normal_bc), "\n")

for (nm in names(subtypes)) {
  lvl <- subtypes[[nm]]
  sub_bc <- unique(substr(
    cd$barcode[cd$paper_BRCA_Subtype_PAM50 %in% lvl & cd$shortLetterCode %in% "TP"], 1, 16))
  cat("\n--- TCGA", nm, "---\n")
  cat("Tumor samples:", length(sub_bc), "\n")

  tcga_sub    <- log2(tpm_dedup[, sub_bc]     + 1)
  tcga_normal <- log2(tpm_dedup[, normal_bc]  + 1)

  expressed   <- rowSums(tcga_sub > 1) >= 0.2 * ncol(tcga_sub)
  tcga_sub    <- tcga_sub[expressed, ]
  tcga_normal <- tcga_normal[rownames(tcga_sub), ]
  cat("Genes after expression filter:", nrow(tcga_sub), "\n")

  saveRDS(tcga_sub,    sprintf("data/tcga_%s_expr.rds",        nm))
  saveRDS(tcga_normal, sprintf("data/tcga_%s_normal_expr.rds", nm))
  write_aracne_matrix(tcga_sub, sprintf("data/tcga_%s_matrix.txt", nm))
}

rm(se, se_pc, tpm_all, df_tpm, tpm_dedup); gc()

}

# ============================================================
# SECTION 2: METABRIC — reuse mtbrc_expr.rds (full matrix) + ComBat per
# subtype (same design as basal_pre_networks.R Section 2)
# ============================================================

message("\nLoading METABRIC clinical data and expression matrix...")
cl   <- readRDS("data/mtbrc_clinical_mrna.rds")
expr <- readRDS("data/mtbrc_expr.rds")

normal_ids <- cl$SAMPLE_ID[cl$CLAUDIN_SUBTYPE == "Normal"]
cat("METABRIC — Normal samples (shared across subtypes):", length(normal_ids), "\n")

for (nm in names(subtypes)) {
  lvl <- subtypes[[nm]]
  sub_ids <- cl$SAMPLE_ID[cl$CLAUDIN_SUBTYPE == lvl]
  cat("\n--- METABRIC", nm, "---\n")
  cat("Tumor samples:", length(sub_ids), "\n")

  mtbrc_sub    <- expr[, colnames(expr) %in% sub_ids]
  mtbrc_normal <- expr[, colnames(expr) %in% normal_ids]

  # rowSums(..., na.rm=TRUE): mtbrc_expr.rds has 16 scattered NAs (microarray);
  # without na.rm, any NA in a row turns the whole rowSums into NA and produces
  # a spurious NA-named row when indexing (seen in LumA: 8 genes affected).
  expressed    <- rowSums(mtbrc_sub > 6, na.rm = TRUE) >= 0.2 * ncol(mtbrc_sub)
  mtbrc_sub    <- mtbrc_sub[expressed, ]
  mtbrc_normal <- mtbrc_normal[rownames(mtbrc_sub), ]
  cat("Genes after expression filter:", nrow(mtbrc_sub), "\n")

  combined     <- cbind(mtbrc_sub, mtbrc_normal)
  group_labels <- c(rep("Tumor", ncol(mtbrc_sub)), rep("Normal", ncol(mtbrc_normal)))
  cohort_map   <- setNames(cl$COHORT, cl$SAMPLE_ID)
  batch_vec    <- cohort_map[colnames(combined)]

  cat("Cohort distribution in", nm, "+ Normal:\n")
  print(table(batch_vec, group_labels))

  mod <- model.matrix(~ group_labels)
  message("Applying ComBat (batch = COHORT, mod = ~group) for ", nm, "...")
  combined_corrected <- ComBat(dat = combined, batch = batch_vec, mod = mod)

  mtbrc_sub    <- combined_corrected[, group_labels == "Tumor"]
  mtbrc_normal <- combined_corrected[, group_labels == "Normal"]

  saveRDS(mtbrc_sub,    sprintf("data/mtbrc_%s_expr.rds",        nm))
  saveRDS(mtbrc_normal, sprintf("data/mtbrc_%s_normal_expr.rds", nm))
  write_aracne_matrix(mtbrc_sub, sprintf("data/mtbrc_%s_matrix.txt", nm))
}

message("\n====  PRE-PROCESSING COMPLETE (LumA / LumB / Her2)  ====")
