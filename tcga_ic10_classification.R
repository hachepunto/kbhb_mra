suppressMessages({
  library(TCGAbiolinks)
  library(iC10)
  library(iC10TrainingData)
})

if (!interactive()) setwd(normalizePath("."))
dir.create("data", showWarnings = FALSE)

# ============================================================
# Apply the published iC10 classifier (Ali et al. 2014, Genome Biology,
# doi:10.1186/s13059-014-0431-1 -- same group as Curtis et al. 2012 Nature,
# doi:10.1038/nature10983) to TCGA-BRCA Basal samples, using GDC-harmonized
# "Masked Copy Number Segment" data (Affymetrix SNP6, hg38). This is a
# trained classifier applied to new data, NOT de novo clustering -- directly
# comparable to the genomic IntClust groups already available for METABRIC
# (data/mtbrc_basal_clinical_extra.tsv, column INTCLUST).
#
# Motivation: in METABRIC, the de novo TMR cluster (data/basal_denovo_clusters.tsv)
# is almost perfectly separated by IntClust (non-canonical genomic subgroups
# concentrate in C1). TCGA has no public IntClust call, so the same check was
# not previously possible there -- this script closes that gap.
#
# This does NOT recompute the de novo cluster, composite score, or any other
# pipeline result -- it only reads data/basal_denovo_clusters.tsv (read-only)
# and adds a new, independent clinical/genomic annotation.
# ============================================================

tcga_basal_expr <- readRDS("data/tcga_basal_expr.rds")
tcga_ids <- colnames(tcga_basal_expr)   # same 195 samples used throughout the pipeline
cat("n TCGA Basal samples:", length(tcga_ids), "\n")

# ============================================================
# 1. Download Masked Copy Number Segment data (Affymetrix SNP6) for exactly
#    these 195 samples
# ============================================================

barcodes16 <- substr(tcga_ids, 1, 16)
q <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Copy Number Variation",
  data.type     = "Masked Copy Number Segment",
  barcode       = barcodes16
)
GDCdownload(q)
cn <- GDCprepare(q)
cat("Segment rows:", nrow(cn), "| unique samples with CN data:", length(unique(cn$Sample)),
    "of", length(tcga_ids), "\n")

# ============================================================
# 2. Build the segmented-data input iC10::matchFeatures() expects and
#    classify
# ============================================================

seg_df <- data.frame(
  ID              = substr(cn$Sample, 1, 16),
  chromosome_name = suppressWarnings(as.integer(cn$Chromosome)),   # drops X/Y -> NA
  loc.start       = cn$Start,
  loc.end         = cn$End,
  seg.mean        = cn$Segment_Mean
)
seg_df <- seg_df[!is.na(seg_df$chromosome_name), ]   # autosomes only (1-22); matches training features
seg_df <- seg_df[seg_df$ID %in% tcga_ids, ]

feat <- matchFeatures(CN = seg_df, CN.by.feat = "probe", ref = "hg38")
cat("Matched CN feature matrix:", nrow(feat$CN), "probes x", ncol(feat$CN), "samples\n")
# No normalizeFeatures() call: per the package docs, "no further normalization
# is needed on the copy number, as log2 ratios are comparable between platforms."

res <- iC10(feat)

ic10_calls <- data.frame(
  sample_id = colnames(feat$CN),
  iC10      = as.character(res$class),
  row.names = NULL
)

missing_ids <- setdiff(tcga_ids, ic10_calls$sample_id)
cat("\n=== Coverage ===\n")
cat("Classified:", nrow(ic10_calls), "of", length(tcga_ids), "TCGA Basal samples\n")
cat("Not classified (no usable CN segment data):", length(missing_ids), "\n")
if (length(missing_ids) > 0) print(missing_ids)

cat("\n=== iC10 group distribution (TCGA Basal cohort) ===\n")
print(table(ic10_calls$iC10, useNA = "always"))

write.table(ic10_calls, "data/tcga_basal_ic10_calls.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/tcga_basal_ic10_calls.tsv saved")

# ============================================================
# 3. Cross-tabulate against the existing de novo TMR cluster
# ============================================================

clust_df <- read.delim("data/basal_denovo_clusters.tsv", stringsAsFactors = FALSE)
cl_named <- setNames(clust_df$cluster, clust_df$sample_id)

cross_df <- merge(ic10_calls, data.frame(sample_id = names(cl_named), cluster = cl_named),
                   by = "sample_id")
cat("\nn samples with both iC10 call and de novo cluster:", nrow(cross_df), "\n")

# -- 3a. Full 10-group table: reference/transparency only. Several groups
#    have n<=5 (IntClust 3, 5, 7 here) -- NOT reliable individually, do not
#    report per-group percentages from this table as a standalone finding.
tt_full <- table(IntClust = cross_df$iC10, DeNovoCluster = cross_df$cluster)
cat("\n=== [reference only, low-n cells present] iC10 group x de novo cluster (TCGA Basal) ===\n")
print(tt_full)
cat("group sizes:\n"); print(table(cross_df$iC10))

full_out <- as.data.frame(tt_full)
names(full_out) <- c("IntClust", "DeNovoCluster", "n")
full_out$group_n <- table(cross_df$iC10)[full_out$IntClust]
full_out$low_n_group <- full_out$group_n <= 5
write.table(full_out, "data/Supplementary_TableS5b_tcga_ic10_denovo_cluster_full.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/Supplementary_TableS5b_tcga_ic10_denovo_cluster_full.tsv saved (", nrow(full_out), " rows; reference table -- flags low_n_group <= 5)")

# -- 3b. Collapsed, well-powered 2x2 table: IntClust {8,9,10} (the common,
#    "canonical" genomic groups, 87% of the cohort here) vs. all other
#    IntClust groups combined. This is the reported result.
cross_df$ic10_group <- ifelse(cross_df$iC10 %in% c("8", "9", "10"), "IntClust{8,9,10}", "Other IntClust")
tt_collapsed <- table(cross_df$ic10_group, cross_df$cluster)
cat("\n=== [reported result] IntClust{8,9,10} vs Other IntClust x de novo cluster ===\n")
print(tt_collapsed)

ft <- fisher.test(tt_collapsed)
ct <- suppressWarnings(chisq.test(tt_collapsed))
cat("\nFisher exact p:", ft$p.value, "\n")
cat("Chi-squared p:", ct$p.value, "(expected counts all > 5:", all(ct$expected > 5), ")\n")

collapsed_out <- as.data.frame(tt_collapsed)
names(collapsed_out) <- c("ic10_group", "DeNovoCluster", "n")
collapsed_out$fisher_p <- ft$p.value
collapsed_out$chisq_p  <- ct$p.value
write.table(collapsed_out, "data/Supplementary_TableS5_tcga_ic10_denovo_cluster.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/Supplementary_TableS5_tcga_ic10_denovo_cluster.tsv saved (", nrow(collapsed_out), " rows)")

message("\n====  tcga_ic10_classification.R COMPLETE  ====")
