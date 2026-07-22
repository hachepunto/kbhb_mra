suppressMessages(library(writexl))

if (!interactive()) setwd(normalizePath("."))
dir.create("data", showWarnings = FALSE)

# ============================================================
# Consolidates the already-generated Supplementary Table .tsv files (data/)
# into a single data/Supplementary_Tables.xlsx, one sheet per table, for
# submission. Purely a packaging step -- reads existing .tsv files, does not
# recompute or modify anything. Re-run the relevant upstream script
# (compare_kbhb_mrs.R, luminal_her2_mra_kbhb.R, figure1_panel_c.R,
# tcga_ic10_classification.R) first if a table needs updating; this script
# only picks up whatever is currently on disk.
#
# Uses writexl, NOT openxlsx: openxlsx 4.2.8.1 (current CRAN release, so
# there is no newer version to update to) unconditionally writes a broken
# <dimension ref="A1"/> tag (regardless of actual sheet size) and a dangling
# worksheet -> drawing/vmlDrawing relationship that points at files it never
# creates (xl/drawings/ is absent from the archive) -- reproduced with a
# zero-styling, zero-column-width minimal workbook, so it is not caused by
# setColWidths()/addStyle() here, it is unconditional on every saveWorkbook()
# call in this openxlsx build. That corrupted every sheet in an earlier
# version of this script's output: Excel/LibreOffice clipped the visible
# print area to the bogus dimension, and strict openpyxl (readonly=FALSE)
# raised KeyError looking for the missing drawing part. write_xlsx() does
# not add per-sheet drawing relationships at all and computes a correct
# <dimension> from the actual data, sidestepping both bugs entirely. Trade-
# off: no column auto-width (write_xlsx has no equivalent); format_headers
# gives a bold header row, which is the only formatting actually required
# here -- correctness over cosmetics.
#
# NOTE on the "Contents" sheet descriptions: these are drafted here from
# this repo's README.md (Result tables section), which is the only
# machine-readable description available in this repository. This script
# does NOT have access to the manuscript's own Supplementary Materials
# section text -- verify these descriptions against that section before
# submission and edit the `contents` table below if wording differs.
# ============================================================

tables <- list(
  S1  = "data/Supplementary_TableS1_ORA_GO_BP.tsv",
  S2  = "data/Supplementary_TableS2_ORA_Reactome.tsv",
  S3  = "data/Supplementary_TableS3_meta_NES_subtypes.tsv",
  S4  = "data/Supplementary_TableS4_denovo_cluster_receptor_status.tsv",
  S4b = "data/Supplementary_TableS4b_denovo_cluster_receptor_status_contingency.tsv",
  S5  = "data/Supplementary_TableS5_tcga_ic10_denovo_cluster.tsv",
  S5b = "data/Supplementary_TableS5b_tcga_ic10_denovo_cluster_full.tsv",
  S6  = "data/Supplementary_TableS6_metabric_ic10_denovo_cluster.tsv",
  S6b = "data/Supplementary_TableS6b_metabric_ic10_denovo_cluster_full.tsv",
  S7  = "data/Supplementary_TableS7_denovo_cluster_mean_NES.tsv"
)
stopifnot(all(file.exists(unlist(tables))))   # fail loudly if an upstream script hasn't been (re-)run yet

contents <- data.frame(
  Sheet       = names(tables),
  Description = c(
    S1  = "Full over-representation analysis (ORA) results: GO Biological Process terms enriched among the 7 Kbhb TMR regulons.",
    S2  = "Full over-representation analysis (ORA) results: Reactome pathways enriched among the 7 Kbhb TMR regulons.",
    S3  = "Meta-analysis NES of the Kbhb TMR panel across PAM50 subtypes (Basal-like, Luminal A, Luminal B, HER2-enriched), TCGA + METABRIC Stouffer combination.",
    S4  = "Editor-requested association: de novo Kbhb-TMR cluster vs. ER/PR/HER2 receptor status, per cohort (Fisher/chi-squared, BH-adjusted). Reported result.",
    S4b = "Supporting contingency counts (cluster x receptor status, per cohort) underlying Table S4.",
    S5  = "TCGA: iC10 genomic subtype (Ali et al. 2014) vs. de novo Kbhb-TMR cluster, collapsed to the well-powered IntClust{8,9,10}-vs-rest 2x2. Reported result.",
    S5b = "TCGA: iC10 genomic subtype vs. de novo Kbhb-TMR cluster, full 10-group table. Reference/transparency only -- several cells have n<=5 (see low_n_group column).",
    S6  = "METABRIC: IntClust vs. de novo Kbhb-TMR cluster, collapsed to the well-powered IntClust10-vs-rest 2x2. Reported result.",
    S6b = "METABRIC: IntClust vs. de novo Kbhb-TMR cluster, full (>=10-group) table. Reference/transparency only -- several cells have n<=5 (see low_n_group column).",
    S7  = "Mean Kbhb-TMR activity (NES) per TMR, per Cohort x de novo cluster block (the same 4 blocks split in Figure 1C)."
  ),
  stringsAsFactors = FALSE
)

sheet_data <- c(
  list(Contents = contents),
  lapply(tables, function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE))
)

write_xlsx(sheet_data, "data/Supplementary_Tables.xlsx", format_headers = TRUE)

cat("\n=== data/Supplementary_Tables.xlsx saved ===\n")
cat("Sheets (in order):", paste(names(sheet_data), collapse = ", "), "\n\n")
row_counts <- sapply(sheet_data, nrow)
print(data.frame(Sheet = names(row_counts), Rows = row_counts), row.names = FALSE)
message("\ndata/Supplementary_Tables.xlsx saved (", length(sheet_data), " sheets)")
