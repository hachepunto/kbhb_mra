library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

if (!interactive()) setwd(normalizePath("."))
dir.create("data",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# Figure 1, panel C — redesigned heatmap (editor's major-revision
# request) + final Figure 1 assembly (A + B + C).
#
# This script does NOT recompute the pipeline. The de novo cluster,
# composite score, and survival/Cox models are computed exactly once,
# in basal_subtype_clinical.R, and saved to data/*.tsv|*.rds — this
# script only reads those outputs, so the two pieces of the pipeline
# can never silently diverge (e.g. from a different seed or package
# version). Panels A and B (NES scatter, meta-NES lollipop) are
# likewise reloaded from the .rds objects saved by compare_kbhb_mrs.R.
#
# What IS new here: fresh clinical annotation tracks (PAM50, IntClust,
# histological type, grade, age, ER/PR/HER2 status) pulled directly
# from source, plus the additional receptor-status-vs-de-novo-cluster
# association analysis requested by the editor.
# ============================================================

# ============================================================
# 0. Reuse already-computed results (read-only)
# ============================================================

mrs_df  <- read.delim("data/kbhb_mrs_comparison.tsv", stringsAsFactors = FALSE)
sig_tfs <- mrs_df$regulator[mrs_df$FDR_meta < 0.05]   # 7 significant TMRs

meta_act <- readRDS("data/kbhb_metaviper_activity.rds")   # 817 TF x 404 samples, z-scored
stopifnot(all(sig_tfs %in% rownames(meta_act)))

tcga_basal_expr  <- readRDS("data/tcga_basal_expr.rds")   # only used for colnames() -> cohort membership
mtbrc_basal_expr <- readRDS("data/mtbrc_basal_expr.rds")
tcga_ids  <- colnames(tcga_basal_expr)
mtbrc_ids <- colnames(mtbrc_basal_expr)
cohort_of <- setNames(ifelse(colnames(meta_act) %in% tcga_ids, "TCGA", "METABRIC"), colnames(meta_act))

clust_df <- read.delim("data/basal_denovo_clusters.tsv", stringsAsFactors = FALSE)
cl_named <- setNames(clust_df$cluster, clust_df$sample_id)
k_best   <- length(unique(clust_df$cluster))

fig1 <- readRDS("data/fig1_nes_scatter.rds")
fig2 <- readRDS("data/fig2_lollipop_metanes.rds")

# ============================================================
# 1. Fresh clinical annotation data (inspected before harmonizing)
# ============================================================

# ---- 1a. METABRIC: everything needed is already in
#      data/mtbrc_basal_clinical_extra.tsv (built by basal_subtype_clinical.R,
#      section 1a-METABRIC). Read-only here as well. ----

mtbrc_clin <- read.delim("data/mtbrc_basal_clinical_extra.tsv", stringsAsFactors = FALSE)
rownames(mtbrc_clin) <- mtbrc_clin$sample_id
mtbrc_clin <- mtbrc_clin[mtbrc_ids, ]

cat("\n================ METABRIC — raw clinical values (pre-harmonization) ================\n")
cat("\n-- INTCLUST --\n");             print(table(mtbrc_clin$INTCLUST, useNA = "always"))
cat("\n-- CLAUDIN_SUBTYPE --\n");      print(table(mtbrc_clin$CLAUDIN_SUBTYPE, useNA = "always"))
cat("\n-- HISTOLOGICAL_SUBTYPE --\n"); print(table(mtbrc_clin$HISTOLOGICAL_SUBTYPE, useNA = "always"))
cat("\n-- GRADE --\n");                print(table(mtbrc_clin$GRADE, useNA = "always"))
cat("\n-- ER_STATUS --\n");            print(table(mtbrc_clin$ER_STATUS, useNA = "always"))
cat("\n-- PR_STATUS --\n");            print(table(mtbrc_clin$PR_STATUS, useNA = "always"))
cat("\n-- HER2_STATUS --\n");          print(table(mtbrc_clin$HER2_STATUS, useNA = "always"))
cat("\n-- AGE_AT_DIAGNOSIS summary --\n"); print(summary(mtbrc_clin$AGE_AT_DIAGNOSIS))
# No "claudin-low" category observed in this Basal-selected subset (100%
# "Basal" in CLAUDIN_SUBTYPE) -> no extra PAM50 category needed for METABRIC.

# ---- 1b. TCGA ----
#
# TCGA clinical data is unavoidably spread across multiple releases (there is
# no single source covering everything we need at full quality) -- so each
# field below is sourced independently, decided strictly by coverage/fitness
# for that specific field, not by convenience. Candidate sources compared:
# (a) `cd`, the colData of data/tcga_brca_rnaseq_se.rds -- already on disk,
# zero new downloads, and literally the table whose PanCanAtlas-derived
# "paper_*" columns defined Basal cohort membership in basal_pre_networks.R
# (line 65); (b) TCGAquery_subtype() ("bs" in basal_subtype_clinical.R) --
# same underlying 2018 PanCanAtlas source as (a), just fetched independently;
# (c) cBioPortal `brca_tcga_pub`, the TCGA Nature 2012 release
# (doi:10.1038/nature11412) -- tried and REVERTED for PAM50, see below;
# (d) the GDC Clinical Supplement (BCR Biotab).
#
# - PAM50 -> (a) cd$paper_BRCA_Subtype_PAM50 -- SAME source/criterion
#   basal_pre_networks.R already used to define Basal cohort membership for
#   ARACNe/msVIPER. An earlier version of this script used brca_tcga_pub
#   (2012) instead, on the reasoning that cd's PAM50 is tautologically 100%
#   "Basal" (zero information) while brca_tcga_pub showed 2 discordant
#   "Luminal A" calls. That trade turned out not to be worth it: brca_tcga_pub
#   only matches 141/195 of our cohort, and -- critically -- 45 of those 141
#   matched patients ALSO have no PAM50 value recorded in that 2012 release
#   (it is not exhaustively annotated even for patients it does include), for
#   99/195 total NA. That is real, uninformative missingness from a source
#   that never had anything to do with how this cohort was actually built --
#   not a meaningful discordance signal. Using cd instead makes this track
#   100% "Basal-like" (mapped from cd's "Basal" label) with zero NA, which is
#   the correct, internally-consistent picture: the figure should show
#   exactly the cohort-definition criterion the study actually used. The "2
#   discordant Luminal A" observation is dropped as a reported finding
#   accordingly -- heterogeneity within the nominal Basal cohort is still
#   evidenced (more robustly) by IntClust (both cohorts) and by the
#   receptor-status association + attenuated Kbhb activity in METABRIC's C1.
# - ER/PR/HER2 status -> (d) GDC BCR Biotab. Neither `cd` nor `bs` (both
#   PanCanAtlas-derived) carry ER/PR/HER2 in any naming variant in the
#   installed TCGAbiolinks version -- confirmed by inspecting colnames(cd)
#   and colnames(bs) directly, there is nothing to fall back on there.
# - Histological type -> (d) GDC BCR Biotab (~194/195, effectively full
#   coverage). cd$paper_BRCA_Pathology covers only 153/195 (78%; the rest are
#   the literal string "NA", not true NA) -- clearly worse.
# - Age -> (a) cd$paper_age_at_initial_pathologic_diagnosis. Same 100%
#   coverage as BCR Biotab's age_at_diagnosis, but already on disk (no extra
#   download), and it is the identical field basal_subtype_clinical.R's Cox
#   analysis already draws from via `bs` -- using it here keeps the age value
#   consistent across both scripts by construction.
# - Grade -> stays NA for TCGA. cd$paper_Tumor_Grade has 195/195 "non-missing"
#   cells, but every one of them is the literal string "NA" (not true NA) --
#   i.e. BRCA histological grade is simply not centrally reported in the
#   PanCanAtlas source for this tumor type.
#
# IntClust NA (grey) for the 4 TCGA samples iC10 couldn't classify (no usable
# copy-number segment data) is the only remaining expected TCGA missingness
# in this heatmap's annotation tracks.

suppressMessages(library(TCGAbiolinks))

cd_tcga <- as.data.frame(colData(readRDS("data/tcga_brca_rnaseq_se.rds")))
m_cd <- match(tcga_ids, cd_tcga$sample)
stopifnot(all(!is.na(m_cd)))   # cd$sample (16-char barcode) covers the full cohort by construction

query_clin <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical",
                        data.type = "Clinical Supplement", data.format = "BCR Biotab")
GDCdownload(query_clin)
gdc_biotab <- GDCprepare(query_clin)[["clinical_patient_brca"]]
gdc_biotab <- gdc_biotab[grepl("^TCGA-", gdc_biotab$bcr_patient_barcode), ]   # drop CDE_ID header remnant rows

m_gdc <- match(substr(tcga_ids, 1, 12), gdc_biotab$bcr_patient_barcode)

# HER2 "final" call: IHC first-line (0/1+ -> Negative, 3+ -> Positive);
# FISH reflex only when IHC was Equivocal/[Not Evaluated]/Indeterminate --
# the standard clinical algorithm, not IHC-only.
her2_ihc  <- trimws(gdc_biotab$her2_status_by_ihc[m_gdc])
her2_fish <- trimws(gdc_biotab$her2_fish_status[m_gdc])
her2_source <- ifelse(her2_ihc %in% c("Positive", "Negative"), "resolved_by_IHC",
                ifelse(her2_fish %in% c("Positive", "Negative"), "resolved_by_FISH_reflex",
                       "unresolved_NA"))
her2_final <- ifelse(her2_ihc == "Positive", "Positive",
               ifelse(her2_ihc == "Negative", "Negative",
               ifelse(her2_fish == "Positive", "Positive",
               ifelse(her2_fish == "Negative", "Negative", NA_character_))))

cat("\n================ TCGA — raw clinical values (pre-harmonization) ================\n")
cat("\n-- cd$paper_BRCA_Subtype_PAM50 (same source as basal_pre_networks.R cohort definition) --\n")
print(table(cd_tcga$paper_BRCA_Subtype_PAM50[m_cd], useNA = "always"))
cat("\nMatched to GDC BCR Biotab (ER/PR/HER2/histology, full cohort):", sum(!is.na(m_gdc)), "of", length(tcga_ids), "\n")
cat("\n-- er_status_by_ihc --\n");   print(table(gdc_biotab$er_status_by_ihc[m_gdc], useNA = "always"))
cat("\n-- pr_status_by_ihc --\n");   print(table(gdc_biotab$pr_status_by_ihc[m_gdc], useNA = "always"))
cat("\n-- her2_status_by_ihc --\n"); print(table(her2_ihc, useNA = "always"))
cat("\n-- her2_fish_status --\n");   print(table(her2_fish, useNA = "always"))
cat("\n-- HER2 resolution path (IHC vs FISH reflex) --\n"); print(table(her2_source, useNA = "always"))
cat("\n-- HER2 final status --\n"); print(table(her2_final, useNA = "always"))
cat("\n-- histological_type --\n");       print(table(gdc_biotab$histological_type[m_gdc], useNA = "always"))
cat("\nMatched to cd (colData, age only, zero new downloads):", sum(!is.na(m_cd)), "of", length(tcga_ids), "\n")
cat("\n-- cd$paper_age_at_initial_pathologic_diagnosis summary --\n")
print(summary(cd_tcga$paper_age_at_initial_pathologic_diagnosis[m_cd]))
cat("\n-- cd$paper_Tumor_Grade (confirms literal \"NA\" string, not true missingness) --\n")
print(table(cd_tcga$paper_Tumor_Grade[m_cd], useNA = "always"))

tcga_clin <- data.frame(
  sample_id         = tcga_ids,
  pam50             = ifelse(cd_tcga$paper_BRCA_Subtype_PAM50[m_cd] %in% "Basal", "Basal-like", NA_character_),
  er_status         = gdc_biotab$er_status_by_ihc[m_gdc],
  pr_status         = gdc_biotab$pr_status_by_ihc[m_gdc],
  her2_status       = her2_final,
  histological_type = gdc_biotab$histological_type[m_gdc],
  age               = cd_tcga$paper_age_at_initial_pathologic_diagnosis[m_cd],
  grade             = NA_character_,
  row.names         = NULL, stringsAsFactors = FALSE
)
write.table(tcga_clin, "data/tcga_basal_clinical_extra.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/tcga_basal_clinical_extra.tsv saved")

# ============================================================
# 2. Harmonize clinical variables across cohorts
# ============================================================

harm_binary_status <- function(x) {
  x <- trimws(x)
  out <- rep(NA_character_, length(x))
  out[x %in% "Positive"] <- "Positive"
  out[x %in% "Negative"] <- "Negative"
  out   # anything else (Indeterminate, Not Performed, Equivocal, [Not Available], NA) -> NA
}

harm_histology <- function(x, idc, ilc) {
  x <- trimws(x)
  out <- rep("Other", length(x))
  out[x %in% idc] <- "IDC"
  out[x %in% ilc] <- "ILC"
  out[is.na(x) | x %in% c("", "[Not Available]")] <- NA_character_
  out
}

tcga_pam50    <- setNames(tcga_clin$pam50, tcga_clin$sample_id)
tcga_er       <- setNames(harm_binary_status(tcga_clin$er_status),   tcga_clin$sample_id)
tcga_pr       <- setNames(harm_binary_status(tcga_clin$pr_status),   tcga_clin$sample_id)
tcga_her2     <- setNames(harm_binary_status(tcga_clin$her2_status), tcga_clin$sample_id)
tcga_hist     <- setNames(
  harm_histology(tcga_clin$histological_type,
                 idc = "Infiltrating Ductal Carcinoma",
                 ilc = "Infiltrating Lobular Carcinoma"),
  tcga_clin$sample_id)
tcga_grade    <- setNames(rep(NA_character_, nrow(tcga_clin)), tcga_clin$sample_id)
tcga_age      <- setNames(tcga_clin$age, tcga_clin$sample_id)

# TCGA has no public IntClust call, but tcga_ic10_classification.R applies
# the published iC10 classifier (Ali et al. 2014) directly to TCGA copy-
# number data, giving an IntClust-equivalent group per sample (191/195
# classified; read-only here, not recomputed). The 4 unclassified samples
# (no usable CN segment data) stay NA/grey, same as any other missing field.
ic10_calls    <- read.delim("data/tcga_basal_ic10_calls.tsv", stringsAsFactors = FALSE)
tcga_intclust <- setNames(rep(NA_character_, nrow(tcga_clin)), tcga_clin$sample_id)
tcga_intclust[ic10_calls$sample_id] <- ic10_calls$iC10
cat("\nTCGA IntClust (via iC10), matched:", sum(!is.na(tcga_intclust)), "of", length(tcga_intclust), "\n")

mtbrc_pam50    <- setNames(ifelse(mtbrc_clin$CLAUDIN_SUBTYPE %in% "Basal", "Basal-like", NA), mtbrc_clin$sample_id)
mtbrc_er       <- setNames(harm_binary_status(mtbrc_clin$ER_STATUS),   mtbrc_clin$sample_id)
mtbrc_pr       <- setNames(harm_binary_status(mtbrc_clin$PR_STATUS),   mtbrc_clin$sample_id)
mtbrc_her2     <- setNames(harm_binary_status(mtbrc_clin$HER2_STATUS), mtbrc_clin$sample_id)
mtbrc_hist     <- setNames(
  harm_histology(mtbrc_clin$HISTOLOGICAL_SUBTYPE, idc = "Ductal/NST", ilc = "Lobular"),
  mtbrc_clin$sample_id)
mtbrc_grade    <- setNames(as.character(mtbrc_clin$GRADE), mtbrc_clin$sample_id)
mtbrc_age      <- setNames(mtbrc_clin$AGE_AT_DIAGNOSIS, mtbrc_clin$sample_id)
mtbrc_intclust <- setNames(mtbrc_clin$INTCLUST, mtbrc_clin$sample_id)

samples <- colnames(meta_act)
pam50_v    <- c(tcga_pam50,    mtbrc_pam50)[samples]
er_v       <- c(tcga_er,       mtbrc_er)[samples]
pr_v       <- c(tcga_pr,       mtbrc_pr)[samples]
her2_v     <- c(tcga_her2,     mtbrc_her2)[samples]
hist_v     <- c(tcga_hist,     mtbrc_hist)[samples]
grade_v    <- c(tcga_grade,    mtbrc_grade)[samples]
age_v      <- c(tcga_age,      mtbrc_age)[samples]
intclust_v <- c(tcga_intclust, mtbrc_intclust)[samples]
cohort_v   <- cohort_of[samples]
denovo_v   <- cl_named[samples]

# ============================================================
# 3. Build HeatmapAnnotation (10 tracks, column = patient/sample)
# ============================================================

# PAM50 / IntClust colors sampled from the legend pixels of the original
# publications, not invented:
#   TCGA Nature 2012 (doi:10.1038/nature11412), Fig. 2c
#   Curtis et al. Nature 2012 (doi:10.1038/nature10983), Fig. 2
pam50_colors <- c(
  "Basal-like"    = "#E21E26",
  "HER2-enriched" = "#F69999",
  "Luminal A"     = "#2178B4",
  "Luminal B"     = "#A6CEE2",
  "Normal-like"   = "#66A943"
)
# IntClust 4 is split by ER status in Curtis et al.; both METABRIC subgroups
# keep the base IntClust-4 hue (#27C2CD) but at different shades so they stay
# visually distinguishable in the legend while signaling they belong to the
# same cluster. TCGA's IntClust comes from applying the iC10 classifier
# directly (see tcga_ic10_classification.R), which does not split group 4 by
# ER status -- it returns plain "4", so that level is mapped onto the same
# base IntClust-4 hue used for "4ER-" (not a new color; same 10-color
# palette shared by both cohorts, as required).
intclust_colors <- c(
  "1" = "#F75822", "2" = "#66C26B", "3" = "#D2327B", "4" = "#27C2CD",
  "4ER-" = "#27C2CD", "4ER+" = "#1B8790", "5" = "#8B171A", "6" = "#FCF346",
  "7" = "#34469F", "8" = "#FDAB19", "9" = "#D18CC0", "10" = "#774DA8"
)
denovo_colors <- setNames(brewer.pal(max(3, k_best), "Dark2")[seq_len(k_best)], sort(unique(denovo_v)))
cohort_colors <- c(TCGA = "#E69F00", METABRIC = "#009E73")
hist_colors   <- c(IDC = "#66C2A5", ILC = "#FC8D62", Other = "#8DA0CB")
# Sequential blues (increasing darkness = higher grade); TCGA is 100% NA for
# this field (grey85, same as every other NA) -- only METABRIC has real
# values here, and they are heavily skewed (Grade 1 n=2, Grade 2 n=17,
# Grade 3 n=187), so Grade 1 is nearly invisible in the heatmap by design,
# not a rendering issue.
grade_colors  <- c("1" = "#C6DBEF", "2" = "#6BAED6", "3" = "#08519C")
er_colors     <- c(Positive = "#D95F02", Negative = "#7570B3")
pr_colors     <- c(Positive = "#E7298A", Negative = "#66A61E")
her2_colors   <- c(Positive = "#E6AB02", Negative = "#A6761D")
age_col_fun   <- colorRamp2(
  seq(min(age_v, na.rm = TRUE), max(age_v, na.rm = TRUE), length.out = 5),
  c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")   # viridis, 5 stops
)

top_anno <- HeatmapAnnotation(
  Cohort            = cohort_v,
  PAM50             = pam50_v,
  IntClust          = intclust_v,
  `De novo cluster` = denovo_v,
  `Histological type` = hist_v,
  Grade             = grade_v,
  Age               = age_v,
  ER                = er_v,
  PR                = pr_v,
  HER2              = her2_v,
  col = list(
    Cohort              = cohort_colors,
    PAM50                = pam50_colors,
    IntClust             = intclust_colors,
    `De novo cluster`    = denovo_colors,
    `Histological type`  = hist_colors,
    Grade                = grade_colors,
    Age                  = age_col_fun,
    ER                   = er_colors,
    PR                   = pr_colors,
    HER2                 = her2_colors
  ),
  na_col              = "grey85",
  annotation_name_side = "left",
  annotation_name_gp  = gpar(fontsize = 9),
  simple_anno_size    = unit(0.6, "cm"),
  annotation_legend_param = list(
    Cohort             = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    PAM50              = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    IntClust           = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    `De novo cluster`  = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    `Histological type`= list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    Grade              = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    Age                = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9)),
    ER                 = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    PR                 = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm")),
    HER2               = list(title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), grid_height = unit(0.45, "cm"))
  ),
  gap                 = unit(1, "mm")
)

# ============================================================
# 4. Main heatmap — same divergent palette/clustering as before, now with
#    columns split into 4 fixed visual blocks (Cohort x De novo cluster) so
#    TCGA and METABRIC sit side by side instead of interleaved by the joint
#    clustering. The split is purely visual: it uses the Cohort and De novo
#    cluster labels already computed and validated (IntClust/receptor tests
#    above), it does not recompute or reassign them. cluster_column_slices =
#    FALSE keeps the 4 blocks fixed in this order (TCGA C1/C2, METABRIC
#    C1/C2) rather than letting ComplexHeatmap reorder the blocks; within
#    each block, columns are still free to cluster locally (same distance/
#    method as before) purely for visual sample ordering. Rows (the 7 TMRs)
#    are NOT split -- a single shared row order/dendrogram across all 4
#    blocks, so co-activity patterns are directly comparable cohort to
#    cohort.
# ============================================================

column_split_df <- data.frame(
  Cohort  = factor(cohort_v[samples], levels = c("TCGA", "METABRIC")),
  Cluster = factor(denovo_v[samples], levels = sort(unique(denovo_v)))
)

ht <- Heatmap(
  meta_act[sig_tfs, samples],
  name                     = "Activity\n(NES)",
  col                      = colorRampPalette(c("#2980B9", "white", "#C0392B"))(100),
  clustering_distance_columns = function(m) as.dist(1 - cor(t(m))),
  clustering_method_columns   = "ward.D2",
  cluster_rows             = TRUE,
  show_column_names         = FALSE,
  top_annotation            = top_anno,
  row_names_gp              = gpar(fontsize = 9),
  column_split             = column_split_df,
  cluster_column_slices    = FALSE,
  column_title_gp          = gpar(fontsize = 10, fontface = "bold"),
  column_gap               = unit(2, "mm"),
  heatmap_legend_param      = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))
)

panel_c_grob <- grid.grabExpr(
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE),
  wrap.grobs = TRUE   # self-contained grob: otherwise it depends on named viewports
                      # that no longer exist once replayed inside patchwork's own layout
)

# ============================================================
# 5. Assemble Figure 1 (A + B + C)
# ============================================================

panel_c <- patchwork::wrap_elements(full = panel_c_grob)

figure1 <- (fig1 | fig2) / panel_c +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(
    caption    = "Huang et al. 2021 gene set | ARACNe-AP + msVIPER | TCGA + METABRIC",
    tag_levels = "A",
    theme      = theme(
      plot.tag = element_text(face = "bold", size = 16)
    )
  )

ggsave("figures/Figure1.pdf", figure1, width = 16, height = 18)
message("figures/Figure1.pdf saved")

# ============================================================
# 6. Additional stats requested by the editor: de novo cluster vs.
#    ER/PR/HER2 status, per cohort (same pattern as the established-
#    subtype association already computed in basal_subtype_clinical.R)
# ============================================================

fisher_or_chisq <- function(tab_x, tab_y, label) {
  ok <- !is.na(tab_x) & !is.na(tab_y) & nzchar(tab_y)
  tt <- table(tab_x[ok], tab_y[ok])
  cat("\n--", label, "-- (n =", sum(ok), ")\n"); print(tt)
  if (nrow(tt) < 2 || ncol(tt) < 2) return(data.frame(test = label, n = sum(ok), p = NA, method = "n/a"))
  ft <- tryCatch(fisher.test(tt, simulate.p.value = TRUE, B = 10000),
                 error = function(e) NULL)
  if (!is.null(ft)) return(data.frame(test = label, n = sum(ok), p = ft$p.value, method = "Fisher (simulated)"))
  ct <- chisq.test(tt)
  data.frame(test = label, n = sum(ok), p = ct$p.value, method = "Chi-squared")
}

contingency_long <- function(tab_x, tab_y, label) {
  ok <- !is.na(tab_x) & !is.na(tab_y) & nzchar(tab_y)
  tt <- table(cluster = tab_x[ok], status = tab_y[ok])
  df <- as.data.frame(tt)
  cbind(test = label, df)
}

cat("\n=== de novo cluster <-> ER/PR/HER2 status: contingency tables ===\n")
receptor_assoc <- rbind(
  fisher_or_chisq(denovo_v[tcga_ids],  er_v[tcga_ids],    "TCGA: de novo cluster vs ER status"),
  fisher_or_chisq(denovo_v[tcga_ids],  pr_v[tcga_ids],    "TCGA: de novo cluster vs PR status"),
  fisher_or_chisq(denovo_v[tcga_ids],  her2_v[tcga_ids],  "TCGA: de novo cluster vs HER2 status"),
  fisher_or_chisq(denovo_v[mtbrc_ids], er_v[mtbrc_ids],   "METABRIC: de novo cluster vs ER status"),
  fisher_or_chisq(denovo_v[mtbrc_ids], pr_v[mtbrc_ids],   "METABRIC: de novo cluster vs PR status"),
  fisher_or_chisq(denovo_v[mtbrc_ids], her2_v[mtbrc_ids], "METABRIC: de novo cluster vs HER2 status")
)
receptor_assoc$FDR <- p.adjust(receptor_assoc$p, method = "BH")
cat("\n=== de novo cluster <-> ER/PR/HER2 status (editor-requested) ===\n")
print(receptor_assoc, row.names = FALSE)
write.table(receptor_assoc, "data/Supplementary_TableS4_denovo_cluster_receptor_status.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/Supplementary_TableS4_denovo_cluster_receptor_status.tsv saved (", nrow(receptor_assoc), " rows)")

receptor_contingency <- rbind(
  contingency_long(denovo_v[tcga_ids],  er_v[tcga_ids],    "TCGA: de novo cluster vs ER status"),
  contingency_long(denovo_v[tcga_ids],  pr_v[tcga_ids],    "TCGA: de novo cluster vs PR status"),
  contingency_long(denovo_v[tcga_ids],  her2_v[tcga_ids],  "TCGA: de novo cluster vs HER2 status"),
  contingency_long(denovo_v[mtbrc_ids], er_v[mtbrc_ids],   "METABRIC: de novo cluster vs ER status"),
  contingency_long(denovo_v[mtbrc_ids], pr_v[mtbrc_ids],   "METABRIC: de novo cluster vs PR status"),
  contingency_long(denovo_v[mtbrc_ids], her2_v[mtbrc_ids], "METABRIC: de novo cluster vs HER2 status")
)
write.table(receptor_contingency, "data/Supplementary_TableS4b_denovo_cluster_receptor_status_contingency.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/Supplementary_TableS4b_denovo_cluster_receptor_status_contingency.tsv saved (", nrow(receptor_contingency), " rows)")

# ============================================================
# 7. METABRIC IntClust <-> de novo cluster (formalizing the by-hand check):
#    full 10(+)-group table (reference; several low-n cells) and the
#    collapsed, well-powered 2x2 -- IntClust10 (the dominant, genomically
#    "pure Basal" group, 143/209 = 68% of the cohort) vs every other IntClust
#    group combined. Same pattern as tcga_ic10_classification.R's
#    IntClust{8,9,10}-vs-rest collapse, but the majority bucket here is just
#    IntClust10 -- METABRIC's IntClust9 splits close to 50/50 between C1/C2
#    (8 vs 8), so it does not belong with the canonical/majority bucket the
#    way TCGA's IntClust9 does.
# ============================================================

mtbrc_intclust_v <- mtbrc_clin$INTCLUST[match(mtbrc_ids, mtbrc_clin$sample_id)]
mtbrc_cluster_v  <- denovo_v[mtbrc_ids]

tt_mtbrc_full <- table(IntClust = mtbrc_intclust_v, DeNovoCluster = mtbrc_cluster_v)
cat("\n=== [reference only, low-n cells present] METABRIC IntClust x de novo cluster ===\n")
print(tt_mtbrc_full)
mtbrc_group_n <- table(mtbrc_intclust_v)
cat("group sizes:\n"); print(mtbrc_group_n)

mtbrc_full_out <- as.data.frame(tt_mtbrc_full)
names(mtbrc_full_out) <- c("IntClust", "DeNovoCluster", "n")
mtbrc_full_out$group_n <- mtbrc_group_n[mtbrc_full_out$IntClust]
mtbrc_full_out$low_n_group <- mtbrc_full_out$group_n <= 5
write.table(mtbrc_full_out, "data/Supplementary_TableS6b_metabric_ic10_denovo_cluster_full.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/Supplementary_TableS6b_metabric_ic10_denovo_cluster_full.tsv saved (", nrow(mtbrc_full_out), " rows)")

mtbrc_ic10_group <- ifelse(mtbrc_intclust_v == "10", "IntClust10", "Other IntClust")
tt_mtbrc_collapsed <- table(mtbrc_ic10_group, mtbrc_cluster_v)
cat("\n=== [reported result] METABRIC IntClust10 vs Other IntClust x de novo cluster ===\n")
print(tt_mtbrc_collapsed)

ft_mtbrc <- fisher.test(tt_mtbrc_collapsed)
ct_mtbrc <- suppressWarnings(chisq.test(tt_mtbrc_collapsed))
cat("\nFisher exact p:", ft_mtbrc$p.value, "\n")
cat("Chi-squared p:", ct_mtbrc$p.value, "(expected counts all > 5:", all(ct_mtbrc$expected > 5), ")\n")

mtbrc_collapsed_out <- as.data.frame(tt_mtbrc_collapsed)
names(mtbrc_collapsed_out) <- c("ic10_group", "DeNovoCluster", "n")
mtbrc_collapsed_out$fisher_p <- ft_mtbrc$p.value
mtbrc_collapsed_out$chisq_p  <- ct_mtbrc$p.value
write.table(mtbrc_collapsed_out, "data/Supplementary_TableS6_metabric_ic10_denovo_cluster.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/Supplementary_TableS6_metabric_ic10_denovo_cluster.tsv saved (", nrow(mtbrc_collapsed_out), " rows)")

# ============================================================
# 8. Mean NES per TMR per block (Cohort x De novo cluster) -- the same 4
#    blocks used to split the heatmap columns (column_split_df above),
#    formalizing the by-hand summary of "does the co-activity pattern look
#    similar across cohorts" into a reportable table.
# ============================================================

mean_nes_block <- do.call(rbind, lapply(sig_tfs, function(tf) {
  vals <- meta_act[tf, samples]
  data.frame(
    TMR      = tf,
    Cohort   = column_split_df$Cohort,
    Cluster  = column_split_df$Cluster,
    NES      = vals
  )
})) %>%
  group_by(TMR, Cohort, Cluster) %>%
  summarise(mean_NES = mean(NES), n = n(), .groups = "drop") %>%
  as.data.frame()

cat("\n=== Mean NES per TMR per block (Cohort x De novo cluster) ===\n")
print(mean_nes_block, row.names = FALSE)
write.table(mean_nes_block, "data/Supplementary_TableS7_denovo_cluster_mean_NES.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/Supplementary_TableS7_denovo_cluster_mean_NES.tsv saved (", nrow(mean_nes_block), " rows)")

# ============================================================
# 9. Supplementary Figure S4 -- % receptor-positive by de novo cluster, by
#    cohort (ggplot2 version of the ER/PR/HER2 rates already computed in
#    Supplementary_TableS4b above; read back in rather than recomputed, to
#    guarantee the figure and the table can never diverge)
# ============================================================

receptor_rates <- receptor_contingency %>%
  tidyr::separate(test, into = c("cohort", "receptor"), sep = ": de novo cluster vs ", remove = FALSE) %>%
  mutate(receptor = sub(" status", "", receptor)) %>%
  group_by(cohort, receptor, cluster) %>%
  mutate(pct = 100 * Freq / sum(Freq)) %>%
  ungroup() %>%
  filter(status == "Positive")

p_receptor_rate <- ggplot(receptor_rates, aes(x = cluster, y = pct, fill = receptor)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3) +
  facet_wrap(~cohort) +
  scale_fill_manual(values = c(ER = "#D95F02", PR = "#E7298A", HER2 = "#E6AB02")) +
  labs(x = "De novo TMR cluster", y = "% receptor-positive", fill = "Receptor") +
  ylim(0, max(receptor_rates$pct) * 1.2) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey90"))
# No embedded caption: panel/figure captions belong in the manuscript figure
# legend, not baked into the image (same convention as the rest of the repo).

ggsave("figures/Supplementary_FigureS4_receptor_rate_by_cluster.pdf", p_receptor_rate, width = 8, height = 5)
message("figures/Supplementary_FigureS4_receptor_rate_by_cluster.pdf saved")

message("\n====  figure1_panel_c.R COMPLETE  ====")
