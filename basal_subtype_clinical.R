library(SummarizedExperiment)
library(survival)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(survminer)
library(patchwork)

if (!interactive()) setwd(normalizePath("."))
dir.create("data",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# TASK 1 — Editor: PAM50 Basal-like is not a single entity.
# Do the TF (TMRs) segregate subtypes within Basal-like? Are they
# associated with worse outcome?
#
# 1a. Correlation with already-established subtypes (TCGA: Lehmann/Bareche
#     TNBC; METABRIC: IntClust)
# 1b. Continuous score (priority) + de novo clustering (secondary)
# ============================================================

mrs_df  <- read.delim("data/kbhb_mrs_comparison.tsv", stringsAsFactors = FALSE)
mrs_sig <- mrs_df[mrs_df$FDR_meta < 0.05, ]
sig_tfs <- mrs_sig$regulator                                      # significant TMRs, derived from the data
tmr_sign <- setNames(sign(mrs_sig$NES_meta), mrs_sig$regulator)    # VEZF1 = -1, rest = +1
cat("Sign (direction) per TMR used for the composite score:\n"); print(tmr_sign)

meta_act <- readRDS("data/kbhb_metaviper_activity.rds")   # 817 TF x 404 samples, z-scored per TF and per cohort
stopifnot(all(sig_tfs %in% rownames(meta_act)))

tcga_basal_expr  <- readRDS("data/tcga_basal_expr.rds")   # log2(TPM+1), 14488 genes x 195 Basal TP samples
mtbrc_basal_expr <- readRDS("data/mtbrc_basal_expr.rds")  # ComBat-corrected, 11976 genes x 209 Basal samples

# ============================================================
# SECTION 1a-TCGA: Lehmann/Bareche TNBC subtype
# ============================================================
# Code and gene signatures vendored from BCTL-Bordet/TNBC_molecularsubtypes
# (Sotiriou lab), which reimplements:
#   Lehmann BD et al. J Clin Invest 2011;121(7):2750-67 (original signatures,
#     genes positively/negatively associated per subtype)
#   Bareche Y et al. Ann Oncol 2018;29(4):895-902 (reassigns "basal_like_2"
#     samples to their second-highest score -> 5 stable subtypes:
#     BL, IM, M, LAR, MSL; drops the "unstable" category from Lehmann's original)

source("external/tnbc_lehmann_bareche_BCTL-Bordet/Functions.R")
sig_load <- local({ e <- new.env(); load("external/tnbc_lehmann_bareche_BCTL-Bordet/lehmann.RData", envir = e); e$sig })

tcga_tnbc_subtype <- TNBCclassif(as.matrix(tcga_basal_expr), version = "bareche", sig = sig_load, coef = FALSE)
cat("\n=== TCGA — TNBC subtypes (Lehmann/Bareche) ===\n")
print(table(tcga_tnbc_subtype, useNA = "always"))

tcga_subtype_df <- data.frame(
  sample_id = names(tcga_tnbc_subtype),
  subtype   = as.character(tcga_tnbc_subtype),
  row.names = NULL
)
write.table(tcga_subtype_df, "data/tcga_basal_tnbc_subtype.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/tcga_basal_tnbc_subtype.tsv saved")

# ============================================================
# SECTION 1a-METABRIC: IntClust + grade/stage (already in the downloaded
# data; GRADE/TUMOR_STAGE live in data_clinical_sample.txt, which was not
# merged into mtbrc_clinical_mrna.rds by brca_tcga_mtbrc.R)
# ============================================================

mtbrc_clin <- readRDS("data/mtbrc_clinical_mrna.rds")
mtbrc_sample_extra <- read.delim("brca_metabric/data_clinical_sample.txt",
                                  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                  comment.char = "#")
mtbrc_clin <- merge(mtbrc_clin,
                     mtbrc_sample_extra[, c("SAMPLE_ID","GRADE","TUMOR_STAGE","TUMOR_SIZE")],
                     by = "SAMPLE_ID", all.x = TRUE)
rownames(mtbrc_clin) <- mtbrc_clin$SAMPLE_ID

mtbrc_clin_basal <- mtbrc_clin[colnames(mtbrc_basal_expr), ]
cat("\n=== METABRIC Basal — IntClust ===\n")
print(table(mtbrc_clin_basal$INTCLUST, useNA = "always"))
cat("\n=== METABRIC Basal — GRADE ===\n")
print(table(mtbrc_clin_basal$GRADE, useNA = "always"))
cat("\n=== METABRIC Basal — TUMOR_STAGE ===\n")
print(table(mtbrc_clin_basal$TUMOR_STAGE, useNA = "always"))

write.table(data.frame(sample_id = rownames(mtbrc_clin_basal), mtbrc_clin_basal, row.names = NULL),
            "data/mtbrc_basal_clinical_extra.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/mtbrc_basal_clinical_extra.tsv saved")

# ============================================================
# SECTION 1a: established subtype <-> TMR activity score correlation
# (per cohort, one TMR at a time; ANOVA + Kruskal-Wallis, BH)
# ============================================================

run_subtype_assoc <- function(activity_mat, sample_ids, subtype_vec, tfs, cohort_label) {
  subtype_vec <- subtype_vec[sample_ids]
  out <- lapply(tfs, function(tf) {
    nes <- activity_mat[tf, sample_ids]
    df  <- data.frame(nes = nes, grp = factor(subtype_vec))
    df  <- df[!is.na(df$grp), ]
    if (nlevels(droplevels(df$grp)) < 2) {
      return(data.frame(cohort = cohort_label, TMR = tf, n = nrow(df),
                         anova_F = NA, anova_p = NA, kw_stat = NA, kw_p = NA))
    }
    df$grp <- droplevels(df$grp)
    fit <- aov(nes ~ grp, data = df)
    a   <- summary(fit)[[1]]
    kw  <- kruskal.test(nes ~ grp, data = df)
    data.frame(cohort = cohort_label, TMR = tf, n = nrow(df),
               anova_F = a["grp", "F value"], anova_p = a["grp", "Pr(>F)"],
               kw_stat = unname(kw$statistic), kw_p = kw$p.value)
  })
  do.call(rbind, out)
}

tcga_ids  <- tcga_subtype_df$sample_id
tcga_subtype_named <- setNames(tcga_subtype_df$subtype, tcga_subtype_df$sample_id)
assoc_tcga <- run_subtype_assoc(meta_act, tcga_ids, tcga_subtype_named, sig_tfs, "TCGA_Lehmann_Bareche")

mtbrc_ids <- rownames(mtbrc_clin_basal)
mtbrc_intclust_named <- setNames(as.character(mtbrc_clin_basal$INTCLUST), rownames(mtbrc_clin_basal))
assoc_mtbrc <- run_subtype_assoc(meta_act, mtbrc_ids, mtbrc_intclust_named, sig_tfs, "METABRIC_IntClust")

assoc_all <- rbind(assoc_tcga, assoc_mtbrc)
assoc_all$anova_FDR <- p.adjust(assoc_all$anova_p, method = "BH")
assoc_all$kw_FDR    <- p.adjust(assoc_all$kw_p,    method = "BH")
cat("\n=== 1a: TMR activity ~ established subtype (ANOVA + Kruskal-Wallis, BH over 14 tests) ===\n")
print(assoc_all, row.names = FALSE)
write.table(assoc_all, "data/basal_subtype_activity_assoc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/basal_subtype_activity_assoc.tsv saved")

# ============================================================
# SECTION 1b (priority): continuous composite score
# ============================================================

flip <- tmr_sign[sig_tfs]
composite <- colMeans(meta_act[sig_tfs, ] * flip)
pc1_input <- t(meta_act[sig_tfs, ] * flip)
pca <- prcomp(pc1_input, scale. = FALSE)
pc1 <- pca$x[, 1]
if (cor(pc1, composite) < 0) pc1 <- -pc1   # orient PC1 to match the composite score
var_pc1 <- summary(pca)$importance["Proportion of Variance", 1]
cat("\nComposite score vs PC1 correlation:", round(cor(composite, pc1), 3),
    "| Variance explained by PC1:", round(100 * var_pc1, 1), "%\n")

cohort_of <- ifelse(colnames(meta_act) %in% colnames(tcga_basal_expr), "TCGA", "METABRIC")
score_df <- data.frame(sample_id = colnames(meta_act), cohort = cohort_of,
                        composite_score = composite, pc1_score = pc1,
                        row.names = NULL)
write.table(score_df, "data/basal_tmr_composite_score.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/basal_tmr_composite_score.tsv saved")

# ── Spearman: composite score vs. clinico-pathological variables ──

suppressMessages(library(TCGAbiolinks))
bs <- TCGAquery_subtype(tumor = "BRCA")
bs$patient12 <- bs$patient
tcga_clin_map <- bs[match(substr(tcga_ids, 1, 12), bs$patient12), ]
tcga_stage_num <- c("Stage_I"=1,"Stage_II"=2,"Stage_III"=3,"Stage_IV"=4)[tcga_clin_map$pathologic_stage]
tcga_age       <- tcga_clin_map$age_at_initial_pathologic_diagnosis

mtbrc_grade_num <- suppressWarnings(as.numeric(mtbrc_clin_basal$GRADE))
mtbrc_stage_num <- suppressWarnings(as.numeric(mtbrc_clin_basal$TUMOR_STAGE))
mtbrc_age       <- mtbrc_clin_basal$AGE_AT_DIAGNOSIS

spearman_test <- function(x, y, label) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 5) return(data.frame(variable = label, n = sum(ok), rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.frame(variable = label, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

comp_tcga  <- composite[tcga_ids]
comp_mtbrc <- composite[mtbrc_ids]

clin_corr <- rbind(
  spearman_test(comp_tcga,  tcga_stage_num,  "TCGA: composite vs pathologic_stage"),
  spearman_test(comp_tcga,  tcga_age,        "TCGA: composite vs age_at_diagnosis"),
  spearman_test(comp_mtbrc, mtbrc_grade_num, "METABRIC: composite vs GRADE"),
  spearman_test(comp_mtbrc, mtbrc_stage_num, "METABRIC: composite vs TUMOR_STAGE"),
  spearman_test(comp_mtbrc, mtbrc_age,       "METABRIC: composite vs AGE_AT_DIAGNOSIS")
)
clin_corr$FDR <- p.adjust(clin_corr$p, method = "BH")
cat("\n=== 1b: Spearman composite score vs. clinico-pathological variables ===\n")
print(clin_corr, row.names = FALSE)
write.table(clin_corr, "data/basal_composite_score_clinical_corr.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/basal_composite_score_clinical_corr.tsv saved")

# Note: TCGA Tumor_Grade is NA for BRCA in TCGAquery_subtype (not centrally
# reported for this tumor type) -> that correlation is skipped.
cat("\nTCGA Tumor_Grade available:", sum(!is.na(tcga_clin_map$Tumor_Grade)), "of", nrow(tcga_clin_map), "\n")

# ── Cox proportional hazards (continuous score) ──

days_na <- function(x) { x <- as.character(x); x[x %in% c("[Not Available]","")] <- NA; as.numeric(x) }
tcga_clin_map$days_to_death_num         <- days_na(tcga_clin_map$days_to_death)
tcga_clin_map$days_to_last_followup_num <- days_na(tcga_clin_map$days_to_last_followup)
tcga_time   <- ifelse(tcga_clin_map$vital_status %in% "Dead",
                       tcga_clin_map$days_to_death_num, tcga_clin_map$days_to_last_followup_num)
tcga_status <- as.integer(tcga_clin_map$vital_status %in% "Dead")

cox_tcga_df <- data.frame(time = tcga_time, status = tcga_status,
                           score = as.numeric(scale(comp_tcga)))
cox_tcga_df <- cox_tcga_df[complete.cases(cox_tcga_df) & cox_tcga_df$time > 0, ]
cox_tcga <- coxph(Surv(time, status) ~ score, data = cox_tcga_df)
cat("\n=== Cox TCGA (z-scaled composite score, OS) ===\n"); print(summary(cox_tcga))

mtbrc_status_bin <- as.integer(grepl("^1", mtbrc_clin_basal$OS_STATUS))
cox_mtbrc_df <- data.frame(time = mtbrc_clin_basal$OS_MONTHS, status = mtbrc_status_bin,
                            score = as.numeric(scale(comp_mtbrc)))
cox_mtbrc_df <- cox_mtbrc_df[complete.cases(cox_mtbrc_df) & cox_mtbrc_df$time > 0, ]
cox_mtbrc <- coxph(Surv(time, status) ~ score, data = cox_mtbrc_df)
cat("\n=== Cox METABRIC (z-scaled composite score, OS) ===\n"); print(summary(cox_mtbrc))

extract_cox <- function(fit, label, n) {
  s <- summary(fit)
  data.frame(cohort = label, n = n, HR = s$coefficients[1, "exp(coef)"],
             lower95 = s$conf.int[1, "lower .95"], upper95 = s$conf.int[1, "upper .95"],
             p = s$coefficients[1, "Pr(>|z|)"])
}
cox_tab <- rbind(
  extract_cox(cox_tcga,  "TCGA",     nrow(cox_tcga_df)),
  extract_cox(cox_mtbrc, "METABRIC", nrow(cox_mtbrc_df))
)
cox_tab$FDR <- p.adjust(cox_tab$p, method = "BH")
cat("\n=== Cox summary (continuous composite score, per cohort) ===\n")
print(cox_tab, row.names = FALSE)
write.table(cox_tab, "data/basal_composite_score_cox.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/basal_composite_score_cox.tsv saved")

# ============================================================
# SECTION 1b (secondary): de novo clustering over the significant TMRs
# ============================================================

clust_input <- t(meta_act[sig_tfs, ])   # samples x TMRs (no sign flip: clustering on the raw pattern)
d  <- as.dist(1 - cor(t(clust_input)))
hc <- hclust(d, method = "ward.D2")

sil_by_k <- sapply(2:4, function(k) {
  cl <- cutree(hc, k = k)
  mean(cluster::silhouette(cl, d)[, "sil_width"])
})
names(sil_by_k) <- paste0("k=", 2:4)
cat("\n=== Average silhouette per k (de novo clustering, significant TMRs) ===\n"); print(sil_by_k)

k_best <- as.integer(sub("k=", "", names(sil_by_k)[which.max(sil_by_k)]))
cl_best <- cutree(hc, k = k_best)
cat("\nChosen k (max silhouette):", k_best, "| average silhouette:", round(max(sil_by_k), 3), "\n")

clust_df <- data.frame(sample_id = names(cl_best), cluster = paste0("C", cl_best),
                        cohort = cohort_of[match(names(cl_best), colnames(meta_act))],
                        row.names = NULL)
write.table(clust_df, "data/basal_denovo_clusters.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/basal_denovo_clusters.tsv saved")

# ── de novo cluster <-> established subtype (chi-squared/Fisher) ──

cl_named <- setNames(clust_df$cluster, clust_df$sample_id)

fisher_or_chisq <- function(tab_x, tab_y, label) {
  ok <- !is.na(tab_x) & !is.na(tab_y) & nzchar(tab_y)
  tt <- table(tab_x[ok], tab_y[ok])
  if (nrow(tt) < 2 || ncol(tt) < 2) return(data.frame(test = label, p = NA, method = "n/a"))
  ft <- tryCatch(fisher.test(tt, simulate.p.value = TRUE, B = 10000),
                 error = function(e) NULL)
  if (!is.null(ft)) return(data.frame(test = label, p = ft$p.value, method = "Fisher (simulated)"))
  ct <- chisq.test(tt)
  data.frame(test = label, p = ct$p.value, method = "Chi-squared")
}

tab1 <- fisher_or_chisq(cl_named[tcga_ids],  tcga_subtype_named,   "TCGA: de novo cluster vs Lehmann/Bareche")
tab2 <- fisher_or_chisq(cl_named[mtbrc_ids], mtbrc_intclust_named, "METABRIC: de novo cluster vs IntClust")
cluster_vs_subtype <- rbind(tab1, tab2)
cluster_vs_subtype$FDR <- p.adjust(cluster_vs_subtype$p, method = "BH")
cat("\n=== de novo cluster <-> established subtype ===\n"); print(cluster_vs_subtype, row.names = FALSE)
write.table(cluster_vs_subtype, "data/basal_denovo_cluster_vs_subtype.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ── de novo cluster <-> survival (KM + log-rank), per cohort ──

km_tcga_df <- data.frame(time = tcga_time, status = tcga_status,
                          cluster = cl_named[tcga_ids])
km_tcga_df <- km_tcga_df[complete.cases(km_tcga_df) & km_tcga_df$time > 0, ]
km_mtbrc_df <- data.frame(time = mtbrc_clin_basal$OS_MONTHS, status = mtbrc_status_bin,
                           cluster = cl_named[mtbrc_ids])
km_mtbrc_df <- km_mtbrc_df[complete.cases(km_mtbrc_df) & km_mtbrc_df$time > 0, ]

logrank_tcga  <- survdiff(Surv(time, status) ~ cluster, data = km_tcga_df)
logrank_mtbrc <- survdiff(Surv(time, status) ~ cluster, data = km_mtbrc_df)
logrank_p <- function(sd) 1 - pchisq(sd$chisq, length(sd$n) - 1)

cat("\nLog-rank TCGA de novo cluster, p =", logrank_p(logrank_tcga), "\n")
cat("Log-rank METABRIC de novo cluster, p =", logrank_p(logrank_mtbrc), "\n")

km_summary <- data.frame(
  cohort = c("TCGA", "METABRIC"),
  n = c(nrow(km_tcga_df), nrow(km_mtbrc_df)),
  k = k_best,
  logrank_p = c(logrank_p(logrank_tcga), logrank_p(logrank_mtbrc))
)
write.table(km_summary, "data/basal_denovo_cluster_survival.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/basal_denovo_cluster_survival.tsv saved")

fit_tcga  <- survfit(Surv(time, status) ~ cluster, data = km_tcga_df)
fit_mtbrc <- survfit(Surv(time, status) ~ cluster, data = km_mtbrc_df)


# ggsurvplot's default pval formatting varies decimal places by significant
# digits (e.g. "0.093" vs "0.81"); pass the log-rank p-value already computed
# above, fixed to 2 decimals, so both panels format consistently.
pval_label <- function(p) paste0("p = ", sprintf("%.2f", p))

p_km_tcga <- ggsurvplot(fit_tcga, data = km_tcga_df,
                         pval = pval_label(logrank_p(logrank_tcga)), risk.table = TRUE,
                         xlab = "Months")
p_km_mtbrc <- ggsurvplot(fit_mtbrc, data = km_mtbrc_df,
                          pval = pval_label(logrank_p(logrank_mtbrc)), risk.table = TRUE,
                          xlab = "Months")

# Single combined figure: TCGA (A, left) | METABRIC (B, right), each panel
# with its survival curve stacked over its risk table. No panel-level titles
# embedded in the image — captions belong in the manuscript figure legend.
panel_tcga  <- wrap_elements(full = p_km_tcga$plot  / p_km_tcga$table  + plot_layout(heights = c(3, 1)))
panel_mtbrc <- wrap_elements(full = p_km_mtbrc$plot / p_km_mtbrc$table + plot_layout(heights = c(3, 1)))

km_combined <- (panel_tcga | panel_mtbrc) +
  plot_annotation(tag_levels = "A")

ggsave("figures/FigureS3.pdf", km_combined, width = 14, height = 7)
message("figures/FigureS3.pdf saved")

# ============================================================
# SECTION: updated heatmap (Fig 1C) with dual annotation
# ============================================================

established <- character(length(colnames(meta_act)))
names(established) <- colnames(meta_act)
established[tcga_ids]  <- tcga_subtype_named[tcga_ids]
established[mtbrc_ids] <- paste0("IntClust_", mtbrc_intclust_named[mtbrc_ids])

col_anno <- data.frame(
  Cohort                 = cohort_of,
  `Established subtype`  = established,
  `De novo cluster`      = cl_named[colnames(meta_act)],
  row.names              = colnames(meta_act),
  check.names            = FALSE
)

# No embedded title: panel/figure captions belong in the manuscript
# figure legend, not baked into the image.
ph <- pheatmap(
  meta_act[sig_tfs, ],
  color                    = colorRampPalette(c("#2980B9","white","#C0392B"))(100),
  clustering_distance_cols = "correlation",
  clustering_method        = "ward.D2",
  cluster_rows             = TRUE,
  show_colnames            = FALSE,
  annotation_col           = col_anno,
  silent                   = TRUE
)

# ============================================================
# SECTION: assemble Figure 1 (submission) — panels A/B built in
# compare_kbhb_mrs.R (NES concordance scatter, meta-NES lollipop),
# reloaded here, combined with the updated panel C above (7 significant
# Kbhb TMRs, established subtype + de novo cluster annotation).
# ============================================================

fig1 <- readRDS("data/fig1_nes_scatter.rds")
fig2 <- readRDS("data/fig2_lollipop_metanes.rds")
panel_c <- patchwork::wrap_elements(full = ph$gtable)

figure1 <- (fig1 | fig2) / panel_c +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    caption    = "Huang et al. 2021 gene set | ARACNe-AP + msVIPER | TCGA + METABRIC",
    tag_levels = "A",
    theme      = theme(
      plot.tag = element_text(face = "bold", size = 16)
    )
  )

ggsave("figures/Figure1.pdf", figure1, width = 16, height = 14)
message("figures/Figure1.pdf saved")

message("\n====  TASK 1 COMPLETE  ====")
