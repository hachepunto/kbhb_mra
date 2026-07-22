library(viper)
library(ggplot2)
library(patchwork)
library(dplyr)

if (!interactive()) setwd(normalizePath("."))
dir.create("data",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ============================================================
# Extension of the MRA to Luminal A/B and HER2-enriched
#
# Same msVIPER + meta-analysis pipeline as mra_kbhb.R / compare_kbhb_mrs.R
# (see the parameter documentation in the results notes for exact values),
# repeated for LumA, LumB, Her2, using the new ARACNe-AP networks
# (data/{tcga,mtbrc}_{subtype}_network.txt) and reusing the same Kbhb gene
# set (data/kbhb_geneset.rds) and the same Normal samples already used.
# ============================================================

subtypes   <- c("LumA", "LumB", "Her2")
kbhb_genes <- readRDS("data/kbhb_geneset.rds")
cat("Kbhb genes (Huang S1):", length(kbhb_genes), "\n")

run_msviper_subtype <- function(subtype) {

  cat("\n============================================================\n")
  cat("SUBTYPE:", subtype, "\n")
  cat("============================================================\n")

  tcga_tumor   <- readRDS(sprintf("data/tcga_%s_expr.rds",         subtype))
  tcga_normal  <- readRDS(sprintf("data/tcga_%s_normal_expr.rds",  subtype))
  mtbrc_tumor  <- readRDS(sprintf("data/mtbrc_%s_expr.rds",        subtype))
  mtbrc_normal <- readRDS(sprintf("data/mtbrc_%s_normal_expr.rds", subtype))

  kbhb_tcga  <- intersect(kbhb_genes, rownames(tcga_tumor))
  kbhb_mtbrc <- intersect(kbhb_genes, rownames(mtbrc_tumor))
  cat("Kbhb genes in TCGA", subtype, ": ", length(kbhb_tcga),  "\n")
  cat("Kbhb genes in METABRIC", subtype, ": ", length(kbhb_mtbrc), "\n")

  message("Building TCGA regulon for ", subtype, "...")
  tcga_regulon <- aracne2regulon(
    afile = sprintf("data/tcga_%s_network.txt", subtype),
    eset  = tcga_tumor
  )
  cat("TCGA regulon", subtype, "— regulators:", length(tcga_regulon), "\n")

  message("Building METABRIC regulon for ", subtype, "...")
  mtbrc_regulon <- aracne2regulon(
    afile = sprintf("data/mtbrc_%s_network.txt", subtype),
    eset  = mtbrc_tumor
  )
  cat("METABRIC regulon", subtype, "— regulators:", length(mtbrc_regulon), "\n")

  message("Kbhb signature, TCGA ", subtype, " (tumor vs Normal)...")
  tcga_sig <- rowTtest(x = tcga_tumor[kbhb_tcga, ],  y = tcga_normal[kbhb_tcga, ])
  message("Kbhb signature, METABRIC ", subtype, "...")
  mtbrc_sig <- rowTtest(x = mtbrc_tumor[kbhb_mtbrc, ], y = mtbrc_normal[kbhb_mtbrc, ])

  message("Null model, TCGA ", subtype, " (1000 permutations)...")
  tcga_null <- ttestNull(
    x = tcga_tumor[kbhb_tcga, ], y = tcga_normal[kbhb_tcga, ],
    per = 1000, repos = TRUE, seed = 1
  )
  message("Null model, METABRIC ", subtype, "...")
  mtbrc_null <- ttestNull(
    x = mtbrc_tumor[kbhb_mtbrc, ], y = mtbrc_normal[kbhb_mtbrc, ],
    per = 1000, repos = TRUE, seed = 1
  )

  message("msVIPER, TCGA ", subtype, "...")
  tcga_mrs <- msviper(ges = tcga_sig$statistic, regulon = tcga_regulon,
                       nullmodel = tcga_null, minsize = 25, verbose = FALSE)
  tcga_mrs <- ledge(tcga_mrs)
  tcga_mrs <- shadow(tcga_mrs, minsize = 25, verbose = FALSE)

  message("msVIPER, METABRIC ", subtype, "...")
  mtbrc_mrs <- msviper(ges = mtbrc_sig$statistic, regulon = mtbrc_regulon,
                        nullmodel = mtbrc_null, minsize = 25, verbose = FALSE)
  mtbrc_mrs <- ledge(mtbrc_mrs)
  mtbrc_mrs <- shadow(mtbrc_mrs, minsize = 25, verbose = FALSE)

  saveRDS(tcga_mrs,  sprintf("data/tcga_%s_mrs.rds",  subtype))
  saveRDS(mtbrc_mrs, sprintf("data/mtbrc_%s_mrs.rds", subtype))
  message("data/tcga_", subtype, "_mrs.rds and data/mtbrc_", subtype, "_mrs.rds saved")

  cat("\n=== Top 15 Master Regulators — TCGA", subtype, "===\n")
  print(summary(tcga_mrs, mrs = 15))
  cat("\n=== Top 15 Master Regulators — METABRIC", subtype, "===\n")
  print(summary(mtbrc_mrs, mrs = 15))

  # ── Stouffer meta-analysis (same method as compare_kbhb_mrs.R) ──

  tcga_nes  <- tcga_mrs$es$nes
  mtbrc_nes <- mtbrc_mrs$es$nes
  common_regs <- intersect(names(tcga_nes), names(mtbrc_nes))
  cat("Regulators in common", subtype, ":", length(common_regs), "\n")

  df <- data.frame(
    regulator  = common_regs,
    NES_TCGA   = tcga_nes[common_regs],
    NES_MTBRC  = mtbrc_nes[common_regs],
    pval_TCGA  = tcga_mrs$es$p.value[common_regs],
    pval_MTBRC = mtbrc_mrs$es$p.value[common_regs],
    row.names  = NULL
  )
  df$NES_meta  <- (df$NES_TCGA + df$NES_MTBRC) / sqrt(2)
  df$pval_meta <- 2 * pnorm(-abs(df$NES_meta))
  df$FDR_meta  <- p.adjust(df$pval_meta, method = "BH")
  df <- df[order(-df$NES_meta), ]
  rownames(df) <- NULL

  write.table(df, sprintf("data/%s_mrs_comparison.tsv", tolower(subtype)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("data/", tolower(subtype), "_mrs_comparison.tsv saved")

  sig_tfs <- df$regulator[df$FDR_meta < 0.05]
  cat("\nSignificant TMRs (meta FDR<0.05) in", subtype, ":",
      paste(sig_tfs, collapse = ", "), "\n")

  df
}

results <- setNames(lapply(subtypes, run_msviper_subtype), subtypes)

# ============================================================
# COMPARISON TABLE: Basal-like vs LumA vs LumB vs HER2
# For the TMRs already significant in Basal-like (FOXM1, CENPA, and the rest)
# ============================================================

basal_df   <- read.delim("data/kbhb_mrs_comparison.tsv", stringsAsFactors = FALSE)
tmrs_basal <- basal_df$regulator[basal_df$FDR_meta < 0.05]   # significant Basal-like TMRs, derived from the data

extract_subtype <- function(df, tfs, label) {
  sub <- df[df$regulator %in% tfs, c("regulator","NES_TCGA","NES_MTBRC","NES_meta","FDR_meta")]
  sub$subtype <- label
  sub
}

comparison_long <- rbind(
  extract_subtype(basal_df,      tmrs_basal, "Basal-like"),
  extract_subtype(results$LumA,  tmrs_basal, "LumA"),
  extract_subtype(results$LumB,  tmrs_basal, "LumB"),
  extract_subtype(results$Her2,  tmrs_basal, "HER2")
)

# TFs absent from the regulon, or filtered out by minsize=25, will not appear;
# they are recorded as explicit NA to leave a record in the table.
full_grid <- expand.grid(regulator = tmrs_basal,
                          subtype   = c("Basal-like","LumA","LumB","HER2"),
                          stringsAsFactors = FALSE)
comparison_full <- merge(full_grid, comparison_long, by = c("regulator","subtype"), all.x = TRUE)
comparison_full$regulator <- factor(comparison_full$regulator, levels = tmrs_basal)
comparison_full$subtype   <- factor(comparison_full$subtype,
                                     levels = c("Basal-like","LumA","LumB","HER2"))
comparison_full <- comparison_full[order(comparison_full$regulator, comparison_full$subtype), ]

cat("\n=== COMPARISON TABLE meta-NES: Basal-like vs LumA vs LumB vs HER2 ===\n")
print(comparison_full, row.names = FALSE)

write.table(comparison_full, "data/Supplementary_TableS3_meta_NES_subtypes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("data/Supplementary_TableS3_meta_NES_subtypes.tsv saved")

# ============================================================
# FIGURE: graded specificity of the Kbhb TMR panel across PAM50 subtypes
# Dot plot (not a plain heatmap / bar chart) — three cell states, and rows
# grouped into 3 specificity tiers by facet (see Results 2.4):
#   Tier 1 - broad/pan-active   : significant in >= 3 of 4 subtypes
#   Tier 2 - partially specific : significant in 2 of 4 subtypes
#   Tier 3 - Basal-like-exclusive: significant in 1 of 4 subtypes (Basal-like only)
# Tiers are derived from the data itself (not hardcoded), so the figure
# stays reproducible if the upstream msVIPER/meta-analysis results change.
# ============================================================

comparison_full$status <- factor(
  ifelse(is.na(comparison_full$NES_meta), "Not evaluable",
         ifelse(comparison_full$FDR_meta < 0.05, "Significant (FDR < 0.05)", "Not significant")),
  levels = c("Significant (FDR < 0.05)", "Not significant", "Not evaluable")
)

# Fill only carries the meta-NES value for significant cells; not-significant
# and not-evaluable cells fall back to na.value (flat grey), so full color
# saturation is reserved for the cells that actually support the claim.
comparison_full$fill_value <- ifelse(comparison_full$status == "Significant (FDR < 0.05)",
                                      comparison_full$NES_meta, NA)

n_sig <- tapply(comparison_full$status == "Significant (FDR < 0.05)",
                 comparison_full$regulator, sum)
tier_lab <- c(
  "Broad\n(≥ 3/4 subtypes)",
  "Partial\n(2/4 subtypes)",
  "Basal-like exclusive\n(1/4 subtypes)"
)
tier_of <- ifelse(n_sig[as.character(comparison_full$regulator)] >= 3, tier_lab[1],
            ifelse(n_sig[as.character(comparison_full$regulator)] == 2, tier_lab[2],
                   tier_lab[3]))
comparison_full$tier <- factor(tier_of, levels = tier_lab)

# Within each tier, order regulators by Basal-like meta-NES (descending),
# consistent with the Fig. 1B / Fig. 2 ranking; reversed because ggplot
# plots the first factor level at the bottom of the y-axis.
basal_rank <- basal_df$regulator[order(-basal_df$NES_meta)]
basal_rank <- basal_rank[basal_rank %in% tmrs_basal]
comparison_full$regulator <- factor(comparison_full$regulator, levels = rev(basal_rank))

cat("\n=== TMR tiers (graded specificity across PAM50 subtypes, n_sig/4) ===\n")
print(data.frame(regulator = names(n_sig), n_sig = as.integer(n_sig))[order(-n_sig), ], row.names = FALSE)

fig_compare <- ggplot(comparison_full, aes(x = subtype, y = regulator)) +
  geom_point(aes(shape = status, size = status, colour = status, stroke = status,
                  fill = fill_value)) +
  facet_grid(tier ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2980B9", mid = "white", high = "#C0392B", midpoint = 0,
    na.value = "grey85", name = "meta-NES\n(Stouffer)"
  ) +
  scale_shape_manual(name = "Status", values = c(
    "Significant (FDR < 0.05)" = 21, "Not significant" = 21, "Not evaluable" = 4),
    # "fill" is a separate continuous scale (meta-NES), not part of this
    # discrete "Status" guide, so ggplot has no fill to show in the legend
    # keys and defaults to blank/white for all three -- misleadingly implying
    # "Not significant" points are pale/white when the actual points are
    # solid grey85 (na.value, set below via fill_value <- NA for that group).
    # override.aes fixes only that one key to match the real point color;
    # "Significant" is left blank on purpose -- its points take a whole
    # range of colors from the continuous scale, so no single fill would be
    # accurate there, and that scale's own colorbar legend covers it.
    guide = guide_legend(override.aes = list(
      fill = c("Significant (FDR < 0.05)" = NA, "Not significant" = "grey85", "Not evaluable" = NA)
    ))) +
  scale_size_manual(name = "Status", values = c(
    "Significant (FDR < 0.05)" = 8, "Not significant" = 3.5, "Not evaluable" = 4.5)) +
  scale_colour_manual(name = "Status", values = c(
    "Significant (FDR < 0.05)" = "black", "Not significant" = "grey55", "Not evaluable" = "grey45")) +
  scale_discrete_manual(aesthetics = "stroke", name = "Status", values = c(
    "Significant (FDR < 0.05)" = 1.3, "Not significant" = 0.6, "Not evaluable" = 1.3)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid        = element_line(color = "grey93", linewidth = 0.3),
    axis.text.x       = element_text(face = "bold"),
    axis.text.y       = element_text(face = "italic"),
    strip.text.y      = element_text(face = "bold", angle = 0, hjust = 0),
    strip.background  = element_rect(fill = "grey92", color = NA),
    panel.spacing.y   = unit(0.6, "lines"),
    legend.position   = "right"
  )

ggsave("figures/Figure5.pdf", fig_compare, width = 9, height = 6.5)
message("figures/Figure5.pdf saved")

message("\n====  MRA + META-ANALYSIS COMPLETE (LumA / LumB / Her2)  ====")
