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

write.table(comparison_full, "data/kbhb_tmr_subtype_comparison.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("data/kbhb_tmr_subtype_comparison.tsv saved")

# ============================================================
# FIGURE: meta-NES of the significant Kbhb TMRs across the 4 subtypes
# (similar to the Fig 2 lollipop, but comparing subtypes instead of cohorts)
# ============================================================

comparison_full$sig <- ifelse(!is.na(comparison_full$FDR_meta) & comparison_full$FDR_meta < 0.05,
                               "FDR < 0.05", "n.s. / not evaluable")

fig_compare <- ggplot(comparison_full,
                       aes(x = subtype, y = NES_meta, fill = sig)) +
  geom_col(width = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  facet_wrap(~ regulator, nrow = 2) +
  scale_fill_manual(values = c("FDR < 0.05" = "#C0392B", "n.s. / not evaluable" = "grey70"),
                     name = NULL) +
  labs(
    title    = "Status of the Kbhb axis (CENPA/FOXM1 + TMRs) by PAM50 subtype",
    subtitle = "meta-NES (Stouffer, TCGA + METABRIC) — Basal-like vs Luminal A/B vs HER2-enriched",
    x = NULL, y = "meta-NES (Stouffer)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold")
  )

ggsave("figures/kbhb_tmr_subtype_comparison.pdf", fig_compare, width = 10, height = 7)
message("figures/kbhb_tmr_subtype_comparison.pdf saved")

message("\n====  MRA + META-ANALYSIS COMPLETE (LumA / LumB / Her2)  ====")
