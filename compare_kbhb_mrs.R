library(viper)
library(ggplot2)
library(ggrepel)
library(ggside)
library(patchwork)
library(scales)
library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# Set your project directory:
# setwd("/path/to/Kbhb")
if (!interactive()) setwd(normalizePath("."))
dir.create("figures", showWarnings = FALSE)

# ============================================================
# Kbhb MRA COMPARISON: TCGA vs METABRIC
#
# 1. Stouffer integration: meta-NES and meta-FDR
# 2. Per-sample VIPER activity (z-scored, combined across cohorts)
# 3. Figures:
#    Fig 1 — NES scatter TCGA vs METABRIC
#    Fig 2 — Lollipop top MRs by meta-NES
#    Fig 3 — VIPER activity heatmap
# 4. ORA of significant TMR regulons (GO-BP + Reactome)
# ============================================================


# ============================================================
# SECTION 1: Load msVIPER results
# ============================================================

tcga_mrs  <- readRDS("data/tcga_basal_mrs.rds")
mtbrc_mrs <- readRDS("data/mtbrc_basal_mrs.rds")

tcga_nes  <- tcga_mrs$es$nes
mtbrc_nes <- mtbrc_mrs$es$nes

cat("TCGA regulators:", length(tcga_nes), "\n")
cat("METABRIC regulators:", length(mtbrc_nes), "\n")


# ============================================================
# SECTION 2: Stouffer integration (meta-NES)
# ============================================================
# msVIPER NES approximates a standard normal. Equal-weight Stouffer:
# Z_meta = (Z_TCGA + Z_METABRIC) / sqrt(2)

common_regs <- intersect(names(tcga_nes), names(mtbrc_nes))
cat("Regulators in common:", length(common_regs), "\n")

df <- data.frame(
  regulator  = common_regs,
  NES_TCGA   = tcga_nes[common_regs],
  NES_MTBRC  = mtbrc_nes[common_regs],
  size_TCGA  = tcga_mrs$es$size[common_regs],
  size_MTBRC = mtbrc_mrs$es$size[common_regs],
  pval_TCGA  = tcga_mrs$es$p.value[common_regs],
  pval_MTBRC = mtbrc_mrs$es$p.value[common_regs],
  row.names  = NULL
)

df$NES_meta  <- (df$NES_TCGA + df$NES_MTBRC) / sqrt(2)
df$pval_meta <- 2 * pnorm(-abs(df$NES_meta))
df$FDR_meta  <- p.adjust(df$pval_meta, method = "BH")

df <- df[order(-df$NES_meta), ]
rownames(df) <- NULL

cat("\nTop 20 MRs by meta-NES:\n")
print(df[1:20, c("regulator","NES_TCGA","NES_MTBRC","NES_meta","FDR_meta")])

saveRDS(df, "data/kbhb_mrs_comparison.rds")
write.table(df, "data/kbhb_mrs_comparison.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/kbhb_mrs_comparison.tsv saved")


# ============================================================
# SECTION 3: Per-sample VIPER activity (combined across cohorts)
# ============================================================
# Run per-sample VIPER on each dataset using its own regulon,
# then integrate scores — manual equivalent of metaVIPER.

# Reconstruct regulons for metaVIPER
tcga_basal  <- readRDS("data/tcga_basal_expr.rds")
mtbrc_basal <- readRDS("data/mtbrc_basal_expr.rds")

tcga_regulon  <- aracne2regulon("data/tcga_basal_network.txt",  eset = tcga_basal)
mtbrc_regulon <- aracne2regulon("data/mtbrc_basal_network.txt", eset = mtbrc_basal)

# Each dataset uses its own (most informative) regulon; scores are
# z-normalized per TF within each dataset before combining.

# Run per-sample VIPER — TCGA
tcga_act <- viper(eset = tcga_basal, regulon = tcga_regulon,
                  minsize = 25, verbose = FALSE)

# Run per-sample VIPER — METABRIC
mtbrc_act <- viper(eset = mtbrc_basal, regulon = mtbrc_regulon,
                   minsize = 25, verbose = FALSE)

# Z-score per TF within each dataset to make activity scales comparable
z_rows     <- function(m) t(scale(t(m)))
tcga_act_z  <- z_rows(tcga_act)
mtbrc_act_z <- z_rows(mtbrc_act)

# Combine: TFs present in both regulons, all samples
common_tfs <- intersect(rownames(tcga_act_z), rownames(mtbrc_act_z))
meta_act   <- cbind(tcga_act_z[common_tfs, ], mtbrc_act_z[common_tfs, ])
cat("TFs in common:", length(common_tfs), "\n")
cat("Combined activity dimensions:", dim(meta_act), "\n")

saveRDS(meta_act, "data/kbhb_metaviper_activity.rds")
message("data/kbhb_metaviper_activity.rds saved")


# ============================================================
# SECTION 4: Figures
# ============================================================

# ── Base theme ───────────────────────────────────────────────

pal_dir <- c(
  "Concordant +"  = "#C0392B",
  "Concordant −" = "#2980B9",
  "Discordant"    = "grey60"
)

theme_kbhb <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 10, color = "grey40"),
    legend.position   = "bottom"
  )

# ── Classify concordance ─────────────────────────────────────

df$concordance <- ifelse(
  df$NES_TCGA > 0 & df$NES_MTBRC > 0, "Concordant +",
  ifelse(df$NES_TCGA < 0 & df$NES_MTBRC < 0, "Concordant −", "Discordant")
)

# Top 15 concordant regulators to label
label_regs <- df$regulator[df$concordance != "Discordant"][
  order(-abs(df$NES_meta[df$concordance != "Discordant"]))[1:15]
]
df$label <- ifelse(df$regulator %in% label_regs, df$regulator, NA)

r_val <- cor(df$NES_TCGA, df$NES_MTBRC)


# ── FIG 1: NES scatter TCGA vs METABRIC ─────────────────────

fig1 <- ggplot(df, aes(x = NES_TCGA, y = NES_MTBRC,
                        color = concordance, label = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(size = -log10(pval_meta + 1e-6)), alpha = 0.75) +
  geom_text_repel(
    na.rm = TRUE, size = 3, fontface = "bold",
    max.overlaps = 20, box.padding = 0.4,
    segment.color = "grey50", segment.size = 0.3
  ) +
  scale_color_manual(values = pal_dir, name = NULL) +
  scale_size_continuous(name = expression(-log[10]~"(p meta)"), range = c(1.5, 5)) +
  annotate("text", x = min(df$NES_TCGA) * 0.9, y = max(df$NES_MTBRC) * 0.9,
           label = paste0("r = ", round(r_val, 2)),
           hjust = 0, size = 4, fontface = "italic") +
  labs(
    title    = "Kbhb MRA — Concordance across cohorts",
    subtitle = "Transcriptional master regulators of the Kbhb program in Basal BRCA",
    x        = "NES (TCGA Basal, n = 195)",
    y        = "NES (METABRIC Basal, n = 209)"
  ) +
  theme_kbhb

ggsave("figures/fig1_nes_scatter.pdf", fig1, width = 7, height = 6)
message("figures/fig1_nes_scatter.pdf saved")


# ── FIG 2: Lollipop top MRs by meta-NES ─────────────────────

top_pos <- head(df[df$NES_meta > 0, ], 10)
top_neg <- tail(df[df$NES_meta < 0, ], 10)
df_lol  <- rbind(top_pos, top_neg)
df_lol$regulator <- factor(df_lol$regulator,
                            levels = df_lol$regulator[order(df_lol$NES_meta)])

df_lol_long <- rbind(
  data.frame(regulator = df_lol$regulator, NES = df_lol$NES_TCGA,
             dataset = "TCGA",     NES_meta = df_lol$NES_meta),
  data.frame(regulator = df_lol$regulator, NES = df_lol$NES_MTBRC,
             dataset = "METABRIC", NES_meta = df_lol$NES_meta)
)

fig2 <- ggplot(df_lol, aes(x = regulator, y = NES_meta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_segment(aes(xend = regulator, yend = 0, color = NES_meta > 0),
               linewidth = 1.2) +
  geom_point(data = df_lol_long,
             aes(y = NES, shape = dataset, fill = dataset),
             size = 3, color = "white", stroke = 0.5) +
  geom_point(aes(color = NES_meta > 0), size = 4) +
  scale_color_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9"),
                     labels = c("TRUE" = "Activated", "FALSE" = "Repressed"),
                     name = "Direction") +
  scale_shape_manual(values = c("TCGA" = 21, "METABRIC" = 24), name = "Cohort") +
  scale_fill_manual(values = c("TCGA" = "#E67E22", "METABRIC" = "#8E44AD"),
                    name = "Cohort") +
  coord_flip() +
  labs(
    title    = "Top Master Regulators of the Kbhb program — meta-NES",
    subtitle = "Stouffer combination (TCGA + METABRIC). Points: individual NES per cohort.",
    x        = NULL,
    y        = "meta-NES (Stouffer)"
  ) +
  theme_kbhb +
  theme(legend.box = "horizontal")

ggsave("figures/fig2_lollipop_metanes.pdf", fig2, width = 8, height = 7)
message("figures/fig2_lollipop_metanes.pdf saved")


# ── FIG 3: VIPER activity heatmap — top MRs per sample ───────

top25    <- df$regulator[order(-abs(df$NES_meta))][1:25]
top25    <- intersect(top25, rownames(meta_act))
act_sub  <- meta_act[top25, ]

sample_meta <- data.frame(
  sample  = colnames(meta_act),
  dataset = ifelse(colnames(meta_act) %in% colnames(tcga_basal),
                   "TCGA", "METABRIC"),
  row.names = colnames(meta_act)
)

pos_regs     <- intersect(df$regulator[df$NES_meta > 0], rownames(act_sub))
if (length(pos_regs) == 0) pos_regs <- rownames(act_sub)[1]
sort_scores  <- colMeans(act_sub[pos_regs, , drop = FALSE])
sample_order <- order(sample_meta$dataset, sort_scores)

df_heat <- as.data.frame(as.table(act_sub[, sample_order]))
colnames(df_heat) <- c("TF", "Sample", "Activity")
df_heat$Dataset <- sample_meta[as.character(df_heat$Sample), "dataset"]
df_heat$TF      <- factor(df_heat$TF,
                           levels = top25[order(df[match(top25, df$regulator), "NES_meta"])])

fig3 <- ggplot(df_heat, aes(x = Sample, y = TF, fill = Activity)) +
  geom_tile() +
  facet_grid(. ~ Dataset, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colors = c("#2980B9", "white", "#C0392B"),
    limits = c(-max(abs(act_sub)), max(abs(act_sub))),
    name   = "VIPER\nActivity"
  ) +
  labs(
    title    = "VIPER activity — per sample",
    subtitle = "Top 25 MRs of the Kbhb program (TCGA + METABRIC)",
    x = NULL, y = NULL
  ) +
  theme_kbhb +
  theme(
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    strip.text    = element_text(face = "bold"),
    panel.spacing = unit(0.15, "cm")
  )

ggsave("figures/fig3_metaviper_heatmap.pdf", fig3, width = 10, height = 6)
message("figures/fig3_metaviper_heatmap.pdf saved")


# ── Combined panel ───────────────────────────────────────────

panel <- (fig1 | fig2) / fig3 +
  plot_annotation(
    title       = "Kbhb transcriptional program in Basal BRCA — MRA",
    caption     = "Huang et al. 2021 gene set | ARACNe-AP + msVIPER | TCGA + METABRIC",
    tag_levels  = "A",
    theme       = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.tag   = element_text(face = "bold", size = 16)
    )
  )

ggsave("figures/panel_kbhb_mra.pdf", panel, width = 16, height = 12)
message("figures/panel_kbhb_mra.pdf saved")


# ============================================================
# SECTION 5: ORA of significant TMR regulons
#            clusterProfiler (GO-BP) + ReactomePA
# ============================================================

sig_tfs <- df$regulator[df$FDR_meta < 0.05]
cat("\nSignificant TMRs for ORA (meta FDR < 0.05):", paste(sig_tfs, collapse = ", "), "\n")

# ── 5.1 Regulon targets from the TCGA regulon ────────────────
# TCGA is the discovery cohort; METABRIC targets give similar results.
# Note: VEZF1 (the sole repressed TMR) lacks significant Reactome enrichment
# because its endothelial/vascular targets are not well-represented in
# Reactome pathways at FDR < 0.05, while GO-BP captures them.

available_tfs <- intersect(sig_tfs, names(tcga_regulon))
missing_tfs   <- setdiff(sig_tfs, names(tcga_regulon))
if (length(missing_tfs) > 0)
  cat("Warning — TMRs absent from TCGA regulon:", paste(missing_tfs, collapse = ", "), "\n")

gene_lists_sym <- lapply(
  setNames(available_tfs, available_tfs),
  function(tf) names(tcga_regulon[[tf]]$tfmode)
)

cat("Regulon sizes (TCGA):\n")
print(sapply(gene_lists_sym, length))

# ── 5.2 Convert to Entrez IDs ────────────────────────────────

sym2entrez <- function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
}

message("Converting gene symbols to Entrez IDs...")
gene_lists_entrez <- lapply(gene_lists_sym, sym2entrez)
universe_entrez   <- sym2entrez(rownames(tcga_basal))

cat("Universe (TCGA genes with Entrez ID):", length(universe_entrez), "\n")

# ── 5.3 ORA — GO Biological Process ──────────────────────────

message("Running ORA GO Biological Process...")
cc_go <- compareCluster(
  geneClusters  = gene_lists_entrez,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  universe      = universe_entrez,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)
cc_go <- clusterProfiler::simplify(cc_go, cutoff = 0.6, by = "p.adjust")

# ── 5.4 ORA — Reactome ───────────────────────────────────────

message("Running ORA Reactome...")
cc_react <- compareCluster(
  geneClusters  = gene_lists_entrez,
  fun           = "enrichPathway",
  organism      = "human",
  universe      = universe_entrez,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# ── 5.5 ORA figure: classic dotplot — top 20 pathways per database ──
# Pathways on Y-axis, TMRs on X-axis, faceted by source (GO-BP | Reactome).
# Top 20 selected by: most TMRs enriching them first, then lowest min p.adjust.

parse_gene_ratio <- function(x) {
  sapply(x, function(r) {
    p <- as.numeric(strsplit(r, "/")[[1]])
    p[1] / p[2]
  })
}

df_go_raw <- as.data.frame(cc_go) %>%
  mutate(source = "GO-BP",    GeneRatioNum = parse_gene_ratio(GeneRatio))
df_rx_raw <- as.data.frame(cc_react) %>%
  mutate(source = "Reactome", GeneRatioNum = parse_gene_ratio(GeneRatio))

df_ora <- bind_rows(df_go_raw, df_rx_raw) %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(Cluster, Description, GeneRatioNum, p.adjust, Count, source) %>%
  mutate(Description = ifelse(source == "GO-BP",
                              paste0("GOBP_", Description),
                              paste0("RCTM_", Description)))

# TMR order: meta-NES descending (left = highest NES on x-axis)
tmr_order <- sig_tfs

# Full pathway ordering and numbered index (reference file, all pathways)
pathway_order <- df_ora %>%
  group_by(source, Description) %>%
  summarise(n_tmr = n_distinct(Cluster), .groups = "drop") %>%
  arrange(source, desc(n_tmr), Description) %>%
  pull(Description)

pathway_table <- data.frame(
  key     = pathway_order,
  PathID  = seq_along(pathway_order),
  source  = ifelse(startsWith(pathway_order, "GOBP_"), "GO-BP", "Reactome"),
  Pathway = sub("^(GOBP|RCTM)_", "", pathway_order),
  stringsAsFactors = FALSE
)
write.table(pathway_table[, c("PathID", "source", "Pathway")],
            "data/ORA_pathway_index.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/ORA_pathway_index.tsv saved (", nrow(pathway_table), " pathways)")

df_ora <- df_ora %>%
  left_join(dplyr::select(pathway_table, key, PathID),
            by = c("Description" = "key")) %>%
  mutate(Cluster = factor(Cluster, levels = tmr_order))

# Helper: truncate long pathway names
trunc_label <- function(x, n = 55) {
  ifelse(nchar(x) > n, paste0(strtrim(x, n - 3), "..."), x)
}

# Select top 20 per database
select_top20 <- function(src) {
  df_ora %>%
    filter(source == src) %>%
    group_by(Description) %>%
    summarise(n_tmr = n_distinct(Cluster), min_padj = min(p.adjust), .groups = "drop") %>%
    arrange(desc(n_tmr), min_padj) %>%
    slice_head(n = 20) %>%
    pull(Description)
}
top20_go <- select_top20("GO-BP")
top20_rx <- select_top20("Reactome")

# Force-include all GO-BP pathways for any TMR absent from the top-20 cut
# (e.g., VEZF1 has only 2 GO-BP and 0 Reactome enrichments)
tmrs_in_gobp <- df_ora %>%
  filter(Description %in% top20_go) %>%
  pull(Cluster) %>% as.character() %>% unique()
missing_tmrs <- setdiff(as.character(tmr_order), tmrs_in_gobp)
if (length(missing_tmrs) > 0) {
  extra_go <- df_ora %>%
    filter(source == "GO-BP", Cluster %in% missing_tmrs) %>%
    pull(Description) %>% unique()
  top20_go <- unique(c(top20_go, extra_go))
  cat("Added", length(extra_go), "extra GO-BP pathways for TMR(s):",
      paste(missing_tmrs, collapse = ", "), "\n")
}

df_top <- df_ora %>%
  filter(Description %in% c(top20_go, top20_rx)) %>%
  mutate(Label = trunc_label(sub("^(GOBP|RCTM)_", "", Description)))

# Within each facet, order pathways by ascending n_tmr (most shared at top)
path_levels <- df_top %>%
  group_by(source, Label) %>%
  summarise(n_tmr = n_distinct(Cluster), .groups = "drop") %>%
  arrange(source, n_tmr) %>%
  pull(Label)

df_top <- df_top %>%
  mutate(Label = factor(Label, levels = unique(path_levels)))

fig_ora <- ggplot(df_top, aes(x = Cluster, y = Label)) +
  geom_point(aes(size = GeneRatioNum, color = p.adjust)) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradientn(
    colors = c("#C0392B", "#8E44AD", "#2980B9"),
    trans  = "log10",
    name   = "p.adjust",
    guide  = guide_colorbar(barheight = 6)
  ) +
  scale_size_continuous(name = "Gene Ratio", range = c(1.5, 6)) +
  scale_x_discrete(drop = FALSE) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey92"),
    strip.background  = element_rect(fill = "grey92", color = NA),
    strip.text        = element_text(face = "bold", size = 10),
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_text(size = 10, face = "bold", angle = 35, hjust = 1)
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "ORA of significant TMR regulons — Top 20 pathways per database"
  )

ggsave("figures/fig_ora_tmrs.pdf", fig_ora, width = 8, height = 12)
message("figures/fig_ora_tmrs.pdf saved")

# ── 5.6 Summary tables ───────────────────────────────────────────────────────

# Table A: pathway counts per TMR broken down by database
tmr_counts <- df_ora %>%
  group_by(Cluster) %>%
  summarise(
    n_GOBP     = n_distinct(PathID[source == "GO-BP"]),
    n_Reactome = n_distinct(PathID[source == "Reactome"]),
    n_total    = n_distinct(PathID),
    .groups    = "drop"
  ) %>%
  arrange(desc(n_total))

cat("\n── ORA pathway counts per TMR ──\n")
print(as.data.frame(tmr_counts))
write.table(tmr_counts, "data/ORA_tmr_counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/ORA_tmr_counts.tsv saved")

# Table B: top 20 pathways most shared across TMRs
top_shared <- df_ora %>%
  group_by(Description, source) %>%
  summarise(
    n_TMRs   = n_distinct(Cluster),
    TMRs     = paste(sort(unique(as.character(Cluster))), collapse = ", "),
    min_padj = min(p.adjust),
    .groups  = "drop"
  ) %>%
  mutate(Pathway = sub("^(GOBP|RCTM)_", "", Description)) %>%
  arrange(desc(n_TMRs), min_padj) %>%
  dplyr::select(source, Pathway, n_TMRs, TMRs, min_padj) %>%
  slice_head(n = 20)

cat("\n── Top 20 most shared pathways across TMRs ──\n")
print(as.data.frame(top_shared))
write.table(top_shared, "data/ORA_top_shared_pathways.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/ORA_top_shared_pathways.tsv saved")

# ── 5.7 Save ORA results ─────────────────────────────────────

saveRDS(list(go = cc_go, reactome = cc_react), "data/kbhb_ora_results.rds")
write.table(as.data.frame(cc_go),    "data/ORA_GO_BP.tsv",    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(cc_react), "data/ORA_Reactome.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/kbhb_ora_results.rds, data/ORA_GO_BP.tsv, data/ORA_Reactome.tsv saved")


message("\n====  COMPARISON COMPLETE  ====")
message("Table: data/kbhb_mrs_comparison.tsv")
message("Figures in figures/:")
message("  fig1_nes_scatter.pdf")
message("  fig2_lollipop_metanes.pdf")
message("  fig3_metaviper_heatmap.pdf")
message("  panel_kbhb_mra.pdf")
message("  fig_ora_tmrs.pdf  (top 20 per DB dotplot)")
message("Data tables in data/:")
message("  ORA_pathway_index.tsv     — all enriched pathways with numeric ID")
message("  ORA_tmr_counts.tsv        — pathways per TMR (GO-BP / Reactome / total)")
message("  ORA_top_shared_pathways.tsv — top 20 pathways shared across most TMRs")
message("Note: shadow network panels (C/D) are part of figures/fig_supp_combined.pdf — see mra_kbhb.R")
