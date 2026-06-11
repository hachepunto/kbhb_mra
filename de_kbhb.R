library(DESeq2)
library(limma)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(ggraph)
library(tidygraph)
library(pheatmap)
library(dplyr)
library(patchwork)

# Set your project directory:
# setwd("/path/to/Kbhb")
if (!interactive()) setwd(normalizePath("."))
dir.create("figures", showWarnings = FALSE)

# ============================================================
# DIFFERENTIAL EXPRESSION OF THE Kbhb PROGRAM — Basal BRCA
#
# Basal Tumor vs Normal in two independent cohorts:
#   TCGA     -> DESeq2 on raw STAR counts (unstranded)
#   METABRIC -> limma on log2 microarray (ComBat-corrected)
#
# Outputs in data/:
#   tcga_de_results.rds   — DESeq2 data.frame (log2FoldChange, padj, ...)
#   mtbrc_de_results.rds  — limma data.frame  (logFC, adj.P.Val, ...)
#
# Figures in figures/:
#   fig_de_volcano_tcga.pdf
#   fig_de_volcano_mtbrc.pdf
#   fig_de_heatmap_kbhb.pdf
# ============================================================


# ============================================================
# SECTION 1: Kbhb gene set and significant MRA TMRs
# ============================================================

kbhb_genes <- readRDS("data/kbhb_geneset.rds")
mrs_df      <- read.delim("data/kbhb_mrs_comparison.tsv", stringsAsFactors = FALSE)
sig_tfs     <- mrs_df$regulator[mrs_df$FDR_meta < 0.05]

cat("Kbhb proteome genes (Huang 2021):", length(kbhb_genes), "\n")
cat("Significant TMRs (meta FDR < 0.05):", length(sig_tfs), "\n")
cat("TMRs:", paste(sig_tfs, collapse = ", "), "\n")


# ============================================================
# SECTION 2: TCGA-BRCA — DESeq2
# ============================================================

# ── 2.1 Load SE and extract raw counts ───────────────────────

message("Loading TCGA SummarizedExperiment...")
se <- readRDS("data/tcga_brca_rnaseq_se.rds")
cd <- as.data.frame(colData(se))
rd <- as.data.frame(rowData(se))

# Restrict to protein-coding genes
pc_idx <- rd$gene_type == "protein_coding"
se_pc  <- se[pc_idx, ]
rd_pc  <- rd[pc_idx, ]

# STAR unstranded counts — required by DESeq2 (raw integers, no normalisation)
counts_raw           <- assay(se_pc, "unstranded")
rownames(counts_raw) <- rd_pc$gene_name

# ── 2.2 Collapse duplicated HGNC symbols ─────────────────────
# Column-wise median across probes/transcripts sharing the same symbol,
# then round to integer (DESeq2 requirement). Mirrors basal_pre_networks.R.

dup_syms <- rownames(counts_raw)[duplicated(rownames(counts_raw))]
message("Duplicated protein_coding symbols: ", length(unique(dup_syms)))

df_counts      <- as.data.frame(counts_raw)
df_counts$.sym <- rownames(counts_raw)
counts_dedup   <- do.call(rbind,
  lapply(
    split(df_counts[, colnames(df_counts) != ".sym"], df_counts$.sym),
    function(chunk) apply(chunk, 2, median, na.rm = TRUE)
  )
)
counts_raw <- round(counts_dedup)
message("Unique genes after symbol collapse: ", nrow(counts_raw))

# ── 2.3 Filter samples: Basal TP vs Normal NT ────────────────

colnames(counts_raw) <- substr(cd$barcode, 1, 16)

basal_bc  <- unique(substr(
  cd$barcode[cd$paper_BRCA_Subtype_PAM50 %in% "Basal" & cd$shortLetterCode %in% "TP"],
  1, 16
))
normal_bc <- unique(substr(
  cd$barcode[cd$shortLetterCode %in% "NT"],
  1, 16
))

basal_bc  <- intersect(basal_bc,  colnames(counts_raw))
normal_bc <- intersect(normal_bc, colnames(counts_raw))

cat("Basal TP samples (TCGA):", length(basal_bc), "\n")
cat("Normal NT samples (TCGA):", length(normal_bc), "\n")

counts_sub <- counts_raw[, c(basal_bc, normal_bc)]
condition  <- factor(
  c(rep("Tumor", length(basal_bc)), rep("Normal", length(normal_bc))),
  levels = c("Normal", "Tumor")
)

# ── 2.4 DESeq2 ───────────────────────────────────────────────

message("Building DESeqDataSet...")
dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData   = data.frame(condition = condition),
  design    = ~ condition
)

# Minimum expression filter: >=10 counts in >=20% of samples (mirrors basal_pre_networks.R)
keep <- rowSums(counts(dds) >= 10) >= 0.2 * ncol(dds)
dds  <- dds[keep, ]
message("Genes after minimum count filter: ", nrow(dds))

message("Running DESeq2 (Basal tumor vs Normal)...")
dds <- DESeq(dds)

res_tcga    <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)
res_tcga_df <- as.data.frame(res_tcga)
res_tcga_df$gene <- rownames(res_tcga_df)
res_tcga_df <- res_tcga_df[!is.na(res_tcga_df$padj), ]

res_tcga_df$is_kbhb <- res_tcga_df$gene %in% kbhb_genes

cat("DE genes TCGA (padj < 0.05):", sum(res_tcga_df$padj < 0.05), "\n")
cat("  Up:  ", sum(res_tcga_df$padj < 0.05 & res_tcga_df$log2FoldChange > 0), "\n")
cat("  Down:", sum(res_tcga_df$padj < 0.05 & res_tcga_df$log2FoldChange < 0), "\n")
cat("  Kbhb genes tested:", sum(res_tcga_df$is_kbhb), "\n")

saveRDS(res_tcga_df, "data/tcga_de_results.rds")
message("data/tcga_de_results.rds saved")


# ============================================================
# SECTION 3: METABRIC — limma
# ============================================================

message("Loading METABRIC matrices (ComBat-corrected)...")
mtbrc_basal  <- readRDS("data/mtbrc_basal_expr.rds")
mtbrc_normal <- readRDS("data/mtbrc_normal_expr.rds")

cat("Basal samples (METABRIC):", ncol(mtbrc_basal), "\n")
cat("Normal samples (METABRIC):", ncol(mtbrc_normal), "\n")

expr_combined <- cbind(mtbrc_basal, mtbrc_normal)
condition_m   <- factor(
  c(rep("Tumor", ncol(mtbrc_basal)), rep("Normal", ncol(mtbrc_normal))),
  levels = c("Normal", "Tumor")
)

design_m     <- model.matrix(~ condition_m)
message("Running limma (eBayes)...")
fit          <- lmFit(expr_combined, design_m)
fit          <- eBayes(fit)
res_mtbrc_df <- topTable(fit, coef = "condition_mTumor", number = Inf, sort.by = "none")
res_mtbrc_df$gene <- rownames(res_mtbrc_df)

res_mtbrc_df$is_kbhb <- res_mtbrc_df$gene %in% kbhb_genes

cat("DE genes METABRIC (adj.P.Val < 0.05):", sum(res_mtbrc_df$adj.P.Val < 0.05), "\n")
cat("  Up:  ", sum(res_mtbrc_df$adj.P.Val < 0.05 & res_mtbrc_df$logFC > 0), "\n")
cat("  Down:", sum(res_mtbrc_df$adj.P.Val < 0.05 & res_mtbrc_df$logFC < 0), "\n")
cat("  Kbhb genes tested:", sum(res_mtbrc_df$is_kbhb), "\n")

saveRDS(res_mtbrc_df, "data/mtbrc_de_results.rds")
message("data/mtbrc_de_results.rds saved")


# ============================================================
# SECTION 4: Figures
# ============================================================

# ── Base theme (consistent with compare_kbhb_mrs.R) ─────────

theme_kbhb <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    legend.position  = "bottom"
  )

# Color palette — sig/ns split for Kbhb genes, TMR orange, rest grey
pal_cat <- c(
  "TMR"        = "#E67E22",
  "Kbhb ↑"    = "#922B21",   # dark red  — significant up
  "Kbhb ↑ ns" = "#F1948A",   # light red — non-significant up
  "Kbhb ↓"    = "#1A5276",   # dark blue — significant down
  "Kbhb ↓ ns" = "#85C1E9",   # light blue — non-significant down
  "Other"      = "grey80"
)

# ── 4.1 Gene classification helper ───────────────────────────
# A Kbhb gene is "significant" only when BOTH padj < fdr_thr AND |LFC| >= fc_thr,
# i.e. it falls in the top-left or top-right quadrant defined by the dashed lines.
# Genes that pass only one criterion stay in the n.s. (light) layer.

classify_genes <- function(df, fdr_col, lfc_col, kbhb_genes, sig_tfs,
                           fdr_thr = 0.05, fc_thr = 0.5) {
  sig  <- df[[fdr_col]] < fdr_thr & abs(df[[lfc_col]]) >= fc_thr
  up   <- df[[lfc_col]] > 0
  kb   <- df$gene %in% kbhb_genes
  df$category <- "Other"
  df$category[ kb &  sig &  up] <- "Kbhb ↑"
  df$category[ kb &  sig & !up] <- "Kbhb ↓"
  df$category[ kb & !sig &  up] <- "Kbhb ↑ ns"
  df$category[ kb & !sig & !up] <- "Kbhb ↓ ns"
  df$category[df$gene %in% sig_tfs] <- "TMR"   # TMR always on top layer
  df$category <- factor(df$category, levels = names(pal_cat))
  df$label    <- ifelse(df$gene %in% sig_tfs, df$gene, NA_character_)
  df
}

# ── 4.2 Generic volcano function ─────────────────────────────

make_volcano <- function(df, lfc_col, fdr_col, title, subtitle,
                         kbhb_genes, sig_tfs, pal, fdr_thr = 0.05, fc_thr = 0.5) {

  df <- classify_genes(df, fdr_col, lfc_col, kbhb_genes, sig_tfs, fdr_thr, fc_thr)
  df$neg_log10p <- -log10(df[[fdr_col]] + 1e-300)

  n_ku    <- sum(df$category == "Kbhb ↑")
  n_kd    <- sum(df$category == "Kbhb ↓")
  n_ku_ns <- sum(df$category == "Kbhb ↑ ns")
  n_kd_ns <- sum(df$category == "Kbhb ↓ ns")

  hline_y <- -log10(fdr_thr)

  ggplot(df, aes(x = .data[[lfc_col]], y = neg_log10p,
                 color = category, label = label)) +
    geom_hline(yintercept = hline_y, linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    geom_point(data = \(d) filter(d, category == "Other"),
               size = 0.4, alpha = 0.25) +
    geom_point(data = \(d) filter(d, category %in% c("Kbhb ↑ ns", "Kbhb ↓ ns")),
               size = 1.2, alpha = 0.55) +
    geom_point(data = \(d) filter(d, category %in% c("Kbhb ↑", "Kbhb ↓")),
               size = 1.6, alpha = 0.85) +
    geom_point(data = \(d) filter(d, category == "TMR"),
               size = 3.5, alpha = 1) +
    geom_text_repel(
      data          = \(d) filter(d, !is.na(label)),
      size          = 3.5, fontface = "bold",
      max.overlaps  = 30, box.padding = 0.5,
      segment.color = "grey40", segment.size = 0.3
    ) +
    scale_color_manual(
      values = pal,
      labels = c(
        "TMR"        = paste0("Sig. MRA TMR (n=", length(sig_tfs), ")"),
        "Kbhb ↑"    = paste0("Kbhb up sig. (n=",    n_ku,    ")"),
        "Kbhb ↑ ns" = paste0("Kbhb up n.s. (n=",    n_ku_ns, ")"),
        "Kbhb ↓"    = paste0("Kbhb down sig. (n=",  n_kd,    ")"),
        "Kbhb ↓ ns" = paste0("Kbhb down n.s. (n=",  n_kd_ns, ")"),
        "Other"      = "Other genes"
      ),
      name = NULL
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                nrow = 2)) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = expression(log[2]~"Fold Change  (Tumor / Normal)"),
      y        = expression(-log[10]~"(FDR)")
    ) +
    theme_kbhb
}

# ── FIG DE-1: Volcano TCGA ───────────────────────────────────

fig_v_tcga <- make_volcano(
  df         = res_tcga_df,
  lfc_col    = "log2FoldChange",
  fdr_col    = "padj",
  title      = "Differential expression — TCGA Basal BRCA",
  subtitle   = paste0("Basal tumor (n=", length(basal_bc),
                      ") vs normal (n=", length(normal_bc),
                      ") · DESeq2 · Kbhb genes and TMRs highlighted"),
  kbhb_genes = kbhb_genes,
  sig_tfs    = sig_tfs,
  pal        = pal_cat
)

ggsave("figures/fig_de_volcano_tcga.pdf", fig_v_tcga, width = 7, height = 6)
message("figures/fig_de_volcano_tcga.pdf saved")


# ── FIG DE-2: Volcano METABRIC ───────────────────────────────

fig_v_mtbrc <- make_volcano(
  df         = res_mtbrc_df,
  lfc_col    = "logFC",
  fdr_col    = "adj.P.Val",
  title      = "Differential expression — METABRIC Basal BRCA",
  subtitle   = paste0("Basal tumor (n=", ncol(mtbrc_basal),
                      ") vs normal (n=", ncol(mtbrc_normal),
                      ") · limma · Kbhb genes and TMRs highlighted"),
  kbhb_genes = kbhb_genes,
  sig_tfs    = sig_tfs,
  pal        = pal_cat
)

ggsave("figures/fig_de_volcano_mtbrc.pdf", fig_v_mtbrc, width = 7, height = 6)
message("figures/fig_de_volcano_mtbrc.pdf saved")



# ── FIG DE-3 & DE-4: Expression heatmaps — DE Kbhb genes ─────
# One heatmap per dataset: rows = DE Kbhb genes, columns = all samples
# (tumor + normal). Values z-scored per gene. Dendrogram on both axes.
# No row labels. Column annotation bar: Tumor vs Normal.

# Shared annotation colours and z-score breaks
# Column named "Condition" (not "Group") to avoid collision with the color scale legend
anno_group_colors <- list(
  Condition = c("Tumor" = "#922B21", "Normal" = "#1A5276")
)
breaks_z <- seq(-2.5, 2.5, length.out = 101)   # cap z-scores for colour saturation

zscore_mat <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.nan(z)] <- 0   # genes with zero variance across samples -> 0
  z
}

# ── FIG DE-3: TCGA ──────────────────────────────────────────

tcga_basal_expr  <- readRDS("data/tcga_basal_expr.rds")
tcga_normal_expr <- readRDS("data/tcga_normal_expr.rds")

kbhb_de_tcga <- res_tcga_df$gene[
  res_tcga_df$gene %in% kbhb_genes & res_tcga_df$padj < 0.05
]
cat("DE Kbhb genes (padj<0.05) TCGA:", length(kbhb_de_tcga), "\n")

tcga_all   <- cbind(tcga_basal_expr, tcga_normal_expr)
tcga_genes <- intersect(kbhb_de_tcga, rownames(tcga_all))
tcga_z     <- zscore_mat(tcga_all[tcga_genes, ])

tcga_col_anno <- data.frame(
  Condition = c(rep("Tumor",  ncol(tcga_basal_expr)),
                rep("Normal", ncol(tcga_normal_expr))),
  row.names = colnames(tcga_all)
)

# Clustering: (1 - Pearson r) distance + Ward D2 linkage.
# Preferred over Euclidean for z-scored expression: groups genes/samples
# by co-expression pattern rather than absolute magnitude.
ph_tcga <- pheatmap(
  tcga_z,
  color                      = colorRampPalette(c("#1A5276", "white", "#922B21"))(100),
  breaks                     = breaks_z,
  cluster_rows               = TRUE,
  cluster_cols               = TRUE,
  clustering_distance_rows   = "correlation",
  clustering_distance_cols   = "correlation",
  clustering_method          = "ward.D2",
  show_rownames              = FALSE,
  show_colnames              = FALSE,
  annotation_col             = tcga_col_anno,
  annotation_colors          = anno_group_colors,
  border_color               = NA,
  main                       = paste0("Kbhb DE genes — TCGA  (n = ", length(tcga_genes), " genes)"),
  silent                     = TRUE
)
pdf("figures/fig_de_heatmap_tcga.pdf", width = 9, height = 7)
grid::grid.draw(ph_tcga$gtable)
dev.off()
message("figures/fig_de_heatmap_tcga.pdf saved")

# ── FIG DE-4: METABRIC ──────────────────────────────────────

kbhb_de_mtbrc <- res_mtbrc_df$gene[
  res_mtbrc_df$gene %in% kbhb_genes & res_mtbrc_df$adj.P.Val < 0.05
]
cat("DE Kbhb genes (adj.P.Val<0.05) METABRIC:", length(kbhb_de_mtbrc), "\n")

# expr_combined = cbind(mtbrc_basal, mtbrc_normal) — already built in Section 3
mtbrc_genes <- intersect(kbhb_de_mtbrc, rownames(expr_combined))
mtbrc_z     <- zscore_mat(expr_combined[mtbrc_genes, ])

mtbrc_col_anno <- data.frame(
  Condition = c(rep("Tumor",  ncol(mtbrc_basal)),
                rep("Normal", ncol(mtbrc_normal))),
  row.names = colnames(expr_combined)
)

ph_mtbrc <- pheatmap(
  mtbrc_z,
  color                      = colorRampPalette(c("#1A5276", "white", "#922B21"))(100),
  breaks                     = breaks_z,
  cluster_rows               = TRUE,
  cluster_cols               = TRUE,
  clustering_distance_rows   = "correlation",
  clustering_distance_cols   = "correlation",
  clustering_method          = "ward.D2",
  show_rownames              = FALSE,
  show_colnames              = FALSE,
  annotation_col             = mtbrc_col_anno,
  annotation_colors          = anno_group_colors,
  border_color               = NA,
  main                       = paste0("Kbhb DE genes — METABRIC  (n = ", length(mtbrc_genes), " genes)"),
  silent                     = TRUE
)
pdf("figures/fig_de_heatmap_mtbrc.pdf", width = 9, height = 7)
grid::grid.draw(ph_mtbrc$gtable)
dev.off()
message("figures/fig_de_heatmap_mtbrc.pdf saved")

# Supplementary figure: 2x2 panel (volcanos top, heatmaps bottom)
# A = volcano TCGA, B = volcano METABRIC, C = heatmap TCGA, D = heatmap METABRIC
fig_supp_de <- (fig_v_tcga | fig_v_mtbrc) /
               (wrap_elements(ph_tcga$gtable) | wrap_elements(ph_mtbrc$gtable)) +
  plot_layout(heights = c(6, 7)) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.title = element_blank())
  )
ggsave("figures/fig_supp_de.pdf", fig_supp_de, width = 14, height = 13)
message("figures/fig_supp_de.pdf saved")


# ============================================================
# SECTION 5: DE summary for the significant TMRs
# ============================================================

cat("\n=== Differential expression of significant TMRs ===\n")

tmr_tcga <- res_tcga_df %>%
  filter(gene %in% sig_tfs) %>%
  select(gene, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  arrange(match(gene, sig_tfs))

tmr_mtbrc <- res_mtbrc_df %>%
  filter(gene %in% sig_tfs) %>%
  select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
  arrange(match(gene, sig_tfs))

cat("\nTCGA (DESeq2):\n")
print(tmr_tcga, row.names = FALSE)
cat("\nMETABRIC (limma):\n")
print(tmr_mtbrc, row.names = FALSE)

write.table(tmr_tcga,  "data/tmr_de_tcga.tsv",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tmr_mtbrc, "data/tmr_de_mtbrc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/tmr_de_tcga.tsv and data/tmr_de_mtbrc.tsv saved")


# ============================================================
# SECTION 6: DE Kbhb gene tables (padj < 0.05 and |LFC| >= 0.5)
# ============================================================

kbhb_sig_tcga <- res_tcga_df %>%
  filter(is_kbhb, padj < 0.05, abs(log2FoldChange) >= 0.5) %>%
  arrange(padj) %>%
  select(gene, log2FoldChange, lfcSE, stat, pvalue, padj, is_kbhb)

kbhb_sig_mtbrc <- res_mtbrc_df %>%
  filter(is_kbhb, adj.P.Val < 0.05, abs(logFC) >= 0.5) %>%
  arrange(adj.P.Val) %>%
  select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, is_kbhb)

cat("Sig. DE Kbhb genes (padj<0.05, |LFC|>=0.5) TCGA:    ", nrow(kbhb_sig_tcga), "\n")
cat("Sig. DE Kbhb genes (adj.P<0.05, |LFC|>=0.5) METABRIC:", nrow(kbhb_sig_mtbrc), "\n")
cat("  In common:", length(intersect(kbhb_sig_tcga$gene, kbhb_sig_mtbrc$gene)), "\n")

write.table(kbhb_sig_tcga,  "data/kbhb_de_tcga.tsv",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(kbhb_sig_mtbrc, "data/kbhb_de_mtbrc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/kbhb_de_tcga.tsv and data/kbhb_de_mtbrc.tsv saved")


# ============================================================
# SECTION 7: Kbhb DE concordance — consensus table, scatter,
#            and UpSet plot
# ============================================================

# ── 7.1 Consensus table ──────────────────────────────────────
# All Kbhb genes present in at least one DE result

kbhb_all <- union(
  res_tcga_df$gene[res_tcga_df$is_kbhb],
  res_mtbrc_df$gene[res_mtbrc_df$is_kbhb]
)

lfc_t  <- setNames(res_tcga_df$log2FoldChange, res_tcga_df$gene)
padj_t <- setNames(res_tcga_df$padj,           res_tcga_df$gene)
lfc_m  <- setNames(res_mtbrc_df$logFC,         res_mtbrc_df$gene)
padj_m <- setNames(res_mtbrc_df$adj.P.Val,     res_mtbrc_df$gene)

consensus <- data.frame(
  gene       = kbhb_all,
  lfc_tcga   = lfc_t[kbhb_all],
  padj_tcga  = padj_t[kbhb_all],
  lfc_mtbrc  = lfc_m[kbhb_all],
  padj_mtbrc = padj_m[kbhb_all],
  row.names  = NULL,
  stringsAsFactors = FALSE
) %>%
  mutate(
    sig_tcga   = !is.na(padj_tcga)  & padj_tcga  < 0.05 & abs(lfc_tcga)  >= 0.5,
    sig_mtbrc  = !is.na(padj_mtbrc) & padj_mtbrc < 0.05 & abs(lfc_mtbrc) >= 0.5,
    concordant = sig_tcga & sig_mtbrc & sign(lfc_tcga) == sign(lfc_mtbrc),
    category   = case_when(
      sig_tcga & sig_mtbrc & lfc_tcga > 0 & lfc_mtbrc > 0 ~ "Concordant up",
      sig_tcga & sig_mtbrc & lfc_tcga < 0 & lfc_mtbrc < 0 ~ "Concordant down",
      sig_tcga & sig_mtbrc                                  ~ "Discordant",
      sig_tcga & !sig_mtbrc                                 ~ "TCGA only",
      !sig_tcga & sig_mtbrc                                 ~ "METABRIC only",
      TRUE                                                   ~ "Not DE"
    )
  )

cat("\nKbhb DE concordance summary:\n")
print(table(consensus$category))

write.table(consensus, "data/kbhb_consensus.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("data/kbhb_consensus.tsv saved")


# ── 7.2 Scatter log2FC TCGA vs METABRIC (Kbhb genes) ────────

pal_conc <- c(
  "Concordant up"   = "#922B21",
  "Concordant down" = "#1A5276",
  "Discordant"      = "#8E44AD",
  "TCGA only"       = "#E67E22",
  "METABRIC only"   = "#27AE60",
  "Not DE"          = "grey80"
)

sc_df <- consensus %>%
  filter(!is.na(lfc_tcga) & !is.na(lfc_mtbrc)) %>%
  mutate(category = factor(category, levels = names(pal_conc)))

# Label strategy: TMRs in the scatter (only those that are Kbhb substrates)
# + top concordant genes to complete n_labels total, ranked by |LFC_TCGA| + |LFC_METABRIC|
n_labels     <- 10
tmr_in_sc    <- sc_df$gene[sc_df$gene %in% sig_tfs]
top_conc     <- sc_df %>%
  filter(concordant, !gene %in% sig_tfs) %>%
  mutate(score = abs(lfc_tcga) + abs(lfc_mtbrc)) %>%
  arrange(desc(score)) %>%
  slice_head(n = n_labels - length(tmr_in_sc)) %>%
  pull(gene)

sc_df <- sc_df %>%
  mutate(label = case_when(
    gene %in% sig_tfs  ~ gene,
    gene %in% top_conc ~ gene,
    TRUE               ~ NA_character_
  ))

r_val <- cor(sc_df$lfc_tcga, sc_df$lfc_mtbrc, use = "complete.obs")

fig_scatter <- ggplot(sc_df, aes(x = lfc_tcga, y = lfc_mtbrc,
                                  color = category, label = label)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55", linewidth = 0.4) +
  geom_point(data = \(d) filter(d, category == "Not DE"),
             size = 0.8, alpha = 0.3) +
  geom_point(data = \(d) filter(d, !category %in%
               c("Not DE", "Concordant up", "Concordant down")),
             size = 1.4, alpha = 0.75) +
  geom_point(data = \(d) filter(d, category %in%
               c("Concordant up", "Concordant down")),
             size = 1.8, alpha = 0.9) +
  geom_point(data = \(d) filter(d, !is.na(label)),
             size = 3.5, alpha = 1) +
  geom_text_repel(                            # top concordant: text = point color
    data          = \(d) filter(d, !is.na(label) & !gene %in% sig_tfs),
    aes(color     = category),
    bg.color      = "white",
    bg.r          = 0.15,
    size          = 3, fontface = "bold",
    max.overlaps  = Inf,
    segment.color = "grey40", segment.size = 0.3,
    box.padding   = 0.5,      point.size   = NA
  ) +
  geom_text_repel(                            # TMR genes: white text
    data          = \(d) filter(d, !is.na(label) & gene %in% sig_tfs),
    color         = "white",
    bg.color      = scales::alpha("grey15", 0.70),
    bg.r          = 0.15,
    size          = 3, fontface = "bold",
    max.overlaps  = Inf,
    segment.color = "grey40", segment.size = 0.3,
    box.padding   = 0.5,      point.size   = NA
  ) +
  scale_color_manual(values = pal_conc, name = NULL) +
  annotate("text",
           x = min(sc_df$lfc_tcga, na.rm = TRUE) * 0.85,
           y = max(sc_df$lfc_mtbrc, na.rm = TRUE) * 0.90,
           label = paste0("r = ", round(r_val, 2)),
           hjust = 0, size = 4, fontface = "italic") +
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 1), nrow = 2)) +
  labs(
    title    = "Kbhb gene concordance — TCGA vs METABRIC",
    subtitle = paste0("log₂FC Basal tumor / normal · Kbhb genes (n=", nrow(sc_df),
                      ") · top ", n_labels, " concordant genes labelled"),
    x        = expression(log[2]~"FC  (TCGA)"),
    y        = expression(log[2]~"FC  (METABRIC)")
  ) +
  theme_kbhb

ggsave("figures/fig_kbhb_scatter.pdf", fig_scatter, width = 7, height = 6)
message("figures/fig_kbhb_scatter.pdf saved")


# ── 7.3 Circular Packing — Kbhb genes by DE category ─────────
# Hierarchy: root → {Not DE, DE container → subcategories}
# Subcategories split TCGA/METABRIC only by direction (up/down).
# Circle area is proportional to gene count.

cp_counts <- consensus %>%
  mutate(cp_cat = case_when(
    category == "Concordant up"                  ~ "Concordant up",
    category == "Concordant down"                ~ "Concordant down",
    category == "TCGA only"     & lfc_tcga  > 0 ~ "TCGA only up",
    category == "TCGA only"                      ~ "TCGA only down",
    category == "METABRIC only" & lfc_mtbrc > 0 ~ "METABRIC only up",
    category == "METABRIC only"                  ~ "METABRIC only down",
    category == "Discordant"                     ~ "Discordant",
    TRUE                                         ~ "Not DE"
  )) %>%
  count(cp_cat, name = "n") %>%
  filter(n > 0)

de_present <- cp_counts$cp_cat[cp_counts$cp_cat != "Not DE"]
not_de_n   <- cp_counts$n[cp_counts$cp_cat == "Not DE"]
de_total   <- sum(cp_counts$n[cp_counts$cp_cat != "Not DE"])

edges_cp <- data.frame(
  from = c("root", "root", rep("DE", length(de_present))),
  to   = c("Not DE", "DE", de_present),
  stringsAsFactors = FALSE
)

nodes_cp <- data.frame(
  name = c("root", "Not DE", "DE", de_present),
  size = c(0L, not_de_n, 0L,
           cp_counts$n[match(de_present, cp_counts$cp_cat)]),
  stringsAsFactors = FALSE
)

pal_cp <- c(
  "root"               = NA,
  "DE"                 = NA,
  "Not DE"             = "grey85",
  "Concordant up"  = "#922B21",
  "Concordant down"  = "#1A5276",
  "TCGA only up"   = "#E59866",
  "TCGA only down"   = "#85C1E9",
  "METABRIC only up" = "#A9DFBF",
  "METABRIC only down" = "#27AE60",
  "Discordant"         = "#8E44AD"
)

graph_cp <- tbl_graph(nodes = nodes_cp, edges = edges_cp, directed = TRUE) %>%
  activate(nodes) %>%
  mutate(
    is_leaf    = node_is_leaf(),
    label      = ifelse(is_leaf & size >= 8,
                        paste0(name, "\n(n=", size, ")"), ""),
    de_label   = ifelse(name == "DE",
                        paste0("Sig. DE\n(n=", de_total, ")"), ""),
    text_col   = "white"
  )

# Pre-compute layout to access node x,y for geom_text_repel
lay_cp <- create_layout(graph_cp, layout = "circlepack", weight = size)

fig_cp <- ggraph(lay_cp) +
  geom_node_circle(                           # DE container: solid outline
    aes(filter = name == "DE"),
    fill = NA, color = "grey40", linetype = "solid", linewidth = 1
  ) +
  geom_node_circle(                           # leaf bubbles: filled
    aes(fill = name, filter = is_leaf),
    color = "white", linewidth = 0.35
  ) +
  ggrepel::geom_text_repel(                  # repelled labels — avoids overlap
    data          = dplyr::filter(lay_cp, is_leaf & size >= 8),
    aes(x = x, y = y, label = label),
    color         = "white",
    bg.color      = scales::alpha("grey15", 0.70),
    bg.r          = 0.15,
    size          = 2.5, fontface = "bold", lineheight = 0.9,
    max.overlaps  = Inf,
    segment.color = "grey50", segment.size  = 0.3,
    box.padding   = 0.4,      point.size    = NA
  ) +
  ggrepel::geom_text_repel(                  # DE container label — same bg trick
    data          = dplyr::filter(lay_cp, name == "DE"),
    aes(x = x, y = y, label = de_label),
    color         = "white",
    bg.color      = scales::alpha("grey15", 0.70),
    bg.r          = 0.15,
    size          = 3, fontface = "bold", lineheight = 0.9,
    force         = 0,           # stay at centroid, do not repel
    segment.color = NA,          # no connector line
    point.size    = NA
  ) +
  scale_fill_manual(
    values   = pal_cp,
    na.value = "transparent",
    breaks   = c("Concordant up", "Concordant down",
                 "TCGA only up",  "TCGA only down",
                 "METABRIC only up", "METABRIC only down",
                 "Discordant", "Not DE"),
    name     = NULL
  ) +
  coord_equal() +
  theme_void() +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3))) +
  labs(
    title    = "Kbhb genes — DE concordance (Circular Packing)",
    subtitle = paste0("n = ", nrow(consensus), " Kbhb genes · padj<0.05 and |LFC|>=0.5")
  )

ggsave("figures/fig_kbhb_circlepack.pdf", fig_cp, width = 8, height = 8)
message("figures/fig_kbhb_circlepack.pdf saved")

# Main concordance panel: A = circular packing (left), B = scatter (right)
fig_concordance <- (fig_cp | fig_scatter) +
  plot_layout(widths = c(0.8, 1.2)) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(plot.title = element_blank())
  )
ggsave("figures/fig_kbhb_concordance.pdf", fig_concordance, width = 15, height = 10)
message("figures/fig_kbhb_concordance.pdf saved")


message("\n====  DIFFERENTIAL EXPRESSION COMPLETE  ====")
message("Tables:")
message("  data/tcga_de_results.rds   — full DESeq2 results  (is_kbhb column)")
message("  data/mtbrc_de_results.rds  — full limma results   (is_kbhb column)")
message("  data/tmr_de_tcga.tsv       — DE stats for sig. TMRs (TCGA)")
message("  data/tmr_de_mtbrc.tsv      — DE stats for sig. TMRs (METABRIC)")
message("  data/kbhb_de_tcga.tsv      — Kbhb DE genes padj<0.05 |LFC|>=0.5 (TCGA)")
message("  data/kbhb_de_mtbrc.tsv     — Kbhb DE genes adj.P<0.05 |LFC|>=0.5 (METABRIC)")
message("  data/kbhb_consensus.tsv    — all Kbhb genes with concordance annotation")
message("Figures (individual):")
message("  figures/fig_de_volcano_tcga.pdf")
message("  figures/fig_de_volcano_mtbrc.pdf")
message("  figures/fig_de_heatmap_tcga.pdf")
message("  figures/fig_de_heatmap_mtbrc.pdf")
message("  figures/fig_kbhb_scatter.pdf")
message("  figures/fig_kbhb_circlepack.pdf")
message("Figures (combined panels):")
message("  figures/fig_supp_de.pdf          — A: volcano TCGA | B: volcano METABRIC | C: heatmap TCGA | D: heatmap METABRIC  (supplementary)")
message("  figures/fig_kbhb_concordance.pdf — A: circlepack | B: scatter  (main)")
