
# circos_tmr_kbhb.R
# Circos: 7 Kbhb TMRs → Kbhb DE genes
#
# Edges: ARACNe TCGA Basal network (TMR → Kbhb DE gene)
# Sectors:
#   – TMR arc  : 7 TMRs (emitters, one color per TMR)
#   – Gene arc : Kbhb DE genes (colored by DE category, ordered by category
#                then descending n_TMRs)
# Tracks:
#   – Outer : labels (all TMRs bold; genes with ≥2 TMRs labeled small)
#   – Inner  : bar = n TMRs targeting each gene (gene sectors only)
# ============================================================


# ============================================================
# SECTION 1: Libraries
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(janitor)
  library(vroom)
  library(circlize)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(grid)
})

# setwd("/path/to/Kbhb")  # set your project directory
dir.create("figures", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)


# ============================================================
# SECTION 2: Load data
# ============================================================

# ── 2.1 Significant TMRs ─────────────────────────────────────
# meta-NES descending order (CENPA highest; VEZF1 sole repressed TMR)
sig_tfs   <- c("CENPA","FOXM1","HMGA1","CHCHD3","E2F7","ZNF232","VEZF1")
tmr_order <- c("CENPA","FOXM1","HMGA1","CHCHD3","E2F7","ZNF232","VEZF1")

# ── 2.2 Kbhb DE consensus ────────────────────────────────────
consensus <- read_tsv("data/kbhb_consensus.tsv", show_col_types = FALSE)

cat("Columns in kbhb_consensus.tsv:", paste(names(consensus), collapse = ", "), "\n")
cat("Category values:\n")
print(table(consensus$category))

# Harmonize: split "TCGA only" / "METABRIC only" into up/down if not already split
if (!"TCGA only up" %in% consensus$category && "TCGA only" %in% consensus$category) {
  lfc_col <- intersect(c("LFC_TCGA","log2FoldChange_TCGA","lfc_tcga"), names(consensus))[1]
  consensus <- consensus %>%
    mutate(category = case_when(
      category == "TCGA only"     & get(lfc_col) > 0 ~ "TCGA only up",
      category == "TCGA only"     & get(lfc_col) < 0 ~ "TCGA only down",
      category == "METABRIC only" & get(lfc_col) > 0 ~ "METABRIC only up",
      category == "METABRIC only" & get(lfc_col) < 0 ~ "METABRIC only down",
      TRUE ~ category
    ))
}

de_categories <- c("Concordant up","Concordant down",
                   "TCGA only up","TCGA only down",
                   "METABRIC only up","METABRIC only down")

kbhb_de <- consensus %>%
  filter(category %in% de_categories) %>%
  dplyr::select(gene, category)

cat("\nKbhb DE genes by category:\n")
print(sort(table(kbhb_de$category), decreasing = TRUE))
cat("Total DE genes:", nrow(kbhb_de), "\n")


# ============================================================
# SECTION 3: ARACNe filter — TCGA Basal network only
# ============================================================
# Retain edges where the regulator is one of the 7 TMRs and
# the target is a differentially expressed Kbhb gene.

aracne_raw <- vroom("data/tcga_basal_network.txt",
                    delim = "\t", col_types = cols(.default = "c"),
                    show_col_types = FALSE) %>%
  clean_names()

# Flexible column detection (ARACNe-AP output header varies by version)
reg_col <- intersect(c("regulator","tf","regulator_gene","source"), names(aracne_raw))[1]
tgt_col <- intersect(c("target","target_gene","destination"),       names(aracne_raw))[1]
mi_col  <- intersect(c("mi","mutual_information"),                  names(aracne_raw))[1]

if (is.na(reg_col) || is.na(tgt_col)) {
  aracne <- aracne_raw %>%
    rename(tf = 1, target = 2, mi = 3) %>%
    mutate(mi = as.numeric(mi))
} else {
  aracne <- aracne_raw %>%
    rename(tf = all_of(reg_col), target = all_of(tgt_col)) %>%
    mutate(mi = as.numeric(.data[[mi_col]]))
}

aracne_filt <- aracne %>%
  filter(tf     %in% sig_tfs,
         target %in% kbhb_de$gene)

cat("\nARACNe edges (TMR → Kbhb DE gene):", nrow(aracne_filt), "\n")
cat("Unique Kbhb genes targeted        :", n_distinct(aracne_filt$target), "\n")
cat("Edges per TMR:\n"); print(sort(table(aracne_filt$tf), decreasing = TRUE))

write_tsv(aracne_filt %>% dplyr::select(tf, target, mi),
          "data/circos_aracne_tmr_kbhb.tsv")
message("data/circos_aracne_tmr_kbhb.tsv saved")


# ============================================================
# SECTION 4: Gene summary
# ============================================================

interactions <- aracne_filt %>%
  rename(gene = target) %>%
  left_join(kbhb_de, by = "gene") %>%
  filter(!is.na(category))

gene_summary <- interactions %>%
  group_by(gene, category) %>%
  summarise(n_tmrs = n_distinct(tf),
            tmrs   = paste(sort(tf), collapse = ","),
            .groups = "drop") %>%
  arrange(category, desc(n_tmrs), gene)

cat("\n── Interactions (ARACNe TCGA Basal) ─────────────\n")
cat("Total edges     :", nrow(interactions), "\n")
cat("Unique genes    :", n_distinct(interactions$gene), "\n")
cat("Per DE category:\n"); print(sort(table(interactions$category), decreasing = TRUE))
cat("\nn_TMRs per gene distribution:\n")
print(table(gene_summary$n_tmrs))

write_tsv(gene_summary, "data/circos_gene_summary.tsv")
write_tsv(interactions, "data/circos_interactions_final.tsv")
message("data/circos_gene_summary.tsv + circos_interactions_final.tsv saved")

# Pre-compute display labels to avoid repeated lookup inside panel.fun
# Gene labels: "SYMBOL (n_tmrs)" — single line, no overlap risk from \n
gene_summary <- gene_summary %>%
  mutate(label = paste0(gene, " (", n_tmrs, ")"))
gene_lbl <- setNames(gene_summary$label, gene_summary$gene)

# TMR labels: "SYMBOL (n_targets)"
tmr_n_map  <- interactions %>%
  group_by(tf) %>% summarise(n = n_distinct(gene), .groups = "drop")
tmr_lbl <- setNames(
  paste0(tmr_n_map$tf, " (", tmr_n_map$n[match(tmr_order, tmr_n_map$tf)], ")"),
  tmr_order
)


# ============================================================
# SECTION 5: Circos plot
# ============================================================

# ── 5.1 Sector order and colors ──────────────────────────────

tmr_plot <- tmr_order   # all 7 TMRs included

cat_order <- c("Concordant up","Concordant down",
               "TCGA only up","TCGA only down",
               "METABRIC only up","METABRIC only down")

gene_order <- gene_summary %>%
  mutate(category = factor(category, levels = cat_order)) %>%
  arrange(category, desc(n_tmrs), gene) %>%
  pull(gene)

cat_col <- c(
  "Concordant up"      = "#C0392B",
  "Concordant down"    = "#2980B9",
  "TCGA only up"       = "#E67E22",
  "TCGA only down"     = "#85C1E9",
  "METABRIC only up"   = "#27AE60",
  "METABRIC only down" = "#82E0AA"
)

n_tmr_plot <- length(tmr_plot)
tmr_col <- setNames(
  colorRampPalette(brewer.pal(min(n_tmr_plot, 9), "Set1"))(n_tmr_plot),
  tmr_plot
)

gene_cat <- setNames(gene_summary$category[match(gene_order, gene_summary$gene)], gene_order)
gene_col <- cat_col[gene_cat]
names(gene_col) <- gene_order

all_sectors <- c(tmr_plot, gene_order)
all_col     <- c(tmr_col, gene_col)

# ── 5.2 Gap vector ────────────────────────────────────────────
# 8° gap at TMR/gene arc boundaries; 3° between DE categories; 0.5° within groups

gap_vec <- setNames(rep(0.5, length(all_sectors)), all_sectors)
gap_vec[tmr_plot[length(tmr_plot)]]     <- 8
gap_vec[gene_order[length(gene_order)]] <- 8

gene_cats_ordered <- gene_cat[gene_order]
cat_change <- which(c(FALSE, gene_cats_ordered[-1] != gene_cats_ordered[-length(gene_cats_ordered)]))
if (length(cat_change) > 0)
  gap_vec[gene_order[cat_change - 1]] <- 3

# ── 5.3 Interaction matrix ────────────────────────────────────
# Binary: 1 = ARACNe edge exists between TMR (row) and gene (col)

mat <- interactions %>%
  distinct(tf, gene) %>%
  mutate(edge = 1L) %>%
  complete(tf = tmr_plot, gene = gene_order, fill = list(edge = 0L)) %>%
  pivot_wider(names_from = gene, values_from = edge, values_fill = 0L) %>%
  column_to_rownames("tf") %>%
  as.matrix()
mat <- mat[tmr_plot, gene_order, drop = FALSE]

col_mat <- matrix(tmr_col[rownames(mat)],
                  nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = dimnames(mat))
col_mat[mat == 0] <- "transparent"

# ── 5.4 Draw ──────────────────────────────────────────────────

bar_ramp <- colorRampPalette(c("grey88","#C0392B"))(n_tmr_plot + 1)

# Circos body only — no legend, no title (added in save blocks via layout)
draw_circos <- function() {
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.after    = gap_vec,
    track.margin = c(0.01, 0.01),
    cell.padding = c(0, 0, 0, 0)
  )

  chordDiagram(
    x                 = mat,
    order             = all_sectors,
    grid.col          = all_col,
    col               = col_mat,
    link.arr.col      = col_mat,
    link.border       = NA,
    directional       = 1,
    direction.type    = "arrows",
    link.arr.length   = 0.04,
    link.arr.width    = 0.04,
    transparency      = 0.35,
    annotationTrack   = "grid",
    preAllocateTracks = list(
      list(track.height = uh(6, "mm")),   # track 1 (outer): labels
      list(track.height = uh(4, "mm"))    # track 2 (inner): n_tmrs bar
    )
  )

  # Track 2 (inner): n_tmrs bar on gene sectors; solid band on TMR sectors
  circos.trackPlotRegion(
    track.index = 2,
    bg.border   = NA,
    panel.fun   = function(x, y) {
      sec <- get.cell.meta.data("sector.index")
      xl  <- get.cell.meta.data("xlim")
      if (sec %in% gene_order) {
        n <- gene_summary$n_tmrs[gene_summary$gene == sec]
        if (length(n) == 0 || is.na(n)) n <- 0L
        circos.rect(xl[1], 0, xl[2], n / n_tmr_plot,
                    col = bar_ramp[n + 1], border = NA)
      } else {
        circos.rect(xl[1], 0, xl[2], 1, col = tmr_col[sec], border = NA)
      }
    }
  )

  # Track 1 (outer): TMRs always (with n_targets); genes if targeted by ≥3 TMRs (with n_tmrs)
  circos.trackPlotRegion(
    track.index = 1,
    bg.border   = NA,
    panel.fun   = function(x, y) {
      sec    <- get.cell.meta.data("sector.index")
      xl     <- get.cell.meta.data("xlim")
      is_tmr <- sec %in% tmr_plot
      n_t    <- if (sec %in% gene_summary$gene)
                  gene_summary$n_tmrs[gene_summary$gene == sec] else 0L
      show_label <- is_tmr || (!is.na(n_t) && n_t >= 3)
      if (show_label) {
        lbl <- if (is_tmr) tmr_lbl[sec] else gene_lbl[sec]
        circos.text(
          x          = mean(xl),
          y          = 0.5,
          labels     = lbl,
          facing     = "clockwise",
          niceFacing = TRUE,
          adj        = c(0, 0.5),
          cex        = if (is_tmr) 0.80 else 0.42,
          font       = if (is_tmr) 2L else 1L
        )
      }
    }
  )
}

# Legend object (drawn in a separate panel)
make_legend <- function() {
  lgd_tmr <- Legend(
    labels    = tmr_plot,
    legend_gp = gpar(fill = tmr_col[tmr_plot], col = NA),
    title     = "TMR",
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  lgd_cat <- Legend(
    labels    = names(cat_col),
    legend_gp = gpar(fill = cat_col, col = NA),
    title     = "DE category",
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  lgd_bar <- Legend(
    labels    = paste0(seq_len(n_tmr_plot), " TMR(s)"),
    legend_gp = gpar(fill = bar_ramp[seq_len(n_tmr_plot) + 1], col = NA),
    title     = "n TMRs targeting gene",
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  packLegend(lgd_tmr, lgd_cat, lgd_bar, direction = "vertical")
}

# ── 5.5 Save ──────────────────────────────────────────────────
# draw() from ComplexHeatmap always uses DEVICE NPC coordinates (root viewport),
# not par(fig) regions — so layout() and par(fig) have no effect on legend position.
#
# On a 20×14" device, circlize auto-widens its canvas in x to keep the circle round
# (height-limited at 14"). The circos circle + label tracks occupy roughly
# x = 0.15 – 0.83 NPC, leaving ~17% empty on the right.
# We simply place the legend at x ≈ 0.85 NPC, inside that empty right margin.

save_circos_figure <- function(file, ...) {
  do.call(if (grepl("\\.pdf$", file)) pdf else png, c(list(file), list(...)))
  par(mar = c(1, 1, 2, 1))
  draw_circos()
  # No embedded title: panel/figure captions belong in the manuscript
  # figure legend, not baked into the image (same convention as elsewhere).
  draw(make_legend(),
       x = unit(0.85, "npc"), y = unit(0.50, "npc"), just = c("left", "center"))
  dev.off()
}

save_circos_figure("figures/fig_circos_tmr_kbhb.pdf", width = 20, height = 14)
message("figures/fig_circos_tmr_kbhb.pdf saved")

# PNG: 20×14 inches at 400 DPI = 8000×5600 px
save_circos_figure("figures/fig_circos_tmr_kbhb.png", width = 8000, height = 5600, res = 400)
message("figures/fig_circos_tmr_kbhb.png saved")


message("\n====  CIRCOS COMPLETE  ====")
message("Figures  : figures/fig_circos_tmr_kbhb.pdf  /  .png")
message("Data tables:")
message("  data/circos_aracne_tmr_kbhb.tsv     — ARACNe-filtered TMR→gene edges (with MI)")
message("  data/circos_gene_summary.tsv         — per-gene: n_tmrs, which TMRs, DE category")
message("  data/circos_interactions_final.tsv   — full edge table with DE category")


# ============================================================
# SECTION 6: Sankey diagram — TMR → DE gene categories
# ============================================================
# Two horizontal axes: TMR (top) → DE gene categories (bottom).
# Flow width proportional to number of ARACNe-supported interactions.
# Colors match the circos palette: tmr_col for TMRs, cat_col for categories.

suppressPackageStartupMessages(library(ggalluvial))

# ── Aggregate interactions by TMR × DE category ──────────────
sankey_df <- interactions %>%
  count(tf, category, name = "n") %>%
  mutate(
    # TMR: meta-NES descending → CENPA at top of the TMR bar
    tf       = factor(tf,       levels = rev(tmr_order)),
    category = factor(category, levels = rev(de_categories))  # rev for coord_flip reading order
  )

# Named color palette covering both TMR names and category names
sankey_pal <- c(tmr_col, cat_col)

# ── Draw ─────────────────────────────────────────────────────
fig_sankey <- ggplot(
  sankey_df,
  aes(y = n, axis1 = category, axis2 = tf)
) +
  geom_alluvium(
    aes(fill = tf),
    width    = 1/5,
    alpha    = 0.80,
    knot.pos = 0.4
  ) +
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width = 1/5, color = "grey30", linewidth = 0.4
  ) +
  geom_label(
    stat  = "stratum",
    aes(label = after_stat(stratum)),
    size  = 3.4,
    label.size    = 0,
    label.padding = unit(0.18, "lines"),
    fill  = "white", alpha = 0.85
  ) +
  scale_fill_manual(values = sankey_pal, guide = "none") +
  scale_x_discrete(
    limits = c("DE category", "TMR"),   # axis1=category at left/bottom, axis2=tf at right/top
    expand = c(0.12, 0.12)
  ) +
  # coord_flip: x (left↔right) ↔ y (top↔bottom)
  # axis1 (TMR)      → top after flip
  # axis2 (category) → bottom after flip
  coord_flip() +
  theme_void(base_size = 12) +
  theme(
    axis.text.y  = element_text(size = 11, face = "bold",
                                margin = margin(r = 6), hjust = 1),
    plot.title   = element_text(hjust = 0.5, size = 12,
                                margin = margin(b = 10)),
    plot.margin  = margin(20, 30, 20, 30)
  ) +
  labs(title = "Kbhb TMR regulons → DE Kbhb gene categories  (ARACNe TCGA Basal)")

ggsave("figures/Supplementary_FigureS5_sankey_tmr_kbhb.pdf", fig_sankey, width = 10, height = 7)
ggsave("figures/Supplementary_FigureS5_sankey_tmr_kbhb.png", fig_sankey, width = 10, height = 7, dpi = 300)
message("figures/Supplementary_FigureS5_sankey_tmr_kbhb.pdf + .png saved")
