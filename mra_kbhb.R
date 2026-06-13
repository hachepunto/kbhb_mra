library(viper)

# Set your project directory:
# setwd("/path/to/Kbhb")
if (!interactive()) setwd(normalizePath("."))

# ============================================================
# MRA DEL PROGRAMA TRANSCRIPCIONAL Kbhb — Basal BRCA
#
# Firma molecular: ~800-1,000 sustratos Kbhb expresados en BRCA
# (Huang et al. 2021, Sci Adv 7:eabe2771)
# Redes: ARACNe-AP Basal TCGA + METABRIC
# Análisis: msVIPER (viper package)
#
# Salidas en data/:
#   kbhb_geneset.rds        — vector de símbolos HGNC del proteoma Kbhb
#   tcga_basal_mrs.rds      — objeto msVIPER TCGA
#   mtbrc_basal_mrs.rds     — objeto msVIPER METABRIC
# ============================================================


# ============================================================
# SECCIÓN 1: Gene set Kbhb — Huang 2021 Tabla S1
# ============================================================

s1 <- read.delim(
  "Huang et al. 2021 Sci Adv/abe2771_table_s1.txt",
  skip             = 17,
  header           = TRUE,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# "Gene names" tiene símbolos separados por ";" — extraer únicos
kbhb_genes <- unique(trimws(unlist(strsplit(
  na.omit(s1[["Gene names"]]), ";"
))))
kbhb_genes <- kbhb_genes[kbhb_genes != ""]

cat("Kbhb genes (Huang S1):", length(kbhb_genes), "\n")
saveRDS(kbhb_genes, "data/kbhb_geneset.rds")


# ============================================================
# SECCIÓN 2: Cargar matrices de expresión
# ============================================================

tcga_basal  <- readRDS("data/tcga_basal_expr.rds")
tcga_normal <- readRDS("data/tcga_normal_expr.rds")
mtbrc_basal  <- readRDS("data/mtbrc_basal_expr.rds")
mtbrc_normal <- readRDS("data/mtbrc_normal_expr.rds")

# Genes Kbhb expresados en cada dataset
kbhb_tcga  <- intersect(kbhb_genes, rownames(tcga_basal))
kbhb_mtbrc <- intersect(kbhb_genes, rownames(mtbrc_basal))

cat("Kbhb genes en TCGA Basal:    ", length(kbhb_tcga),  "\n")
cat("Kbhb genes en METABRIC Basal:", length(kbhb_mtbrc), "\n")


# ============================================================
# SECCIÓN 3: Construir regulones desde redes ARACNe-AP
# ============================================================
# La red tiene 4 columnas: Regulator, Target, MI, pvalue
# aracne2regulon con format="3col" lee las primeras 3 (ignora pvalue)
# Requiere la matriz de expresión con la que se construyó la red

message("Construyendo regulon TCGA...")
tcga_regulon <- aracne2regulon(
  afile  = "data/tcga_basal_network.txt",
  eset   = tcga_basal
)
cat("Regulon TCGA — reguladores:", length(tcga_regulon), "\n")

message("Construyendo regulon METABRIC...")
mtbrc_regulon <- aracne2regulon(
  afile  = "data/mtbrc_basal_network.txt",
  eset   = mtbrc_basal
)
cat("Regulon METABRIC — reguladores:", length(mtbrc_regulon), "\n")


# ============================================================
# SECCIÓN 4: Firma Kbhb restringida (t-estadísticos tumor vs normal)
# ============================================================
# Se usan SOLO los sustratos Kbhb como firma para msVIPER.
# Esto identifica los TFs que regulan preferentemente el programa
# Kbhb, no los reguladores generales del tumor.

message("Calculando firma Kbhb TCGA (tumor Basal vs Normal)...")
tcga_sig  <- rowTtest(
  x = tcga_basal[kbhb_tcga, ],
  y = tcga_normal[kbhb_tcga, ]
)

message("Calculando firma Kbhb METABRIC...")
mtbrc_sig <- rowTtest(
  x = mtbrc_basal[kbhb_mtbrc, ],
  y = mtbrc_normal[kbhb_mtbrc, ]
)


# ============================================================
# SECCIÓN 5: Modelo nulo por permutación
# ============================================================

message("Generando modelo nulo TCGA (1000 permutaciones)...")
tcga_null <- ttestNull(
  x     = tcga_basal[kbhb_tcga, ],
  y     = tcga_normal[kbhb_tcga, ],
  per   = 1000,
  repos = TRUE,
  seed  = 1
)

message("Generando modelo nulo METABRIC...")
mtbrc_null <- ttestNull(
  x     = mtbrc_basal[kbhb_mtbrc, ],
  y     = mtbrc_normal[kbhb_mtbrc, ],
  per   = 1000,
  repos = TRUE,
  seed  = 1
)


# ============================================================
# SECCIÓN 6: msVIPER
# ============================================================

message("Corriendo msVIPER TCGA...")
tcga_mrs <- msviper(
  ges       = tcga_sig$statistic,
  regulon   = tcga_regulon,
  nullmodel = tcga_null,
  minsize   = 25,
  verbose   = FALSE
)
tcga_mrs <- ledge(tcga_mrs)
tcga_mrs <- shadow(tcga_mrs, minsize = 25, verbose = FALSE)

message("Corriendo msVIPER METABRIC...")
mtbrc_mrs <- msviper(
  ges       = mtbrc_sig$statistic,
  regulon   = mtbrc_regulon,
  nullmodel = mtbrc_null,
  minsize   = 25,
  verbose   = FALSE
)
mtbrc_mrs <- ledge(mtbrc_mrs)
mtbrc_mrs <- shadow(mtbrc_mrs, minsize = 25, verbose = FALSE)


# ============================================================
# SECCIÓN 7: Guardar resultados
# ============================================================

saveRDS(tcga_mrs,  "data/tcga_basal_mrs.rds")
saveRDS(mtbrc_mrs, "data/mtbrc_basal_mrs.rds")
message("data/tcga_basal_mrs.rds y data/mtbrc_basal_mrs.rds guardados")


# ── Vista rápida de resultados ───────────────────────────────

cat("\n=== Top 20 Master Regulators — TCGA Basal ===\n")
print(summary(tcga_mrs, mrs = 20))

cat("\n=== Top 20 Master Regulators — METABRIC Basal ===\n")
print(summary(mtbrc_mrs, mrs = 20))


# ============================================================
# SECCIÓN 8: Supplementary Figure S1 — msVIPER + Shadow network
#            Panel A: TCGA barplot    · Panel B: METABRIC barplot
#            Panel C: TCGA shadow     · Panel D: METABRIC shadow
# ============================================================

library(ggplot2)
library(patchwork)
library(ggraph)
library(tidygraph)
library(dplyr)

# Load results if running this section standalone
if (!exists("tcga_mrs"))  tcga_mrs  <- readRDS("data/tcga_basal_mrs.rds")
if (!exists("mtbrc_mrs")) mtbrc_mrs <- readRDS("data/mtbrc_basal_mrs.rds")

dir.create("figures", showWarnings = FALSE)

# ── A/B: msVIPER barplots ─────────────────────────────────────

msviper_to_df <- function(mrs_obj, top_n = 25) {
  nes  <- mrs_obj$es$nes
  pval <- mrs_obj$es$p.value
  size <- mrs_obj$es$size
  df_out <- data.frame(
    regulator = names(nes),
    NES       = nes,
    pval      = pval,
    FDR       = p.adjust(pval, method = "BH"),
    size      = size[names(nes)],
    row.names = NULL
  )
  head(df_out[order(-abs(df_out$NES)), ], top_n)
}

plot_msviper_bar <- function(df, title, subtitle) {
  df$regulator <- factor(df$regulator, levels = df$regulator[order(df$NES)])
  df$sig <- ifelse(df$FDR < 0.001, "***",
             ifelse(df$FDR < 0.01,  "**",
             ifelse(df$FDR < 0.05,  "*", "")))
  df$label_x <- ifelse(df$NES >= 0,
                        df$NES + max(abs(df$NES)) * 0.03,
                       -max(abs(df$NES)) * 0.03 + df$NES)
  ggplot(df, aes(x = NES, y = regulator, fill = NES > 0)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(aes(x = label_x, label = sig, hjust = ifelse(NES >= 0, 0, 1)),
              size = 3.2, vjust = 0.75, color = "grey20") +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey30") +
    scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9")) +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    labs(title = title, subtitle = subtitle,
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 9, color = "grey45"),
      axis.text.y        = element_text(size = 9)
    )
}

df_tcga  <- msviper_to_df(tcga_mrs,  top_n = 25)
df_mtbrc <- msviper_to_df(mtbrc_mrs, top_n = 25)

panel_A <- plot_msviper_bar(
  df_tcga,
  title    = "Panel A — TCGA Basal BRCA",
  subtitle = "msVIPER · Kbhb signature (1,322 genes) · n = 195 tumor vs 113 normal"
)
panel_B <- plot_msviper_bar(
  df_mtbrc,
  title    = "Panel B — METABRIC Basal BRCA",
  subtitle = "msVIPER · Kbhb signature (1,213 genes) · n = 209 tumor vs 148 normal"
)

# ── C/D: Shadow network ───────────────────────────────────────

# Cross-platform TMR classification from Stouffer integration
# (produced by compare_kbhb_mrs.R; falls back to empty if unavailable)
if (file.exists("data/kbhb_mrs_comparison.rds")) {
  df_cmp     <- readRDS("data/kbhb_mrs_comparison.rds")
  sig_tfs_sh <- df_cmp$regulator[df_cmp$FDR_meta < 0.05]
  tmr_rep_sh <- sig_tfs_sh[df_cmp$NES_meta[match(sig_tfs_sh, df_cmp$regulator)] < 0]
} else {
  warning("kbhb_mrs_comparison.rds not found — TMR node coloring will be skipped")
  sig_tfs_sh <- character(0)
  tmr_rep_sh <- character(0)
}

# Extract shadow edges from msviper $shadow slot.
# Arrow direction: from = master (explains), to = shadow (explained).
# Handles Format A (named list) and Format B (matrix / data.frame).
.extract_shadow_edges <- function(mrs, cohort_label) {
  sh <- mrs$shadow
  if (is.null(sh))
    return(data.frame(from = character(0), to = character(0),
                      cohort = character(0), stringsAsFactors = FALSE))
  if (is.list(sh) && !is.data.frame(sh)) {
    pairs <- do.call(rbind, lapply(names(sh), function(tgt) {
      m <- as.character(sh[[tgt]])
      if (length(m) == 0) return(NULL)
      data.frame(from = m, to = tgt, cohort = cohort_label,
                 stringsAsFactors = FALSE)
    }))
    return(if (is.null(pairs))
      data.frame(from = character(0), to = character(0),
                 cohort = character(0), stringsAsFactors = FALSE)
      else pairs)
  }
  df_sh <- as.data.frame(sh, stringsAsFactors = FALSE)
  if (nrow(df_sh) == 0)
    return(data.frame(from = character(0), to = character(0),
                      cohort = character(0), stringsAsFactors = FALSE))
  cn <- tolower(colnames(df_sh))
  to_col   <- if ("shadow" %in% cn) which(cn == "shadow")[1] else 2L
  from_col <- if (any(cn %in% c("regulator", "master", "targets")))
                which(cn %in% c("regulator", "master", "targets"))[1] else 1L
  data.frame(
    from   = as.character(df_sh[[from_col]]),
    to     = as.character(df_sh[[to_col]]),
    cohort = cohort_label,
    stringsAsFactors = FALSE
  )
}

.classify_shadow_nodes <- function(edges, nes_vec, cohort_label,
                                   tmr_set, repressed_set) {
  if (is.null(edges) || nrow(edges) == 0) return(NULL)
  all_nodes <- unique(c(edges$from, edges$to))
  nes_vals  <- nes_vec[all_nodes]
  type <- dplyr::case_when(
    all_nodes %in% tmr_set                                      ~ "TMR",
    all_nodes %in% edges$from & !(all_nodes %in% edges$to)      ~ "Regulator",
    TRUE                                                         ~ "Shadow target"
  )
  fill <- dplyr::case_when(
    type == "TMR" & all_nodes %in% repressed_set ~ "#2980B9",
    type == "TMR"                                 ~ "#C0392B",
    type == "Regulator"                           ~ "#E67E22",
    TRUE                                          ~ "grey65"
  )
  data.frame(
    name      = all_nodes,
    cohort    = cohort_label,
    nes       = nes_vals,
    type      = type,
    fill      = fill,
    nes_label = ifelse(!is.na(nes_vals), sprintf("NES\n%+.2f", nes_vals), ""),
    stringsAsFactors = FALSE
  )
}

.make_shadow_plot <- function(edges, nodes, title_label) {
  if (is.null(nodes) || nrow(nodes) == 0 || nrow(edges) == 0)
    return(ggplot() +
             labs(title = title_label, subtitle = "No shadow pairs") +
             theme_bw(base_size = 11))
  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE) %>%
    activate(nodes) %>%
    mutate(fill = fill, nes_label = nes_label)
  ggraph(g, layout = "sugiyama") +
    geom_edge_arc(
      aes(),
      arrow     = arrow(length = unit(3, "mm"), type = "closed"),
      end_cap   = circle(7, "mm"),
      start_cap = circle(7, "mm"),
      color     = "grey40",
      linewidth = 0.7,
      strength  = 0.25
    ) +
    geom_node_point(aes(fill = I(fill)),
                    shape = 21, size = 14, color = "white", stroke = 1.2) +
    geom_node_text(aes(label = name),
                   fontface = "bold", size = 3.2, color = "white") +
    geom_node_text(aes(label = nes_label, y = y - 0.35),
                   size = 2.4, color = "grey20", vjust = 0) +
    labs(title    = title_label,
         subtitle = "Arrow: A → B = A partially explains the inferred activity of B") +
    theme_void(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.margin   = margin(15, 15, 15, 15)
    )
}

edges_tcga_sh  <- .extract_shadow_edges(tcga_mrs,  "TCGA")
edges_mtbrc_sh <- .extract_shadow_edges(mtbrc_mrs, "METABRIC")

nodes_tcga_sh  <- .classify_shadow_nodes(edges_tcga_sh,  tcga_mrs$es$nes,
                                          "TCGA",     sig_tfs_sh, tmr_rep_sh)
nodes_mtbrc_sh <- .classify_shadow_nodes(edges_mtbrc_sh, mtbrc_mrs$es$nes,
                                          "METABRIC", sig_tfs_sh, tmr_rep_sh)

panel_C <- .make_shadow_plot(edges_tcga_sh,  nodes_tcga_sh,
                              "Panel C — TCGA Basal: Shadow pairs")
panel_D <- .make_shadow_plot(edges_mtbrc_sh, nodes_mtbrc_sh,
                              "Panel D — METABRIC Basal: Shadow pairs")

# ── Legend ───────────────────────────────────────────────────

legend_df <- data.frame(
  x     = c(1,   1,   5,   5),
  y     = c(2,   1,   2,   1),
  fill  = c("#C0392B", "#2980B9", "#E67E22", "grey65"),
  label = c("TMR (activated, meta-FDR < 0.05)",
            "TMR (repressed, meta-FDR < 0.05)",
            "Cohort-significant regulator",
            "Shadow target (explained node)"),
  stringsAsFactors = FALSE
)

p_legend <- ggplot(legend_df, aes(x, y)) +
  geom_point(aes(fill = I(fill)), shape = 21, size = 5,
             color = "white", stroke = 1) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.18, size = 3) +
  xlim(0.7, 9.5) + ylim(0.5, 2.8) +
  theme_void() +
  labs(title = "Node type") +
  theme(
    plot.title  = element_text(face = "bold", size = 9, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

# ── Assemble and save ────────────────────────────────────────

fig_supp_combined <- (panel_A | panel_B) / (panel_C | panel_D) / p_legend +
  plot_layout(heights = c(1, 1, 0.15)) +
  plot_annotation(
    title   = "Supplementary Figure S1 — msVIPER results and regulatory hierarchy per cohort",
    caption = paste0(
      "Top 25 MRs by |NES| · FDR per cohort (BH) · * <0.05  ** <0.01  *** <0.001\n",
      "Shadow pairs from msviper $shadow · NES from per-cohort msVIPER analysis"
    ),
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.caption = element_text(size = 8, color = "grey50", hjust = 0)
    )
  )

ggsave("figures/fig_supp_combined.pdf", fig_supp_combined, width = 14, height = 16)
ggsave("figures/fig_supp_combined.png", fig_supp_combined,
       width = 14, height = 16, dpi = 300)
message("figures/fig_supp_combined.pdf/.png saved")

message("\n====  MRA KBHB COMPLETO  ====")
message("Objetos msVIPER: data/tcga_basal_mrs.rds · data/mtbrc_basal_mrs.rds")
message("Figura suplementaria: figures/fig_supp_combined.pdf/.png")
