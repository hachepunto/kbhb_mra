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

message("\n====  MRA KBHB COMPLETO  ====")
