library(TCGAbiolinks)
library(SummarizedExperiment)

# Set your project directory:
# setwd("/path/to/Kbhb")
if (!interactive()) setwd(normalizePath("."))

message(format(Sys.time()), " — Iniciando descarga TCGA-BRCA RNA-seq")

# ============================================================
# Query: todas las muestras TCGA-BRCA, RNA-seq STAR Counts
# ============================================================
query <- GDCquery(
  project          = "TCGA-BRCA",
  data.category    = "Transcriptome Profiling",
  data.type        = "Gene Expression Quantification",
  workflow.type    = "STAR - Counts"
)

message(format(Sys.time()), " — Query OK: ", nrow(getResults(query)), " archivos")

# ============================================================
# Descarga (chunks de 100 para no saturar el GDC)
# ============================================================
GDCdownload(query, files.per.chunk = 100)

message(format(Sys.time()), " — Descarga completa. Preparando SummarizedExperiment...")

# ============================================================
# Preparar objeto de expresión
# ============================================================
brca_se <- GDCprepare(query)

message(format(Sys.time()), " — SE creado: ", nrow(brca_se), " genes x ", ncol(brca_se), " muestras")

# ============================================================
# Guardar
# ============================================================
saveRDS(brca_se, "tcga_brca_rnaseq_se.rds")

message(format(Sys.time()), " — Guardado en tcga_brca_rnaseq_se.rds")
message("LISTO")
