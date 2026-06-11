
# tmr_motifs.R
# Motif-Integrated Network of TMRs for paper:
# "Two Cohorts, One Network: Consensus Master Regulators 
# Orchestrating Papillary Thyroid Carcinoma"
# Hugo Tovar, National Institute of Genomic Medicine, Mexico 
# hatovar@inmegen.gob.mx

#############################################
# 1) Load required libraries and set up folders
#############################################

# Motif databases and tools
suppressPackageStartupMessages({
  library(JASPAR2022)
  library(TFBSTools)
  
  # Genomics utilities
  library(biomaRt)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  
  # Data wrangling & I/O (load only what's needed)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(vroom)
  library(janitor)
  
  # Visualization
  library(circlize)
  library(RColorBrewer)
  library(metap)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(grid)
})

# Create output directories (idempotent)
outputsFolder <- file.path(getwd(), "TF_motives_results/")
plotsFolder   <- file.path(outputsFolder, "plots/")
dir.create(outputsFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(plotsFolder,   showWarnings = FALSE, recursive = TRUE)

# Load helper functions
source("helpers.R")

# utilities
`%||%` <- function(a, b) if (is.null(a)) b else a

#############################################
# 2) Obtaining motifs
#############################################

# Load the list of significant TMRs (adjusted p < 0.05)
# Expect a TSV with at least columns: 'tf' and 'meta_padj'
meta_tmrs <- vroom::vroom("meta_results/meta_mrs_results.tsv", .name_repair = janitor::make_clean_names) 

sign_tmrs <- meta_tmrs %>%
 filter(meta_padj < 0.05) %>% 
 pull(tf)
 
# Helper: fetch the latest PFMatrix version for a given TF name from JASPAR
get_pfm_for_tf <- function(tf_name) {
  ms <- TFBSTools::getMatrixSet(
    JASPAR2022,
    opts = list(all_versions = TRUE, name = tf_name)
  )
  if (length(ms) == 0) {
    warning("No motif found in JASPAR for: ", tf_name)
    return(NULL)
  }
  # Choose the latest version by ID suffix (e.g., MA0593.1 -> 1)
  ids <- vapply(ms, TFBSTools::ID, character(1))
  vers <- suppressWarnings(as.numeric(sub(".*\\.(\\d+)$", "\\1", ids)))
  ms[[which.max(vers)]]
}

# Build a PFM list for TMRs present in JASPAR (by TF name)
pfm_list <- purrr::map(sign_tmrs, get_pfm_for_tf) %>%
  rlang::set_names(sign_tmrs) %>%
  purrr::compact()

found_tfs   <- names(pfm_list)
missing_tfs <- setdiff(sign_tmrs, found_tfs)



# Search for exact match in tf_name (case insensitive)
find_rows_exact <- function(tbl, tf_symbol) {
  if (is.null(tbl)) return(tibble())
  tbl %>%
    filter(str_to_upper(tf_name) == str_to_upper(tf_symbol))
}

# Fallback "contains" (useful for micro-variants; avoids very loose matches)
# Uses word borders and allows optional -/_ (ETV-5 vs ETV5).
find_rows_fuzzy <- function(tbl, tf_symbol) {
  if (is.null(tbl)) return(tibble())
  pat <- str_replace_all(tf_symbol, "[-_ ]", "[-_ ]?")
  pat <- paste0("\\b", pat, "\\b")
  tbl %>%
    filter(str_detect(tf_name, regex(pat, ignore_case = TRUE)))
}

# Motifs priority
motif_type_priority <- c("ChIP-seq","ChIP-seq/ChIP-exo","SELEX","PBM","B1H","EMSA","DNase","PWM","Other")
prio <- function(x) match(x, motif_type_priority, nomatch = length(motif_type_priority)+1L)

choose_best_motif <- function(df_tfrows) {
  if (nrow(df_tfrows) == 0) return(NULL)
  df_tfrows %>%
    mutate(
      tf_status     = tf_status %||% NA,
      motif_type    = motif_type %||% NA,
      m_source_year = suppressWarnings(as.integer(m_source_year)),
      pr_tf_status  = case_when(tf_status == "D" ~ 0L, tf_status == "I" ~ 1L, TRUE ~ 2L),
      pr_type       = prio(motif_type),
      pr_year       = -m_source_year
    ) %>%
    arrange(pr_tf_status, pr_type, pr_year, motif_id) %>%
    slice(1)
}

# Given a TF, first search for exact in primary; if not, exact in all/plus; if not, fuzzy in all three
get_cisbp_motif_meta <- function(tf_symbol, cis_primary, cis_all = NULL, cis_plus = NULL) {
  cand <- bind_rows(
    find_rows_exact(cis_primary, tf_symbol),
    find_rows_exact(cis_all,     tf_symbol),
    find_rows_exact(cis_plus,    tf_symbol)
  ) %>% distinct()

  if (nrow(cand) == 0) {
    cand <- bind_rows(
      find_rows_fuzzy(cis_primary, tf_symbol),
      find_rows_fuzzy(cis_all,     tf_symbol),
      find_rows_fuzzy(cis_plus,    tf_symbol)
    ) %>% distinct()
  }

  if (nrow(cand) == 0) return(NULL)
  choose_best_motif(cand)
}

try_read_one <- function(tf) {
  meta <- get_cisbp_motif_meta(tf, cis_primary, cis_all, cis_plus)
  if (is.null(meta) || is.null(meta$motif_id) ||
      is.na(meta$motif_id) || meta$motif_id %in% c(".", "")) {
    return(tibble::tibble(tf=tf, motif_id=meta$motif_id %||% NA_character_, ok=FALSE, reason="bad_id"))
  }
  pfm <- try(read_cisbp_pfm(meta$motif_id, pwms_dir, verbose = TRUE), silent = TRUE)
  if (inherits(pfm, "try-error") || is.null(pfm)) {
    return(tibble::tibble(tf=tf, motif_id=meta$motif_id, ok=FALSE, reason="read_fail"))
  }
  tibble::tibble(tf=tf, motif_id=meta$motif_id, ok=TRUE, reason="ok")
}


# ---- CIS-BP configuration ----
unzip("extra_data/Homo_sapiens_2025_08_23_3_48_am.zip")
cisbp_dir <- "extra_data/Homo_sapiens_2025_08_23_3_48_am" 
pwms_dir  <- file.path(cisbp_dir, "pwms_all_motifs")

cis_primary <- vroom::vroom(paste0(cisbp_dir, "/TF_Information.txt"), delim="\t", col_types=cols(.default="c")) %>% janitor::clean_names()
cis_all     <- vroom::vroom(paste0(cisbp_dir, "/TF_Information_all_motifs.txt"),  delim="\t", col_types=cols(.default="c")) %>% janitor::clean_names()
cis_plus    <- vroom::vroom(paste0(cisbp_dir, "/TF_Information_all_motifs_plus.txt"), delim="\t", col_types=cols(.default="c")) %>% janitor::clean_names()

cisbp_pwm_list <- purrr::map(missing_tfs, function(tf) {
  meta <- get_cisbp_motif_meta(tf, cis_primary, cis_all, cis_plus)
  if (is.null(meta) || is.null(meta$motif_id) || meta$motif_id %in% c(".", "")) {
    warning("CIS-BP has no motifs (neither exact nor fuzzy) to: ", tf)
    return(NULL)
  }
  pfm <- read_cisbp_pfm(meta$motif_id, pwms_dir, verbose = FALSE)
  if (is.null(pfm)) {
    warning("I couldn't read PFM for ", tf, " (motif_id=", meta$motif_id, ")")
    return(NULL)
  }
  pwm <- TFBSTools::toPWM(pfm, pseudocounts = 0.2)
  pwm@tags$source        <- "CIS-BP"
  pwm@tags$tf_status     <- meta$tf_status
  pwm@tags$motif_type    <- meta$motif_type
  pwm@tags$m_source_year <- meta$m_source_year
  pwm
}) %>%
  rlang::set_names(missing_tfs) %>%
  purrr::compact()

audit <- purrr::map_dfr(missing_tfs, try_read_one)
audit %>% count(ok, reason)
audit %>% filter(!ok)


data.frame(
  tf = missing_tfs,
  ok = missing_tfs %in% names(cisbp_pwm_list)
) |> dplyr::count(ok)


# 1) Loading of the manualy recovered
recovered <- vroom::vroom("extra_data/tfs_found_jaspar.txt", show_col_types = FALSE)

manual_pfms <- recovered %>%
  filter(!is.na(id_jaspar)) %>%
  distinct(tf, id_jaspar) %>%
  mutate(
    pfm = map(id_jaspar, ~ {
      ms <- TFBSTools::getMatrixSet(JASPAR2022, opts = list(collection="CORE", tax_group="vertebrates", ID = .x))
      ms[[.x]]
    })
  ) %>%
  select(tf, pfm) %>%
  deframe()

# 2) Ensures source labels for everyone
tag_source <- function(x, src) {
  if (is.null(x@tags$source)) x@tags$source <- src
  x
}
pfm_list         <- imap(pfm_list,         ~ tag_source(.x, "JASPAR"))                # auto JASPAR
manual_pfms      <- imap(manual_pfms,      ~ tag_source(.x, "JASPAR(manual)"))        # manual JASPAR
cisbp_pwm_list   <- imap(cisbp_pwm_list,   ~ { .x@tags$source <- "CIS-BP"; .x })      # ya son PWMs; OK

# 3) Convert everything to PWM (if anything is still PFM)
to_pwm_if_needed <- function(x) {
  if (inherits(x, "PWMatrix")) return(x)
  if (inherits(x, "PFMatrix")) return(TFBSTools::toPWM(x, pseudocounts = 0.2))
  stop("Object of unrecognized motif:", class(x)[1])
}
pwm_jaspar_auto   <- imap(pfm_list,    ~ to_pwm_if_needed(.x))
pwm_jaspar_manual <- imap(manual_pfms, ~ to_pwm_if_needed(.x))

# 4) Priority unification: JASPAR manual > JASPAR auto > CIS-BP
#    - To avoid overwriting, we build in reverse order and let the last one win.
pwm_list_extended <- list()
# base: CIS-BP (lowest priority)
pwm_list_extended <- modifyList(pwm_list_extended, cisbp_pwm_list, keep.null = FALSE)
# then JASPAR auto (overwrites if exists in CIS-BP)
pwm_list_extended <- modifyList(pwm_list_extended, pwm_jaspar_auto, keep.null = FALSE)
# at the end JASPAR manual (highest priority)
pwm_list_extended <- modifyList(pwm_list_extended, pwm_jaspar_manual, keep.null = FALSE)

# 5) Coverage Summary
found_tfs_all   <- sort(intersect(sign_tmrs, names(pwm_list_extended)))
missing_after   <- setdiff(sign_tmrs, names(pwm_list_extended))

coverage_tbl <- tibble(
  total_sign_tmrs = length(sign_tmrs),
  n_with_pwm      = length(found_tfs_all),
  n_missing       = length(missing_after)
)

message("Founds: ", coverage_tbl$n_with_pwm, "/", coverage_tbl$total_sign_tmrs,
        " (missing ", coverage_tbl$n_missing, ")")

# 6) effective source used
motif_source <- function(pwm) pwm@tags$source %||% NA_character_
source_tbl <- tibble(
  tf = found_tfs_all,
  source = vapply(pwm_list_extended[found_tfs_all], motif_source, character(1))
) %>%
  mutate(source = ifelse(is.na(source), "UNKNOWN", source)) %>%
  count(source, name = "n")

# 7) Which of the original 'missing_tfs' did CIS-BP resolve?
resolved_by_cisbp <- intersect(names(cisbp_pwm_list), setdiff(sign_tmrs, names(pfm_list)))
resolved_tbl <- tibble(
  resolved_by_cisbp = length(resolved_by_cisbp),
  still_missing     = length(missing_after)
)

# 8) Metadata table for reproducibility
metadata_table <- purrr::imap_dfr(pwm_list_extended, function(pwm, tf) {
  get_tag_chr <- function(x) {
    if (is.null(x)) return(NA_character_)
    as.character(x)
  }
  tibble::tibble(
    tf         = tf,
    source     = get_tag_chr(pwm@tags$source),
    id         = get_tag_chr(tryCatch(TFBSTools::ID(pwm),   error = function(e) NA_character_)),
    name       = get_tag_chr(tryCatch(TFBSTools::name(pwm), error = function(e) NA_character_)),
    motif_type = get_tag_chr(pwm@tags$motif_type),
    tf_status  = get_tag_chr(pwm@tags$tf_status),
    year       = get_tag_chr(pwm@tags$m_source_year),
    file       = get_tag_chr(pwm@tags$file),
    data_type  = get_tag_chr(pwm@tags$data_type)
  )
}) %>% dplyr::arrange(tf)

# 9) Persist results
dir.create(file.path(outputsFolder, "motifs"), showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(metadata_table, file.path(outputsFolder, "motifs", "motif_metadata_combined.tsv"))

# Save the combo list to RDS
saveRDS(pwm_list_extended, file.path(outputsFolder, "motifs", "pwm_list_extended_list.rds"))


#############################################
# 3) Obtaining promoter regions
#############################################

# Query Ensembl for TSS data
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

tss_data <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "chromosome_name",
    "strand",
    "transcription_start_site"
  ),
  filters = "hgnc_symbol",
  values = sign_tmrs,
  mart = ensembl
)

# Report genes not found
not_found <- setdiff(sign_tmrs, unique(tss_data$hgnc_symbol))
if(length(not_found)) {
  warning("These genes were not found in Ensembl:\n", paste(not_found, collapse=", "))
}

# Select one representative TSS per gene
# Rule: for + strand take the minimal TSS (most upstream); for - strand take the maximal TSS.
# If a gene appears on both strands (rare), the first strand encountered is used and a warning is issued.
tss_upstream <- tss_data %>%
  group_by(hgnc_symbol) %>%
  dplyr::slice(if (unique(strand) == 1) which.min(transcription_start_site)
        else                 which.max(transcription_start_site)) %>%
  ungroup() %>%
  mutate(
    tss = transcription_start_site,
    tss_promoter_start = if_else(strand == 1, tss - 2000, tss - 200),
    tss_promoter_end   = if_else(strand == 1, tss + 200, tss + 2000)
  )

# Build GRanges for promoters (UCSC hg38 uses "chr" prefix)
promoter_gr <- GRanges(
  seqnames = paste0("chr", tss_upstream$chromosome_name),
  ranges   = IRanges(start = pmax(tss_upstream$tss_promoter_start, 1),
                     end   = tss_upstream$tss_promoter_end),
  strand   = ifelse(tss_upstream$strand == 1, "+", "-"),
  gene     = tss_upstream$hgnc_symbol
)

# Extract promoter sequences from BSgenome (hg38)
promoter_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, promoter_gr)

# Name sequences as TF_chr:start-end(strand)
names(promoter_seqs) <- paste0(
  tss_upstream$hgnc_symbol, "_",
  seqnames(promoter_gr), ":",
  start(promoter_gr), "-", end(promoter_gr),
  strand(promoter_gr)
)

# Save promoters as FASTA for external use
# Biostrings::writeXStringSet(promoter_seqs, filepath = file.path(outputsFolder, "promoters_50tmrs_hg38.fa"))

#############################################
# 4) Motif hit search
#############################################

# Build all TF–target combinations
combos <- expand_grid(
  TF_from = names(pwm_list_extended),
  TF_to   = names(promoter_seqs)
)

# Search helper: run PWM scan for a (TF_from, TF_target) pair
search_pair <- function(TF_from, TF_to) {
  pwm <- pwm_list_extended[[TF_from]]
  seq <- promoter_seqs[[TF_to]]
  
  hits <- searchSeq(pwm, seq, min.score = "85%", strand = "*")
  if (length(hits) == 0) return(NULL)
  
  df <- as.data.frame(hits)
  df$TF_from <- TF_from
  df$TF_to   <- TF_to
  df
}

# 3. Execute scans over all pairs and bind results
binding_hits <- combos %>%
  pmap_dfr(~ search_pair(..1, ..2))

tss_up2 <- tss_upstream
colnames(tss_up2)[ colnames(tss_up2) == "strand" ] <- "strand_tss"

binding_hits2 <- binding_hits %>%
  mutate(
    TF_target = str_remove(TF_to, "_chr.*")
  )

# annotate genomic coordinates
binding_hits_annot <- binding_hits2 %>%
  left_join(
    tss_up2 %>% dplyr::select(hgnc_symbol, chromosome_name, strand_tss,
                       tss_promoter_start, tss_promoter_end),
    by = c("TF_target" = "hgnc_symbol")
  ) %>%
  mutate(
    genomic_start = if_else(strand_tss == 1L,
                            tss_promoter_start + start,
                            tss_promoter_end   - end),
    genomic_end   = if_else(strand_tss == 1L,
                            tss_promoter_start + end,
                            tss_promoter_end   - start),
    chromosome    = paste0("chr", chromosome_name)
  ) %>%
  dplyr::select(TF_from, TF_target, chromosome, genomic_start, genomic_end, everything())

# Save tabular outputs
readr::write_tsv(binding_hits_annot, paste0(outputsFolder,"binding_all50x37_only_motifs.tsv"))

#############################################
# 5) Summary table
#############################################

# 50 TFs list
tf_list <- tibble(tf = unique(c(binding_hits_annot$TF_from, binding_hits_annot$TF_target)))

# For each TF as emitter (TF_from): total hits and number of distinct targets
tf_from_summary <- binding_hits_annot %>%
  group_by(tf = TF_from) %>%
  summarise(
    n_hits    = n(),
    n_targets = n_distinct(TF_target),
    .groups   = "drop"
  )

# For each TF as target (TF_target): total incoming hits
tf_to_summary <- binding_hits_annot %>%
  group_by(tf = TF_target) %>%
  summarise(
    n_strikes = n(),
    .groups   = "drop"
  )

# "Beaters": number of distinct emitters that target each TF
beaters_summary <- binding_hits_annot %>%
  group_by(tf = TF_target) %>%
  summarise(
    n_beaters = n_distinct(TF_from),
    .groups   = "drop"
  )

tf_summary <- tf_list %>%
  left_join(tf_from_summary,   by = "tf") %>%
  left_join(tf_to_summary,     by = "tf") %>%
  left_join(beaters_summary,   by = "tf") %>%
  replace_na(list(
    n_hits    = 0,
    n_targets = 0,
    n_strikes = 0,
    n_beaters = 0
  ))

# Assemble final summary table
tf_summary <- tf_summary %>%
  left_join(
    tss_upstream %>%
      dplyr::select(
        hgnc_symbol,
        chromosome      = chromosome_name,
        genomic_start  = tss_promoter_start,
        genomic_end    = tss_promoter_end,
        strand          = strand
      ),
    by = c("tf" = "hgnc_symbol")
  ) %>%
  dplyr::select(
    tf,
    chromosome,
    genomic_start,
    genomic_end,
    strand,
    n_hits,
    n_targets,
    n_strikes,
    n_beaters
  ) %>% 
  arrange(desc(n_targets), desc(n_hits))

# Save and preview
print(tf_summary, n = 50)
readr::write_tsv(tf_summary, paste0(outputsFolder, "interactions_summary_only_motifs.tsv"))

#############################################
# 6) Only motif hits the circos plot
#############################################

# Build interaction matrix in the exact TF order used in the summary
tf_order <- tf_summary$tf

interaction_matrix <- binding_hits_annot %>%
  count(TF_from, TF_target, name = "n") %>%
  complete(TF_from = tf_order, TF_target = tf_order, fill = list(n = 0)) %>%
  pivot_wider(names_from = TF_target, values_from = n) %>%
  tibble::column_to_rownames("TF_from") %>%
  as.matrix()

# Ensure rows/cols are in the same order
interaction_matrix <- interaction_matrix[tf_order, tf_order, drop = FALSE]
stopifnot(identical(rownames(interaction_matrix), tf_order))
stopifnot(identical(colnames(interaction_matrix), tf_order))

## 2) Colors and sector labels

base_pal   <- colorRampPalette(brewer.pal(11, "Spectral"))(length(tf_order))
set.seed(42)  # reproducible shuffle
tf_colors  <- setNames(sample(base_pal), tf_order)
names(tf_colors) <- tf_order

sector_labels <- paste0(tf_summary$tf, " (", tf_summary$n_targets, "/", tf_summary$n_beaters, ")")
names(sector_labels) <- tf_summary$tf   

## Per-link color: inherit emitter color; transparent if no link
color_mat <- matrix(tf_colors[rownames(interaction_matrix)],
                    nrow = nrow(interaction_matrix), ncol = ncol(interaction_matrix),
                    dimnames = dimnames(interaction_matrix))
color_mat[interaction_matrix == 0] <- "transparent"

## Drawing helper (used for on-screen and file outputs)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))

chordDiagram(
  x               = interaction_matrix,
  order           = tf_order,
  grid.col        = tf_colors,
  col             = color_mat,
  link.arr.col    = color_mat,
  link.border     = NA,
  directional     = TRUE,
  direction.type  = "arrows",
  link.arr.length = 0.05,
  link.arr.width  = 0.05,
  transparency    = 0.25,
  annotationTrack = "grid"
)

# Safety check
stopifnot(identical(get.all.sector.index(), tf_order))

## Add outer labels (outside the ring)
circos.trackPlotRegion(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x       = mean(xlim),
      y       = ylim[2] + mm_y(2.5),
      labels  = sector_labels[sec],
      facing  = "clockwise",
      niceFacing = FALSE,
      adj     = c(0, 0.5),
      cex     = 0.75
    )
  }
)
title("Circos plot of directional interactions among 50 thyroid Transcriptional Master Regulators (Outgoing vs. Incoming TFs)")


#############################################
# 7) Expresión analysis
#############################################

deg_tcga <- vroom::vroom("results/DEGs_limma_TCGA.tsv", .name_repair = janitor::make_clean_names)
deg_geo <- vroom::vroom("results/DEGs_limma_GEO.tsv", .name_repair = janitor::make_clean_names)

tmr_df <- deg_tcga %>%
  dplyr::select(gene_name, t_tcga = t, log_fc_tcga = log_fc, p_tcga = p_value) %>%
  inner_join(
    deg_geo %>% select(gene_name, t_geo = t, log_fc_geo = log_fc, p_geo = p_value),
    by = "gene_name"
  ) %>%
  filter(gene_name %in% tf_order) %>%
  mutate(
    direction_tcga = sign(log_fc_tcga),
    direction_geo  = sign(log_fc_geo),
    same_direction = direction_tcga == direction_geo
  )

table(tmr_df$same_direction)
# TRUE
#   50

tmr_meta_fisher <- tmr_df %>%
  rowwise() %>%
  mutate(p_meta_fisher = sumlog(c(p_tcga, p_geo))$p) %>%
  ungroup() %>%
  mutate(p_meta_fisher_adj = p.adjust(p_meta_fisher, "BH"))

tmr_meta <- tmr_df %>%
  mutate(
    se_tcga = abs(log_fc_tcga) / pmax(abs(t_tcga), .Machine$double.eps),
    se_geo  = abs(log_fc_geo)  / pmax(abs(t_geo),  .Machine$double.eps),
    w_tcga  = 1 / (se_tcga^2),
    w_geo   = 1 / (se_geo^2),

    lfc_meta = (w_tcga*log_fc_tcga + w_geo*log_fc_geo) / (w_tcga + w_geo),
    se_meta  = sqrt(1/(w_tcga + w_geo)),
    z_meta_iv = lfc_meta / se_meta,
    p_meta_iv = 2*pnorm(abs(z_meta_iv), lower.tail = FALSE),
    p_meta_iv_adj = p.adjust(p_meta_iv, "BH"),
    dir = if_else(lfc_meta > 0, "up", "down")
  )


meta_sig <- tmr_meta %>%
  transmute(gene_name, dir, sig = p_meta_iv_adj < 0.05)

tf_mode <- binding_hits_annot %>%
  semi_join(meta_sig, by = c("TF_from"   = "gene_name")) %>%
  semi_join(meta_sig, by = c("TF_target" = "gene_name")) %>%
  left_join(rename(meta_sig, dir_from = dir, sig_from = sig),   by = c("TF_from"   = "gene_name")) %>%
  left_join(rename(meta_sig, dir_to   = dir, sig_to   = sig),   by = c("TF_target" = "gene_name")) %>%
  mutate(
    mode = case_when(
      sig_from & sig_to & dir_from == dir_to ~ "possible_activation",
      sig_from & sig_to & dir_from != dir_to ~ "possible_repression",
      TRUE                                   ~ "uncertain"
    )
  )



# ========= 1) META direction by TF =========
meta_dir <- tmr_meta %>%
  transmute(tf = gene_name, dir, sig = p_meta_iv_adj < 0.05)

stopifnot(all(tf_order %in% meta_dir$tf))

# ========= 2) Classify links: activation/repression =========
edge_df <- binding_hits_annot %>%
  count(TF_from, TF_target, name = "n_hits") %>%
  left_join(rename(meta_dir, dir_from = dir, sig_from = sig), by = c("TF_from" = "tf")) %>%
  left_join(rename(meta_dir, dir_to   = dir, sig_to   = sig), by = c("TF_target" = "tf")) %>%
  mutate(
    mode = case_when(
      sig_from & sig_to & dir_from == dir_to ~ "possible_activation",
      sig_from & sig_to & dir_from != dir_to ~ "possible_repression",
      TRUE                                   ~ "uncertain"
    )
  )

#############################################
# 8) ARACNe network support
#############################################

aracne_tcga <- read_tsv("tcga_tumor_network.txt") %>% 
      rename(tf = Regulator, target = Target, mi = MI)
aracne_geo <- read_tsv("geo_tumor_network.txt") %>% 
      rename(tf = Regulator, target = Target, mi = MI)


aracne_union <- bind_rows(
  aracne_tcga %>% mutate(ds = "tcga"),
  aracne_geo  %>% mutate(ds = "geo")
) %>%
  count(tf, target, name = "aracne_support")

edge_df2 <- edge_df %>%
  left_join(aracne_union, by = c("TF_from" = "tf", "TF_target" = "target")) %>%
  mutate(
    aracne_support = replace_na(aracne_support, 0L)
  )

# Filter by minimum support
min_support <- 1
edge_plot <- edge_df2 %>%
  filter(aracne_support >= min_support)

# ========= 4) Matrix for chord and colors =========

edge_clean <- edge_plot %>%
  transmute(TF_from, TF_target, weight = n_hits * aracne_support) %>%
  filter(!is.na(TF_from), !is.na(TF_target))

# “observed” matrix (only for present pairs)
mat_obs <- xtabs(weight ~ TF_from + TF_target,
                data = edge_clean) %>%
            as.matrix()

full_mat <- matrix(0, nrow = length(tf_order), ncol = length(tf_order),
                   dimnames = list(tf_order, tf_order))
common_r <- intersect(rownames(mat_obs), tf_order)
common_c <- intersect(colnames(mat_obs), tf_order)
full_mat[common_r, common_c] <- mat_obs[common_r, common_c]

interaction_mat <- full_mat

#############################################
# 9) Circus plot ARACNe support + Expresión analysis
#############################################

base_pal <- colorRampPalette(brewer.pal(11, "Spectral"))(length(tf_order))
set.seed(42)
tf_colors <- setNames(sample(base_pal), tf_order)

# Link color matrix inheriting from the emitter
color_mat <- matrix(tf_colors[rownames(interaction_mat)],
                    nrow = nrow(interaction_mat), ncol = ncol(interaction_mat),
                    dimnames = dimnames(interaction_mat))
color_mat[interaction_mat == 0] <- "transparent"

grid_cols <- tf_colors


################################################################################
# ==== FIGURE 2 ====
################################################################################

# Save to PDF and PNG
pdf_file <- paste0(plotsFolder, "TF_binding_circos.pdf")
png_file <- paste0(plotsFolder, "TF_binding_circos.png")

# Dimensions in inches or pixels
pdf_wh  <- 16
png_wh  <- 12400
png_res    <- 600 
dir_cols <- c(up = "#E41A1C", down = "#377EB8")  # rojo=up, azul=down
names(dir_cols) <- names(dir_cols)

# 1) PDF save
pdf(pdf_file, width = pdf_wh+0.5, height = pdf_wh)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
chordDiagram(
  x               = interaction_mat,
  order           = tf_order,
  grid.col        = tf_colors,
  col             = color_mat,
  link.arr.col    = color_mat,
  link.border     = NA,
  directional     = TRUE,
  direction.type  = "arrows",
  link.arr.length = 0.05,
  link.arr.width  = 0.05,
  transparency    = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = uh(4, "mm")),
    list(track.height = uh(3, "mm"))
  ))
old_pad <- circos.par("cell.padding")
circos.par(cell.padding = c(0, 0, 0, 0)) 
circos.trackPlotRegion(
  track.index = 2,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    dir  <- meta_dir$dir[match(sec, meta_dir$tf)]
    col  <- ifelse(is.na(dir), "grey85", dir_cols[dir])
    xl   <- get.cell.meta.data("xlim")
    circos.rect(xl[1], 0, xl[2], 1, col = col, border = NA)
  })
circos.par(cell.padding = old_pad)
circos.trackPlotRegion(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xl   <- get.cell.meta.data("xlim")
    circos.text(
      x       = mean(xl),
      y       = 0.5,
      labels  = sector_labels[sec],
      facing  = "clockwise",
      niceFacing = FALSE,
      adj     = c(0, 0.5),
      cex     = 0.75
    )
  })
draw(
  lg_pack,
  x    = unit(0.06, "npc"), 
  y    = unit(0.06, "npc"), 
  just = c("left", "bottom"))
title("Circos plot of directional interactions among 50 thyroid Transcriptional Master Regulators (Outgoing vs. Incoming TFs)")
dev.off()

# 2) PNG save
png(png_file, width = png_wh, height = png_wh, res = png_res)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
chordDiagram(
  x               = interaction_mat,
  order           = tf_order,
  grid.col        = tf_colors,
  col             = color_mat,
  link.arr.col    = color_mat,
  link.border     = NA,
  directional     = TRUE,
  direction.type  = "arrows",
  link.arr.length = 0.05,
  link.arr.width  = 0.05,
  transparency    = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = uh(4, "mm")),
    list(track.height = uh(3, "mm"))
  ))
old_pad <- circos.par("cell.padding")
circos.par(cell.padding = c(0, 0, 0, 0))
circos.trackPlotRegion(
  track.index = 2,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    dir  <- meta_dir$dir[match(sec, meta_dir$tf)]
    col  <- ifelse(is.na(dir), "grey85", dir_cols[dir])
    xl   <- get.cell.meta.data("xlim")
    circos.rect(xl[1], 0, xl[2], 1, col = col, border = NA)
  })
circos.par(cell.padding = old_pad)
circos.trackPlotRegion(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xl   <- get.cell.meta.data("xlim")
    circos.text(
      x       = mean(xl),
      y       = 0.5, 
      labels  = sector_labels[sec], 
      facing  = "clockwise",
      niceFacing = FALSE,
      adj     = c(0, 0.5),
      cex     = 0.75
    )
  })
draw(
  lg_pack,
  x    = unit(0.08, "npc"),  
  y    = unit(0.08, "npc"),   
  just = c("left", "bottom"))
title("Circos plot of directional interactions among 50 thyroid Transcriptional Master Regulators (Outgoing vs. Incoming TFs)")
dev.off()


# === EDGES =======================================================
# Asegura que partimos sólo de columnas "limpias"
edge_base <- edge_plot %>%
  select(TF_from, TF_target, n_hits, aracne_support, mode) %>%
  distinct()

# Meta-info desde tmr_meta con nombres ya diferenciados
meta_from <- tmr_meta %>%
  transmute(tf = gene_name,
            lfc_from = lfc_meta,
            p_from   = p_meta_iv_adj,
            dir_from = dir)

meta_to <- tmr_meta %>%
  transmute(tf = gene_name,
            lfc_to = lfc_meta,
            p_to   = p_meta_iv_adj,
            dir_to = dir)

edges_cyto <- edge_base %>%
  left_join(meta_from, by = c("TF_from"   = "tf")) %>%
  left_join(meta_to,   by = c("TF_target" = "tf")) %>%
  mutate(weight = n_hits * aracne_support) %>%
  filter(weight > 0) %>%
  transmute(
    source      = TF_from,
    target      = TF_target,
    interaction = "regulates",
    directed    = TRUE,
    weight,                 # = n_hits * soporte ARACNe
    n_hits,
    aracne_support,         # 1..2
    mode,                   # possible_activation / possible_repression / uncertain
    dir_from, dir_to,       # up / down
    lfc_from, p_from,
    lfc_to,   p_to
  ) %>%
  arrange(desc(weight))

# === NODES =======================================================
# Usa los TF presentes en la red final (más compacto para Cytoscape)
tf_in_network <- sort(unique(c(edges_cyto$source, edges_cyto$target)))

# Grados coherentes con la red final (si ya tienes interaction_mat, úsalo; si no, calcúlalo aquí)
# Construimos una matriz  TF_from x TF_target desde edges_cyto (peso > 0)
mat_obs <- xtabs(weight ~ source + target, data = edges_cyto) %>% as.matrix()

# Out-/In-degree por existencia de arista (>0) y strengths por suma de pesos
node_degrees <- tibble(tf = tf_in_network) %>%
  mutate(
    n_targets_mi = purrr::map_int(tf, ~ if (.x %in% rownames(mat_obs)) sum(mat_obs[.x, ] > 0) else 0L),
    n_beaters_mi = purrr::map_int(tf, ~ if (.x %in% colnames(mat_obs)) sum(mat_obs[, .x] > 0) else 0L),
    out_strength = purrr::map_dbl(tf, ~ if (.x %in% rownames(mat_obs)) sum(mat_obs[.x, ], na.rm = TRUE) else 0),
    in_strength  = purrr::map_dbl(tf, ~ if (.x %in% colnames(mat_obs)) sum(mat_obs[, .x], na.rm = TRUE) else 0)
  )

# Posiciones promotoras (ajusta nombres si difieren)
pos_info <- tss_upstream %>%
  select(hgnc_symbol, chromosome_name, tss_promoter_start, tss_promoter_end, strand) %>%
  mutate(chromosome = paste0("chr", chromosome_name)) %>%
  transmute(
    tf = hgnc_symbol,
    chromosome,
    promoter_start = tss_promoter_start,
    promoter_end   = tss_promoter_end,
    strand
  ) %>%
  distinct(tf, .keep_all = TRUE)

nodes_cyto <- tibble(tf = tf_in_network) %>%
  left_join(tmr_meta %>% select(gene_name, lfc_meta, p_meta_iv_adj, dir),
            by = c("tf" = "gene_name")) %>%
  left_join(node_degrees, by = "tf") %>%
  left_join(pos_info,     by = "tf") %>%
  transmute(
    id = tf, label = tf,
    lfc_meta, p_meta_iv_adj, dir,
    n_targets_mi, n_beaters_mi,
    out_strength, in_strength,
    chromosome, promoter_start, promoter_end, strand
  )

# === EXPORT =====================================================
dir.create(paste0(outputsFolder, "cytoscape"), showWarnings = FALSE, recursive = TRUE)
write_tsv(edges_cyto, paste0(outputsFolder, "cytoscape/edges_aracne_expr.tsv"))
write_tsv(nodes_cyto, paste0(outputsFolder, "cytoscape/nodes_tf_meta_coords.tsv"))
