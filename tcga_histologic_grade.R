library(readxl)
library(dplyr)

if (!interactive()) setwd(normalizePath("."))
dir.create("data", showWarnings = FALSE)

# ============================================================
# Derive TCGA-BRCA histologic grade (Nottingham) from Thennavan et al. 2021
# Data S2 (Cell Genomics 1:100067, PMID 35465400), as requested by the IJMS
# editor in the second Minor Revision round. TCGA's central clinical data
# (PanCanAtlas paper_Tumor_Grade) does not report grade for this tumor type
# (see figure1_panel_c.R, section 1b) -- Thennavan's pathologist panel
# re-scored the three Nottingham components directly from histology slides,
# so grade can be reconstructed from there instead.
#
# Grade is NOT provided pre-computed in Data S2: the three Nottingham
# components (tubule formation, nuclear pleomorphism, mitotic count) are
# scored 1-3 each (Heng et al. 2017, PMID 27861902); grade is derived here
# by summing them and applying the standard Elston-Ellis cutoff.
# ============================================================

# ---- 1. Read Data S2 (sheet has a trailing space in its name -- use
#      positional index to avoid depending on that literal string) ----

thenn <- read_excel(
  "PMID_35465400/1-s2.0-S2666979X21000835-mmc3.xlsx",
  sheet = 1
)
stopifnot(nrow(thenn) == 1063, ncol(thenn) == 15)
stopifnot(!anyDuplicated(thenn[["Sample CLID"]]))   # one row per patient

# ---- 2. Parse the three Nottingham components to integer score (1-3) ----
# Cells are free text, e.g. "(score = 3) <10%" or the literal string "NA"
# (not a true missing value) when the component wasn't scored.

parse_score <- function(x) {
  x <- trimws(x)
  out <- rep(NA_integer_, length(x))
  hit <- regexpr("score = [0-9]", x)
  ok <- hit > 0
  out[ok] <- as.integer(sub(".*score = ([0-9]).*", "\\1", x[ok]))
  out
}

tubule_formation     <- parse_score(thenn[["Epithelial tubule formation"]])
nuclear_pleomorphism <- parse_score(thenn[["Nuclear pleomorphism"]])
mitotic_count        <- parse_score(thenn[["Mitosis"]])

cat("\n================ Thennavan Data S2 -- parsed component scores (1063 patients) ================\n")
cat("\n-- tubule_formation --\n");     print(table(tubule_formation, useNA = "always"))
cat("\n-- nuclear_pleomorphism --\n"); print(table(nuclear_pleomorphism, useNA = "always"))
cat("\n-- mitotic_count --\n");        print(table(mitotic_count, useNA = "always"))

# ---- 3. Nottingham score / grade. NA if any component is missing -- no
#      imputation. ----

nottingham_score <- tubule_formation + nuclear_pleomorphism + mitotic_count
grade <- cut(
  nottingham_score,
  breaks = c(-Inf, 5, 7, Inf),
  labels = c("1", "2", "3")
)
grade <- as.character(grade)

thenn_grade <- data.frame(
  patient_id           = thenn[["Sample CLID"]],
  tubule_formation     = tubule_formation,
  nuclear_pleomorphism = nuclear_pleomorphism,
  mitotic_count        = mitotic_count,
  nottingham_score     = nottingham_score,
  grade                = grade,
  stringsAsFactors     = FALSE
)

# ---- 4. Join to our 195 Basal-like TCGA samples, by 12-char patient
#      barcode (Thennavan is patient-level; our sample_id is 16-char,
#      aliquot-level). NOT joined on Thennavan's own PAM50/CLOW subtype
#      label -- see tareas_minor_revision2.md for why (179 "Basal" there
#      != our 195, different subtyping criterion and cohort coverage). ----

tcga_clin <- read.delim("data/tcga_basal_clinical_extra.tsv", stringsAsFactors = FALSE)
stopifnot(nrow(tcga_clin) == 195)

patient_id <- substr(tcga_clin$sample_id, 1, 12)
n_dup_patients <- sum(duplicated(patient_id))
if (n_dup_patients > 0) {
  cat("\nNOTE:", n_dup_patients, "of the 195 samples share a patient barcode with",
      "another sample in the cohort (TCGA replicate aliquots/portions, e.g. 01A/01B",
      "of the same tumor) -- both receive the same Thennavan grade by construction,",
      "the join is many(samples)-to-one(patient), not 1:1:\n")
  dup_ids <- patient_id[duplicated(patient_id) | duplicated(patient_id, fromLast = TRUE)]
  print(tcga_clin$sample_id[patient_id %in% unique(dup_ids)])
}

out <- data.frame(sample_id = tcga_clin$sample_id, patient_id = patient_id,
                   stringsAsFactors = FALSE) %>%
  left_join(thenn_grade, by = "patient_id")

stopifnot(nrow(out) == 195)   # join must not fan out

# ---- 5. Write output ----

write.table(out, "data/tcga_basal_grade_thennavan.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("data/tcga_basal_grade_thennavan.tsv saved")

# ---- 6. Log the numbers that go into the manuscript ----

n_with_grade    <- sum(!is.na(out$grade))
n_without_grade <- sum(is.na(out$grade))
n_not_in_thennavan <- sum(is.na(match(out$patient_id, thenn_grade$patient_id)))

cat("\n================ TCGA Basal-like (n=195) -- histologic grade from Thennavan et al. 2021 ================\n")
cat("With grade:   ", n_with_grade,    "\n")
cat("Without grade:", n_without_grade, "\n")
cat("  of which, patient not present at all in Thennavan's 1,063-patient cohort:",
    n_not_in_thennavan, "\n")
cat("  of which, patient present but >=1 Nottingham component unscored:",
    n_without_grade - n_not_in_thennavan, "\n")
cat("\nGrade distribution (TCGA Basal-like, n =", n_with_grade, "with grade available):\n")
print(table(out$grade, useNA = "no"))
