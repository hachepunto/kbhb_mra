# kbhb_mra

[![DOI](https://zenodo.org/badge/1266709469.svg)](https://doi.org/10.5281/zenodo.20768838)

**Master Regulators of the Kbhb Transcriptional Program in Basal Breast Cancer**

Computational pipeline to identify and characterize transcriptional master regulators (TMRs) of the β-hydroxybutyrylation (Kbhb) proteome in PAM50 Basal-like breast cancer, using two independent patient cohorts (TCGA-BRCA and METABRIC).

---

## Overview

1. Download and normalize TCGA-BRCA RNA-seq and METABRIC microarray data
2. Infer cohort-specific transcriptional regulatory networks with ARACNe-AP
3. Identify Kbhb-program TMRs via msVIPER (restricted signature)
4. Cross-cohort meta-analysis (Stouffer method)
5. Differential expression of Kbhb genes in Basal vs. Normal
6. Visualization: ORA dotplot, circos plot, Sankey diagram

---

## Scripts

| Script | Description |
|--------|-------------|
| `brca_tcga_mtbrc.R` | Download TCGA-BRCA RNA-seq (TCGAbiolinks) and METABRIC (cBioPortal); pre-process; batch correction (ComBat); filter protein-coding genes |
| `basal_pre_networks.R` | Prepare expression matrices for ARACNe-AP network inference |
| `mra_kbhb.R` | msVIPER MRA with Kbhb signature; shadow analysis; Stouffer meta-analysis |
| `de_kbhb.R` | Differential expression (DESeq2 / limma) and cross-cohort concordance classification |
| `compare_kbhb_mrs.R` | ORA of TMR regulons (clusterProfiler + ReactomePA); figures; summary tables |
| `basal_subtype_clinical.R` | Basal-subtype associations (Lehmann/Bareche, IntClust), continuous composite score, de novo clustering, survival (KM/Cox); editor-requested de novo cluster vs. histologic grade association |
| `tcga_histologic_grade.R` | Derives TCGA-BRCA histologic (Nottingham) grade from Thennavan et al. 2021 Data S2 (tubule formation, nuclear pleomorphism, mitotic count), joined to the Basal cohort by 12-char patient barcode — TCGA's central clinical data does not report grade for this tumor type |
| `figure1_panel_c.R` | Redesigned Figure 1 panel C (ComplexHeatmap: PAM50/IntClust/de novo cluster + clinical annotations) and final Figure 1 assembly; editor-requested ER/PR/HER2 vs. de novo cluster association. TCGA clinical fields are consolidated per-field by coverage across 4 candidate sources (`cd` colData, `TCGAquery_subtype()`, cBioPortal `brca_tcga_pub` 2012, GDC BCR Biotab) — see in-script comments for the field-by-field justification |
| `tcga_ic10_classification.R` | Applies the published iC10 classifier (Ali et al. 2014, Genome Biology) to TCGA-BRCA Basal samples via GDC copy-number segments, to obtain an IntClust-equivalent genomic annotation for TCGA (not available from any public TCGA source otherwise); cross-tabulates against the de novo TMR cluster |
| `circos_tmr_kbhb.R` | Circos plot (TMR → DE Kbhb genes, ARACNe support) + Sankey diagram |

---

## Dependencies

### R packages
```r
# Data retrieval & processing
TCGAbiolinks, SummarizedExperiment, sva, DESeq2, limma

# Network inference client
# ARACNe-AP runs separately on a Linux server — see basal_pre_networks.R

# MRA
viper

# Visualization
circlize, ComplexHeatmap, ggplot2, ggrepel, ggalluvial,
clusterProfiler, ReactomePA, org.Hs.eg.db,
RColorBrewer, patchwork, pheatmap, ggraph, tidygraph

# Genomic subtype classification
iC10, iC10TrainingData (CRAN; iC10 also requires Bioconductor package `impute`)

# Utilities
dplyr, tidyr, purrr, readr, tibble, stringr, vroom, janitor, jsonlite
```

### External tool
ARACNe-AP (Lachmann et al. 2016) — run independently on cluster before `mra_kbhb.R`.

---

## Input data

Place the following files in `data/` before running:

| File | Source |
|------|--------|
| `kbhb_genes.txt` | Huang et al. 2021 *Sci Adv* — Kbhb proteome gene symbols |
| `tcga_basal_network.txt` | ARACNe-AP consolidated network (TCGA Basal) |
| `mtbrc_basal_network.txt` | ARACNe-AP consolidated network (METABRIC Basal) |

TCGA and METABRIC expression data are downloaded automatically by `download_tcga_rnaseq.R` and `brca_tcga_mtbrc.R`.

---

## Execution order

```bash
Rscript brca_tcga_mtbrc.R        # 1. Download + pre-process TCGA and METABRIC
Rscript basal_pre_networks.R     # 2. Export matrices for ARACNe-AP
# → Run ARACNe-AP on cluster (external, manual — see "Manual/external steps" below)
Rscript mra_kbhb.R               # 3. MRA + meta-analysis
Rscript de_kbhb.R                # 4. Differential expression
Rscript compare_kbhb_mrs.R       # 5. ORA and summary figures
Rscript basal_subtype_clinical.R # 5a. Subtype associations, composite score, de novo clustering, survival
Rscript tcga_histologic_grade.R  # 5b. TCGA histologic grade from Thennavan et al. 2021 Data S2
Rscript figure1_panel_c.R        # 5c. Panel C (heatmap) + Figure 1 assembly (reads 5a and 5b's outputs)
Rscript tcga_ic10_classification.R # 5d. TCGA IntClust-equivalent genomic classification (iC10)
Rscript circos_tmr_kbhb.R        # 6. Circos + Sankey visualization
Rscript build_supplementary_tables_xlsx.R # 7. Consolidate Supplementary_Table*.tsv into .xlsx
```

`tcga_histologic_grade.R` must run before `figure1_panel_c.R`: the latter reads
`data/tcga_basal_grade_thennavan.tsv` for the Figure 1C Grade track and errors out
(`stopifnot`) if it's missing. Likewise, `build_supplementary_tables_xlsx.R` only
packages whatever `Supplementary_Table*.tsv` files are already on disk — run the
scripts that produce them first if a table is missing or stale.

There is no single `run_all.R`: several steps above are not scripted R at all
(ARACNe-AP, a separate compute-cluster job) or depend on external services that can
fail, stall, or require credentials/interaction (GDC, cBioPortal) — a one-shot
driver would misrepresent this pipeline as more push-button than it is.

### Manual/external steps

- **ARACNe-AP** (Lachmann et al. 2016) runs outside R, on a Linux compute cluster,
  between `basal_pre_networks.R` and `mra_kbhb.R` — see `run_aracne_ap_luminal_her2.sh`
  for the shell-script pattern used (the Basal-cohort run predates that script and
  isn't itself checked in). Its output (`data/tcga_basal_network.txt`,
  `data/mtbrc_basal_network.txt`) is already versioned in this repo, so this step
  only needs to be repeated if those networks are intentionally regenerated —
  **do not re-run it** otherwise (hours of compute; see the note in `CLAUDE.md`).
- **GDC downloads** happen automatically via `TCGAbiolinks::GDCquery()`/`GDCdownload()`
  inside `brca_tcga_mtbrc.R` (RNA-seq + clinical), `tcga_ic10_classification.R`
  (copy-number segments), and `figure1_panel_c.R` (BCR Biotab clinical supplement) —
  no manual step required, but each needs network access to the GDC API and can be
  slow or fail transiently.
- **METABRIC** raw files are fetched from cBioPortal by `brca_tcga_mtbrc.R`
  (`download.file()`, skipped if `brca_metabric/` is already present) — same
  external-network caveat as GDC.
- **Thennavan et al. 2021 Data S2** (`PMID_35465400/*.xlsx`) is third-party
  supplementary material already included in this repository (see `CLAUDE.md`) —
  nothing to download, `tcga_histologic_grade.R` reads it in place.

---

## Key outputs

Final manuscript figures and result tables are versioned in this repository under
`figures/` and `data/`. Figure numbers correspond to the submitted manuscript; the
script that generates the underlying plot/table is noted alongside.

### Figures

| File | Description | Generated by |
|------|-------------|--------------|
| `figures/Figure1.pdf` | 3-panel: (A) NES scatter TCGA vs. METABRIC, (B) meta-NES lollipop, (C) activity heatmap of the 7 significant Kbhb TMRs annotated by PAM50, IntClust, de novo cluster, and clinical variables (histology, grade, age, ER/PR/HER2) | `compare_kbhb_mrs.R` (panels A/B) + `basal_subtype_clinical.R` (stats, de novo cluster) + `figure1_panel_c.R` (panel C + final assembly) |
| `figures/Figure2.pdf` | Lollipop plot, meta-analysis NES | `compare_kbhb_mrs.R` |
| `figures/Figure3.pdf` | ORA dotplot of TMR regulons (GO-BP + Reactome) | `compare_kbhb_mrs.R` |
| `figures/Figure4.pdf` | Circos plot: TMR → DE Kbhb genes (ARACNe support) | `circos_tmr_kbhb.R` |
| `figures/Figure5.pdf` | Graded specificity of the Kbhb TMR panel across PAM50 subtypes (dot plot, 3 tiers: broad / partial / Basal-like-exclusive) | `luminal_her2_mra_kbhb.R` |
| `figures/Figure6.pdf` | Graphical summary of the BHB–Kbhb–TMR axis | manual |
| `figures/FigureS1.pdf` | msVIPER results and shadow-pair regulatory hierarchy per cohort | `mra_kbhb.R` |
| `figures/FigureS2.pdf` | DE volcano plots and heatmaps, TCGA and METABRIC | `de_kbhb.R` |
| `figures/FigureS3.pdf` | Kaplan-Meier survival, de novo TMR clusters: (A) TCGA, (B) METABRIC | `basal_subtype_clinical.R` |
| `figures/Supplementary_FigureS4_receptor_rate_by_cluster.pdf` | % ER/PR/HER2-positive by de novo TMR cluster, per cohort (ggplot2; read back from Table S4b, cannot diverge from it) | `figure1_panel_c.R` |
| `figures/Supplementary_FigureS5_sankey_tmr_kbhb.pdf` | Sankey diagram: TMR regulons → DE Kbhb gene categories | `circos_tmr_kbhb.R` |

### Result tables

| File | Description | Generated by |
|------|-------------|--------------|
| `data/kbhb_mrs_comparison.tsv` | Cross-cohort Stouffer meta-analysis of TMR NES | `compare_kbhb_mrs.R` |
| `data/kbhb_consensus.tsv` | Cross-cohort DE concordance table (Kbhb genes) | `de_kbhb.R` |
| `data/kbhb_de_tcga.tsv`, `data/kbhb_de_mtbrc.tsv` | Full DE results, Kbhb genes (TCGA / METABRIC) | `de_kbhb.R` |
| `data/tmr_de_tcga.tsv`, `data/tmr_de_mtbrc.tsv` | DE results restricted to the 7 TMRs | `de_kbhb.R` |
| `data/ORA_tmr_counts.tsv` | Pathway counts per TMR regulon | `compare_kbhb_mrs.R` |
| `data/ORA_top_shared_pathways.tsv` | Top 20 most shared pathways across TMRs | `compare_kbhb_mrs.R` |
| `data/Supplementary_TableS1_ORA_GO_BP.tsv`, `data/Supplementary_TableS2_ORA_Reactome.tsv`, `data/ORA_pathway_index.tsv` | Full ORA results (GO-BP, Reactome) and pathway index | `compare_kbhb_mrs.R` |
| `data/circos_*.tsv` | TMR–gene interaction data underlying the circos plot | `circos_tmr_kbhb.R` |
| `data/tcga_basal_network.txt`, `data/mtbrc_basal_network.txt` | ARACNe-AP consolidated regulatory networks (TCGA / METABRIC Basal) | ARACNe-AP (external) |
| `data/tcga_basal_clinical_extra.tsv` | TCGA Basal clinical annotation (PAM50, ER/PR/HER2, histology, age) as consolidated for Figure 1C — see per-field source notes in the Scripts table above | `figure1_panel_c.R` |
| `data/tcga_basal_ic10_calls.tsv` | iC10 genomic subtype call per TCGA Basal sample | `tcga_ic10_classification.R` |
| `data/Supplementary_TableS4_denovo_cluster_receptor_status.tsv`, `data/Supplementary_TableS4b_..._contingency.tsv` | Editor-requested association: de novo TMR cluster vs. ER/PR/HER2 status, per cohort (Fisher/chi-squared + BH; S4b = supporting contingency counts) | `figure1_panel_c.R` |
| `data/Supplementary_TableS5_tcga_ic10_denovo_cluster.tsv`, `data/Supplementary_TableS5b_..._full.tsv` | iC10 group vs. de novo TMR cluster (TCGA): collapsed, well-powered IntClust{8,9,10}-vs-rest 2×2 (S5, reported result) and the full 10-group table (S5b, reference; several low-n cells) | `tcga_ic10_classification.R` |
| `data/Supplementary_TableS6_metabric_ic10_denovo_cluster.tsv`, `data/Supplementary_TableS6b_..._full.tsv` | IntClust vs. de novo TMR cluster (METABRIC): collapsed IntClust10-vs-rest 2×2 (S6, reported result) and the full ≥10-group table (S6b, reference; several low-n cells) | `figure1_panel_c.R` |
| `data/Supplementary_TableS7_denovo_cluster_mean_NES.tsv` | Mean Kbhb-TMR NES per TMR, per Cohort × De novo cluster block (the same 4 blocks split in Figure 1C) | `figure1_panel_c.R` |

### Not included in this repository

Processed and raw TCGA-BRCA / METABRIC expression matrices are **not** versioned
here. They are third-party data, already publicly available and cited via GDC
(TCGA-BRCA) and cBioPortal (METABRIC), and are regenerated locally by
`brca_tcga_mtbrc.R` and `basal_pre_networks.R`. Large intermediate `.rds` objects
(full msVIPER null models, raw `SummarizedExperiment` downloads) are likewise
excluded and are regenerable from the scripts.

---

## Citation

Manuscript in preparation — IJMS (MDPI).

Code/repository: [10.5281/zenodo.20768838](https://doi.org/10.5281/zenodo.20768838)

Kbhb proteome: Huang et al. (2021) *Sci Adv* 7:eabe2771  
MRA framework: Alvarez et al. (2016) *Nat Genet* 48:838  
ARACNe-AP: Lachmann et al. (2016) *Bioinformatics* 32:2233  
TF catalog: Lambert et al. (2018) *Cell* 172:458
