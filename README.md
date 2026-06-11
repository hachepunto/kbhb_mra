# kbhb_mra

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
| `download_tcga_rnaseq.R` | Download TCGA-BRCA RNA-seq from GDC (TCGAbiolinks) |
| `brca_tcga_mtbrc.R` | Pre-process TCGA and METABRIC; batch correction (ComBat); filter protein-coding genes |
| `basal_pre_networks.R` | Prepare expression matrices for ARACNe-AP network inference |
| `mra_kbhb.R` | msVIPER MRA with Kbhb signature; shadow analysis; Stouffer meta-analysis |
| `de_kbhb.R` | Differential expression (DESeq2 / limma) and cross-cohort concordance classification |
| `compare_kbhb_mrs.R` | ORA of TMR regulons (clusterProfiler + ReactomePA); figures; summary tables |
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
| `metabric_basal_network.txt` | ARACNe-AP consolidated network (METABRIC Basal) |

TCGA and METABRIC expression data are downloaded automatically by `download_tcga_rnaseq.R` and `brca_tcga_mtbrc.R`.

---

## Execution order

```bash
Rscript download_tcga_rnaseq.R   # 1. Download TCGA
Rscript brca_tcga_mtbrc.R        # 2. Pre-process both cohorts
Rscript basal_pre_networks.R     # 3. Export matrices for ARACNe-AP
# → Run ARACNe-AP on cluster
Rscript mra_kbhb.R               # 4. MRA + meta-analysis
Rscript de_kbhb.R                # 5. Differential expression
Rscript compare_kbhb_mrs.R       # 6. ORA and summary figures
Rscript circos_tmr_kbhb.R        # 7. Circos + Sankey visualization
```

---

## Key outputs

| File | Description |
|------|-------------|
| `figures/fig_circos_tmr_kbhb.pdf` | Main figure: circos plot TMR → DE Kbhb genes |
| `figures/fig_sankey_tmr_kbhb.pdf` | Supp. figure: Sankey TMR → DE categories |
| `figures/fig_ora_tmrs.pdf` | ORA dotplot (GO-BP + Reactome) |
| `data/kbhb_consensus.tsv` | Cross-cohort DE concordance table |
| `data/ORA_tmr_counts.tsv` | Pathway counts per TMR regulon |
| `data/ORA_top_shared_pathways.tsv` | Top 20 most shared pathways across TMRs |

---

## Citation

Manuscript in preparation — IJMS (MDPI).

Kbhb proteome: Huang et al. (2021) *Sci Adv* 7:eabe2771  
MRA framework: Alvarez et al. (2016) *Nat Genet* 48:838  
ARACNe-AP: Lachmann et al. (2016) *Bioinformatics* 32:2233  
TF catalog: Lambert et al. (2018) *Cell* 172:458
