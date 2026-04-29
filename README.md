# Nolt et al.   
# Analysis code and processed data

Analysis scripts, processed data, and rendered Quarto reports backing the bulk-RNA-seq and GSEA results in:

> Nolt G.L., MacLean S.M., Golden L.R., Thorpe S.P., Stephens I.O., Funnell J.L., Brock C.R., Lucido C.C., Smith L., Hernandez G., Arbones-Mainar J.M., Wood C., Morganti J.M., Johnson L.A. *Aged APOE4 mice show limited immunometabolic adaptation to APOE2 microglial replacement.* Neurobiology of Aging (in revision). DOI: TBD.

## Scope

This repository contains the code and processed data needed to reproduce the panels of the manuscript that were generated in R / Quarto:

- **Figure 4B** — Volcano plot of differentially expressed genes (4s2M vs 4s2 LPS, whole brain).
- **Figure 4C** — GSEA bar plot (GO Biological Process, collapsed pathways).

Other figures in the manuscript (Figs 1–3, 5, 6, all supplementary figures) come from immunofluorescence, qPCR, Luminex cytokine, lipidomics, and metabolomics analyses performed outside this codebase. See *Data availability* below for where to find them.

## Data availability

| Data | Location |
|---|---|
| Raw RNA-seq FASTQs | NCBI GEO accession (_will be activated upon submission_) |
| Salmon gene-level counts (the input to DESeq2) | `nf-core_results/salmon_counts/salmon.merged.gene_counts.tsv` |
| Sample metadata (`mouseID` ↔ `altID` ↔ `genotype`) | `20250905_Bulk seq IDs.xlsx` |
| Differential expression result table | `analysis/20250909-DEGs.csv` |
| Lipidomics & metabolomics scripts and data | _PENDING_ |
| Luminex cytokine, immunofluorescence / HALO, and qPCR data (Figs 1–3, 5, 6) | Source-data tables of the published article; statistics performed in GraphPad Prism (see Methods). |

## Reproducing the figures

All notebooks are self-contained Quarto documents (`embed-resources: true`) and render to a single HTML file alongside the `.qmd`.


## Upstream pipeline (not run locally)

FASTQs were processed with [`nf-core/rnaseq` v3.19.0](https://nf-co.re/rnaseq/3.19.0) on the University of Kentucky Morgan HPC cluster, executed by [`Nextflow` v24.10.4](https://www.nextflow.io) and containerized with Singularity. Reads were aligned to the *Mus musculus* GRCm39 primary assembly (Ensembl release 114) using STAR, and transcript abundance was quantified with Salmon. The exact submission script is `fastq/20250730-nextflow_rnaseq.sh`, and the sample sheet handed to the pipeline is `20250907-sample_sheet.csv` (built by `20250907-sample_sheet_script.R`).

Reviewers can re-run the full upstream pipeline by retrieving the FASTQs from GEO (accession above) and invoking the submission script after adjusting the cluster-specific paths.

## Software environment

R 4.5.x with the following packages (versions current at time of analysis):

- `tidyverse` 2.0.0
- `DESeq2` 1.50.2
- `EnhancedVolcano` 1.24.0
- `fgsea`, `msigdbr` (mouse `M5/BP` sets)
- `kableExtra`, `patchwork`

## Caveats

- **Sample `M2` is excluded from all differential-expression analyses** as E2 expression is incompatible with a microglia-specific switch and is most likely a genotyping error. 



## Contact

- Lance A. Johnson — corresponding author, johnson.lance@uky.edu

