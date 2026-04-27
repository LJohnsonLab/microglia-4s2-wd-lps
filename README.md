# Nolt et al. — analysis code and processed data

Analysis scripts, processed data, and rendered Quarto reports backing the bulk-RNA-seq, GSEA, and Luminex cytokine results in:

> Nolt G.L., MacLean S.M., Golden L.R., Thorpe S.P., Stephens I.O., Funnell J.L., Brock C.R., Lucido C.C., Smith L., Hernandez G., Arbones-Mainar J.M., Wood C., Morganti J.M., Johnson L.A. *Aged APOE4 mice show limited immunometabolic adaptation to APOE2 microglial replacement.* Neurobiology of Aging (in revision). DOI: TBD.

## Scope

This repository contains the code and processed data needed to reproduce the panels of the manuscript that were generated in R / Quarto:

- **Figure 4B** — Volcano plot of differentially expressed genes (4s2M vs 4s2 LPS, whole brain).
- **Figure 4C** — GSEA bar plot (GO Biological Process, collapsed pathways).
- **Figure 5B–I** — Luminex cytokine summary tables (plasma pre / post LPS and brain).

Other figures in the manuscript (Figs 1–3, 5A, 6, all supplementary figures) come from immunofluorescence, qPCR, lipidomics, and metabolomics analyses performed outside this codebase. See *Data availability* below for where to find them.

## Data availability

| Data | Location |
|---|---|
| Raw RNA-seq FASTQs | NCBI GEO accession (_will be activated upon submission_) |
| Salmon gene-level counts (the input to DESeq2) | `nf-core_results/salmon_counts/salmon.merged.gene_counts.tsv` |
| Sample metadata (`mouseID` ↔ `altID` ↔ `genotype`) | `20250905_Bulk seq IDs.xlsx` |
| Luminex cytokine measurements | `Raw Data Johnson_MacLean Cytokine Results 2023 0612.csv` |
| Differential expression result table | `analysis/20250909-DEGs.csv` |
| Lipidomics & metabolomics scripts and data | _PENDING_ |
| Immunofluorescence / HALO and qPCR data (Figs 1–3, 5A, 6) | Source-data tables of the published article; statistics performed in GraphPad Prism (see Methods). |

## Reproducing the figures

All notebooks are self-contained Quarto documents (`embed-resources: true`) and render to a single HTML file alongside the `.qmd`.

### Figure 4B — DESeq2 differential expression
```bash
quarto render analysis/20250909-switch_lps.qmd
```
- **Inputs**: `nf-core_results/salmon_counts/salmon.merged.gene_counts.tsv`, `20250905_Bulk seq IDs.xlsx`.
- **Outputs**: `analysis/20250909-switch_lps.html` (volcano plot, PCA, sample table) and `analysis/20250909-DEGs.csv` (full DE table, used by Fig. 4C).
- DE design: `~ genotype`, reference level `+/+` (so positive log2FC = up in `MC/+` / 4s2M).
- Filtering: ≥10 counts in ≥5 samples, then median count > 10.
- Sample `M2` is dropped as a QC outlier (sentinel `select(-M2)` in the notebook).

### Figure 4C — GSEA against GO Biological Process
```bash
quarto render analysis/20251101-GSEA.qmd
```
- **Input**: `analysis/20250909-DEGs.csv` plus mouse GO:BP gene sets retrieved at render time via `msigdbr` (MSigDB collection `M5`, subcollection `BP`).
- **Output**: `analysis/20251101-GSEA.html` — full ranked-pathway table and bar plot of normalized enrichment scores.
- Method: `fgsea` with `minSize = 30`, `maxSize = 500`, FDR < 0.01, then `collapsePathways` to remove redundant gene-set overlap. Reproduces the 188 → 58 pathway counts cited in Methods §2.10.3.

### Figure 5B–I — Luminex cytokines
```bash
quarto render analysis/20250912-switch_lps_bulk+cytoq.qmd
```
- **Inputs**: `Raw Data Johnson_MacLean Cytokine Results 2023 0612.csv`, `20250905_Bulk seq IDs.xlsx`.
- **Output**: `analysis/20250912-switch_lps_bulk+cytoq.html` — within-genotype longitudinal tables, between-genotype comparisons stratified by condition, and a Benjamini–Hochberg-adjusted summary across the 15 (5 cytokines × 3 conditions) genotype tests.

## Upstream pipeline (not run locally)

FASTQs were processed with [`nf-core/rnaseq` v3.19.0](https://nf-co.re/rnaseq/3.19.0) on the University of Kentucky Morgan HPC cluster, executed by [`Nextflow` v24.10.4](https://www.nextflow.io) and containerized with Singularity. Reads were aligned to the *Mus musculus* GRCm39 primary assembly (Ensembl release 114) using STAR, and transcript abundance was quantified with Salmon. The exact submission script is `fastq/20250730-nextflow_rnaseq.sh`, and the sample sheet handed to the pipeline is `20250907-sample_sheet.csv` (built by `20250907-sample_sheet_script.R`).

Reviewers can re-run the full upstream pipeline by retrieving the FASTQs from GEO (accession above) and invoking the submission script after adjusting the cluster-specific paths.

## Software environment

R 4.5.x with the following packages (versions current at time of analysis):

- `tidyverse` 2.0.0
- `DESeq2` 1.50.2
- `EnhancedVolcano` 1.24.0
- `fgsea`, `msigdbr` (mouse `M5/BP` sets)
- `compareGroups` 4.10.2
- `kableExtra`, `patchwork`

## Caveats

- **Sample `M2` is excluded from all differential-expression analyses.** as E2 expression is incompatible with a microglia-specific switch and is most likely a sample-tracking / genotyping error or off-target Cre recombination. 



## Contact

- Lance A. Johnson — corresponding author, johnson.lance@uky.edu

