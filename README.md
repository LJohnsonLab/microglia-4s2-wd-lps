# Nolt et al.   
# Analysis code and processed data

Analysis scripts, processed data, and rendered Quarto reports backing the bulk-RNA-seq and GSEA results in:

> Nolt G.L., MacLean S.M., Golden L.R., Thorpe S.P., Stephens I.O., Funnell J.L., Brock C.R., Lucido C.C., Smith L., Hernandez G., Arbones-Mainar J.M., Wood C., Morganti J.M., Johnson L.A. *Aged APOE4 mice show limited immunometabolic adaptation to APOE2 microglial replacement.* Neurobiology of Aging (in revision). DOI: TBD.

## Scope

This repository contains the code and processed data needed to reproduce the panels of the manuscript that were generated in R / Quarto:

- **Figure 4B** — Volcano plot of differentially expressed genes (4s2M vs 4s2 LPS, whole brain).
- **Figure 4C** — GSEA bar plot (GO Biological Process, collapsed pathways).
- **Figures 1–3 (HFD cohort) — lipidomics & metabolomics** — full per-feature two-way Sex × Genotype ANOVA tables and emmeans post hoc contrasts, in `lipidomics-metabolomics-public/`.

Other figures in the manuscript (Fig 5 cytokines, Fig 6, supplementary figures) come from immunofluorescence, qPCR, and Luminex analyses performed outside this codebase. See *Data availability* below for where to find them.

## Data availability

| Data | Location |
|---|---|
| Raw RNA-seq FASTQs | NCBI GEO accession (_will be activated upon submission_) |
| Salmon gene-level counts (the input to DESeq2) | `nf-core_results/salmon_counts/salmon.merged.gene_counts.tsv` |
| Sample metadata (`mouseID` ↔ `altID` ↔ `genotype`) | `20250905_Bulk seq IDs.xlsx` |
| Differential expression result table | `analysis/20250909-DEGs.csv` |
| Lipidomics & metabolomics — scripts, raw measurements, and full ANOVA / emmeans output tables (HFD cohort, Figs 1–3) | `lipidomics-metabolomics-public/` |
| Luminex cytokine, immunofluorescence / HALO, and qPCR data (Figs 5, 6) | Source-data tables of the published article; statistics performed in GraphPad Prism (see Methods). |

## Reproducing the figures

All notebooks are self-contained Quarto documents (`embed-resources: true`) and render to a single HTML file alongside the `.qmd`.

### Lipidomics & metabolomics (Figs 1–3, HFD cohort)

The HFD cohort lipidomics and metabolomics analyses are independent of the
LPS-cohort RNA-seq pipeline and live in `lipidomics-metabolomics-public/`.
Each analysis is a single self-contained R script that consumes the raw IOS
measurements and the matching metadata, performs tissue-weight normalization
and log2 transform, fits a two-way Sex × Genotype linear model, and writes
both the per-feature ANOVA table and the BH-adjusted emmeans post hoc
contrasts:

```r
# Lipidomics — per-lipid and per-subclass
setwd("lipidomics-metabolomics-public/Georgia HFD/Lipids")
source("../../GN_HFD_Lipidomics_ANOVA.R")
# → TwoWay_SexGen_withWeightNorm_outputs_v7b/

# Metabolomics — per-metabolite, with global BH across all post hoc tests
setwd("lipidomics-metabolomics-public/Georgia HFD/Metabolites")
source("../../GN_HFD_Metabolomics_ANOVA.R")
# → TwoWay_Metabolites_SexGen_outputs_v2_globalBH/
```

The repository ships both the input CSVs and the full output tables, so
the results can either be inspected directly or reproduced by rerunning the
scripts.

**Multiple-testing correction.** Every table uses Benjamini–Hochberg, but
the *family* over which BH is applied differs by table. Each output column
is named so the family is unambiguous:

- **Per-feature 2-way ANOVA tables** (`per_lipid_2way_ANOVA.csv`,
  `per_subclass_2way_ANOVA.csv`, `per_metabolite_2way_ANOVA.csv`).
  Column `p_adj_BH` — BH applied across **all features** within each
  ANOVA term (Sex, Genotype, Sex × Genotype).
- **Lipidomics emmeans post hoc tables** (`per_lipid_emmeans_*`,
  `per_subclass_emmeans_*`). Two adjusted columns are reported side by
  side:
  - `p_adj` — BH applied **within each individual lipid**, across the
    pairwise contrasts of that single feature (the value the manuscript
    text quotes).
  - `p_adj_BH_global` — BH applied **across every feature × every
    pairwise contrast** within the post hoc family (ALL_pairs,
    Sex_within_Genotype, or Genotype_within_Sex). This matches the
    correction used for metabolomics and is the more conservative test.
- **Metabolomics emmeans post hoc tables**
  (`per_metabolite_emmeans_*_GLOBALBH.csv`). Column `p_adj_BH_global`
  only — same family-wide BH as the lipid `p_adj_BH_global` column.

Both omics tables therefore expose the same family-wide correction
(`p_adj_BH_global`), and the lipid tables additionally retain the
within-feature `p_adj` for traceability with the manuscript text.


## Upstream pipeline (not run locally)

FASTQs were processed with [`nf-core/rnaseq` v3.19.0](https://nf-co.re/rnaseq/3.19.0) on the University of Kentucky Morgan HPC cluster, executed by [`Nextflow` v24.10.4](https://www.nextflow.io) and containerized with Singularity. Reads were aligned to the *Mus musculus* GRCm39 primary assembly (Ensembl release 114) using STAR, and transcript abundance was quantified with Salmon. The exact submission script is `fastq/20250730-nextflow_rnaseq.sh`, and the sample sheet handed to the pipeline is `20250907-sample_sheet.csv` (built by `20250907-sample_sheet_script.R`).

The full upstream pipeline can be re-run by retrieving the FASTQs from GEO (accession above) and invoking the submission script after adjusting the cluster-specific paths.

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

