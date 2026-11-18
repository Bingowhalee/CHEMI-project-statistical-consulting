# CHEMI project – gut microbiota, microbial metabolites and gut health markers

R code for a statistical consulting project (Statistical Consulting course, Leiden University, 2026) carried out for the Department of Human Nutrition and Health, Wageningen University & Research.

**Authors:** Jue Wang and Simar Yadav

## Research questions

In 60 patients with stage I–III colorectal cancer:

1. Which gut bacterial genera are associated with plasma microbial metabolites (short-chain fatty acids, branched-chain fatty acids, tryptophan-derived indoles, TMAO) and with markers of gut health (faecal calprotectin, serum LBP), adjusting for cancer stage and renal function (creatinine)?
2. Which of these metabolites and markers are most strongly related to the overall microbiota composition?

## Data

The patient data are **not included** in this repository because they are confidential and subject to a data agreement with the client. To reproduce the results, place the data file provided by the client in the same folder as the script:

```
Data_CHEMI_statistics_LU_20260929.csv
```

Expected layout: columns 1–18 contain cancer stage, plasma creatinine and the 16 outcomes (metabolites and markers); all remaining columns contain the relative abundances of the bacterial genera.

## Requirements

- R 4.6.1 (other recent versions should work)
- CRAN packages: `vegan`, `glmnet`, `ggplot2`, `tidyr`, `pheatmap`, `patchwork`
- Bioconductor packages: `mixOmics`, `globaltest`

```r
install.packages(c("vegan", "glmnet", "ggplot2", "tidyr", "pheatmap", "patchwork"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("mixOmics", "globaltest"))
```

## How to run

1. Put `R_data.R` and the data file in the same folder.
2. Set this folder as the working directory, e.g. `setwd("path/to/folder")`.
3. Run the script from top to bottom (`source("R_data.R")`), or section by section in RStudio.

All figures and tables are written to the subfolder `Final/`, which is created automatically. A single random seed (`4562550`) is set before every step that involves randomness, so the results are fully reproducible.

Run time: a few minutes in total. The elastic net with leave-one-out cross-validation (section 6a) and the sparse PLS stability analysis (section 6b) take the longest.

## Structure of the script

| Section | Content | Main output in `Final/` |
|---|---|---|
| 0. Setup | Packages, seed, data import | – |
| 1. Prepare data | Prevalence filter (≥ 10% and ≥ 30% of patients), Hellinger transformation of genera, log2 transformation and standardisation of outcomes, readable genus names | – |
| 2. Describe data | Key numbers, Table 1, prevalence of genera | `table1.csv`, `prevalence_histogram.png` |
| 3. Exploration | PCoA, PERMANOVA and dispersion test for cancer stage, metabolite arrows (envfit) | `pcoa_arrows.png` |
| 4. Research question 1 | Genus-by-genus linear regression adjusted for stage and creatinine, Benjamini–Hochberg FDR | `regression_results.csv`, `heatmap_regression.png` |
| 5. Research question 2 | Global test per outcome (all genera together), ranked by FDR-adjusted p-value | `hierarchy_globaltest.png` |
| 6a. Elastic net | Penalised regression with all genera, leave-one-out cross-validation, predictive R² versus covariates only, genus selection frequency | `elastic_net_stable.png`, `elastic_net_final.png` |
| 6b. Sparse PLS | Fixed sparsity (20 genera per component) with stability selection over 200 subsamples of 80% of the patients; heat maps of stable genera | `spls_stability.csv`, `spls_stability.png`, `spls_heatmap.png`, `spls_heatmap30.png` |
| 7. Sensitivity analyses | Adjusting for stage only; genera present in ≥ 30% of patients | `sensitivity_table.csv`, `sensitivity_scatter.png`, `heatmap_regression_30.png`, `heatmap_regression_stage.png` |
| 8. Diagnostics and key outcomes | Residual plots, influential patients and Spearman correlations for key associations; acetate and LBP; comparison of regression and elastic net | `acetate_lbp_genera.png`, diagnostic plots, `comparison_regression_elasticnet.csv` |

## Report

The methods and results are described in the accompanying consulting report: *Gut Microbiota Composition in Relation to Plasma Microbial Metabolites and Markers of Gut Health in Patients with Colorectal Cancer: A Statistical Consulting Report* (Wang & Yadav, 2026).
