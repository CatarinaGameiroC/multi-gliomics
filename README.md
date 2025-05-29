# Uncovering Glioma Heterogeneity through Sparse Multi-Omics Analysis  

This repository contains code for analyzing glioma heterogeneity using two sparse dimensionality reduction methods: **MOFA** (Multi-Omics Factor Analysis) and **sGCCA** (Sparse Generalized Canonical Correlation Analysis).

## Structure
- Scripts for each method are organized in their respective folders: `MOFA/` and `sGCCA/`.
- All preprocessing and analysis steps are provided in `.Rmd` files.

## Datasets
- Only the **mutations** dataset is included directly in the repository. Only the mutations dataset and the corresponding ground-truth labels for each patient are included directly in this repository.
- The ground-truth labels are based on the classification described in [this study](https://www.biorxiv.org/content/10.1101/2023.02.19.529134v3.full.pdf).
- Other omics datasets (e.g., mRNA, miRNA, methylation) can be downloaded using the code in the provided `.Rmd` notebooks.
- **Note:** The **DNA methylation** dataset is ~8 GB. Due to its size, we recommend using the GDC Data Transfer Tool:
  1. Obtain the manifest file using the code in the script.
  2. Download the data via terminal using the [GDC client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool).


## Notes
- This project is part of a Master’s thesis.




