These files produce the results shown in figures 1-3, focusing on cHSPC inter-individual heterogeneity among healthy individuals.
To reproduce the results, you should:
- Download the code from this repository.
- Download dependencies for the code: pheatmap, zoo, dplyr, metacell, openxlsx, tgstat, Matrix, anndata, umap, glmnet, RColorBrewer.
- Download required (large) files from our AWS bucket.
- Download required (small) files from this repository.
- Edit `params.R` as follows:
  - set `MODEL.DIR` to the path in which the AWS files are located (they should all be in the same directory).
  - set `DATA.DIR` to the root directory in which the small data files (downloaded from this repository) are located. The directory structure for the small files should match the directory structure in the repository.
  - set `BASE.FIG.DIR` to a directory where figures will be saved.
  - set `SUPP.TABLE.DIR` to a directory where tables will be saved.
- Execute the `fig1` function in `analyze_pbmc_model.R` to create Figure 1 and its corresponding Extended Data Figures and tables.
- Execute the `fig2` function in `analyze_indivs.R` to create Figures 2 and 3, and their corresponding Extended Data Figures and tables.

