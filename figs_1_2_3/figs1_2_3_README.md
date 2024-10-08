Files in the directory produce the results shown in figures 1-3, focusing on cHSPC inter-individual heterogeneity among healthy individuals.
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


The code files in this directory are:
- `analyze_pbmc_model.R` - creates Figure 1 and its corresponding Extended Data Figures and tables.
- `analyze_indivs.R` - creates Figures 2-3 and their corresponding Extended Data Figures and tables.
- `utils.R` - some utility functions.
- `params.R` - parameters file.

The files needed to run the code are listed below.

AWS files:
- 148_indiv_ref_cells.h5ad - UMI matrix and cell metadata for the cells in the reference model containig 148 individuals.
- 148_indiv_ref_metacells.h5ad - metacell model for the cells from 148_indiv_ref_cells.h5ad.
- 148_indiv_ref_metacells_unfiltered.h5ad - metacell model for cells belonging to the 148 individuals before cell filtering.
- bm_hca_ref_metacells.h5ad - metacell model for the bone marrow data from the Human Cell Atlas.
- bm_our_ref_metacells.h5ad - metacell model for CD34-enriched bone marrow we collected for this study.
- pbmc_160k_ref_cells.h5ad - UMI matrix and cell metadata for PBMC cells from Zheng et al. 2017.
- pbmc_160k_ref_metacells.h5ad - metacells model for the cells from pbmc_160k_ref_cells.h5ad.
- pbmc_garvan_ref_cells.h5ad - UMI matrix and cell metadata for PBMC cells from Yazar et al. 2022 from the Garvan Institute.
- pbmc_garvan_ref_metacells.h5ad - metacells model for the cells from pbmc_garvan_ref_cells.h5ad.
- bm_palantir_ref_cells.h5ad - UMI matrix and cell metadata for CD34 bone marrow cells from Setty et al. 2019 (Palantir method paper).
- bm_palantir_ref_metacells.h5ad - metacells model for the cells from bm_palantir_ref_cells.h5ad.

Files found in Github (under figs_1_2_3/input in this repository):
- sex_genes - list of genes that are strongly differentially expressed between the sexes.
- tfs.txt - list of transcription factors.
- bm_colors.csv - colors used to display cell types in the bone marrow.
- gene_batch_effect_pvals - contains p-values for a test testing for each gene whether its expression suffers from batch effects.
- n122_got.Rdata - contains data from the GoT experiment performed on individual N122.
- indiv_id_map_df - a table mapping between sample ids to individual id (needed since some individuals have two samples).
- age_sex_original_cohort.csv, age_sex_new_cohort.csv - the age and sex of individuals (divided into two files based on when data from these individuals was collected).
- cbc_original_cohort.csv, cbc_longitudinal_original_cohort.csv, cbc_longitudinal_new_cohort.csv - blood counts for the individuals. The longitudinal files are separated based on when data from these individuals was collected. cbc_original_cohort.csv contains only one count per individual, and is used in case no longitudinal data was available.
- rdw_vs_mut/cohort.csv - contains metadata on the cohort used to test association between RDW and mutation status.
- rdw_vs_mut/arch.xlsx - contains mutation data for the cohort used to test association between RDW and mutation status.
- mutation_files/arch3.csv, mutation_files/arch4.csv - contains mutation data for the 148 individuals in the scRNA reference.
