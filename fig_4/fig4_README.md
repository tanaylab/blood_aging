These files produce the results shown in figure 4, focusing on applications of the cHSPC reference model (described in figures 1-3) to diagnostics.
<br/><br/>

Setup and run order:
1. Create a conda environment using `blood_aging/fig_4/blood_aging_conda_spec_file.txt`.
2. pip install libraries, see below.
3. copy large files from our AWS bucket as described below.
4. generate_donor_features_and_estimate_karyotype.ipynb
5. generate_combined_illumina_ultima_feature_data.ipynb
6. mds_xgboost_classification.ipynb
7. generate_blood_aging_paper_figs_and_tables.ipynb
<br/><br/>

pip install command:
```console
$ pip install metacells rpy2 numpy scipy scikit-learn statsmodels seaborn biopython openpyxl ipykernel notebook scikit-learn-intelex threadpoolctl ipympl adjustText ruptures xgboost shap ipython
```
<br/><br/>
For the code to run, you have to download some large files from our AWS bucket and move them to the expected paths, e.g. using:
```console
$ cp blood_aging_aws_bucket/fig4/ult_mds_cells.h5ad blood_aging/output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/cells_with_metacell_attrs.h5ad
$ cp blood_aging_aws_bucket/fig4/ult_mds_metacells.h5ad blood_aging/output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/metacells_with_projection.h5ad
$ cp blood_aging_aws_bucket/fig4/illu_mds_cells.h5ad blood_aging/output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/cells_with_metacell_attrs.h5ad
$ cp blood_aging_aws_bucket/fig4/illu_mds_metacells.h5ad blood_aging/output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/metacells_with_projection.h5ad
$ cp blood_aging_aws_bucket/fig4/illu_ref_cells.h5ad blood_aging/output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_c.h5ad
$ cp blood_aging_aws_bucket/fig4/illu_ref_metacells.h5ad blood_aging/output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_mc.h5ad
$ cp blood_aging_aws_bucket/fig4/genes.gtf blood_aging/input/human_genome_files/genes.gtf
$ cp blood_aging_aws_bucket/fig4/genome.fa blood_aging/input/human_genome_files/genome.fa
```


