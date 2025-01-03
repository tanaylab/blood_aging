These files produce the results shown in figure 4, focusing on applications of the cHSPC reference model (described in figures 1-3) to diagnostics.
<br/><br/>

Setup and run order (we used a CentOS Linux 7 machine):
1. Download this repository into a location with at least 100G free storage space, e.g., using `git clone https://github.com/tanaylab/blood_aging.git`.
2. Download required (large) files from our AWS bucket at https://blood-aging.s3.ap-south-1.amazonaws.com/.
3. Create and activate a conda environment using `blood_aging/fig_4/blood_aging_conda_spec_file.txt`.
4. `pip install -r blood_aging/fig_4/blood_aging_pip_requirements.txt`
5. `cd blood_aging/fig_4` (i.e., fig_4 directory in the local copy of this repository)
6. `python setup.py path_to_downloaded_aws_bucket_fig_4` (this does some quick and dirty path fixing, and also copies the downloaded large files (required for Figure 4 analysis) into the expected locations)
7. `code/generate_donor_features_and_estimate_karyotype.ipynb`
8. `code/generate_combined_illumina_ultima_feature_data.ipynb`
9. `code/mds_xgboost_classification.ipynb`
10. `code/generate_blood_aging_paper_figs_and_tables.ipynb`

