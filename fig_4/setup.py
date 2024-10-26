import os
import shutil
import sys
import pathlib

curr_dir_path = os.getcwd()
assert curr_dir_path.endswith('/fig_4'), f'curr_dir_path: {curr_dir_path}'

print(f'quick and dirty fixing of a path in each ipynb')
for ipynb_path in [
    'code/generate_donor_features_and_estimate_karyotype.ipynb',
    'code/generate_combined_illumina_ultima_feature_data.ipynb',
    'code/mds_xgboost_classification.ipynb',
    'code/generate_blood_aging_paper_figs_and_tables.ipynb',
]:
    with open(ipynb_path) as f:
        ipynb_contents = f.read()
    fixed_ipynb_contents = ipynb_contents.replace("os.chdir('dummy_to_replace_with_blood_aging_dir_path/code')", f"os.chdir('{curr_dir_path}/code')")

    with open(ipynb_path, 'w') as f:
        f.write(fixed_ipynb_contents)


print(f'quick and dirty copying of files from local copy of aws bucket into where the code expects these files to be (sorry for this ugliness!)')
aws_bucket_fig4_dir_path = sys.argv[1]

for source_dir_path, dest_dir_path in [
    (os.path.join(aws_bucket_fig4_dir_path, 'N211_ad_files', 'final_only_N211_only_demux_16_01_22_1'), 'output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_only_N211_only_demux_16_01_22_1'),
    (os.path.join(aws_bucket_fig4_dir_path, 'N211_ad_files', 'final_only_N211_only_demux_26_03_23_1'), 'output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_only_N211_only_demux_26_03_23_1'),
    (os.path.join(aws_bucket_fig4_dir_path, 'N180_ad_files', 'final_only_N200_only_demux_26_03_23_1'), 'output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_only_N200_only_demux_26_03_23_1'),
    (os.path.join(aws_bucket_fig4_dir_path, 'N180_ad_files', 'final_only_N200_only_demux_11_04_21_2'), 'output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_only_N200_only_demux_11_04_21_2'),
    (os.path.join(aws_bucket_fig4_dir_path, 'N180_ad_files', 'final_only_N200_only_demux_10_01_22_1'), 'output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_only_N200_only_demux_10_01_22_1'),
]:
    print(f'start copying {source_dir_path}')
    dest_containing_dir_path = os.path.dirname(dest_dir_path)
    pathlib.Path(dest_containing_dir_path).mkdir(parents=True, exist_ok=True)
    shutil.copytree(source_dir_path, dest_dir_path)

for source_path, dest_path in [
    (os.path.join(aws_bucket_fig4_dir_path, 'ult_mds_cells.h5ad'), 'output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/cells_with_metacell_attrs.h5ad'),
    (os.path.join(aws_bucket_fig4_dir_path, 'ult_mds_metacells.h5ad'), 'output_and_given_intermediate_output/240623_pb_ult_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/metacells_with_projection.h5ad'),
    (os.path.join(aws_bucket_fig4_dir_path, 'illu_mds_cells.h5ad'), 'output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/cells_with_metacell_attrs.h5ad'),
    (os.path.join(aws_bucket_fig4_dir_path, 'illu_mds_metacells.h5ad'), 'output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_mds_cyto_normal_excluding_atlas/metacells_with_projection.h5ad'),
    (os.path.join(aws_bucket_fig4_dir_path, 'illu_ref_cells.h5ad'), 'output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_c.h5ad'),
    (os.path.join(aws_bucket_fig4_dir_path, 'illu_ref_metacells.h5ad'), 'output_and_given_intermediate_output/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_mc.h5ad'),
    (os.path.join(aws_bucket_fig4_dir_path, 'genes.gtf'), 'input/human_genome_files/genes.gtf'),
    (os.path.join(aws_bucket_fig4_dir_path, 'genome.fa'), 'input/human_genome_files/genome.fa'),
]:
    print(f'start copying {source_path}')
    dest_dir_path = os.path.dirname(dest_path)
    pathlib.Path(dest_dir_path).mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_path, dest_path)

