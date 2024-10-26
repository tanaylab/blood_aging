import numpy as np
import pandas as pd
import itertools
import re
import os
import collections
import datetime
import pickle
import random

from generic import generic_utils
from mds import mds_in_out_dir_paths

# https://stackoverflow.com/questions/16424493/pandas-setting-no-of-max-rows/16433953#16433953
# pd.describe_option('display')
pd.set_option('display.large_repr', 'truncate')
pd.set_option('display.max_rows', 40)
pd.set_option('display.min_rows', 40)


CELL_RANGER_OUT_MTX_RELATIVE_DIR_PATH = 'outs/filtered_feature_bc_matrix'
CELL_RANGER_UMI_COUNT_OF_BARCODES_DF_CSV_RELATIVE_FILE_PATH = os.path.join(CELL_RANGER_OUT_MTX_RELATIVE_DIR_PATH, 'cell_ranger_umi_count_of_barcodes_df.csv')
CELL_RANGER_OUT_SUMMARY_RELATIVE_FILE_PATH = 'outs/web_summary.html'
CELL_RANGER_BAM_RELATIVE_FILE_PATH = 'outs/possorted_genome_bam.bam'

DATE_PART_IN_EXP_NAME_REGEX = r'[0-9]{2}_[0-9]{2}_[0-9]{2}'
DONOR_ID_REGEX = r'(N|NS|T|G)[0-9]{1,3}'

SPECIFIC_EXPS_TO_EXCLUDE_NAMES_EXCEPT_FOR_NIMROD = [
    
    'demux_08_02_23_1', # contaminated by at least two other wells, and only old healthy.
    # 'N276_bm_05_02_23_1', # contaminated by at least two other wells, but we don't have another option for it.
    'demux_03_07_22_1_ultima', # much more than expected common barcodes with multiple experiments, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). tech_reps_common_barcode_analysis is inconclusive...
    # 'demux_03_07_22_2_ultima', # much more than expected common barcodes with multiple experiments, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). according to tech_reps_common_barcode_analysis, seems like demux_03_07_22_2_ultima only contaminated others, and wasn't contaminated by others, and so we can use it, i think.
    'demux_03_07_22_1', # giving up on it because it has some horrible droplets (clog?), and a better tech replicate.
    # 'demux_03_07_22_2', # high fraction of unassigned droplets suggests pretty high ambient noise. 2 healthy, one mastocytosis, one PV.
    # 'N334_bm_15_01_23_1', # much more than expected common barcodes with multiple experiments, and probably contaminated by at least one other well, but we don't have another option for it.
    'demux_09_06_22_1_ultima', # contaminated by at least one other well, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). also, according to tech_reps_common_barcode_analysis, contaminated by other wells.
    'demux_10_07_22_1_ultima', # contaminated by at least one other well, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). also, according to tech_reps_common_barcode_analysis, contaminated by other wells.
    'demux_13_11_22_2', # contaminated by at least one other well, and demux_13_11_22_1, even though has much more than expected common barcodes with multiple experiments, might be pretty clean (i.e., contaminated others but wasn't contaminated itself?), and presumably with sufficient UMIs per cell.
    # 'demux_13_11_22_1', # much more than expected common barcodes with at least one experiments, but might be pretty clean (i.e., contaminated others but wasn't contaminated itself?). we don't have another option for it. 
    'demux_08_12_22_1', # contaminated by at least two other wells, and only healthy.
    # 'N279_bm_26_10_22_1', # much more than expected common barcodes with multiple experiments, and was also sequenced by illumina (hopefully with sufficient UMIs per cell). much better UMIs per cell than N279_bm_26_10_22_1_illumina, and according to tech_reps_common_barcode_analysis, seems like N279_bm_26_10_22_1 only contaminated others, and wasn't contaminated by others, and so we can use it, i think.
    # 'N280_bm_06_11_22_1', # contaminated by at least two other wells, but we don't have another option for it. also, UMIs per cell is not great.
    # 'N251_04_12_22_1', # contaminated by at least one other well, but we don't have another option for it. clearly contaminated according to according to tech_reps_common_barcode_analysis.
    # 'demux_04_12_22_1', # much more than expected common barcodes with multiple experiments, but might be pretty clean (i.e., contaminated others but wasn't contaminated itself?), and presumably with sufficient UMIs per cell.
    # 'N251_bm_04_12_22_1', # contaminated by at least two other wells, but we don't have another option for it. also, UMIs per cell is not great. also, according to tech_reps_common_barcode_analysis, contaminated by other wells. # TODO: uncomment after we resequence.
    # 'demux_15_01_23_1', # much more than expected common barcodes with multiple experiments, but might be pretty clean (i.e., contaminated others but wasn't contaminated itself?), and presumably with sufficient UMIs per cell.
    'demux_02_05_22_1_ultima', # contaminated by at least one other well, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). also, according to tech_reps_common_barcode_analysis, contaminated by other wells.
    # 'demux_02_05_22_2_ultima', # was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). TODO: consider uncommenting because we already have the same cells in demux_02_05_22_2? 
    'demux_19_06_22_1_ultima', # much more than expected common barcodes with at least one experiment, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). tech_reps_common_barcode_analysis is inconclusive, and i suspect it is contaminated a bit, so don't use demux_19_06_22_1_ultima.
    'demux_13_06_22_1_ultima', # contaminated by at least two other wells, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). also, according to tech_reps_common_barcode_analysis, contaminated by other wells.
    'N279_26_10_22_1', # much more than expected common barcodes with at least one experiment, and was also sequenced by illumina (hopefully with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). # tech_reps_common_barcode_analysis is inconclusive, but comparing the umi count per barcode hists of N279_26_10_22_1 and N279_26_10_22_1_illumina suggests N279_26_10_22_1 was indeed contaminated (maybe by a donor from out of the plate??).
    # 'demux_27_11_22_1', # contaminated by at least one other well, but we don't have another option for it.
    # 'demux_27_02_23_1', # this is actually 28.02.23. contaminated by at least one other well. 3 healthy and 1 mastocytosis. We don't have another option for it
    # 'demux_06_06_22_1_ultima', # much more than expected common barcodes with at least one experiment, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). according to tech_reps_common_barcode_analysis, seems like demux_06_06_22_1_ultima only contaminated others, and wasn't contaminated by others, and so we can use it, i think.
    # 'demux_06_06_22_2_ultima', # much more than expected common barcodes with at least one experiment, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). according to tech_reps_common_barcode_analysis, seems like demux_06_06_22_1_ultima only contaminated others, and wasn't contaminated by others, and so we can use it, i think.
    'demux_19_06_22_2_ultima', # contaminated by at least two other wells, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). also, according to tech_reps_common_barcode_analysis, contaminated by other wells.
    'demux_10_07_22_2_ultima', # contaminated by at least one other well, and was also sequenced by illumina (presumably with sufficient UMIs per cell, or at least higher than the UMIs per cell we got from ultima). also, according to tech_reps_common_barcode_analysis, contaminated by other wells.
    # 'demux_14_11_22_1', # much more than expected common barcodes with at least one experiments, but we don't have another option for it.
    # 'demux_09_06_22_2_ultima', # was also sequenced by illumina (presumably with sufficient UMIs per cell). TODO: consider uncommenting because we already have the same cells in demux_09_06_22_2?
    # 'N12_bm_04_12_22_1', # more than expected common barcodes with multiple experiments, but we don't have another option for it. 230829: according to tech_reps_common_barcode_analysis (specifically, miniseq_230816_N311_BM_4_12_22), seems like N12_bm_04_12_22_1 only contaminated others and wasn't contaminated by others, so we can use the ultima sample, i think.
    # 'demux_06_11_22_1', # more than expected common barcodes with multiple experiments, but we don't have another option for it. according to tech_reps_common_barcode_analysis, seems like it only contaminated others, and wasn't contaminated by others, and so we can use it, i think.
    # TODO: consider discarding experiments by ultima/illumina that were sequenced by both (even when the ultima was not involved in contamination), to not have in the same model the same cell more than once.
    
    
    # NOTE: and so, the libs that we probably want to resquence (not considering that some experiments are less important due to less interesting patients - this might mean we can give up on some of the following libraries):
    # healthy:
    # demux_08_12_22_1
    # demux_27_02_23_1 (the one we are missing is actually the pooled sample from 28.02.23 - due to swapped labels)
    # demux_08_02_23_1

    
    # containing unhealthy:
    # N251_bm_04_12_22_1 (the one we are missing is actually the N12 (N310) BM sample from 04.12.22 - due to swapped labels) (unexpected donors identified here) (contains many N12 cells which we can't tell whether they are from her BM or PB)
    # demux_04_12_22_1 (contains a low number of N12 cells which we can't tell whether they are from her BM or PB) (much more than expected common barcodes with multiple experiments, but might be pretty clean)
    
    # N276_bm_05_02_23_1 (unexpected donors identified here)
    # N334_bm_15_01_23_1 (seems to contain other donors (which i couldn't identify))
    # N280_bm_06_11_22_1 (unexpected donors identified here)
    # N251_04_12_22_1 (seems to contain other donors (which i couldn't identify))
    # demux_27_11_22_1 (unexpected donors identified here)
    
    # demux_14_11_22_1 (much more than expected common barcodes with another experiment)
    # demux_06_11_22_1 (more than expected common barcodes with multiple experiments)

    # N12_bm_04_12_22_1 (the one we are missing is actually the N251 (N311) BM sample from 04.12.22 - due to swapped labels) (more than expected common barcodes with multiple experiments, but might be pretty clean)
    # demux_15_01_23_1 (much more than expected common barcodes with multiple experiments, but might be pretty clean)
    
    
    # non-contaminated or also sequenced by illumina, but maybe with insufficient coverage:
    # demux_29_01_23_1 (demux_29_01_23_1 seq saturation: 16.3%)
    # N279_bm_26_10_22_1_illumina (N279_bm_26_10_22_1_illumina seq saturation: 5.9%) (N279_bm_26_10_22_1 has much more than expected common barcodes with multiple experiments)
    # N279_26_10_22_1_illumina (N279_26_10_22_1_illumina seq saturation: 45.5%) (N279_26_10_22_1 much more than expected common barcodes with multiple experiments)
    # N257_bm_06_10_22_1 (N257_bm_06_10_22_1 seq saturation: 8.8%) (N257_bm_06_10_22_1_illumina has lower median UMI count)
    # demux_09_02_22_1 (demux_09_02_22_1 seq saturation: 15.6%)


    # weirdly, it seems like the following two illumina libraries somehow contaminated each other a bit. see tech_reps_common_barcode_analysis, which seems to me to suggest that the contamination from demux_23_11_20_1_illumina to demux_30_11_20_1_illumina was larger than the contamination of the opposite direction (which wasn't negligible).
    'demux_23_11_20_1_illumina', # also sequenced by ultima
    'demux_30_11_20_1_illumina', # also sequenced by ultima. # NOTE: 240523: my mistake, i think. seems like it really should have been demux_30_11_20_2_illumina, but it's too late now and we will keep it this way.

    # problematic.
    # 'demux_28_11_21_1',
    # 'demux_07_02_22_1',
]

SC_RNA_SEQ_PREPROCESSING_PARAMS = {
    'constant_random_seed_for_reproducibility': 1234,

    'format1_mip_genotype_csv_file_paths': [
        # '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/july22/total_varscan_July_2022.tsv',
        # '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/july22_deep/total_varscan.tsv',
        # '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/november22/total_varscan.tsv',
        # '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/january23/total_varscan.tsv',
        
        # IMPORTANT NOTE: After filtering to keep only single base substitutions, I fixed the data from the freq files such that I replace the depth column with the depth of the major+minor in this row, and then recalculated the VAF, and then for each exp,donor,SNP, i keep only the row with the highest depth.
        # '/net/mraid14/export/data/tgdata/users/orenmil/PBMC_CD34/genotyping_data/freq_files/230301_combined_freq_files.csv'


        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/230301_combined_freq_files.csv',
        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/230607_combined_freq_files.csv',
        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/230705_combined_freq_files.csv',
        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/230706_combined_freq_files.csv',
        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/230917_combined_freq_files.csv',
    ],
    'format2_mip_genotype_csv_file_paths': [
        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/240311_combined_freq_files.csv',
        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/240521_combined_freq_files.csv',
        '/dummy/dummy/dummy/raid/PBMC_CD34/genotyping_data/freq_files/240603_combined_freq_files.csv',
    ],

    

    'blood_sample_id_df_csv_file_path': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'blood_sample_id_df.csv'),
    'scrna_exp_id_df_csv_file_path': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'scrna_exp_id_df.csv'),
    'cytopenia_follow_up_xlsx_file_path': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/241009_cytopenia_follow_up_filtered_fixed.xlsx'),
    'experiments_that_orenmil_processed_metadata_df_csv_file_path': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiments_that_orenmil_processed_metadata.csv'),
    # 'exp_name_and_donor_ids_agg_csv_file_path': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/230521_combined_adjusted_experiment_donors_agg_df.csv'),


    'default_extra_num_of_donors_for_vireo_and_soup': 2,

    'bleeding_dates_with_only_targeted_dna_sequencing_data': {
        '25.10.22', # MDS_data_28-May_2023_open.xlsx says "ARCH only".
        '06.02.23', 
        '28.11.22', 
    },

    'umi_count_per_barcode_exp_types_to_ignore_when_analyzing_composition': {
        'give_up_on',
        'gave_up_but_could_truncate_and_not_bias_composition',
        'conditional_on_more_data_can_truncate_and_not_bias_composition',
        'after_force_cells_borderline_can_truncate_and_not_bias_composition',
        'cant_avoid_biased_composition',
    },
    'umi_count_per_barcode_exp_types_that_can_be_used_when_analyzing_composition': {
        'can_truncate_and_not_bias_composition',
        'probably_can_truncate_and_not_bias_composition',
        'can_truncate_and_presumably_slightly_bias_composition', # NOTE: not sure, but i think better to have it here.
    },
    
    # NOTE: ugh. the name is not completely accurate for all. sometimes the lower threshold is also to minimize FPs, if those are horrible such as in N279_26_10_22_1_illumina.
    'type_to_exp_name_to_forced_num_of_filtered_barcodes_and_min_umi_count': { # fn and fp - false negatives and false positives, respectively
        'give_up_on': {
            # ok that we give up on them
            'N328_15_12_22_1': (None, None), # giving up on them because if irrelevancy to this project.
            'NS20_29_09_21_1': (None, None), # giving up on them because if irrelevancy to this project.
            'T6_25_05_22_1': (None, None), # giving up on them because if irrelevancy to this project.
            'demux_13_11_22_2': (None, None), # (29409, None), # really horrible droplet umi count hist. no idea where to put a threshold. also, giving up on it because it is contaminated. (but that's ok, as we have demux_13_11_22_1, which the miniseq would hopefully prove to be clean) 
            'N279_26_10_22_1': (None, None), # giving up on it because it is probably contaminated. (but that's ok, as we have a better tech rep)
            'demux_03_07_22_1_ultima': (None, None), # giving up on it because it is probably contaminated. (but that's ok, as we have a better tech rep)
            'demux_03_07_22_1': (None, None), # giving up on it because it has some horrible droplets (clog?), and a better tech replicate.
            'demux_02_05_22_1_ultima': (None, None), # giving up on it because it is contaminated. (but that's ok, as we have a better tech rep)
            'demux_08_12_22_1': (None, None), # giving up on it because it is contaminated. only healthy.
            'demux_09_06_22_1_ultima': (None, None), # giving up on it because it is contaminated. (but that's ok, as we have a better tech rep)
            'demux_13_06_22_1_ultima': (None, None), # giving up on it because it is contaminated. (but that's ok, as we have a better tech rep)
        },
        'gave_up_but_could_truncate_and_not_bias_composition': {
            'demux_08_02_23_1': (None, 2**9.6), # giving up on it because it is contaminated. only healthy. # but it looks fine.
            'demux_10_07_22_1_ultima': (None, 2**10.6), # giving up on it because it is contaminated. (but that's ok, as we have a better tech rep)
            'demux_10_07_22_2_ultima': (None, 2**11), # giving up on it because it is contaminated. (but that's ok, as we have a better tech rep)
            'demux_19_06_22_1_ultima': (None, 2**11), # giving up on it because it might be contaminated, (but that's ok, as we have a better tech rep)
            'demux_19_06_22_2_ultima': (None, 2**10.5), # giving up on it because it is contaminated. (but that's ok, as we have a better tech rep)
            'demux_23_11_20_1_illumina': (None, 2**11), # giving up on it because it is probably contaminated. (but that's ok, as we have a better tech rep)
            'demux_30_11_20_1_illumina': (None, 2**10.5), # giving up on it because it is probably contaminated. (but that's ok, as we have a better tech rep)
        },

        

        'can_truncate_and_not_bias_composition': {
            # seemingly ok experiments, but with cell ranger automatically choosing a too low threshold
            'demux_01_01_23_1': (None, 2**10),
            'demux_03_11_21_1': (12557, 2**10.1),
            'demux_03_03_22_1': (None, 2**10),
            'demux_06_06_22_1': (None, 2**11), # the noise level seems somewhat high, though.
            'demux_02_05_22_1': (None, 2**9.7), # cellranger "main" threshold is above, but with some cells below. anyway it is borderline enough that we don't have to change this, i think.
            'demux_06_06_22_1_ultima': (None, 2**11), # the noise level seems somewhat high, though.
            'demux_06_06_22_2_ultima': (None, 2**10),
            'demux_06_06_22_2': (None, 2**10),
            'demux_09_01_22_1': (None, 2**9.5), # the noise level seems somewhat high, though.
            'demux_09_06_22_1': (None, 2**9.6), # though it has a weird peak in the middle...
            'demux_09_06_22_2': (None, 2**9.6),
            'demux_09_06_22_2_ultima': (None, 2**10),
            'demux_10_01_22_1': (None, 2**10),
            'demux_10_07_22_1': (None, 2**10.2),
            'demux_10_07_22_2': (None, 2**10.2),
            'demux_12_06_22_1': (None, 2**10.7),
            'demux_12_06_22_1_ultima': (None, 2**10.5),
            'demux_13_11_22_1': (None, 2**11),
            'demux_16_01_22_1': (None, 2**9.5), # the noise level seems somewhat high, though.
            'demux_17_08_20_1': (None, 2**11),
            'demux_19_06_22_1': (None, 2**10.5),
            'demux_19_06_22_2': (None, 2**10.5),
            'demux_22_06_23_2': (None, 2**10),
            'demux_23_11_20_1': (None, 2**10.5),
            'demux_30_11_20_1': (None, 2**10),
            'demux_26_03_23_1': (None, 2**10),
            'demux_28_02_23_1': (None, 2**9.8), # cellranger "main" threshold is above, but with some cells below. anyway it is borderline enough that we don't have to change this, i think.
            'N188_bm_15_01_23_1': (None, 2**10),
            'N242_bm_12_03_23_1': (None, 2**10),
            'demux_14_02_22_1': (None, 2**9.2),
            'demux_02_03_22_1': (None, 2**10),
            'G8_13_07_23_1': (None, 2**9.8),
            'demux_g_03_07_23_1': (None, 2**9.5),
            'demux_22_02_21_1': (None, 2**9.3),
            'demux_27_11_22_1': (33532, 2**9.5), 
            'demux_n_13_08_23_1': (None, 2**10), 
            'demux_n_15_08_23_2': (None, 2**9.8), 
            'demux_n_15_08_23_3': (None, 2**10.5), 
            'demux_07_11_21_1': (None, 2**9.5), 
            'demux_07_02_22_1': (None, 2**9.2), 
            'demux_04_12_22_1': (None, 2**9.1), 
            'demux_02_01_22_1': (None, 2**10), 
            'demux_n_14_02_24_1': (None, 2**11), 
            'demux_n_01_01_24_1': (None, 2**12), 
            'demux_g_22_02_24_1': (None, 2**10), 
            'demux_g_11_01_24_1': (None, 2**9.8), 
            'demux_g_08_02_24_1': (None, 2**10.3), 
            'demux_g_08_01_24_1': (None, 2**11), 
            'demux_g_01_02_24_1': (None, 2**10.5), 
            'demux_04_01_21_2': (None, 2**9.5), 
            'demux_01_02_21_1': (None, 2**9.5), 
            'demux_07_03_21_2': (None, 2**9), 
            'demux_07_12_20_1': (None, 2**10), 
            'demux_11_01_21_1': (None, 2**9.5), 
            'demux_11_01_21_2': (None, 2**9.5), 
            'demux_11_04_21_1': (None, 2**9), 
            'demux_13_12_22_2': (None, 2**9.3),
            'demux_21_01_21_1': (None, 2**10),
            'demux_21_01_21_2': (None, 2**9.5),
            'demux_21_02_21_2': (None, 2**9),
            'demux_21_12_20_1': (None, 2**9),
            'demux_21_12_20_2': (None, 2**9.3),
            'demux_28_12_20_1': (None, 2**9),
            'demux_28_12_20_2': (None, 2**9),

       
            'demux_27_02_23_1': (None, 2**10), 
            'N334_bm_15_01_23_1': (None, 2**9), # TODO: this might be really bad, if the peak we see is a clog. it would mean

            # maybe ok, though probably with fps we can't avoid
            'demux_01_02_21_2': (None, 2**10.2), # not sure there are really many fps here, because it seems like souporcell and vireo get rid of most fps.
            'demux_01_02_21_2_illumina': (None, 2**8.7), # demux_01_02_21_2 seems much better, but when comparing compositions of some of the tech reps of donors here, they looked pretty much identical in demux_01_02_21_2 and demux_01_02_21_2_illumina, so i guess the composition actually isn't biased in demux_01_02_21_2_illumina.
            'demux_19_12_21_1': (16728, 2**10.8), # not sure there are really many fps here, because it seems like souporcell and vireo get rid of most fps.
            'demux_14_11_22_1': (74238, 2**9.1), # not sure there are really many fps here, because it seems like souporcell and vireo get rid of most fps.
            'demux_11_04_21_2': (None, 2**9.6), # not sure there are really many fps here, because it seems like souporcell and vireo get rid of most fps.
            'N48_bm_12_03_23_1': (47350, 2**9.5),
            'N365_bm_12_03_23_1': (46165, 2**9.3),
            'demux_15_01_23_1': (21185, 2**9),
            'demux_n_30_11_23_1': (None, 2**9.5),
        },
        
        'probably_can_truncate_and_not_bias_composition': {
            'demux_13_06_22_1': (None, 2**9.2),
            'demux_28_02_22_1': (None, 2**9), # cellranger "main" threshold is above, but with some cells below. anyway it is borderline enough that we don't have to change this, i think.
            'demux_22_02_21_2': (None, 2**9.7), # a lot of fps, so i can't really know where to put the threshold. hopefully discarding unassigned/doublets would be sufficient.
            'demux_28_11_21_1': (None, 2**9),
            'demux_30_01_22_1': (None, 2**9.1), # the noise level seems somewhat high, though.
            'demux_09_02_22_1': (None, 2**9), # the noise level seems somewhat high, though. borderline.
            'demux_06_11_22_1': (None, 2**9.1), 
            'N12_bm_04_12_22_1': (64709, 2**9), # biased, but hopefully only very slightly.
            'demux_21_02_21_1': (None, 2**9), # maybe biased, but hopefully only very slightly.
            'demux_13_12_22_1': (None, 2**10), # maybe biased, but hopefully only very slightly. could force-cells to maybe make it slightly better, but doesn't sound like it is worth it. also, all healthy.

            
            'N279_bm_26_10_22_1': (73783, 2**9), # the weird peak (around 2**9.9) is not enriched for genotype_doublet_enriched_mcs.
            'N196_bm_01_01_23_1': (63849, 2**9), # NOTE: ugh. maybe cant_avoid_biased_composition.
            'demux_22_06_23_1': (66181, 2**9.1), # this is ok as we discard 3,donor3 (N365)
            'demux_n_bm_15_08_23_1': (None, 2**9), # NOTE: ugh. maybe cant_avoid_biased_composition.
            'demux_n_bm_19_09_23_1': (None, 2**9.5), 
            'demux_n_19_09_23_1': (None, 2**9), 
            'demux_n_04_12_23_2': (None, 2**9.5), 
            'demux_n_18_12_23_1': (None, 2**9), 
            'demux_n_05_12_23_1': (None, 2**10), 
            'demux_n_04_12_23_1': (None, 2**9), # NOTE: ugh. i really don't like the umi count hist for N191.
            'demux_g_28_09_23_1': (None, 2**9.5),
            'demux_g_07_12_23_1': (None, 2**9), 
            'demux_n_12_02_24_1': (None, 2**11), 
            'demux_g_24_07_23_1': (None, 2**9), 
            'demux_g_12_09_23_1': (None, 2**10.5), 

            'demux_07_03_21_1': (None, 2**9.5),
            
            'demux_03_07_22_2': (None, 2**11.2),
            'demux_03_07_22_2_ultima': (None, 2**12),
            'demux_01_03_21_1': (None, 2**9),
            'demux_03_03_22_2': (None, 2**10), # a bit borderline, but good enough, i think
            'demux_07_12_20_2': (None, 2**9),
            'demux_14_02_22_2': (None, 2**9.5),
            'demux_23_11_22_1': (None, 2**10),
            'demux_27_02_23_2': (None, 2**10), # confusingly, the better tech rep is demux_28_02_23_1 (until (and if) we resequence)
            'demux_g_28_03_24_1': (None, 2**9.5),
        },

        'can_truncate_and_presumably_slightly_bias_composition': {
            'N251_04_12_22_1': (None, 2**9),
            'N251_bm_04_12_22_1': (None, 2**9),
            'demux_n_24_07_23_1': (None, 2**8.2),
            'demux_g_26_06_23_1': (None, 2**8),
            'demux_05_02_23_1': (None, 2**8.3),
            'demux_29_01_23_1': (None, 2**7.5),
            'demux_n_bm_24_07_23_1': (40761, 2**9),
            'N276_bm_05_02_23_1': (108298, 2**9),
            'demux_g_28_12_23_1': (None, 2**9.7), # the bias for the donor with very few cells might not be so slight, but their umi count distribution is suspicious...
            'demux_04_01_21_1': (None, 2**8.8),
            'demux_07_03_21_1_illumina': (None, 2**8.9),
            'demux_11_04_21_2_illumina': (None, 2**8.9),
            'demux_22_02_21_2_illumina': (None, 2**8.9),
            'demux_28_02_22_2': (None, 2**9),

            'demux_02_05_22_2': (None, 2**9), # demux_02_05_22_1 is presumably better
            'demux_02_05_22_2_ultima': (None, 2**9), # demux_02_05_22_1 is presumably better
            'N279_26_10_22_1_illumina': (None, 2**12.8), # ugh. the annoying peak is enriched for genotype_doublet_enriched_mcs. ugh. so i should treat it as made of FPs. 240629: comparing this one to demux_n_30_11_23_1, the HSPC composition seems pretty much the same. so i guess the composition is at most slightly biased (sounds unlikely that the difference between bio reps was exactly cancelled by the bias).
        },

        'conditional_on_more_data_can_truncate_and_not_bias_composition': {
            # conditional on more MIP genotyping data

        },

        'after_force_cells_borderline_can_truncate_and_not_bias_composition': {
        },

        'cant_avoid_biased_composition': {
            'demux_n_20_07_23_1': (None, 2**9),
            'demux_n_20_07_23_2': (None, 2**9),

            'N280_bm_06_11_22_1': (None, 2**9),
            'N257_bm_06_10_22_1': (None, 2**9), # too low... also a lot of fps, so i can't really know where to put the threshold.
            'demux_g_26_09_23_1': (None, 2**9), 
            'demux_n_18_12_23_2': (None, 2**9),
            'demux_n_bm_12_02_24_1': (None, 2**9),
            
            # a better tech rep exists
            'N257_bm_06_10_22_1_illumina': (None, 2**8), # N257_bm_06_10_22_1 is slightly better, but in both can't avoid biased composition, i guess.
            'N279_bm_26_10_22_1_illumina': (None, 2**6.7), # N279_bm_26_10_22_1 is much better
            
        },
    },
    'exp_name_to_forced_num_of_donors_for_vireo_and_soup': {
        # use this temporarily in order to not constantly change the metadata file.
        'demux_08_12_22_1': 10,
        'demux_08_02_23_1': 9,
        'N276_bm_05_02_23_1': 8,
        'N280_bm_06_11_22_1': 10, 
        'N251_bm_04_12_22_1': 10, 
        'demux_28_02_23_1': 8,

        'N251_04_12_22_1': 7,
        # 'demux_g_26_06_23_1': 5, # vireo created an empty GT_donors.vireo.vcf for 6 donors. i hoped this won't happen for 5 donors. it did happen.
        'demux_n_18_12_23_2': 3, # vireo created an empty GT_donors.vireo.vcf for 4 donors, and also for 3 donors. souporcell failed for 4 donors but succeeded for 3! yay.
        'demux_n_19_09_23_1': 3, # souporcell failed for 5 and for 4 donors but succeeded for 3! yay.
        'demux_n_04_12_23_2': 2, # souporcell failed for 4 and for 3 donors but succeeded for 2! yay.
        
        'demux_n_bm_12_02_24_1': 3, # souporcell failed for 4 donors but succeeded for 3! yay.
        'demux_04_01_21_2': 4, # souporcell failed for 6 donors. also vireo for 6 donors kind of failed to identify N220..
        'demux_01_02_21_1': 5, # souporcell failed for 6 donors. 
    },
    'exp_name_to_forced_vireo_seed': {
        'demux_n_18_12_23_2': 3117, # vireo created an empty GT_donors.vireo.vcf for the usual seed. changing to this seed didn't help.
    },
    'exp_name_to_souporcell_extra_args': {
        'demux_g_26_09_23_1': ['--ignore', 'True'],
    },
    'exp_names_to_skip_souporcell': [
        # common souporcell doublets.err upon error:
        # <num> loaded 0 counts, is this a problem?
        # thread 'main' panicked at 'called `Option::unwrap()` on a `None` value', src/main.rs:63:102
        # note: run with `RUST_BACKTRACE=1` environment variable to display a backtrace
        
        'NS20_29_09_21_1', # failed with an error.
        # for each of the following - failed with an error, i guess. (see doublets.err (it is (pretty much?) the same as for NS20_29_09_21_1))
        'N48_bm_12_03_23_1',
        'N242_bm_12_03_23_1',
        'N196_bm_01_01_23_1',
        'N188_bm_15_01_23_1',
        'N279_26_10_22_1_illumina',
        'N257_bm_06_10_22_1_illumina',
        'G8_13_07_23_1',
        'demux_n_24_07_23_1',
        'N12_bm_04_12_22_1',
        'demux_15_01_23_1', # ugh. not a single donor. ugh.
        'demux_g_24_07_23_1', # ugh. not a single donor. ugh.
        'demux_g_07_12_23_1', # ugh. not a single donor. ugh.
        'demux_n_18_12_23_1', # ugh. not a single donor. ugh.
        # 'demux_n_19_09_23_1', souporcell failed for 5 and for 4 donors but succeeded for 3! yay.
        # 'demux_n_04_12_23_2', # souporcell failed for 4 and for 3 donors but succeeded for 2! yay.
        # 'demux_n_18_12_23_2', # souporcell failed for 4 donors but succeeded for 3! yay.
        # 'demux_n_bm_12_02_24_1', # souporcell failed for 4 donors but succeeded for 3! yay.
    ],
    'exp_names_to_skip_souporcell_which_contain_multiple_donors': [
        'demux_15_01_23_1',
        'demux_g_24_07_23_1',
        'demux_g_07_12_23_1',
        'demux_n_18_12_23_1',
    ],
    
    'donor_table_paths': {
        'nili_clinical_data_xlsx': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/240722_MDS_clinical_data.xlsx'),
        'nili_minimal_clinical_df_csv': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/240804_minimal_clinical_df.csv'),
        'arch4_donor_id_and_exp_date_df_csv': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/240802_arch4_donor_id_and_exp_date_df.csv'),
        'minimal_cbc_df_csv': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/240802_minimal_cbc_df.csv'),
        'n1_donor_id_df_csv': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/240802_n1_donor_id_df.csv'),
        # 'nili_clinical_data_xlsx': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/old/240707_MDS_clinical_data.xlsx'),
        'nili_treatment_xlsx': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/240719_Final donors with treatment status.xlsx'),
        'gal_dadi_clinical_data_xlsx': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/240703_MPN_clinical_data.xlsx'),
        'technical_replicates_xlsx': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/230522_Tech_Reps_Clinical.xlsx'),
        'prev_donor_id_xlsx': os.path.join(mds_in_out_dir_paths.EXPERIMENT_METADATA_DIR_PATH, 'experiment_donor_tables/230618_10x_all_previous_instances.xlsx'),
    },

    'exps_to_skip_names': [
        'T6_25_05_22_1', # ICF, 2 years old
        'B_cells_demux_08_02_23_1', # B cells
        'NS20_29_09_21_1', # almost only T cells
        'N328_15_12_22_1', # JMML, 2 years old

        
        # 'demux_11_04_21_1', # partially processed in order to compare to demux_11_04_21_1_elembio
        # 'demux_07_03_21_2', # partially processed in order to compare to demux_07_03_21_2_elembio

        'demux_11_04_21_1_elembio',
        'demux_07_03_21_2_elembio',

        *SPECIFIC_EXPS_TO_EXCLUDE_NAMES_EXCEPT_FOR_NIMROD,
    ],

    'exp_name_to_mouse_human_doublet_barcode_file_path': {
        'demux_g_26_09_23_1': '/dummy/dummy/dummy/raid/mds/from_dror/dror_cellranger_out_demux_g_26_09_23_1_orenmil_out/not_only_human_barcodes.txt',
    },
    'exp_name_to_barcode_whitelist_file_paths': {
        'N251_bm_04_12_22_1': [
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230816_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230816_N310_BM_4_12_22/outs/filtered_feature_bc_matrix/barcodes.tsv', # better without it due to potential false positives.
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_041222_N310BM/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'N251_04_12_22_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_041222_N311PB/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'N276_bm_05_02_23_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_050223_N340_BM/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'N280_bm_06_11_22_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_061122_N280BM/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'demux_13_11_22_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_131122_1/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'demux_15_01_23_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_150123_pool/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'N334_bm_15_01_23_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_271122_150123_N334_BM/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'demux_27_11_22_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_271122_150123_N334_BM/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
        'demux_27_02_23_1': [
            '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/misc/230907_miniseq_quick_dirty_analysis/cell_ranger_out/miniseq_230907_280223/outs/filtered_feature_bc_matrix/barcodes.tsv',
        ],
    },

    'exp_name_to_other_exp_date_part_presumably_accidentally_mixed_and_donor_ids': {
        # IMPORTANT NOTE: for 10 (some of the experiments that were also sequenced by illumina, and all of those that were only sequenced by ultima.) of the following experiments with concatenated fastq files, looked at single fastq files, and all of them were contaminated, so it isn't that the problem was in a specific ultima run, and so sounds more probable that the wells themselves were contaminated, and not something downstream...

        'demux_02_05_22_1_ultima': {('19_06_22', ('N265','N266','N267','N268'))}, # Nili_02.05.22_1 # NOTE: concatenated fastq file (022619, 020138). no healthy in orig exp.
        'demux_10_07_22_1_ultima': {('09_06_22', ('N250','N251','N252','N255'))}, # Nili_10.07.22_1 # NOTE: concatenated fastq file (022679, 020128). 2/4 healthy in orig exp.
        'demux_10_07_22_2_ultima': {('14_11_22', ('N286','N287','N288','N289','N290'))}, # Nili_10.07.22_2 # NOTE: concatenated fastq file (022679, 020128). 2/4 healthy in orig exp.
        'demux_09_06_22_1_ultima': {('10_07_22', ('N273','N274','N275','N276'))}, # Nili_09.06.22_1 # NOTE: concatenated fastq file (022619, 020138). 1/4 healthy in orig exp.
        'demux_19_06_22_2_ultima': {('06_06_22', ('N249','N253','N254','N257')), ('28_02_22', ('N235',))}, # Nili_19.06.22_2 # NOTE: concatenated fastq file (022679, 020128). 3/4 healthy in orig exp.
        

        'demux_13_06_22_1_ultima': {('02_05_22', ('N245', 'N246', 'N247', 'N248')), ('29_09_21', ('NS20',))}, # TODO: should verify this when we have all MIP genotyping, i.e., including NS20. # NOTE: concatenated fastq file (022679, 020128). 2/4 healthy in orig exp.
        'demux_08_02_23_1': {('15_12_22', ('N328',)), ('05_02_23', ('N276',)), ('Elia_exp', ())}, # Nili_15.12.22 (020149)
        'N251_04_12_22_1': {('unknown', ())}, # 230709: given current available genotyping data, seems like a contamination of N12??? but why is N234 there?? unclear. not important enough, i guess. # NOTE: concatenated fastq file (022684, 020169)
        'demux_08_12_22_1': {('13_11_22', ('N282', 'N283', 'N284', 'N285')), ('03_07_22', ('N269', 'N270', 'N271', 'N272'))}, # thought that maybe we also have N12 here, but I think it is more probable it really is N282. # NOTE: concatenated fastq file (022684, 020169). 5/5 healthy in orig exp.
        'demux_27_11_22_1': {('26_10_22', ('N279',))}, # NOTE: concatenated fastq file (022684, 020169). 1/5 healthy in orig exp.
        'demux_27_02_23_1': {('08_02_23', ('N350','N351','N352','N353'))}, # contaminated by B_cells_demux_08_02_23_1 (according to common barcode analysis).
        'demux_13_11_22_2': {('04_12_22', ('N251', 'N12')), ('03_07_22', ('N271',))}, # NOTE: 230709: N251 and N12 are not here anymore. i think they disappeared after i rerun cell-ranger while forcing the number of cells to be 29409 (as previously it automatically chose an extremely low min_umi_count threshold). (020152)
        'N276_bm_05_02_23_1': {('Elia_exp', ('N307', 'N264')), ('15_12_22', ('N328',)), ('08_02_23', ('N350','N351','N352','N353'))}, # (020149)
        'N280_bm_06_11_22_1': {('04_12_22', ('N12',)), ('26_10_22', ('N279',)), ('02_05_22', ('N245', 'N246', 'N247', 'N248')), ('13_06_22', ('N261', 'N262', 'N263', 'N264'))}, # (020152)
        
        'N251_bm_04_12_22_1': {('04_12_22', ('N306', 'N307', 'N308', 'N309', 'N251')), ('15_01_23', ('N334',))}, # NOTE: as it seems like N12 is the most dominant here, i guess this experiment is actually N12_bm_04_12_22_1 (wells accidentally swapped), and also contaminated. # NOTE: concatenated fastq file (022684, 020169) # 230712: really unsure about N306,N307,N308,N309 here. definitely we have N334 here. i guess also N251. but much less sure about the others...
        
        # NOTE: looking at common barcodes, seems like they haven't contaminated each other! but they were contaminated by others (or contaminated others, though this sounds less likely, as they didn't contaminate each other), which might be the cause for the pattern that made me think they contaminated each other..
        # 'demux_03_07_22_1_ultima': {('03_07_22', ())}, # Nili_03.07.22_1 # this is a guess according to weird bimodal distribution. see /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/out_figs/with_doublets_metacells/per_exp_norm_umi_count_hists and /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/out_figs/with_doublets_metacells/suspicious_demux_03_07_22_ultima. # NOTE: concatenated fastq file (022679, 020128). 2/4 healthy in orig exp.
        # 'demux_03_07_22_2_ultima': {('03_07_22', ())}, # Nili_03.07.22_2 # dito demux_03_07_22_1_ultima # NOTE: concatenated fastq file (022679, 020128). 2/4 healthy in orig exp.
    },

    'list_of_donor_id_and_bleeding_date_and_numbered_donor_id': [
        # NOTE: please make sure the numbered_donor_id here does not contain '__' (because we add __ in get_df_with_donor_id_replaced_by_numbered_donor_id, and later assume it is the only occurrence of '__')
        ('N276', '10.07.22', 'N276_1'),
        ('N276', '05.02.23', 'N276_2'),
        ('N325', '01.01.23', 'N325_1'),
        ('N325', '15.01.23', 'N325_2'),
        ('N251', '09.06.22', 'N251_1'),
        ('N251', '04.12.22', 'N251_2'),
        ('N260', '12.06.22', 'N260_1'),
        ('N260', '27.02.23', 'N260_2'),
        ('N200', '11.04.21', 'N200_1'),
        ('N200', '10.01.22', 'N200_2'),
        ('N200', '26.03.23', 'N200_3'),
        ('N198', '02.01.22', 'N198_1'),
        ('N198', '05.02.23', 'N198_2'),
        ('N196', '19.12.21', 'N196_1'),
        ('N196', '01.01.23', 'N196_2'),
        ('N199', '17.08.20', 'N199_1'),
        ('N199', '02.01.22', 'N199_2'),
        ('N235', '17.08.20', 'N235_1'),
        ('N235', '28.02.22', 'N235_2'),
        ('N215', '22.02.21', 'N215_1'),
        ('N215', '30.01.22', 'N215_2'),
        ('N211', '16.01.22', 'N211_1'), # NOTE: 230625: fixed swapping between N211_1 and N211_2 now.
        ('N211', '26.03.23', 'N211_2'),
        ('N307', '04.12.22', 'N307_1'),
        ('N307', '28.02.23', 'N307_2'),
        ('N242', '03.03.22', 'N242_1'),
        ('N242', '12.03.23', 'N242_2'),
        # NOTE: please make sure the numbered_donor_id here does not contain '__' (because we add __ in get_df_with_donor_id_replaced_by_numbered_donor_id, and later assume it is the only occurrence of '__')
        ('N208', '10.01.22', 'N208_1'),
        ('N208', '26.03.23', 'N208_2'),
        ('N240', '02.03.22', 'N240_1'),
        ('N240', '01.01.23', 'N240_2'),
        ('N188', '07.11.21', 'N188_1'),
        ('N188', '15.01.23', 'N188_2'),
        ('N12', '03.11.21', 'N12_1'),
        ('N12', '04.12.22', 'N12_2'),
        ('N48', '22.02.21', 'N48_1'), 
        ('N48', '03.11.21', 'N48_2'),
        ('N48', '12.03.23', 'N48_3'),
        ('N217', '30.11.20', 'N217_1'),
        ('N217', '30.01.22', 'N217_2'),
        ('N224', '30.11.20', 'N224_1'),
        ('N224', '07.02.22', 'N224_2'),
        ('N232', '01.02.21', 'N232_1'),
        # NOTE: please make sure the numbered_donor_id here does not contain '__' (because we add __ in get_df_with_donor_id_replaced_by_numbered_donor_id, and later assume it is the only occurrence of '__')
        ('N232', '14.02.22', 'N232_2'),
        ('N233', '22.02.21', 'N233_1'),
        ('N233', '14.02.22', 'N233_2'),
        ('N210', '11.04.21', 'N210_1'),
        ('N210', '10.01.22', 'N210_2'),
        ('N1', '01.01.23', 'N1_1'),
        ('N1', '14.08.23', 'N1_2_24h'),
        ('N1', '15.08.23', 'N1_3'),
        ('N264', '22.06.23', 'N264_1_24h'),
        ('N264', '14.08.23', 'N264_2_24h'),
        ('N264', '15.08.23', 'N264_3'),
        ('N275', '10.07.22', 'N275_1'),
        ('N275', '22.06.23', 'N275_2_36h'),
        ('N285', '13.11.22', 'N285_1'),
        # NOTE: please make sure the numbered_donor_id here does not contain '__' (because we add __ in get_df_with_donor_id_replaced_by_numbered_donor_id, and later assume it is the only occurrence of '__')
        ('N285', '22.06.23', 'N285_2_48h'),
        ('N257', '06.06.22', 'N257_1'),
        ('N257', '06.10.22', 'N257_2'),
        ('N365', '12.03.23', 'N365_1'),
        ('N365', '22.06.23', 'N365_2'),
        # NOTE: please make sure the numbered_donor_id here does not contain '__' (because we add __ in get_df_with_donor_id_replaced_by_numbered_donor_id, and later assume it is the only occurrence of '__')
        ('G8', '13.07.23', 'G8_1'),
        ('G8', '24.07.23', 'G8_2_11d_frozen'),
        ('N279', '26.10.22', 'N279_1'),
        ('N279', '30.11.23', 'N279_2'),
        ('N191', '28.11.21', 'N191_1'),
        ('N191', '04.12.23', 'N191_2'),
        ('N192', '28.11.21', 'N192_1'),
        ('N192', '05.12.23', 'N192_2'),
        ('N280', '06.11.22', 'N280_1'),
        ('N280', '15.08.23', 'N280_2'),
        ('N281', '06.11.22', 'N281_1'),
        ('N281', '30.11.23', 'N281_2'),
        ('N259', '12.06.22', 'N259_1'),
        ('N259', '07.12.23', 'N259_2'),
        
        
        ('N238', '04.01.21', 'N238_1'),
        ('N238', '02.03.22', 'N238_2'),
        ('N254', '11.01.21', 'N254_1'),
        ('N254', '06.06.22', 'N254_2'),
        ('N223', '01.02.21', 'N223_1'),
        ('N223', '07.02.22', 'N223_2'),
        ('N273', '21.02.21', 'N273_1'),
        ('N273', '10.07.22', 'N273_2'),
        ('N225', '01.03.21', 'N225_1'),
        ('N225', '07.02.22', 'N225_2'),
        ('N226', '01.03.21', 'N226_1'),
        ('N226', '07.02.22', 'N226_2'),
        ('N231', '07.03.21', 'N231_1'),
        ('N231', '14.02.22', 'N231_2'),
        ('N244', '11.04.21', 'N244_1'),
        ('N244', '03.03.22', 'N244_2'),
        ('N269', '07.12.20', 'N269_1'),
        ('N269', '03.07.22', 'N269_2'),
        ('N213', '07.12.20', 'N213_1'),
        ('N213', '16.01.22', 'N213_2'),
        ('N220', '04.01.21', 'N220_1'),
        ('N220', '02.03.22', 'N220_2'),
        ('N268', '21.12.20', 'N268_1'),
        ('N268', '19.06.22', 'N268_2'),
        ('N234', '28.12.20', 'N234_1'),
        ('N234', '14.02.22', 'N234_2'),
        ('N204', '09.01.22', 'N204_1'),
        ('N204', '18.12.23', 'N204_2'),
        
        ('N157_nimrod', '22.02.21', 'N157_nimrod_1'),
        ('N157_nimrod', '30.01.22', 'N157_nimrod_2'),
        ('N159_nimrod', '22.02.21', 'N159_nimrod_1'),
        ('N159_nimrod', '14.02.22', 'N159_nimrod_2'),
        ('N180_nimrod', '11.04.21', 'N180_nimrod_1'),
        ('N180_nimrod', '10.01.22', 'N180_nimrod_2'),
        ('N181_nimrod', '11.04.21', 'N181_nimrod_1'),
        ('N181_nimrod', '10.01.22', 'N181_nimrod_2'),
        ('N78_nimrod', '30.11.20', 'N78_nimrod_1'),
        ('N78_nimrod', '07.02.22', 'N78_nimrod_2'),
        ('N80_nimrod', '30.11.20', 'N80_nimrod_1'),
        ('N80_nimrod', '30.01.22', 'N80_nimrod_2'),
        ('N91_nimrod', '01.02.21', 'N91_nimrod_1'),
        ('N91_nimrod', '14.02.22', 'N91_nimrod_2'),
        ('N184_nimrod', '28.12.20', 'N184_nimrod_1'),
        ('N184_nimrod', '14.02.22', 'N184_nimrod_2'),
        ('N260_nimrod', '12.06.22', 'N260_nimrod_1'),
        ('N260_nimrod', '27.02.23', 'N260_nimrod_2'),
        ('N307_nimrod', '04.12.22', 'N307_nimrod_1'),
        ('N307_nimrod', '28.02.23', 'N307_nimrod_2'),
    ],
    'delayed_sample_numbered_donor_ids': ['N264_1_24h', 'N275_2_36h', 'N285_2_48h', 'N264_2_24h', 'N1_2_24h'],


    'clinical_table_donor_id_and_exp_date_to_forced_diagnosis_and_sex': {
        # NOTE: use this as a fix to the clinical table.
        ('T6', '25.05.22'): ('ICF', 'female'), # Immunodeficiency, centromeric region instability, facial anomalies syndrome (ICF)
        ('NS20', '29.09.21'): ('not_MDS', 'male'),
        ('N224', '30.11.20'): ('normal', 'male'), 
    },



    # comment, manual_comment (for searching...)
    # 'N74': 'we have a sample from 23_11_20, while 2 years later, he had acute leukemia, presumably without MDS symptoms before the acute leukemia',
    
    
    
    
    
    
    
    
    # 'N244': 'while skimming over donors with highest numbers of cells (such that we truncate to 3*median), noticed N244 there. indeed, he dominates both demux_03_03_22_1 and demux_03_03_22_2. he is also in the healthy atlas as 179 - in demux_11_04_21_1 he has 6074 cells while the other 3 donors had 8923, 1244, 1227 (so also had high concentration of CD34+ cells in PB here, if we assume the two donors with 1.2k cells had normal concentrations). somewhat weird given that he is healthy (though 68 years old).',

    
    
    
    
    
    # 'N275': 'N275 is the son of N56 and N295 (aka N55).',

    
    
    
    
    
    
    # 'N336': '240228: genotyping_by_sex_DEGs.ipynb shows that N336 is an outlier. he has more cells than you would expect with high xist_rps4y1_log_ratio. yet looking at vireo/souporcell and matching to genotype, seems like no detectable problem at that stage.',
    # 'N376': '240404: died very recently.',
    # 'N400': '240404: died very recently.',



    'matching_mip_to_vireo_data': {
        # 'rho_type': 'homozygous_only',
        # 'rho_type': 'homozygous_only_weighted',
        # 'rho_type': 'binom_pval',
        'rho_type': 'log2_ratio',
        # 'rho_type': 'raw',
        # 'rho_type': 'discretized',
        'max_aa_vaf': 0.15,
        'min_bb_vaf': 0.85,
        # 'min_vireo_depth': 5,
        'min_vireo_depth': 10,
        'min_mip_depth': 20, # 230618: this seemed to be too high for some donors, e.g., N74. so tried 15 instead. it didn't seem to help (i guess due to the lower signal to noise ratio), so went back to 20.
        # 'min_mip_depth': 5, 
        # 'min_mip_depth': 15,
        # 'min_mip_depth': 10,
        # 'min_mip_depth': 50,
        # 'min_mip_depth': 100,
        'min_num_of_snps_used': 5,
        'max_suspicious_diff_in_rho_between_best_to_second_best_expected_donor': 0.1,
        'max_suspicious_best_expected_donor_num_of_snps_used': 10,
    },
    'matching_mip_to_vireo_manual_examination_comments': {
        'high_confidence_matching': {
            'elimination, vireo donor matches all other expected MIP donors badly',
            'convincing VAF scatter', 
            'somewhat convincing VAF scatter;elimination, vireo donor matches all other MIP donors badly',
            'somewhat convincing VAF scatter;virtually complete demultiplexing agreement between souporcell and vireo, suggesting different vireo donors are indeed different donors, together with elimination - high confidence in matching other vireo donors to other MIP donors',
        },
        'low_confidence_matching': {
            'somewhat convincing VAF scatter;elimination, vireo donor matches all other expected MIP donors badly',
            'somewhat convincing VAF scatter',
        },
        'not_matching': {
            'matches badly, vireo donor matches another MIP donor better',
            'matches badly',
        },
        'unclear': {
            'unclear',
        },
    },
    'cell_ranger_web_summary_alerts': {
        # TODO: do something about each type of alert?
        'demux_03_07_22_1': [
            ('Low Fraction Reads in Cells', '55.2%'),
        ],
        'demux_03_07_22_2': [
            ('Low Fraction Reads in Cells', '62.8%'),
        ],
        'demux_08_02_23_1': [
            ('Low Fraction Reads in Cells', '59.7%'),
        ],
        'demux_03_07_22_1_ultima': [
            ('Low Fraction Reads in Cells', '55.9%'),
        ],
        'demux_03_07_22_2_ultima': [
            ('Low Fraction Reads in Cells', '66.5%'),
        ],
        'demux_10_07_22_1_ultima': [
            ('Low Fraction Reads in Cells', '68.3%'),
        ],
        'demux_13_11_22_2': [
            ('Low Fraction Reads in Cells', '30.3%'),
        ],
        'demux_08_12_22_1': [
            ('Low Fraction Reads in Cells', '22.0%'),
        ],
        'N279_26_10_22_1': [
            ('Low Fraction Reads in Cells', '64.4%'),
        ],
        'demux_n_24_07_23_1': [
            ('Low Fraction Reads in Cells', '44.6%'),
        ],
        'N251_04_12_22_1': [
            ('Low Fraction Reads in Cells', '66.8%'),
        ],
        'demux_g_28_03_24_1': [
            ('Low Fraction Reads in Cells', '60.7%'),
        ],
        'demux_n_18_12_23_1': [
            ('Low Fraction Reads Confidently Mapped To Transcriptome', '26.3%'),
        ],
        'demux_g_26_09_23_1': [
            ('Low Fraction Reads Confidently Mapped To Transcriptome', '8.6%'), 
        ],
    },
    'valid_freq_file_donor_ids_not_in_usual_format': [
        'NS20',
    ],
    'other_donor_id_to_donor_id_i_use': {
        # these are donors with more than one bleeding days.
        # by the way, all DNA is from peripheral blood (because BM samples are too small and precious)
        # i.e., prepare_and_combine_mip_df.prepare_and_combine_mip_df_csv replaces each other donor_id with the donor_id i use, which is simply the first id i got for each donor (that has more than a single donor_id).
	    'N310': 'N12',
	    'N332': 'N188',
	    'N326': 'N196',
	    'N341': 'N198',
	    'N327': 'N240',
	    'N311': 'N251',
	    'N364': 'N48',
	    'N340': 'N276',
	    'N330': 'N325',
	    'N180': 'N200',
	    'N401': 'N200',
	    'N372': 'N200',
	    'N78': 'N224',
	    'N181': 'N210',
	    'N88': 'N213',
	    'N157': 'N215',
	    'N98': 'N216',
	    'N80': 'N217',
	    'N116': 'N220',
	    'N131': 'N223',
	    'N164': 'N225',
	    'N165': 'N226',
	    'N175': 'N231',
	    'N91': 'N232',
	    'N159': 'N233',
	    'N99': 'N234',
	    'N150': 'N219',
	    'N109': 'N238',
	    'N179': 'N244',
	    'N122': 'N254',
	    'N97': 'N268',
	    'N79': 'N269',
	    'N147': 'N273',
	    'N358': 'N260',
	    'N361': 'N307',
	    'N366': 'N242',
	    'N370': 'N211',
	    'N371': 'N208',
        
        'N374': 'N365',

        
        'N38': 'N235',
        'N6': 'N315',
        'N10': 'N264',
        'T34': 'N264',
        'N31': 'N275',
        'N25': 'N323',
        'N36': 'N199',
        'N62': 'N236',
        'N55': 'N295',
        'N324': 'N1',
        'N395': 'N1',
        'N397': 'N1',
	    
        'N378': 'N285',
        'N379': 'N264',
        'N396': 'N264',
        'N398': 'N264',
        'N380': 'N275',

        'N389': 'N338',
        'NS23': 'N382',
        'N393': 'N280', 
        'N399': 'N252',
        
        'T35': 'N320',

        'N406': 'N279',
        'N424': 'N251',
        'N408': 'N367',
        'N404': 'N281',
        'N420': 'N339',
        'N401': 'N200',
        'N411': 'N192',
        'N407': 'N191',
        'N409': 'N382',
        'N412': 'N204',
        'G21': 'N259', 
        
        'N222': 'N162',
        'N221': 'N163',
        'N400': 'G13', 
        
        # 'NS2': 'N209',
        # 'NS6': 'N185',
        # 'NS12': 'N190',
        # 'NS1': 'N194',
        # 'NS14': 'N205',
        # 'NS24': 'N249',
        # 'NS8': 'N263',
    },

    '240723_donor_id_i_use_to_nimrod_donor_id': {
        'N200': 'N180',
        'N210': 'N181',
        'N213': 'N88',
        'N215': 'N157',
        'N216': 'N98',
        'N217': 'N80',
        'N219': 'N150',
        'N220': 'N116',
        'N223': 'N131',
        'N224': 'N78',
        'N225': 'N164',
        'N226': 'N165',
        'N231': 'N175',
        'N232': 'N91',
        'N233': 'N159',
        'N234': 'N99',
        'N238': 'N109',
        'N244': 'N179',
        'N254': 'N122',
        'N268': 'N97',
        'N269': 'N79',
        'N273': 'N147',
        'N275': 'N31',
        'N295': 'N55',
        'N323': 'N25',
    },
    'matching_mip_to_vireo_manual_examination': {
        # NOTE: 230306: only looking at log2 ratio when the result is the same as before (now with freqs, previously with varscan), and only looking at the VAF scatter when the result is different etc. changing to 'convincing VAF scatter' if the log2 ratios are convincing (as they are based on the VAF scatters).

        'demux_09_02_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N230', 'donor0'): 50, # but not enough data.
                ('N229', 'donor1'): 65, # but not enough data.
                ('N227', 'donor2'): 40, # but not enough data.
                ('N228', 'donor3'): 95,

                # alternatives
                ('N227', 'donor0'): np.nan, # not enough data.
                ('N227', 'donor1'): 10, # but not enough data
                
                ('N229', 'donor0'): 10, # out of 100
                ('N229', 'donor2'): -30,
                
                ('N230', 'donor1'): 10, # out of 100
                ('N230', 'donor2'): 10, # out of 100

                # hard ones:
                # (('N227', 'donor0'), ('N230', 'donor0'))
            },
            'pair_with_match_confidence_manual_comment': [
                ('N230', 'donor0', 'convincing VAF scatter'),
                ('N229', 'donor1', 'convincing VAF scatter'),
                ('N227', 'donor2', 'convincing VAF scatter'),
                ('N228', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_14_02_22_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N231', 'donor0', 'convincing VAF scatter'),
                ('N233', 'donor1', 'convincing VAF scatter'),
                ('N232', 'donor2', 'convincing VAF scatter'),
                ('N234', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_14_02_22_2': {
            # around 11k doublets here, but we already have more than 3k barcodes for each donor, so it seems like trying to rescue these suspicious potentially doublet barcodes is not a good idea.
            
            'pair_with_match_confidence_manual_comment': [
                ('N234', 'donor0', 'convincing VAF scatter'),
                ('N232', 'donor1', 'convincing VAF scatter'),
                ('N231', 'donor2', 'convincing VAF scatter'),
                ('N233', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },

        'demux_28_02_22_1': {
            # N235 is not myelofibrotic, but his cells dominated the experiment like often happens with cells of donors with myelofibrosis.
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N235', 'donor0'): 95,
                ('N235', 'donor1'): 90,
                ('N235', 'donor2'): 90,
                ('N235', 'donor3'): 95,

                # alternatives
                ('N219', 'donor0'): -10,
                ('N219', 'donor1'): 0,
                ('N219', 'donor2'): -10,
                ('N219', 'donor3'): -5,
                
                ('N236', 'donor0'): -10,
                ('N236', 'donor1'): -20,
                ('N236', 'donor2'): -15,
                ('N236', 'donor3'): -15,
                
                ('N237', 'donor0'): np.nan, # not enough data.
                ('N237', 'donor1'): np.nan, # not enough data.
                ('N237', 'donor2'): np.nan, # not enough data.
                ('N237', 'donor3'): np.nan, # not enough data.
                # ['N219', 'N235', 'N236', 'N237']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N235', 'donor0', 'convincing VAF scatter'),
                ('N235', 'donor1', 'convincing VAF scatter'),
                ('N235', 'donor2', 'convincing VAF scatter'),
                ('N235', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_28_02_22_2': {
            # 'N219', 'N235', 'N236', 'N237'
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N235', 'donor0'): 90,
                ('N235', 'donor1'): 98,
                ('N235', 'donor2'): 90,
                ('N235', 'donor3'): 95,

                ('N237', 'donor0'): np.nan, # not enough data.
                ('N237', 'donor1'): np.nan, # not enough data.
                ('N237', 'donor2'): np.nan, # not enough data.
                ('N237', 'donor3'): np.nan, # not enough data.

                # alternatives - all of the following match badly.
                # ('N219', 'donor0'): None, 
                # ('N219', 'donor1'): None, 
                # ('N219', 'donor2'): None, 
                # ('N219', 'donor3'): None, 
            
                # ('N236', 'donor0'): None, 
                # ('N236', 'donor1'): None, 
                # ('N236', 'donor2'): None, 
                # ('N236', 'donor3'): None, 
                
            },
            'pair_with_match_confidence_manual_comment': [
                ('N235', 'donor0', 'convincing VAF scatter'),
                ('N235', 'donor1', 'convincing VAF scatter'),
                ('N235', 'donor2', 'convincing VAF scatter'),
                ('N235', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        
        
        'demux_03_11_21_1': {
            
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N48', 'donor0'): 70,
                ('N60', 'donor1'): 90,
                ('N12', 'donor2'): 95,
                ('N209', 'donor3'): 55, # but maybe not enough data

                # alternatives - all alternatives for expected donors match badly.
                # ['N209', 'N12', 'N48', 'N60']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N48', 'donor0', 'convincing VAF scatter'),
                ('N60', 'donor1', 'convincing VAF scatter'),
                ('N12', 'donor2', 'convincing VAF scatter'),
                ('N209', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_07_11_21_1': {
            
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N187', 'donor0'): 90, 
                ('N188', 'donor1'): 90,
                ('N185', 'donor2'): 85,
                ('N186', 'donor3'): 99, 

                # alternatives - all alternatives for expected donors match badly.
                # ['N185', 'N186', 'N187', 'N188']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N187', 'donor0', 'convincing VAF scatter'),
                ('N188', 'donor1', 'convincing VAF scatter'),
                ('N185', 'donor2', 'convincing VAF scatter'),
                ('N186', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_28_11_21_1': {
            # relatively high (0.445) 10x sequencing saturation (https://kb.10xgenomics.com/hc/en-us/articles/115003646912).
            
            # high ribosomal protein genes. batch? or by chance three different people with high ribosomal expression? 
            
            # but what about N190? i would guess it is batchy, but could test later.

            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N193', 'donor0'): 90, 
                ('N190', 'donor1'): 90,
                ('N192', 'donor2'): 90,
                ('N191', 'donor3'): 95, 

                # alternatives - all alternatives for expected donors match badly.
                # ['N190', 'N191', 'N192', 'N193']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N193', 'donor0', 'convincing VAF scatter'),
                ('N190', 'donor1', 'convincing VAF scatter'),
                ('N192', 'donor2', 'convincing VAF scatter'),
                ('N191', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        
        
        'demux_19_12_21_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N189', 'donor0', 'convincing VAF scatter'),
                ('N197', 'donor1', 'convincing VAF scatter'),
                ('N194', 'donor3', 'convincing VAF scatter'),
                ('N196', 'donor5', 'convincing VAF scatter'),
                ('N195', 'donor6', 'convincing VAF scatter'),
                

                ('N189', 'donor2', 'unclear'),
                ('N194', 'donor2', 'unclear'),
                ('N195', 'donor2', 'unclear'),
                ('N196', 'donor2', 'unclear'),
                ('N197', 'donor2', 'unclear'),
                

                # donor4 is irrelevant. vireo assigned only 210 cells to them.
            ],
            
            # 'rematch_anyway': True,
        },
        
        'demux_02_01_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N199', 'donor0'): 60, # but not enough data
                (np.nan, 'donor1'): None, # donor1 is irrelevant. vireo assigned only 2 cells to donor1
                ('N198', 'donor2'): 75, 
                ('N201', 'donor3'): 90, 


                ('N198', 'donor3'): 15, # but not enough data

                # all other alternatives for expected donors match badly.
                # ['N198', 'N199', 'N201', 'N202']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N199', 'donor0', 'convincing VAF scatter'),
                ('N198', 'donor2', 'convincing VAF scatter'),
                ('N201', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_09_01_22_1': {
            # ['N203', 'N204', 'N205', 'N206']
            'pair_with_match_confidence_manual_comment': [
                ('N204', 'donor0', 'somewhat convincing VAF scatter;elimination, vireo donor matches all other MIP donors badly'),

                # donor1 is irrelevant. vireo assigned only 1 cell to them.
                
                ('N205', 'donor2', 'convincing VAF scatter'), 
                ('N206', 'donor3', 'convincing VAF scatter'), 
                ('N203', 'donor4', 'convincing VAF scatter'), 

                # donor5 is irrelevant. vireo assigned only 2 cells to them.
            ],

            # 'rematch_anyway': True,
        },
        
        'demux_10_01_22_1': {
            

            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N208', 'donor0'): 90, 
                ('N200', 'donor1'): 85, 
                ('N207', 'donor2'): 0, 
                ('N207', 'donor3'): 80, 
                
                ('N210', 'donor2'): 65, # but maybe not enough data

                # all other alternatives for expected donors match badly.
                # ['N200', 'N207', 'N208', 'N210']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N208', 'donor0', 'convincing VAF scatter'),
                ('N200', 'donor1', 'convincing VAF scatter'),
                ('N210', 'donor2', 'convincing VAF scatter'), 
                ('N207', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_16_01_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N212', 'donor0'): 70, # but maybe not enough data
                ('N213', 'donor1'): 0, 
                ('N213', 'donor2'): 70, # but maybe not enough data
                ('N214', 'donor3'): 95, 

                ('N213', 'donor0'): 10, # but not enough data

                ('N211', 'donor1'): 75, # but not enough data
                
                # all other alternatives for expected donors match badly.
                # ['N211', 'N212', 'N213', 'N214']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N212', 'donor0', 'convincing VAF scatter'), 
                ('N211', 'donor1', 'convincing VAF scatter'), 
                ('N213', 'donor2', 'convincing VAF scatter'), 
                ('N214', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_30_01_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                (np.nan, 'donor0'): None, # donor0 is irrelevant. vireo assigned only 30 cells to donor0
                ('N218', 'donor1'): 99, 
                ('N215', 'donor2'): 95, 
                ('N217', 'donor3'): 70, # but maybe not enough data

                # all other alternatives for expected donors match badly.
                # ['N215', 'N216', 'N217', 'N218']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N218', 'donor1', 'convincing VAF scatter'),
                ('N215', 'donor2', 'convincing VAF scatter'),
                ('N217', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_07_02_22_1': {
            # ['N223', 'N224', 'N225', 'N226']
            'pair_with_match_confidence_manual_comment': [
                ('N225', 'donor0', 'convincing VAF scatter'), 
                ('N224', 'donor1', 'convincing VAF scatter'),
                ('N223', 'donor2', 'convincing VAF scatter'),
                ('N226', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_21_02_21_1': {
            # ['N144', 'N148', 'N149', 'N151']
            'pair_with_match_confidence_manual_comment': [
                # donor1 and donor2 are irrelevant. each has less than 15 cells.
                ('N151', 'donor0', 'convincing VAF scatter'), 
                ('N144', 'donor3', 'convincing VAF scatter'),
                ('N149', 'donor4', 'convincing VAF scatter'),
                ('N148', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_02_03_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N238', 'donor0'): 85, # but maybe not enough data
                ('N240', 'donor1'): 80, # but maybe not enough data
                ('N239', 'donor2'): 80, # but maybe not enough data
                ('N220', 'donor3'): 80, # but maybe not enough data
                
                # all other alternatives for expected donors match badly.
                # ['N220', 'N238', 'N239', 'N240']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N238', 'donor0', 'convincing VAF scatter'), 
                ('N240', 'donor1', 'convincing VAF scatter'), 
                ('N239', 'donor2', 'convincing VAF scatter'), 
                ('N220', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_03_03_22_1': {
            
            'pair_with_match_confidence_manual_comment': [
                ('N242', 'donor0', 'convincing VAF scatter'), 
                ('N241', 'donor1', 'convincing VAF scatter'), 
                ('N244', 'donor2', 'convincing VAF scatter'), 
                ('N243', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_03_03_22_2': {
            
            'pair_with_match_confidence_manual_comment': [
                ('N241', 'donor0', 'convincing VAF scatter'),
                ('N242', 'donor1', 'convincing VAF scatter'),
                ('N243', 'donor2', 'convincing VAF scatter'),
                ('N244', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_02_05_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N245', 'donor0'): 5, # but not enough data
                ('N245', 'donor1'): 95, 
                (np.nan, 'donor2'): None, # donor2 is irrelevant. vireo assigned only 85 cells to donor2
                (np.nan, 'donor3'): None, # donor3 is irrelevant. vireo assigned only 71 cells to donor3
                
                ('N248', 'donor0'): 55, # but maybe not enough data
                ('N247', 'donor0'): -5, # but not enough data
                
                # all other alternatives for expected donors match badly. ???
                # ['N245', 'N246', 'N247', 'N248']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N248', 'donor0', 'somewhat convincing VAF scatter'), 
                ('N245', 'donor1', 'convincing VAF scatter'), 
            ],
            # 'donor_ids_for_which_extra_mip_data_would_be_nice_to_have': {'N248'}, # but not necessary. commented out because i guess the limiting factor here is scRNA-seq data, rather than MIP data.
            # 'rematch_anyway': True,
        },
        
        'demux_02_05_22_2': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N246', 'donor0'): -10, # but maybe not enough data
                ('N245', 'donor1'): 95,
                (np.nan, 'donor2'): None, # donor2 is irrelevant. vireo assigned only 84 cells to donor2
                (np.nan, 'donor3'): None, # donor3 is irrelevant. vireo assigned only 104 cells to donor3
                
                ('N248', 'donor0'): 50, # but not enough data
                ('N245', 'donor0'): 10, # but not enough data
                
                # all other alternatives for expected donors match badly.
                # ['N245', 'N246', 'N247', 'N248']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N248', 'donor0', 'convincing VAF scatter'), 
                ('N245', 'donor1', 'convincing VAF scatter'), 
            ],
            # 'donor_ids_for_which_extra_mip_data_would_be_nice_to_have': {'N248'}, # but not necessary. commented out because i guess the limiting factor here is scRNA-seq data, rather than MIP data.
            # 'rematch_anyway': True,
        },
        'demux_06_06_22_1': {
            # NOTE: in 230501_MDS_10x_age_and_gender.xlsx there is a mistake about N257. it says N257 is BM only, but indeed N257 PB is in demux_06_06_22_1. we also have a BM from N257 on 06.10.22 (in a batch of its own).
            
            # ['N249', 'N253', 'N254', 'N257']
            'pair_with_match_confidence_manual_comment': [
                ('N253', 'donor0', 'convincing VAF scatter'), 
                ('N249', 'donor1', 'convincing VAF scatter'), 
                ('N257', 'donor2', 'convincing VAF scatter'), 
                ('N254', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_06_06_22_2': {
            # ['N249', 'N253', 'N254', 'N257']
            'pair_with_match_confidence_manual_comment': [
                ('N249', 'donor0', 'convincing VAF scatter'), 
                ('N253', 'donor1', 'convincing VAF scatter'), 
                ('N254', 'donor2', 'convincing VAF scatter'), 
                ('N257', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_17_08_20_1': {
            # ['N35', 'N199', 'N37', 'N235']
            'pair_with_match_confidence_manual_comment': [
                ('N235', 'donor1', 'convincing VAF scatter'), 
                ('N35', 'donor3', 'convincing VAF scatter'), 
                ('N37', 'donor4', 'convincing VAF scatter'), 
                ('N199', 'donor5', 'convincing VAF scatter'), 

                # donor0 and donor2 are irrelevant. each has less than 8 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_26_06_23_1': {
            # ['G1', 'G2', 'G3', 'G4']

            # NOTE: vireo produced an empty GT_donors.vireo.vcf. also when i ran it again with 5 donors instead of 6.

            'pair_with_match_confidence_manual_comment': [
                ('dummy', 'donor1', 'unclear'), 
                # ('', 'donor3', 'convincing VAF scatter'), 
                # ('', 'donor4', 'convincing VAF scatter'), 
                # ('', 'donor5', 'convincing VAF scatter'), 

                # donor3 and donor4 are irrelevant. each has less than 9 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_03_07_23_1': {
            # ['G5', 'G6', 'G7']
            'pair_with_match_confidence_manual_comment': [
                ('dummy', 'donor1', 'unclear'), 
                # ('', 'donor3', 'convincing VAF scatter'), 
                # ('', 'donor4', 'convincing VAF scatter'), 
                # ('', 'donor5', 'convincing VAF scatter'), 

                # donor0 and donor4 are irrelevant. each has less than 50 cells.
            ],
            # 'rematch_anyway': True,
        },
        'G8_13_07_23_1': {
            # ['G8']
            'pair_with_match_confidence_manual_comment': [
                ('G8', 'donor0', 'convincing VAF scatter'), 
                ('G8', 'donor2', 'convincing VAF scatter'), 

                # donor1 is irrelevant - 0 cells.
            ],
            'donor_ids_for_which_extra_mip_data_is_required': {'G8'},
            # 'rematch_anyway': True,
        },
        'demux_n_20_07_23_1': {
            # ['N381', 'N382', 'N383', 'N384', 'N385', 'N386']
            'pair_with_match_confidence_manual_comment': [
                ('N383', 'donor1', 'convincing VAF scatter'), 
                # each of the other donors has 644 cells, and none with a clear match to MIP genotyping data.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_20_07_23_2': {
            # ['N381', 'N382', 'N383', 'N384', 'N385', 'N386']
            'pair_with_match_confidence_manual_comment': [
                ('N383', 'donor0', 'convincing VAF scatter'), 
                # each of the other donors has 928 cells, and none with a clear match to MIP genotyping data.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_24_07_23_1': {
            # ['N387']
            'pair_with_match_confidence_manual_comment': [
                ('N387', 'donor0', 'unclear'), 
                ('N387', 'donor1', 'unclear'), 
                ('N387', 'donor2', 'unclear'), 

            ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_24_07_23_1': {
            # ['N387']
            'pair_with_match_confidence_manual_comment': [
                ('N387', 'donor0', 'convincing VAF scatter'), 
                ('N387', 'donor1', 'convincing VAF scatter'), 
                ('N387', 'donor2', 'convincing VAF scatter'), 

            ],
            # 'rematch_anyway': True,
        },
        
        'demux_09_06_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                (np.nan, 'donor0'): None, # donor0 is irrelevant. vireo assigned only 45 cells to donor0
                (np.nan, 'donor2'): None, # donor2 is irrelevant. vireo assigned only 29 cells to donor2
                ('N251', 'donor1'): 95, 
                ('N250', 'donor3'): 85, # though it looks a bit weird. the vireo vals aren't clustered enough to my liking (at 0, 0.5, 1), maybe?
                
                # all other alternatives for expected donors match badly.
                # ['N250', 'N251', 'N252', 'N255']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N251', 'donor1', 'convincing VAF scatter'),
                ('N250', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_09_06_22_2': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N250', 'donor0'): 90, # though it looks a bit weird. the vireo vals aren't clustered enough to my liking (at 0, 0.5, 1), maybe?
                (np.nan, 'donor2'): None, # donor2 is irrelevant. vireo assigned only 87 cells to donor2
                (np.nan, 'donor3'): None, # donor3 is irrelevant. vireo assigned only 38 cells to donor3
                ('N251', 'donor1'): 98, 
                
                # all other alternatives for expected donors match badly.
                # ['N250', 'N251', 'N252', 'N255']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N250', 'donor0', 'convincing VAF scatter'), 
                ('N251', 'donor1', 'convincing VAF scatter'), # are N251 and N311 the same person? yes.
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_12_06_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N259', 'donor0'): 90,
                ('N260', 'donor1'): 95, 
                ('N260', 'donor2'): 0, 
                ('N258', 'donor3'): 80, # but maybe not enough data

                ('N256', 'donor2'): 85, # but maybe not enough data

                
                # all other alternatives for expected donors match badly.
                # ['N256', 'N258', 'N259', 'N260']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N259', 'donor0', 'convincing VAF scatter'), 
                ('N260', 'donor1', 'convincing VAF scatter'), 
                ('N256', 'donor2', 'convincing VAF scatter'), 
                ('N258', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },

        'demux_13_06_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                (np.nan, 'donor1'): None, # donor1 is irrelevant. vireo assigned only 175 cells to donor1
                (np.nan, 'donor3'): None, # donor3 is irrelevant. vireo assigned only 239 cells to donor3
                ('N262', 'donor0'): 95, 
                ('N264', 'donor2'): -10, # but not enough data
                
                ('N261', 'donor2'): 65, # but not enough data

                
                
                # all other alternatives for expected donors match badly.
                # ['N261', 'N262', 'N263', 'N264']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N262', 'donor0', 'convincing VAF scatter'), 
                ('N261', 'donor2', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_19_06_22_1': {
            'pair_to_manual_score_only_considering_vaf_scatter': {
                # rho predicted matches
                ('N266', 'donor0'): 95,
                ('N267', 'donor1'): 98,
                ('N268', 'donor2'): 0,
                ('N268', 'donor3'): 98,

                ('N265', 'donor2'): 60, # but not enough data

                # all other alternatives for expected donors match badly.
                # ['N265', 'N266', 'N267', 'N268']
            },
            'pair_with_match_confidence_manual_comment': [
                ('N266', 'donor0', 'convincing VAF scatter'),
                ('N267', 'donor1', 'convincing VAF scatter'),
                ('N265', 'donor2', 'convincing VAF scatter'), 
                ('N268', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },

        'demux_19_06_22_2': {
            # ['N265', 'N266', 'N267', 'N268']
            'pair_with_match_confidence_manual_comment': [
                ('N265', 'donor0', 'convincing VAF scatter'), 
                ('N268', 'donor1', 'convincing VAF scatter'),
                ('N266', 'donor2', 'convincing VAF scatter'),
                ('N267', 'donor3', 'convincing VAF scatter'),

            ],
            # 'rematch_anyway': True,
        },
        
        'demux_03_07_22_1': {
            # ['N269', 'N270', 'N271', 'N272']
            'pair_with_match_confidence_manual_comment': [
                ('N269', 'donor2', 'convincing VAF scatter'), 
                ('N270', 'donor3', 'convincing VAF scatter'), 
                ('N272', 'donor4', 'convincing VAF scatter'),
                ('N271', 'donor5', 'convincing VAF scatter'), 
                
                # donor0 and donor1 are irrelevant. each has less than 311 cells.
            ],
            # 'rematch_anyway': True,
        },

        'demux_03_07_22_2': {
            # ['N269', 'N270', 'N271', 'N272']
            'pair_with_match_confidence_manual_comment': [
                ('N272', 'donor0', 'convincing VAF scatter'), 
                ('N270', 'donor1', 'convincing VAF scatter'), 
                ('N269', 'donor2', 'convincing VAF scatter'),
                ('N271', 'donor5', 'convincing VAF scatter'), 

                # donor3 and donor4 are irrelevant. each has less than 243 cells.
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_10_07_22_1': {
            # ['N273', 'N274', 'N275', 'N276']
            'pair_with_match_confidence_manual_comment': [
                ('N274', 'donor0', 'convincing VAF scatter'), 
                ('N276', 'donor1', 'convincing VAF scatter'), 
                ('N273', 'donor2', 'convincing VAF scatter'), 
                ('N275', 'donor3', 'convincing VAF scatter'), 

                # NOTE: N275 is the son of N56 and N295 (aka N55).
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_10_07_22_2': {
            # ['N273', 'N274', 'N275', 'N276']
            'pair_with_match_confidence_manual_comment': [
                ('N275', 'donor0', 'convincing VAF scatter'), 
                ('N274', 'donor1', 'convincing VAF scatter'), 
                ('N273', 'donor2', 'convincing VAF scatter'), 
                ('N276', 'donor3', 'convincing VAF scatter'), 
                
                # NOTE: N275 is the son of N56 and N295 (aka N55).
            ],
            # 'rematch_anyway': True,
        },

        'demux_02_05_22_1_ultima': {
            # ['N245', 'N246', 'N247', 'N248']
            # now that we know about the accidental mixing:
            # 'N245', 'N246', 'N247', 'N248', 'N265', 'N266', 'N267', 'N268'
            'pair_with_match_confidence_manual_comment': [
                ('N266', 'donor4', 'convincing VAF scatter'), 
                ('N267', 'donor5', 'convincing VAF scatter'), 
                ('N245', 'donor6', 'convincing VAF scatter'), 
                ('N268', 'donor9', 'convincing VAF scatter'), 
                
                # donor3 has less than ~800 cells. and their cells don't match any donor much better than others. (N112 has the best match)

                # donor0, donor1, donor2, donor7 and donor8 are irrelevant. each has less than ~500 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_02_05_22_2_ultima': {
            # ['N245', 'N246', 'N247', 'N248']
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor1, donor4 and donor5 are irrelevant. each has less than ~100 cells.

                ('N248', 'donor2', 'convincing VAF scatter'), 
                ('N245', 'donor3', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_06_06_22_2_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor2 and donor5 each has less than 31 cells.
                ('N249', 'donor0', 'convincing VAF scatter'), 
                ('N257', 'donor1', 'convincing VAF scatter'), 
                ('N254', 'donor3', 'convincing VAF scatter'), 
                ('N253', 'donor4', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_09_06_22_1_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor1, donor2, donor6 and donor7 each has less than ~300 cells.
                
                ('N251', 'donor3', 'convincing VAF scatter'),
                ('N273', 'donor4', 'convincing VAF scatter'),
                ('N274', 'donor5', 'convincing VAF scatter'),
                ('N250', 'donor8', 'convincing VAF scatter'), # this actually seems to contain two donors - also N276...
                ('N275', 'donor9', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_09_06_22_2_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor2, donor3 and donor4 each has less than ~150 cells.


                # donor0 is irrelevant. vireo assigned only 5 cells to them.
                # donor1 is irrelevant. vireo assigned only 3 cells to them.

                ('N251', 'donor1', 'convincing VAF scatter'), 
                ('N250', 'donor5', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_12_06_22_1_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor2 and donor3 each has less than ~5 cells.
                ('N259', 'donor0', 'convincing VAF scatter'), 
                ('N256', 'donor1', 'convincing VAF scatter'), 
                ('N258', 'donor4', 'convincing VAF scatter'), 
                ('N260', 'donor5', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_13_06_22_1_ultima': {
            # 'N261', 'N262', 'N263', 'N264'
            'pair_with_match_confidence_manual_comment': [
                # donor2, donor3 and donor5 each has less than ~370 cells.
                
                ('N262', 'donor0', 'convincing VAF scatter'), 
                
                ('N261', 'donor4', 'unclear'), 
                # ('N245', 'donor4', 'unclear'), 
                
                ('N261', 'donor1', 'matches badly'), 
                ('N262', 'donor1', 'matches badly'), 
                ('N263', 'donor1', 'matches badly'), 
                ('N264', 'donor1', 'matches badly'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_19_06_22_1_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor0 and donor4 each has less than ~20 cells.

                ('N267', 'donor1', 'convincing VAF scatter'), 
                ('N265', 'donor2', 'convincing VAF scatter'), 
                ('N268', 'donor3', 'convincing VAF scatter'), 
                ('N266', 'donor5', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_19_06_22_2_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor2 and donor7 each has less than ~400 cells.
                
                ('N265', 'donor1', 'convincing VAF scatter'),
                ('N253', 'donor3', 'convincing VAF scatter'),
                ('N257', 'donor4', 'convincing VAF scatter'),
                ('N267', 'donor5', 'convincing VAF scatter'),
                ('N268', 'donor6', 'convincing VAF scatter'),
                ('N249', 'donor8', 'convincing VAF scatter'),
                ('N266', 'donor9', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_03_07_22_1_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor1 and donor5 are irrelevant. each has less than ~400 cells.
                
                ('N269', 'donor0', 'convincing VAF scatter'),
                ('N270', 'donor2', 'convincing VAF scatter'),
                ('N271', 'donor3', 'convincing VAF scatter'),
                ('N272', 'donor4', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_03_07_22_2_ultima': {
            'pair_with_match_confidence_manual_comment': [
                # donor2 and donor4 are irrelevant. each has less than ~500 cells.
                ('N272', 'donor0', 'convincing VAF scatter'), 
                ('N270', 'donor1', 'convincing VAF scatter'), 
                ('N271', 'donor3', 'convincing VAF scatter'), 
                ('N269', 'donor5', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_10_07_22_1_ultima': {
            # 'N250', 'N251', 'N252', 'N255', 'N273', 'N274', 'N275', 'N276'
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor2, donor4 and donor5 each has less than ~500 cells.
                ('N251', 'donor1', 'convincing VAF scatter'),
                ('N274', 'donor3', 'convincing VAF scatter'),
                ('N276', 'donor6', 'convincing VAF scatter'),
                ('N250', 'donor7', 'unclear'),
                ('N275', 'donor8', 'convincing VAF scatter'),
                ('N273', 'donor9', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_10_07_22_2_ultima': {
            # 'N273', 'N274', 'N275', 'N276', 'N286', 'N287', 'N288', 'N289', 'N290'
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor3, donor6 and donor8 each has less than ~600 cells.
                ('N286', 'donor1', 'convincing VAF scatter'),
                ('N274', 'donor2', 'convincing VAF scatter'),
                ('N273', 'donor4', 'convincing VAF scatter'),
                ('N276', 'donor5', 'convincing VAF scatter'),
                ('N275', 'donor7', 'convincing VAF scatter'),
                ('N288', 'donor9', 'convincing VAF scatter'),
                ('N286', 'donor10', 'unclear'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_01_02_21_2': {
            # NOTE: for some reason, a huge fraction of genotype doublets here have very low log_norm_cell_ranger_umi_count

            # 'N130', 'N133', 'N135', 'N232'
            'pair_with_match_confidence_manual_comment': [
                # donor2 and donor3 each has less than ~400 cells.

                ('N232', 'donor0', 'convincing VAF scatter'), 
                ('N135', 'donor1', 'convincing VAF scatter'), 
                ('N130', 'donor4', 'convincing VAF scatter'), 
                ('N133', 'donor5', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_01_02_21_2_illumina': {
            # 'N130', 'N133', 'N135', 'N232'
            'pair_with_match_confidence_manual_comment': [
                ('N133', 'donor0', 'convincing VAF scatter'), 
                ('N130', 'donor1', 'convincing VAF scatter'), 
                ('N232', 'donor4', 'convincing VAF scatter'), 
                ('N135', 'donor5', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_22_02_21_1': {
            # 'N48', 'N154', 'N158', 'N160', 'N215'
            
            'pair_with_match_confidence_manual_comment': [
                # donor3 and donor6 each has less than 172 cells.
                ('N160', 'donor0', 'convincing VAF scatter'), 
                ('N158', 'donor1', 'convincing VAF scatter'), 
                ('N48', 'donor2', 'convincing VAF scatter'), 
                ('N154', 'donor4', 'convincing VAF scatter'), 
                ('N215', 'donor5', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_22_02_21_2': {
            # 'N152', 'N153', 'N155', 'N156', 'N233'
            'pair_with_match_confidence_manual_comment': [
                # donor1 and donor4 each has less than ~350 cells.
                
                ('N156', 'donor0', 'convincing VAF scatter'), 
                ('N155', 'donor2', 'convincing VAF scatter'), 
                ('N152', 'donor3', 'convincing VAF scatter'), 
                ('N233', 'donor5', 'convincing VAF scatter'), 
                ('N153', 'donor6', 'convincing VAF scatter'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_22_02_21_2_illumina': {
            # 'N152', 'N153', 'N155', 'N156', 'N233'
            'pair_with_match_confidence_manual_comment': [
                ('N233', 'donor0', 'convincing VAF scatter'), 
                ('N156', 'donor3', 'convincing VAF scatter'), 
                ('N155', 'donor4', 'convincing VAF scatter'), 
                ('N152', 'donor5', 'convincing VAF scatter'), 
                # ('N153', 'donor6', 'somewhat convincing VAF scatter'), # soup split this one evenly between a donor, unassigned and doublet... and anyway it wasn't very convincing, and only 536 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_07_03_21_1': {
            # 'N166', 'N167', 'N168', 'N169', 'N171', 'N176'
            'pair_with_match_confidence_manual_comment': [
                # donor3 and donor6 each has less than ~300 cells.

                ('N176', 'donor0', 'convincing VAF scatter'),
                ('N169', 'donor1', 'convincing VAF scatter'),
                ('N168', 'donor2', 'convincing VAF scatter'),
                ('N166', 'donor4', 'convincing VAF scatter'),
                ('N171', 'donor5', 'convincing VAF scatter'),
                ('N167', 'donor7', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_07_03_21_1_illumina': {
            # 'N166', 'N167', 'N168', 'N169', 'N171', 'N176'
            'pair_with_match_confidence_manual_comment': [
                ('N176', 'donor0', 'convincing VAF scatter'),
                ('N169', 'donor3', 'convincing VAF scatter'),
                ('N166', 'donor4', 'convincing VAF scatter'),
                ('N167', 'donor7', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_11_04_21_2': {
            # NOTE: for some reason, a huge fraction of genotype doublets here have very low log_norm_cell_ranger_umi_count
            
            # 'N177', 'N178', 'N184', 'N200', 'N210'
            'pair_with_match_confidence_manual_comment': [
                # ('N184', 'donor0', 'convincing VAF scatter'), # less than 300 cells
                ('N178', 'donor1', 'convincing VAF scatter'),
                ('N184', 'donor2', 'convincing VAF scatter'),
                ('N210', 'donor3', 'convincing VAF scatter'),
                ('N200', 'donor4', 'convincing VAF scatter'),
                # ('N184', 'donor5', 'somewhat convincing VAF scatter'), # the problem is that N177 also matches donor5 very well. ~600 cells 
                ('N177', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_11_04_21_2_illumina': {
            # 'N177', 'N178', 'N184', 'N200', 'N210'
            'pair_with_match_confidence_manual_comment': [
                ('N210', 'donor0', 'convincing VAF scatter'),
                ('N178', 'donor1', 'convincing VAF scatter'),
                ('N177', 'donor2', 'convincing VAF scatter'),
                ('N200', 'donor3', 'convincing VAF scatter'),
                ('N184', 'donor4', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_23_11_20_1': {
            # 'N74', 'N75', 'N76', 'N77'
            'pair_with_match_confidence_manual_comment': [
                # donor3 and donor4 each has less than ~10 cells.
                ('N77', 'donor0', 'convincing VAF scatter'),
                ('N76', 'donor1', 'convincing VAF scatter'),
                ('N74', 'donor2', 'convincing VAF scatter'),
                ('N75', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_23_11_20_1_illumina': {
            # 'N74', 'N75', 'N76', 'N77'
            'pair_with_match_confidence_manual_comment': [
                ('N74','donor2', 'convincing VAF scatter'),
                ('N75','donor3', 'convincing VAF scatter'),
                ('N76','donor4', 'convincing VAF scatter'),
                ('N77','donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_30_11_20_1': { # NOTE: 240523: my mistake, i think. seems like it really should have been demux_30_11_20_2, but it's too late now and we will keep it this way.
            # 'N217', 'N224', 'N81', 'N82'
            'pair_with_match_confidence_manual_comment': [
                ('N82', 'donor0', 'convincing VAF scatter'),
                ('N81', 'donor1', 'convincing VAF scatter'),
                ('N224', 'donor2', 'convincing VAF scatter'),
                ('N217', 'donor5', 'convincing VAF scatter'),
                # donor3 is irrelevant. vireo assigned only ~20 cells to them.
                # donor4 is irrelevant. vireo assigned only ~60 cells to them.
            ],
            # 'rematch_anyway': True,
        },
        'demux_30_11_20_1_illumina': { # NOTE: 240523: my mistake, i think. seems like it really should have been demux_30_11_20_2_illumina, but it's too late now and we will keep it this way.
            # 'N217', 'N224', 'N81', 'N82'
            'pair_with_match_confidence_manual_comment': [
                ('N82', 'donor1', 'convincing VAF scatter'),
                ('N224', 'donor2', 'convincing VAF scatter'),
                ('N81', 'donor3', 'convincing VAF scatter'),
                ('N217', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'T6_25_05_22_1': {
            'pair_with_match_confidence_manual_comment': [
                (np.nan, 'donor3', 'matches badly'), # this is just so that the df won't be empty. we are not going to have the genotype of T6, i guess.
                
                # TODO: ugh. not sure what to do about this experiment. we don't have MIP genotyping data for T6, IIRC.
            ],
            # 'rematch_anyway': True,
        },
        'NS20_29_09_21_1': {
            'pair_with_match_confidence_manual_comment': [
                ('NS20', 'donor0', 'unclear'),
                ('NS20', 'donor1', 'unclear'),
                ('NS20', 'donor2', 'unclear'),
                
                # TODO: fix when NS20 MIP data is available
            ],
            # 'rematch_anyway': True,
        },
        'N328_15_12_22_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N328', 'donor0', 'convincing VAF scatter'),
                ('N328', 'donor1', 'convincing VAF scatter'),
                ('N328', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_01_01_23_1': {
            # 'N196', 'N240', 'N323', 'N1', 'N325'

            'pair_with_match_confidence_manual_comment': [
                # donor1 and donor3 are irrelevant. vireo assigned less than ~10 cells to them.
                
                ('N1', 'donor0', 'convincing VAF scatter'),
                ('N325', 'donor2', 'convincing VAF scatter'),
                ('N240', 'donor4', 'convincing VAF scatter'),
                ('N323', 'donor5', 'convincing VAF scatter'),
                ('N196', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_15_01_23_1': {
            # 'N188', 'N325', 'N329', 'N331', 'N333', 'N334'
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor1 and donor7 each has less than 8 cells.
                ('N325', 'donor2', 'convincing VAF scatter'),
                ('N329', 'donor3', 'convincing VAF scatter'),
                ('N331', 'donor4', 'convincing VAF scatter'),
                ('N188', 'donor5', 'convincing VAF scatter'),
                ('N334', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_29_01_23_1': {
            # 'N335', 'N336', 'N337', 'N338', 'N339'
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor1, donor3 and donor5 each has less than ~200 cells.

                ('N336', 'donor2', 'convincing VAF scatter'),
                ('N337', 'donor4', 'convincing VAF scatter'),
                ('N335', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_05_02_23_1': {
            # 'N276','N198','N342','N343','N344'
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor1, donor2, donor3, donor4 and donor5 each has less than ~350 cells.
                ('N198', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_08_02_23_1': {
            # 'N350','N351','N352','N353'
            'pair_with_match_confidence_manual_comment': [
                

                ('N352', 'donor0', 'convincing VAF scatter'),
                ('N350', 'donor1', 'convincing VAF scatter'),
                ('N351', 'donor2', 'convincing VAF scatter'),
                ('N353', 'donor7', 'convincing VAF scatter'),
                ('N328', 'donor8', 'convincing VAF scatter'),
                # donor4 and donor6 each has less than 614 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_27_02_23_1': {
            # NOTE: my guess: this is not really demux_27_02_23_1, but demux_28_02_23_1

            # demux_28_02_23_1: 'N307', 'N360', 'N362', 'N363'
            # demux_27_02_23_1: 'N354', 'N355', 'N356', 'N357', 'N260', 'N359'
            'pair_with_match_confidence_manual_comment': [
                ('N351', 'donor0', 'convincing VAF scatter'),
                ('N360', 'donor1', 'convincing VAF scatter'),
                ('N363', 'donor2', 'convincing VAF scatter'),
                ('N353', 'donor4', 'convincing VAF scatter'),
                ('N352', 'donor6', 'convincing VAF scatter'),
                ('N350', 'donor7', 'convincing VAF scatter'),
                ('N362', 'donor8', 'convincing VAF scatter'),
                ('N307', 'donor9', 'convincing VAF scatter'),
                # donor3 and donor5 each has less than 79 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_27_02_23_2': {
            # 'N354', 'N355', 'N356', 'N357', 'N260', 'N359'
            'pair_with_match_confidence_manual_comment': [
                ('N356', 'donor1', 'convincing VAF scatter'),
                ('N355', 'donor2', 'convincing VAF scatter'),
                ('N357', 'donor3', 'convincing VAF scatter'),
                ('N354', 'donor5', 'convincing VAF scatter'),
                ('N359', 'donor6', 'convincing VAF scatter'),
                ('N260', 'donor7', 'convincing VAF scatter'),
                # donor0 and donor4 each has less than 29 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_28_02_23_1': {
            # NOTE: my guess: this is not really demux_28_02_23_1, but demux_27_02_23_1

            # demux_27_02_23_1: 'N354', 'N355', 'N356', 'N357', 'N260', 'N359'

            # 'N360','N307','N362','N363'
            'pair_with_match_confidence_manual_comment': [
                # donor0 and donor3 each has less than 17 cells.
                ('N355', 'donor1', 'convincing VAF scatter'),
                ('N359', 'donor2', 'convincing VAF scatter'),
                ('N357', 'donor4', 'convincing VAF scatter'),
                ('N354', 'donor5', 'convincing VAF scatter'),
                ('N260', 'donor6', 'convincing VAF scatter'),
                ('N356', 'donor7', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_26_03_23_1': {
            # 'N368', 'N369', 'N211', 'N208', 'N200'
            'pair_with_match_confidence_manual_comment': [
                ('N369', 'donor2', 'convincing VAF scatter'),
                ('N368', 'donor3', 'convincing VAF scatter'),
                ('N208', 'donor4', 'convincing VAF scatter'),
                ('N211', 'donor5', 'convincing VAF scatter'),
                ('N200', 'donor6', 'convincing VAF scatter'),
                # donor0 and donor1 each has less than 93 cells.
            ],
            # 'rematch_anyway': True,
        },
        'N279_26_10_22_1': {
            # 'N279'
            'pair_with_match_confidence_manual_comment': [
                ('N279', 'donor0', 'convincing VAF scatter'),
                ('N279', 'donor1', 'convincing VAF scatter'),
                ('N279', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_06_11_22_1': {
            # 'N277', 'N280', 'N281'
            'pair_with_match_confidence_manual_comment': [
                ('N281', 'donor0', 'convincing VAF scatter'),
                ('N277', 'donor2', 'convincing VAF scatter'),
                ('N280', 'donor3', 'convincing VAF scatter'), 
                # donor1 is 264 cells.
                # donor4 is 10 cells.
            ],
            # 'rematch_anyway': True,
        },
        'N251_04_12_22_1': {
            # 'N251'
            'pair_with_match_confidence_manual_comment': [
                # donor1, donor2, donor3, donor4, donor6 and donor4 each has less than 777 cells.
                ('N251', 'donor5', 'convincing VAF scatter'),
                # ('?', 'donor0', 'convincing VAF scatter'),
                
            ],
            # 'rematch_anyway': True,
        },
        'demux_13_11_22_1': {
            # 'N282', 'N283', 'N284', 'N285'
            'pair_with_match_confidence_manual_comment': [
                ('N283', 'donor1', 'convincing VAF scatter'),
                ('N284', 'donor3', 'convincing VAF scatter'),
                ('N285', 'donor4', 'convincing VAF scatter'),
                ('N282', 'donor5', 'convincing VAF scatter'),
                # donor0 and donor2 each has less than 21 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_13_11_22_2': {
            # 'N282', 'N283', 'N284', 'N285'
            'pair_with_match_confidence_manual_comment': [
                ('N285', 'donor1', 'convincing VAF scatter'),
                ('N282', 'donor2', 'convincing VAF scatter'),
                ('N284', 'donor5', 'convincing VAF scatter'),
                ('N283', 'donor7', 'convincing VAF scatter'),
                # each of donor0, donor3, donor4, donor6 has less than 581 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_14_11_22_1': {
            # 'N286', 'N287', 'N288', 'N289', 'N290'
            'pair_with_match_confidence_manual_comment': [
                # donor1, donor3 each has less than 439 cells.
                ('N289', 'donor0', 'convincing VAF scatter'), # oh my. before force-cells we had no cluster for this donor!
                ('N288', 'donor2', 'convincing VAF scatter'),
                ('N290', 'donor4', 'convincing VAF scatter'),
                ('N286', 'donor5', 'convincing VAF scatter'),
                ('N287', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_23_11_22_1': {
            # 'N291', 'N292', 'N293', 'N294', 'N295'
            'pair_with_match_confidence_manual_comment': [
                # donor0 and donor3 each has less than 59 cells.
                ('N294', 'donor1', 'convincing VAF scatter'),
                ('N291', 'donor2', 'somewhat convincing VAF scatter'),
                ('N295', 'donor4', 'convincing VAF scatter'),
                ('N292', 'donor5', 'convincing VAF scatter'),
                ('N293', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_04_12_22_1': {
            # 'N12', 'N306', 'N307', 'N308', 'N309'
            
            # NOTE: N100 and N234 are the parents of N306
            
            'pair_with_match_confidence_manual_comment': [
                # donor0, donor1 and donor2 each has less than ~40 cells.

                ('N308', 'donor3', 'convincing VAF scatter'),
                ('N309', 'donor4', 'convincing VAF scatter'),
                ('N306', 'donor5', 'convincing VAF scatter'),
                ('N12', 'donor6', 'convincing VAF scatter'), # NOTE: this is bad. we can't know whether these cells came from BM or PB, because of well contamination (or whatever accident happened).
                ('N307', 'donor7', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_08_12_22_1': {
            # 'N12', 'N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316'
            'pair_with_match_confidence_manual_comment': [
                # ('N12', 'donor0', 'convincing VAF scatter'), # 6.918863237274595
                # ('N282', 'donor0', 'somewhat convincing VAF scatter'),
                ('N283', 'donor1', 'convincing VAF scatter'), # 6.832890014164741
                ('N269', 'donor2', 'convincing VAF scatter'), # 5.129283016944966
                ('N271', 'donor3', 'convincing VAF scatter'), # 6.339850002884624
                ('N270', 'donor5', 'convincing VAF scatter'), # 5.954196310386875
                ('N315', 'donor6', 'convincing VAF scatter'), # 7.011227255423254 # less than ~300 cells. but clear.
                ('N285', 'donor7', 'convincing VAF scatter'), # 8.174925682500678
                ('N272', 'donor8', 'convincing VAF scatter'), # 7.03342300153745
                # ('N12', 'donor9', 'convincing VAF scatter'), # 6.906890595608519
                # ('N282', 'donor9', 'somewhat convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_06_06_22_1_ultima': {
            # 'N249', 'N253', 'N254', 'N257'
            'pair_with_match_confidence_manual_comment': [
                # donor0 and donor4 each has less than ~20 cells.

                ('N253', 'donor1', 'convincing VAF scatter'),
                ('N249', 'donor2', 'convincing VAF scatter'),
                ('N254', 'donor3', 'convincing VAF scatter'),
                ('N257', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_13_12_22_1': {
            # 'N317', 'N318', 'N319', 'N320', 'N321', 'N322'
            'pair_with_match_confidence_manual_comment': [
                # donor3 and donor7 each has less than ~150 cells.
                ('N317', 'donor0', 'convincing VAF scatter'),
                ('N322', 'donor1', 'convincing VAF scatter'),
                ('N318', 'donor2', 'convincing VAF scatter'),
                ('N321', 'donor4', 'convincing VAF scatter'),
                ('N319', 'donor5', 'convincing VAF scatter'),
                ('N320', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_13_12_22_2': {
            # 'N317', 'N318', 'N319', 'N320', 'N321', 'N322'
            'pair_with_match_confidence_manual_comment': [
                # donor0 and donor6 each has less than ~2 cells.
                ('N317', 'donor1', 'convincing VAF scatter'),
                ('N322', 'donor2', 'convincing VAF scatter'),
                ('N321', 'donor3', 'convincing VAF scatter'),
                ('N318', 'donor4', 'convincing VAF scatter'),
                ('N319', 'donor5', 'convincing VAF scatter'),
                ('N320', 'donor7', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_27_11_22_1': {
            # 'N279', 'N296', 'N297', 'N298', 'N299', 'N300'
            'pair_with_match_confidence_manual_comment': [
                # donor0 and donor4 each has less than 16 cells.
                ('N296', 'donor1', 'convincing VAF scatter'),
                ('N299', 'donor2', 'convincing VAF scatter'),
                ('N279', 'donor3', 'convincing VAF scatter'),
                ('N297', 'donor5', 'convincing VAF scatter'),
                ('N298', 'donor6', 'convincing VAF scatter'),
                ('N300', 'donor7', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N276_bm_05_02_23_1': {
            # 'N276'
            'pair_with_match_confidence_manual_comment': [
                ('N276', 'donor0', 'convincing VAF scatter'), 
                ('dummy', 'donor1', 'unclear'), 
                ('N307', 'donor3', 'convincing VAF scatter'), 
                ('N353', 'donor4', 'convincing VAF scatter'), 
                ('dummy', 'donor5', 'unclear'), 
                ('N328', 'donor6', 'convincing VAF scatter'),
                ('N264', 'donor7', 'convincing VAF scatter'), 
                # donor2 has less than 422 cells.
            ],
            # 'rematch_anyway': True,
        },
        'N48_bm_12_03_23_1': {
            # 'N48'
            'pair_with_match_confidence_manual_comment': [
                # donor5 has less than ~10 cells.

                ('N48', 'donor0', 'convincing VAF scatter'),
                ('N48', 'donor1', 'convincing VAF scatter'),
                ('N48', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N365_bm_12_03_23_1': {
            # 'N365'
            'pair_with_match_confidence_manual_comment': [
                ('N365', 'donor0', 'convincing VAF scatter'),
                ('N365', 'donor1', 'convincing VAF scatter'),
                ('N365', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N242_bm_12_03_23_1': {
            # 'N365'
            'pair_with_match_confidence_manual_comment': [
                ('N242', 'donor0', 'convincing VAF scatter'),
                ('N242', 'donor1', 'convincing VAF scatter'),
                ('N242', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N279_bm_26_10_22_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N279', 'donor0', 'convincing VAF scatter'),
                ('N279', 'donor1', 'convincing VAF scatter'),
                ('N279', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N280_bm_06_11_22_1': {
            'pair_with_match_confidence_manual_comment': [
                # all vireo donors other than donor0, donor1 and donor8 look very noisy and unreliable.
                ('N279', 'donor0', 'convincing VAF scatter'), 
                ('N280', 'donor1', 'convincing VAF scatter'), # NOTE: this is bad. we can't know whether these cells came from BM or PB, because of well contamination (or whatever accident happened).
                ('N12', 'donor8', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N12_bm_04_12_22_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N251', 'donor0', 'convincing VAF scatter'),
                ('N251', 'donor1', 'convincing VAF scatter'),
                ('N251', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N188_bm_15_01_23_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N188', 'donor0', 'convincing VAF scatter'),
                ('N188', 'donor1', 'convincing VAF scatter'),
                ('N188', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N257_bm_06_10_22_1_illumina': {
            'pair_with_match_confidence_manual_comment': [
                ('N257', 'donor0', 'convincing VAF scatter'),
                ('N257', 'donor1', 'convincing VAF scatter'),
                ('N257', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N279_26_10_22_1_illumina': {
            'pair_with_match_confidence_manual_comment': [
                ('N279', 'donor0', 'convincing VAF scatter'),
                ('N279', 'donor1', 'convincing VAF scatter'),
                ('N279', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N279_bm_26_10_22_1_illumina': {
            'pair_with_match_confidence_manual_comment': [
                ('N279', 'donor0', 'convincing VAF scatter'),
                ('N279', 'donor1', 'convincing VAF scatter'),
                ('N279', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N196_bm_01_01_23_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N196', 'donor0', 'convincing VAF scatter'),
                ('N196', 'donor1', 'convincing VAF scatter'),
                ('N196', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N251_bm_04_12_22_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N334', 'donor1', 'convincing VAF scatter'),
                ('N12', 'donor2', 'convincing VAF scatter'), # NOTE: this is bad. we can't know whether these cells came from BM or PB, because of well contamination (or whatever accident happened).

                # all donors except donor1 and donor2 have less than ~1000 cells (except donor3 and donor6), and are spread across soup donors.
            ],
            # 'rematch_anyway': True,
        },
        'N257_bm_06_10_22_1': {
            # NOTE: vireo produced an empty GT_donors.vireo.vcf.

            'pair_with_match_confidence_manual_comment': [
                ('N257', 'donor0', 'unclear'),
                # ('N196', 'donor1', 'convincing VAF scatter'),
                # ('N196', 'donor2', 'somewhat convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'N334_bm_15_01_23_1': {
            'pair_with_match_confidence_manual_comment': [
                ('N334', 'donor0', 'convincing VAF scatter'),
                # ('N269', 'donor1', 'convincing VAF scatter'), # ugh. doesn't make sense as N269 is in an experiment with 3 others, and doesn't dominate the experiment. also, when comparing to the vireo-soup cluster of that experiment, it doesn't look like N269
                ('N334', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_22_06_23_2': {
            # this is the 24,36,48 hour experiment.

            # ['N377', 'N285', 'N264', 'N275']
            'pair_with_match_confidence_manual_comment': [
                ('N377', 'donor0', 'convincing VAF scatter'),
                ('N264', 'donor2', 'convincing VAF scatter'),
                ('N275', 'donor4', 'convincing VAF scatter'),
                ('N285', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_22_06_23_1': {
            # ['N373', 'N375', 'N376', 'N365']
            # one of them is N365, i guess.
            'pair_with_match_confidence_manual_comment': [
                ('N376', 'donor1', 'convincing VAF scatter'),
                ('N365', 'donor3', 'convincing VAF scatter'),
                ('N375', 'donor4', 'convincing VAF scatter'),
                # donor0, donor2 and donor5 each has less than 89 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_13_08_23_1': {
            # ['N338', 'N388', 'N390', 'N391']
            'pair_with_match_confidence_manual_comment': [
                ('N391', 'donor0', 'convincing VAF scatter'),
                ('N388', 'donor1', 'convincing VAF scatter'),
                ('N390', 'donor3', 'convincing VAF scatter'),
                ('N338', 'donor5', 'convincing VAF scatter'),
                # donor2 and donor4 each has less than 17 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_15_08_23_1': {
            # ['N392', 'N280']
            'pair_with_match_confidence_manual_comment': [
                ('N280', 'donor0', 'convincing VAF scatter'),
                ('N280', 'donor2', 'convincing VAF scatter'),
                # donor1 and donor3 each has less than 36 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_15_08_23_2': {
            # ['N392', 'N280', 'N264', 'N1']

            'pair_with_match_confidence_manual_comment': [
                ('N280', 'donor1', 'convincing VAF scatter'),
                ('N1', 'donor2', 'convincing VAF scatter'),
                ('N264', 'donor3', 'convincing VAF scatter'),
                ('N392', 'donor4', 'convincing VAF scatter'),
                # donor0 and donor5 each has less than 71 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_15_08_23_3': {
            # ['N264', 'N1', 'N252']

            'pair_with_match_confidence_manual_comment': [
                ('N1', 'donor0', 'convincing VAF scatter'),
                ('N264', 'donor1', 'convincing VAF scatter'),
                ('N252', 'donor2', 'convincing VAF scatter'),
                # donor3 and donor4 each has less than 27 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_24_07_23_1': {
            # ['G8', 'G9', 'G10', 'G11']

            'pair_with_match_confidence_manual_comment': [
                ('G9', 'donor0', 'convincing VAF scatter'),
                ('G9', 'donor1', 'convincing VAF scatter'),
                ('G9', 'donor2', 'convincing VAF scatter'),
                ('G9', 'donor3', 'convincing VAF scatter'),
                ('G9', 'donor4', 'convincing VAF scatter'),
                ('G9', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_12_09_23_1': {
            # ['G12', 'G13', 'G14']

            'pair_with_match_confidence_manual_comment': [
                ('N320', 'donor2', 'convincing VAF scatter'), 
                ('G12', 'donor1', 'convincing VAF scatter'), 
                # ('N1', 'donor0', 'convincing VAF scatter'),
                
                # donor3 and donor4 each has less than 92 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_28_09_23_1': {
            # ['G16', 'G17', 'G18', 'G19', 'G20']
            'pair_with_match_confidence_manual_comment': [
                ('G20', 'donor1', 'convincing VAF scatter'),
                ('G19', 'donor2', 'convincing VAF scatter'),
                ('G17', 'donor3', 'convincing VAF scatter'),
                
                # donor0, donor5 and donor6 each has less than 78 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_07_12_23_1': {
            # ['N259', 'G22', 'G23']

            'pair_with_match_confidence_manual_comment': [
                ('G23', 'donor1', 'somewhat convincing VAF scatter'),
                ('G22', 'donor3', 'convincing VAF scatter'),
                ('N259', 'donor4', 'convincing VAF scatter'),
                
                # donor0 and donor2 each has less than 39 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_18_12_23_1': {
            # ['N204', 'N413']

            'pair_with_match_confidence_manual_comment': [
                ('N413', 'donor0', 'convincing VAF scatter'), 
                ('N204', 'donor1', 'convincing VAF scatter'), # though not as good as i would hope.
                
                # donor2 and donor3 each has less than 94 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_18_12_23_2': {
            # vireo created an empty GT_donors.vireo.vcf. ugh. changing from 4 to 3 donors didn't help. also changing the seed.
            # ['N281', 'N413']

            'pair_with_match_confidence_manual_comment': [
                ('dummy', 'donor0', 'unclear'), 
                # ('N259', 'donor3', 'convincing VAF scatter'),
                
                # donor1 and donor2 each has less than 333 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_19_09_23_1': {
            # ['N200']

            'pair_with_match_confidence_manual_comment': [
                ('N200', 'donor0', 'convincing VAF scatter'),
                ('N200', 'donor1', 'convincing VAF scatter'),
                ('N200', 'donor2', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_19_09_23_1': {
            # ['N200', 'N402', 'N403']

            'pair_with_match_confidence_manual_comment': [
                ('N403', 'donor0', 'convincing VAF scatter'),
                ('N403', 'donor1', 'convincing VAF scatter'),
                
                # donor2, donor3 and donor4 each has less than 59 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_30_11_23_1': {
            # ['N279', 'N281', 'N405']

            'pair_with_match_confidence_manual_comment': [
                ('N281', 'donor1', 'convincing VAF scatter'),
                ('N279', 'donor4', 'convincing VAF scatter'),
                ('N405', 'donor0', 'convincing VAF scatter'),
                # donor2 and donor3 each has less than 6 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_04_12_23_1': {
            # ['N191', 'N367']

            'pair_with_match_confidence_manual_comment': [
                ('N191', 'donor0', 'convincing VAF scatter'),
                ('N367', 'donor1', 'convincing VAF scatter'),
                
                # donor2 and donor3 each has less than 35 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_04_12_23_2': {
            # ['N191', 'N367']

            'pair_with_match_confidence_manual_comment': [
                ('N367', 'donor0', 'convincing VAF scatter'),
                ('N367', 'donor2', 'convincing VAF scatter'),
                ('N367', 'donor3', 'convincing VAF scatter'),
                
                # donor1 has less than 3 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_05_12_23_1': {
            # ['N192', 'N382', 'N410']

            'pair_with_match_confidence_manual_comment': [
                ('N410', 'donor0', 'convincing VAF scatter'),
                ('N192', 'donor1', 'convincing VAF scatter'),
                ('N382', 'donor2', 'convincing VAF scatter'),
                
                # donor3 has less than 51 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_26_09_23_1': {
            # ['G15']

            'pair_with_match_confidence_manual_comment': [
                ('G15', 'donor1', 'convincing VAF scatter'),
                
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_01_01_24_1': {
            # ['N414', 'N415', 'N416']

            'pair_with_match_confidence_manual_comment': [
                ('N415', 'donor1', 'somewhat convincing VAF scatter'), 
                ('N414', 'donor2', 'convincing VAF scatter'),
                ('N416', 'donor4', 'convincing VAF scatter'),
                
                # donor0 and donor3 each has less than 22 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_12_02_24_1': {
            # ['N418', 'N419']

            'pair_with_match_confidence_manual_comment': [
                ('N419', 'donor0', 'convincing VAF scatter'),
                ('N418', 'donor2', 'convincing VAF scatter'),

                # donor1 and donor3 each has less than 40 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_12_02_24_1': {
            # ['N418', 'N419', 'N339', 'N421']

            'pair_with_match_confidence_manual_comment': [
                ('N339', 'donor0', 'unclear'),
                ('N421', 'donor4', 'convincing VAF scatter'),
                
                # donor1, donor2, donor3 and donor5 each has less than 77 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_14_02_24_1': {
            # ['N422', 'N423']

            'pair_with_match_confidence_manual_comment': [
                ('N422', 'donor0', 'convincing VAF scatter'),
                ('N423', 'donor2', 'convincing VAF scatter'),
                
                # donor1 and donor3 each has less than 43 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_28_12_23_1': {
            # ['G27', 'G28', 'G29', 'G30', 'G31']
            'pair_with_match_confidence_manual_comment': [
                ('G31', 'donor0', 'convincing VAF scatter'),
                ('G31', 'donor6', 'convincing VAF scatter'),
                
                # donor1, donor3 and donor5 each has less than 193 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_08_01_24_1': {
            # ['G32', 'G33', 'G34']
            'pair_with_match_confidence_manual_comment': [
                ('G32', 'donor0', 'convincing VAF scatter'),
                ('G33', 'donor2', 'convincing VAF scatter'),
                ('G34', 'donor4', 'convincing VAF scatter'),
                
                # donor1 and donor3 each has less than 4 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_11_01_24_1': {
            # ['G35', 'G36', 'N417']

            'pair_with_match_confidence_manual_comment': [
                ('N417', 'donor1', 'convincing VAF scatter'),
                ('G36', 'donor2', 'convincing VAF scatter'),
                
                # donor0 and donor4 each has less than 35 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_01_02_24_1': {
            # ['G37', 'G38', 'G39']
            'pair_with_match_confidence_manual_comment': [
                ('G39', 'donor0', 'convincing VAF scatter'),
                ('G37', 'donor2', 'somewhat convincing VAF scatter'),
                ('G38', 'donor4', 'convincing VAF scatter'),
                
                # donor1 and donor3 each has less than 5 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_08_02_24_1': {
            # ['G40', 'G41', 'G42']
            'pair_with_match_confidence_manual_comment': [
                ('G41', 'donor0', 'convincing VAF scatter'),
                ('G42', 'donor1', 'convincing VAF scatter'),
                ('G40', 'donor3', 'convincing VAF scatter'),
                
                # donor2 and donor4 each has less than 41 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_22_02_24_1': {
            # ['G43', 'G44', 'G45', 'G46', 'G47', 'G48']
            'pair_with_match_confidence_manual_comment': [
                ('G46', 'donor2', 'convincing VAF scatter'),
                ('G47', 'donor5', 'convincing VAF scatter'),
                ('G45', 'donor6', 'convincing VAF scatter'),
                ('G48', 'donor7', 'convincing VAF scatter'),
                
                # donor0 and donor1 each has less than 7 cells.
            ],
            # 'rematch_anyway': True,
        },
        'demux_04_01_21_2': {
            # ['N108', 'N110', 'N220', 'N111']

            'pair_with_match_confidence_manual_comment': [
                ('N111', 'donor0', 'convincing VAF scatter'),
                ('N108', 'donor1', 'convincing VAF scatter'),
                ('N110', 'donor2', 'convincing VAF scatter'),
                ('N220', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        
        'demux_01_02_21_1': {
            # ['N113', 'N132', 'N134', 'N223']
            'pair_with_match_confidence_manual_comment': [
                ('N134', 'donor0', 'convincing VAF scatter'),
                ('N223', 'donor3', 'convincing VAF scatter'),
                ('N132', 'donor4', 'convincing VAF scatter'),
                ('N113', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_01_03_21_1': {
            # ['N161', 'N162', 'N163', 'N225', 'N226']
            'pair_with_match_confidence_manual_comment': [
                ('N162', 'donor0', 'convincing VAF scatter'),
                # NOTE: N158 is the mother of N161
                ('N161', 'donor1', 'convincing VAF scatter'),
                ('N163', 'donor2', 'convincing VAF scatter'),
                ('N226', 'donor4', 'somewhat convincing VAF scatter'), 
                ('N225', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_04_01_21_1': {
            # ['N107', 'N112', 'N114', 'N115', 'N238']
            'pair_with_match_confidence_manual_comment': [
                ('N107', 'donor0', 'convincing VAF scatter'),
                ('N114', 'donor1', 'convincing VAF scatter'),
                ('N238', 'donor3', 'convincing VAF scatter'),
                ('N115', 'donor5', 'convincing VAF scatter'),
                ('N112', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_07_12_20_1': {
            # ['N269', 'N84', 'N86', 'N89']
            'pair_with_match_confidence_manual_comment': [
                ('N269', 'donor0', 'convincing VAF scatter'),
                ('N89', 'donor1', 'convincing VAF scatter'),
                ('N86', 'donor2', 'convincing VAF scatter'),
                ('N84', 'donor4', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_07_12_20_2': {
            # ['N213', 'N59', 'N85', 'N87']
            'pair_with_match_confidence_manual_comment': [
                ('N85', 'donor0', 'convincing VAF scatter'),
                ('N87', 'donor2', 'convincing VAF scatter'),
                ('N213', 'donor3', 'convincing VAF scatter'),
                ('N59', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_11_01_21_1': {
            # ['N117', 'N119', 'N121', 'N18']
            'pair_with_match_confidence_manual_comment': [
                ('N119','donor0', 'convincing VAF scatter'),
                ('N117','donor1', 'convincing VAF scatter'),
                ('N18','donor2', 'convincing VAF scatter'),
                ('N121','donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_11_01_21_2': {
            # ['N118', 'N120', 'N123', 'N254']
            'pair_with_match_confidence_manual_comment': [
                ('N120', 'donor0', 'convincing VAF scatter'),
                ('N254', 'donor1', 'convincing VAF scatter'),
                ('N123', 'donor3', 'convincing VAF scatter'),
                ('N118', 'donor4', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_21_01_21_1': {
            # ['N126', 'N127', 'N16', 'N17']
            'pair_with_match_confidence_manual_comment': [
                ('N126', 'donor1', 'convincing VAF scatter'),
                ('N16', 'donor2', 'convincing VAF scatter'),
                ('N17', 'donor3', 'convincing VAF scatter'),
                ('N127', 'donor4', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_21_01_21_2': {
            # ['N124', 'N125', 'N128', 'N129']
            'pair_with_match_confidence_manual_comment': [
                ('N129', 'donor0', 'convincing VAF scatter'),
                ('N128', 'donor1', 'convincing VAF scatter'),
                ('N124', 'donor4', 'convincing VAF scatter'),
                ('N125', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_21_02_21_2': {
            # ['N143', 'N145', 'N146', 'N219', 'N273']
            'pair_with_match_confidence_manual_comment': [
                ('N143', 'donor0', 'convincing VAF scatter'),
                ('N273', 'donor1', 'convincing VAF scatter'),
                ('N219', 'donor2', 'convincing VAF scatter'),
                ('N145', 'donor4', 'convincing VAF scatter'),
                ('N146', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_21_12_20_1': { # NOTE: surprisingly, N92 is also here...
            # ['N268', 'N90', 'N93', 'N94']
            'pair_with_match_confidence_manual_comment': [
                ('N90', 'donor0', 'convincing VAF scatter'),
                ('N93', 'donor1', 'convincing VAF scatter'),
                ('N94', 'donor3', 'convincing VAF scatter'),
                ('N268', 'donor4', 'convincing VAF scatter'),
                ('N92', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_21_12_20_2': {
            # ['N216', 'N92', 'N95', 'N96']
            'pair_with_match_confidence_manual_comment': [
                ('N216', 'donor1', 'convincing VAF scatter'),
                ('N96', 'donor3', 'convincing VAF scatter'),
                ('N95', 'donor4', 'convincing VAF scatter'),
                ('N92', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_28_12_20_1': {
            # ['N100', 'N102', 'N103', 'N234']
            'pair_with_match_confidence_manual_comment': [
                ('N100', 'donor0', 'convincing VAF scatter'),
                ('N102', 'donor3', 'convincing VAF scatter'),
                ('N234', 'donor4', 'convincing VAF scatter'),
                ('N103', 'donor5', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_28_12_20_2': {
            # ['N101', 'N104', 'N105', 'N106']
            'pair_with_match_confidence_manual_comment': [
                ('N104', 'donor0', 'convincing VAF scatter'),
                ('N106', 'donor1', 'convincing VAF scatter'),
                ('N105', 'donor2', 'convincing VAF scatter'),
                ('N101', 'donor4', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_28_03_24_1': {
            # ['G49', 'G50']
            'pair_with_match_confidence_manual_comment': [
                ('G50', 'donor0', 'convincing VAF scatter'),
                ('G49', 'donor3', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_11_04_21_1': {
            # ['N83', 'N182', 'N183', 'N244']
            'pair_with_match_confidence_manual_comment': [
                ('N183', 'donor0', 'convincing VAF scatter'),
                ('N83', 'donor1', 'convincing VAF scatter'),
                ('N244', 'donor2', 'convincing VAF scatter'),
                ('N182', 'donor4', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_07_03_21_2': {
            # ['N174', 'N173', 'N172', 'N170', 'N231']
            'pair_with_match_confidence_manual_comment': [
                ('N173', 'donor0', 'convincing VAF scatter'),
                ('N170', 'donor2', 'convincing VAF scatter'),
                ('N231', 'donor3', 'convincing VAF scatter'),
                ('N174', 'donor4', 'convincing VAF scatter'),
                ('N172', 'donor6', 'convincing VAF scatter'),
            ],
            # 'rematch_anyway': True,
        },
    },

    'matching_vireo_to_souporcell': {
        'min_min_agreement_proportion': 0.75, # say we have 0,1,2 in souporcell and donor0,donor1,donor2 in vireo. for '1' and 'donor2', we calculate the proportion of '1' and 'donor2' assigned cells out of all '1' assigned cells, as well as out of all 'donor2' assigned cells. if '1' and 'donor2' is the best out of all '1' and also out of all 'donor2', than min_agreement_proportion is the lower proportion of the two proportion we calculated. otherwise, it is 0.
        'min_cell_count_soup_and_vireo_agree_on': 500, # used for auto identification
        'min_cell_count_for_auto_identified_bad_soup_vireo_pairs': 100,
        'default_min_num_of_cells_per_vireo_soup_donor': 500, # used to determine whether to automatically analyze vireo-soup donor (e.g., log(ratio) and VAF scatter).
    },
    'matching_vireo_to_souporcell_manual_examination': {
        'demux_28_02_22_1': {
            
            
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment_before_230707/demux_03_07_22_2_ultima/vaf_df_csvs/0_donor5_vaf_df.csv' # N269
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_07_02_22_1/vaf_df_csvs/1_donor0_vaf_df.csv' # N225
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_2/vaf_df_csvs/1_unassigned_vaf_df.csv' # N219? more probably a mix of N219, N235 and N236. ugh.
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_1/vaf_df_csvs/0_unassigned_vaf_df.csv' # N219? more probably a mix of N219, N235 and N236. ugh.
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_1/vaf_df_csvs/1_donor1_vaf_df.csv' # N235
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_1/vaf_df_csvs/1_unassigned_vaf_df.csv' # N235

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('0', 'donor0'), ('0', 'donor1'), ('0', 'donor2'), ('0', 'donor3')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            #     # ({('1', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            #     # ({('2', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            #     # ({('3', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                *list(itertools.product('0123', [f'donor{x}' for x in '0123'])), # high confidence for MIP to vireo matches, and it looks like for each donor, vireo VAFs are clustered at 0, 0.5 and 1, as expected. so just ignoring souporcell.
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('0', 'unassigned'),
            ],
            'manually_added_soup_vireo_mip_triplets': [
                # looked at the VAF scatters. and they all looked very supportive of N235, though not perfect (the ones on which soup and vireo disagreed looked better, though i guess that's not too bad, because we give low confidence to these ones).
                ('1', 'unassigned', 'N235'), # homozygous matches,mismatches: 161,12 (with freqs: 603,13)
                ('2', 'unassigned', 'N235'), # homozygous matches,mismatches: 156,6 (with freqs: 594,7)
                ('3', 'unassigned', 'N235'), # homozygous matches,mismatches: 158,10 (with freqs: 591,10)
            ],
            # 'rematch_anyway': True,
        },
        'demux_28_02_22_2': {
            

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('1', 'donor0'), ('1', 'donor1'), ('1', 'donor2'), ('1', 'donor3')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # #     ({('0', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            #     ({('1', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            # #     ({('2', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            # #     ({('3', 'unassigned')}, ['N219', 'N235', 'N236', 'N237']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                *list(itertools.product('0123', [f'donor{x}' for x in '0123'])), # ditto demux_28_02_22_1 argument.
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('1', 'unassigned'),
            ],
            'manually_added_soup_vireo_mip_triplets': [
                # looked at the VAF scatters. and 3/4 looked very supportive of N235, though not perfect. (the ones on which soup and vireo disagreed looked better, though i guess that's not too bad, because we give low confidence to these ones).
                ('0', 'unassigned', 'N235'), # homozygous matches,mismatches: 140,6 (with freqs: 533,8)
                # ('1', 'unassigned', 'N235'), # homozygous matches,mismatches: 17,4 - definitely not good enough. but this is expected - only 1013 cells. (with freqs: ugh. N235: 85,5, N219: 50,3, N236: 36,3, N237: 3,0)
                ('2', 'unassigned', 'N235'), # homozygous matches,mismatches: 142,6 (with freqs: 539,9)
                ('3', 'unassigned', 'N235'), # homozygous matches,mismatches: 144,7 (with freqs: 542,8)
            ],
            # 'rematch_anyway': True,
        },
        'demux_28_11_21_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_07_02_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_07_11_21_1': {
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor1')}, ['N185', 'N186', 'N187', 'N188']),
            #     ({('0', 'donor2')}, ['N185', 'N186', 'N187', 'N188']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # seems like souporcell decided there are only 3 donors (even though we specified 4 in the cmd argument it received). souporcell thought that vireo's donor1 and donor2 are the same person. but donor1 matches well N188, and clearly doesnt match the MIP data of the three other donors. similarly, for donor2 and N185. so the MIP data indicates it couldn't be that the cells of donor1 and donor2 are of the same person, so we just ignore souporcell here.
                ('0', 'donor1'), # with freqs: N188: 72,3. others were bad.
                ('0', 'donor2'), # with freqs: N185: 266,3. others were bad.
            ],
            # 'rematch_anyway': True,
        },

        'demux_19_12_21_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('3', 'donor1')},
            #     # {('6', 'donor1')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('0', 'donor3')}, ['N189', 'N194', 'N195', 'N196', 'N197']),
            #     # ({('0', 'donor4')}, ['N189', 'N194', 'N195', 'N196', 'N197']),
            #     # ({('0', 'unassigned')}, ['N189', 'N194', 'N195', 'N196', 'N197']), # ugh. can't tell with freqs. multiple look okay-ish.
            #     # ({('1', 'donor4')}, ['N189', 'N194', 'N195', 'N196', 'N197']),
            #     # ({('0', 'doublet')}, ['N189', 'N194', 'N195', 'N196', 'N197']), # with freqs, impressively matches somewhat okay-ish all donors. definitely discarding these.
            # ],
            
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     # ('4', 'donor0'),
            #     # ('5', 'donor2'),
            #     # ('2', 'donor3'),
            #     # ('1', 'donor4'),
            #     # ('0', 'donor5'),
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('5', 'unassigned'), 
                ('5', 'donor2'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_02_01_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('2', 'donor0')}, ['N198', 'N199', 'N201', 'N202']),
            #     # ({('2', 'donor3')}, ['N198', 'N199', 'N201', 'N202']),
                
            #     # ({('1', 'donor2')}, ['N198', 'N199', 'N201', 'N202']),
            #     # ({('3', 'donor2')}, ['N198', 'N199', 'N201', 'N202']),
            #     ({('unassigned', 'donor2')}, ['N198', 'N199', 'N201', 'N202']), 
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # NOTE: TL;DR: the 8319 ('unassigned', 'donor2') barcodes are probably N198, but we already have 20k barcodes of N198, for which our confidence is higher, so we give up on these. the rest seems to clearly match vireo's assignments.

                # souporcell thinks donor0 and donor3 are the same donor (2), and that donor2 is actually two donors (1 and 3). but if we assume donor0 doesnt contain multiple donors, than it seems likely that it is N199, as it doesn't match any other donor. donor3 matches N201 very nicely, which shouldn't be the case if donor3 is N199. also, donor0 clearly doesn't match N201. 
                # With regard to donor2, it matches N198 nicely, and it looks like its vireo VAFs are clustered at 0, 0.5 and 1, as expected. if donor2 were two donors, then we should have seen a different pattern (if we believe souporcell, then 1/3 of the cells belong to 3, 1/3 to 1, and 1/3 are unassigned or doublets). also, both ('1', 'donor2') and ('3', 'donor2') match N198 very nicely. the problem is that both also match somewhat nicely N199. still, it seems reasonable to assume they belong to N198.
                # 230307: With current data freqs, the distinction seems clear.
                ('2', 'donor0'),
                ('2', 'donor3'),
                ('1', 'donor2'),
                ('3', 'donor2'), # N198's log ratio is much higher than all others, so seems good enough to me.
                # ('unassigned', 'donor2'), 
            ],

            # 'rematch_anyway': True,
        },
        'demux_09_01_22_1': {
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_09_01_22_1/vaf_df_csvs/2_donor2_vaf_df.csv' # N205?
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_09_01_22_1/vaf_df_csvs/4_donor2_vaf_df.csv' # N205
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_09_01_22_1/vaf_df_csvs/unassigned_donor2_vaf_df.csv' # N205
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_09_01_22_1/vaf_df_csvs/0_donor3_vaf_df.csv' # N206

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('0', 'donor2')}, ['N203', 'N204', 'N205', 'N206']),
            #     # ({('0', 'donor3')}, ['N203', 'N204', 'N205', 'N206']),
            #     # ({('1', 'donor0')}, ['N203', 'N204', 'N205', 'N206']),
            #     # ({('3', 'donor0')}, ['N203', 'N204', 'N205', 'N206']),

            #     ({('unassigned', 'donor0')}, ['N203', 'N204', 'N205', 'N206']), # 5541 cells which i hope to rescue. Nope. no VAF scatter looked good enough. With freqs: N205: 17,0. as we already have almost 25k barcodes of N205 with less uncertainty, we can give up on those.
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('5', 'donor0'),
                ('unassigned', 'donor2'),
                ('4', 'donor2'),
                ('0', 'donor3'),
                ('1', 'donor4'),
                ('2', 'donor2'),
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_03_03_22_1': {
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor1')}, ['N241', 'N242', 'N243', 'N244']),
            #     ({('0', 'donor3')}, ['N241', 'N242', 'N243', 'N244']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # souporcell thinks donor1 and donor3 are the same donor. As N243 doesn't match donor1, and N241 doesn't match donor3, it seems unlikely that both are the same donor, so ignoring souporcell about that.
                ('0', 'donor1'), # with freqs, N241: 63,2
                ('0', 'donor3'), # with freqs, N243: 72,2
            ],
        },
        'demux_03_03_22_2': {
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('3', 'donor0')}, ['N241', 'N242', 'N243', 'N244']),
            #     ({('3', 'donor2')}, ['N241', 'N242', 'N243', 'N244']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # souporcell thinks donor0 and donor2 are the same donor. As N243 doesn't match donor0, and N241 doesn't match donor2, it seems unlikely that both are the same donor, so ignoring souporcell about that.
                ('3', 'donor0'), # with freqs, N241: 61,0
                ('3', 'donor2'), # with freqs, N243: 76,1
            ],
        },
        'demux_02_05_22_1': {
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('2', 'donor0')}, ['N245', 'N246', 'N247', 'N248']), # not really ambigious, but just in case. # with freqs: N248: 44,2. not great, but good enough.

            #     ({('0', 'donor1')}, ['N245', 'N246', 'N247', 'N248']),
            #     ({('1', 'donor1')}, ['N245', 'N246', 'N247', 'N248']),
            #     ({('3', 'donor1')}, ['N245', 'N246', 'N247', 'N248']),
                
            #     ({('unassigned', 'donor1')}, ['N245', 'N246', 'N247', 'N248']), # 8319 cells which i hope to rescue.
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # souporcell thinks donor1 is actually three donors (0,1,3). but donor1 matches N245 extremely well, and it looks like its vireo VAFs are clustered at 0, 0.5 and 1, as expected. also, each of ('0', 'donor1'), ('1', 'donor1'), ('3', 'donor1') matches N245 extremely well, and doesn't match other donors. thus, we ignore souporcell here.
                ('0', 'donor1'), # with freqs: N245: 713,7
                ('1', 'donor1'), # with freqs: N245: 695,3
                ('3', 'donor1'), # with freqs: N245: 708,6
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('unassigned', 'donor1', 'N245'), # homozygous matches,mismatches: 118,1. excellent. with freqs: N245: 608,1.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_02_05_22_2': {
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor0')}, ['N245', 'N246', 'N247', 'N248']),

            #     ({('1', 'donor1')}, ['N245', 'N246', 'N247', 'N248']),
            #     ({('2', 'donor1')}, ['N245', 'N246', 'N247', 'N248']),
            #     ({('3', 'donor1')}, ['N245', 'N246', 'N247', 'N248']),
                
            #     ({('unassigned', 'donor1')}, ['N245', 'N246', 'N247', 'N248']), # 11071 cells which i hope to rescue.
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # souporcell thinks donor1 is actually three donors (1,2,3). but donor1 matches N245 extremely well, and it looks like its vireo VAFs are clustered at 0, 0.5 and 1, as expected. also, each of ('1', 'donor1'), ('2', 'donor1'), ('3', 'donor1') matches N245 extremely well, and doesn't match other donors. thus, we ignore souporcell here.
                ('1', 'donor1'), # with freqs: N245: 698,7
                ('2', 'donor1'), # with freqs: N245: 711,15 (15 seems a bit high, but for other donors the match is extremely bad)
                ('3', 'donor1'), # with freqs: N245: 697,6
                ('0', 'donor0'), # ugh. with freqs: N248: 31,0 and N245: 27,1. yet in all other cases checked here, which pretty clearly belong to N245, it was very clear that N248 does not match. so it seems safe to guess this is indeed N248 here, otherwise, the match would have been very bad (e.g., 31,10).
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('unassigned', 'donor1', 'N245'), # homozygous matches,mismatches: 133,1. excellent. with freqs: N245: 665,4.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_09_06_22_1': {
            # ('N251', 'donor1', 'convincing VAF scatter'), 
            # ('N250', 'donor3', 'convincing VAF scatter'), 

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor3')}, ['N250', 'N251', 'N252', 'N255']),
            #     ({('1', 'donor3')}, ['N250', 'N251', 'N252', 'N255']), # ugh. with freqs: N252: 16,2, N255: 15,2
            #     ({('2', 'donor1')}, ['N250', 'N251', 'N252', 'N255']),
            #     ({('3', 'donor1')}, ['N250', 'N251', 'N252', 'N255']),
            #     # ({('1', 'donor0'), ('1', 'donor2'), ('1', 'donor3'), ('1', 'unassigned')}, ['N250', 'N251', 'N252', 'N255']),
                
            #     ({('unassigned', 'donor1')}, ['N250', 'N251', 'N252', 'N255']), # 4639 cells which i hope to rescue.
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # souporcell thinks donor3 is actually two donors (0,1). 
                ('0', 'donor3'), # ('0', 'donor3') matches very nicely N250, and the clustering around 0, 0.5, 1 looks much better than it looked for all donor3 cells, I think. 
                # souporcell thinks donor1 is actually two donors (2,3). 
                ('2', 'donor1'), # ('2', 'donor1') matches very nicely N251. with freqs: N251: 1175,9
                ('3', 'donor1'), # ('3', 'donor1') matches very nicely N251. with freqs: N251: 1182,6
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('unassigned', 'donor1', 'N251'), # homozygous matches,mismatches: 159,1. excellent. with freqs: N251: 772,4
            ],
            
            # ('1', 'donor3'), # ('1', 'donor3') seems to not match N250! i think it is the first time i see souporcell being better than vireo. ('1', 'donor3') also seems to not match N251 (though the number of SNPs is small). ('1', 'donor3') maybe maybe matches N252, and maybe maybe maybe N255, but the number of SNPs is very small.
            # {('1', 'donor0'), ('1', 'donor2'), ('1', 'donor3'), ('1', 'unassigned')} maybe matches N252, so maybe indeed this is N252. we shall see when we have more MIP data.
            # 'rematch_anyway': True,
        },
        'demux_17_08_20_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor3'),
                ('2', 'donor5'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_26_06_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('1', 'donor3'),
            #     ('2', 'donor5'),
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                # all 4 look convincing in soup_vireo_donor_match_summary.png. my code considers them as low confidence because vireo failed for demux_g_26_06_23_1...
                ('0', 'donor5', 'G4'),
                ('1', 'donor0', 'G2'),
                ('5', 'donor2', 'G3'),
                ('2', 'donor1', 'G1'),
            ],
            # 'rematch_anyway': True,
        },
        'G8_13_07_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/G8_13_07_23_1/vaf_df_csvs/0_donor0_vaf_df.csv', # G8? looks the same as 0,donor2 and 0,unassigned (pretty perfect for min_total_ref_and_alt_depth=20, but not perfect for min_total_ref_and_alt_depth=10)
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/G8_13_07_23_1/vaf_df_csvs/0_donor2_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/G8_13_07_23_1/vaf_df_csvs/0_unassigned_vaf_df.csv',
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],

            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'),
                ('0', 'donor2'),
            ],
            'manually_added_soup_vireo_mip_triplets': [ # TODO: 240521: ugh. still not enough. update when we have MIP genotyping data for G8
                ('0', 'unassigned', 'G8'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_24_07_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],

            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('1', 'donor3'),
            #     ('2', 'donor5'),
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'donor0', 'N387'),
                ('0', 'donor1', 'N387'),
                ('0', 'donor2', 'N387'),
                ('0', 'unassigned', 'N387'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_24_07_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],

            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'),
                ('0', 'donor1'),
                ('0', 'donor2'),
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N387'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_20_07_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/0_donor1_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/0_unassigned_vaf_df.csv', # i guess this is the same as 0,donor1, which i guess is N385, the only myelofibrosis patient in this experiment.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/7_unassigned_vaf_df.csv', # i guess this is not the same as 0,donor1.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/unassigned_donor1_vaf_df.csv', # i suspect this one is some mix of donors.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('1', 'donor3'),
            #     ('2', 'donor5'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_n_20_07_23_2': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/0_unassigned_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/0_donor1_vaf_df.csv', 
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/7_unassigned_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_1/vaf_df_csvs/unassigned_donor1_vaf_df.csv', 
            
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_2/vaf_df_csvs/4_donor0_vaf_df.csv', # i guess this is N385, the only myelofibrosis patient in the experiment.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_2/vaf_df_csvs/4_unassigned_vaf_df.csv', # i guess this is the same as 4,donor0
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_2/vaf_df_csvs/6_unassigned_vaf_df.csv', # seems the same as demux_n_20_07_23_1,7,unassigned
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_20_07_23_2/vaf_df_csvs/unassigned_donor0_vaf_df.csv', # similar to N385, but i guess some mix??

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('1', 'donor3'),
            #     ('2', 'donor5'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_09_06_22_2': {
            # ('N250', 'donor0', 'convincing VAF scatter'), 
            # ('N251', 'donor1', 'convincing VAF scatter'), 

            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('2', 'donor0')}, ['N250', 'N251', 'N252', 'N255']),
                
            #     ({('0', 'donor1')}, ['N250', 'N251', 'N252', 'N255']),
            #     ({('1', 'donor1')}, ['N250', 'N251', 'N252', 'N255']),

            #     ({('unassigned', 'donor1')}, ['N250', 'N251', 'N252', 'N255']), # 4412 cells which i hope to rescue.
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # souporcell thinks donor1 is actually two donors (0,1). 
                ('0', 'donor1'), # ('0', 'donor1') matches very nicely N251. with freqs: N251: 1145,8.
                ('1', 'donor1'), # ('1', 'donor1') matches very nicely N251. with freqs: N251: 1130,5.
                ('2', 'donor0'), # with freqs: N250: 643,6
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('unassigned', 'donor1', 'N251'), # homozygous matches,mismatches: 175,7. not excellent, but good enough. with freqs: N251: 902,4.
            ],
        },
        'demux_13_06_22_1': {
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('2', 'donor1'), ('2', 'donor2'), ('2', 'donor3')}, ['N261', 'N262', 'N263', 'N264']),
            #     # ({('2', 'donor1'), ('2', 'donor3')}, ['N261', 'N262', 'N263', 'N264']), # not enough cells here
            #     ({('2', 'donor2')}, ['N261', 'N262', 'N263', 'N264']),
            #     ({('0', 'donor0')}, ['N261', 'N262', 'N263', 'N264']),
            #     ({('1', 'donor0')}, ['N261', 'N262', 'N263', 'N264']),
            #     ({('3', 'donor0')}, ['N261', 'N262', 'N263', 'N264']),
              
            #     ({('unassigned', 'donor0')}, ['N261', 'N262', 'N263', 'N264']), # 7993 cells which i hope to rescue.
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # souporcell thinks donor1, donor2, and donor3 are the same donor (2).
                # {('2', 'donor1'), ('2', 'donor2'), ('2', 'donor3')} matched N261 worse, compared to how donor2 matched N261. specifically, less SNPs were considered, even though there were a bit more cells for {('2', 'donor1'), ('2', 'donor2'), ('2', 'donor3')}. i guess this is because adding ('2', 'donor1') and ('2', 'donor3') "diluted" the VAFs of donor2, and made them less likely to seem homozygous (<=0.15 or >=0.85). in other words, this result seems to me to suggest donor1, donor2, and donor3 are not the same donor (i.e., souporcell was wrong, again).
                ('2', 'donor2'), # with freqs: N261: 86,0
                # souporcell thinks that donor0 is actually three donors (0,1,3). however, donor0 vireo VAFs are clustered very nicely around 0, 0.5 and 1 (can see this when looking also at donors donor0 doesnt match).
                ('0', 'donor0'), # ('0', 'donor0') matches N262 nicely, and doesn't match all other donors. also, has nice clusters around 0,0.5,1, as expected. with freqs: N262: 714,13
                ('1', 'donor0'), # ('1', 'donor0') matches N262 very nicely, and doesn't match all other donors. also, has nice clusters around 0,0.5,1, as expected. with freqs: N262: 708,5
                ('3', 'donor0'), # ('3', 'donor0') matches N262 very nicely, and doesn't match all other donors. also, has nice clusters around 0,0.5,1, as expected. with freqs: N262: 716,4
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('unassigned', 'donor0', 'N262'), # homozygous matches,mismatches: 106,0. perfect. with freqs: N262: 587,1
            # ],
            # 'rematch_anyway': True,
        },
        'demux_03_07_22_1': {
            # almost 7k doublets here, but around 5.4k are doublet according to both soup and vireo, or doubelt according to one and unassigned according to the other. i.e., there is no reasonable soup,vireo cluster that one thinks is a doublet and the other thinks is a single donor.
    
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1/vaf_df_csvs/3_donor4_vaf_df.csv', # N272
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1/vaf_df_csvs/3_unassigned_vaf_df.csv', # N272? seems like yes
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1/vaf_df_csvs/0_unassigned_vaf_df.csv',
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'unassigned')}, ['N269', 'N270', 'N271', 'N272']), # 567 cells which i hope to rescue. Nope. no VAF scatter looked good enough. (not enough data) 
            #     ({('1', 'unassigned')}, ['N269', 'N270', 'N271', 'N272']), # 337 cells which i hope to rescue. Nope. no VAF scatter looked good enough. (not enough data)
            #     ({('2', 'unassigned')}, ['N269', 'N270', 'N271', 'N272']), # 410 cells which i hope to rescue. Nope. no VAF scatter looked good enough. (not enough data)
            #     ({('3', 'unassigned')}, ['N269', 'N270', 'N271', 'N272']), # 462 cells which i hope to rescue. Nope. no VAF scatter looked good enough. (not enough data)
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor5'), # it was considered bad by the auto examination due to low min_agreement_proportion. but it looks fine.
                ('1', 'donor2'), # it was considered bad by the auto examination due to low min_agreement_proportion. but it looks fine.
                ('2', 'donor3'), # it was considered bad by the auto examination due to low min_agreement_proportion. but it looks fine.
                ('3', 'donor4'), # it was considered bad by the auto examination due to low min_agreement_proportion. but it looks fine.
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('3', 'unassigned', 'N272'),
            ],
        },
        'demux_03_07_22_2': {
            # around 5.2k doublets here, but around 4.5k are doublet according to both soup and vireo, or doubelt according to one and unassigned according to the other. i.e., there is no reasonable soup,vireo cluster that one thinks is a doublet and the other thinks is a single donor.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor2'), # it was considered bad by the auto examination due to low min_agreement_proportion. but it looks fine.
                ('2', 'donor1'), # it was considered bad by the auto examination due to low min_agreement_proportion. but it looks fine.
            ],
            # 'rematch_anyway': True,
        },
        'demux_02_05_22_1_ultima': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('2', 'donor0')}, ['N245', 'N246', 'N247', 'N248', 'N265', 'N266', 'N267', 'N268']), 
            #     ({('7', 'donor3')}, ['N245', 'N246', 'N247', 'N248', 'N265', 'N266', 'N267', 'N268']), 
            #     ({('3', 'donor6')}, ['N245', 'N246', 'N247', 'N248', 'N265', 'N266', 'N267', 'N268']), 
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('unassigned', 'donor7'), 
            #     ('0', 'donor7'), 
            #     ('1', 'donor7'), 
            #     ('4', 'donor7'), 
            #     ('6', 'donor7'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('5', 'unassigned'),
            #     ('2', 'donor0'), # reject pairing to any expected donor.
            #     ('7', 'donor3'), # reject pairing to any expected donor.
            #     ('3', 'donor6'), # reject pairing to any expected donor.
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('2', 'donor0', 'N266'), # what.
            #     ('7', 'donor3', 'N267'), # what.
            #     ('3', 'donor6', 'N268'), # what.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_02_05_22_2_ultima': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor3')}, ['N245', 'N246', 'N247', 'N248']), 
            #     ({('1', 'donor3')}, ['N245', 'N246', 'N247', 'N248']), 
            #     ({('2', 'donor3')}, ['N245', 'N246', 'N247', 'N248']), 
            #     ({('unassigned', 'donor3')}, ['N245', 'N246', 'N247', 'N248']), 
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('2', 'donor3'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_09_06_22_1_ultima': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # #     {('4', 'donor8')}
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('8', 'donor8')}, ['N250', 'N251', 'N252', 'N255', 'N273', 'N274', 'N275', 'N276']), 
            #     ({('4', 'donor8')}, ['N250', 'N251', 'N252', 'N255', 'N273', 'N274', 'N275', 'N276']), 
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('8', 'donor8'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     # ('0', 'donor4'), # this is pretty much the same as donor4. the log ratio is much lower than we usually see, and the VAF scatter doesnt look good.
            #     ('5', 'unassigned'), 
            #     ('3', 'donor3'), # reject pairing to any expected donor.
            #     ('2', 'donor4'), # reject pairing to any expected donor.
            #     ('4', 'donor7'), # reject pairing to any expected donor.
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('4', 'donor8', 'N276'), # ugh. seems like vireo merged two different donors (N276 and N250)
            ],
            # 'rematch_anyway': True,
        },
        'demux_09_06_22_2_ultima': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('3', 'donor5')}
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('3', 'donor5')}, ['N250', 'N251', 'N252', 'N255']), 
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('3', 'donor3'), # this is pretty much the same as donor3 - only 82 cells less (out of ~2893).
            #     ('1', 'donor2'), 
            #     ('2', 'donor2'), 
            #     ('unassigned', 'donor2'), 
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     # ('3', 'donor5', 'N252'), # too borderline. only 328 cells, so giving up on them.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_13_06_22_1_ultima': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_13_06_22_1_ultima/vaf_df_csvs/1_donor1_vaf_df.csv', # NS20? definitely the same as /dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/NS20_29_09_21_1/vaf_df_csvs/0_unassigned_vaf_df.csv, which is very probably NS20.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/NS20_29_09_21_1/vaf_df_csvs/0_unassigned_vaf_df.csv', # NS20?

            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_13_06_22_1_ultima/vaf_df_csvs/1_donor1_vaf_df.csv' # NS20?
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/NS20_29_09_21_1/vaf_df_csvs/0_unassigned_vaf_df.csv' # NS20?

            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_13_06_22_1_ultima/vaf_df_csvs/0_unassigned_vaf_df.csv' # N245?
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_13_06_22_1_ultima/vaf_df_csvs/3_donor4_vaf_df.csv' # N261
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_13_06_22_1_ultima/vaf_df_csvs/0_donor4_vaf_df.csv' # N245
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor4')}, ['N261', 'N262', 'N263', 'N264', 'N245', 'N246', 'N247', 'N248']), 
            #     ({('0', 'unassigned')}, ['N261', 'N262', 'N263', 'N264', 'N245', 'N246', 'N247', 'N248']), 
            #     ({('3', 'donor4')}, ['N261', 'N262', 'N263', 'N264', 'N245', 'N246', 'N247', 'N248']), 
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('0', 'donor1'), # this is pretty much the same as donor1
            #     ('2', 'donor3'),
            #     ('3', 'donor3'),
            #     ('unassigned', 'donor3'),
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('0', 'donor4'), # reject only to add in a triplet.
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'donor4', 'N245'), # what.
                # ('0', 'unassigned', 'N245'), # what. was about unsure about it, but compared to '0', 'donor4' and they seem the same, so ok.
                ('3', 'donor4', 'N261'), 
                # ('1', 'donor1', 'NS20'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_19_06_22_2_ultima': {
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_2/vaf_df_csvs/0_unassigned_vaf_df.csv' # N235
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_19_06_22_2_ultima/vaf_df_csvs/7_unassigned_vaf_df.csv' # N235?? seems like clearly yes
            # # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_19_06_22_2_ultima/vaf_df_csvs/6_donor3_vaf_df.csv' # not N235
            # # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_2/vaf_df_csvs/2_unassigned_vaf_df.csv' # N235

            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_19_06_22_2_ultima/vaf_df_csvs/3_unassigned_vaf_df.csv', # N254? i guess yes, but too noisy to my liking.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_06_06_22_1_ultima/vaf_df_csvs/1_donor3_vaf_df.csv', # N254
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('7', 'unassigned')}, ['N235', 'N249', 'N253', 'N254', 'N257', 'N265', 'N266', 'N267', 'N268']), 
            #     # ({('4', 'donor2')}, ['N265', 'N266', 'N267', 'N268', 'N249', 'N253', 'N254', 'N257']), 
            #     # ({('0', 'donor2')}, ['N265', 'N266', 'N267', 'N268', 'N249', 'N253', 'N254', 'N257']), 
            #     # ({('5', 'donor3')}, ['N265', 'N266', 'N267', 'N268', 'N249', 'N253', 'N254', 'N257']), 

            #     # ({('3', 'donor0')}, ['N265', 'N266', 'N267', 'N268']), 
            #     # ({('3', 'unassigned')}, ['N265', 'N266', 'N267', 'N268']), 
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('9', 'donor1'), 
                ('6', 'donor3'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('7', 'donor1'), # reject pairing to any expected donor.
            #     ('4', 'donor2'), # reject pairing to any expected donor.
            #     ('5', 'donor3'), # reject pairing to any expected donor.
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('7', 'unassigned', 'N235'), # what.
                # ('3', 'unassigned', 'N254'), # when comparing vaf_dfs, it was too noisy to my liking.
            ],
            # 'rematch_anyway': True,
        },
        'demux_03_07_22_1_ultima': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1_ultima/vaf_df_csvs/2_donor0_vaf_df.csv', # N269
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1_ultima/vaf_df_csvs/1_unassigned_vaf_df.csv', # N269? very bad.
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # #     # {('2', 'unassigned'), ('2', 'donor1')},
            # #     # {('2', 'unassigned'), ('2', 'donor1'), ('2', 'donor3')},
            # #     # {('2', 'donor3')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor3')}, ['N269', 'N270', 'N271', 'N272']), 
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('4', 'donor4'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('0', 'unassigned'),
            #     ('3', 'unassigned'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_03_07_22_2_ultima': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1_ultima/vaf_df_csvs/2_donor0_vaf_df.csv', # N269
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1_ultima/vaf_df_csvs/1_unassigned_vaf_df.csv', # N269? very bad.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_2_ultima/vaf_df_csvs/0_donor5_vaf_df.csv', # N269
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_2_ultima/vaf_df_csvs/4_unassigned_vaf_df.csv', # N269? very bad. but why so similar to /dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1_ultima/vaf_df_csvs/1_unassigned_vaf_df.csv.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # #     # {('2', 'unassigned'), ('2', 'donor1')},
            # #     # {('2', 'unassigned'), ('2', 'donor1'), ('2', 'donor3')},
            # #     # {('2', 'donor3')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'donor5')}, ['N269', 'N270', 'N271', 'N272']), 
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor5'),
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('4', 'unassigned'), # 
            # ],
            # 'rematch_anyway': True,
        },
        'demux_10_07_22_2_ultima': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('7', 'donor10')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('4', 'unassigned')}, ['N275', 'N276', 'N273', 'N274', 'N286', 'N287', 'N288', 'N289', 'N290']), 
            #     ({('0', 'unassigned')}, ['N275', 'N276', 'N273', 'N274', 'N286', 'N287', 'N288', 'N289', 'N290']), 
            #     ({('7', 'unassigned')}, ['N275', 'N276', 'N273', 'N274', 'N286', 'N287', 'N288', 'N289', 'N290']), 
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor9'),
                ('7', 'donor1'), 
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('7', 'donor10'), # too few cells, so can't be sure
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('4', 'unassigned', 'N290'),
                ('0', 'unassigned', 'N288'),
                ('7', 'unassigned', 'N286'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_10_07_22_1_ultima': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # #     {('7', 'donor2')},
            # #     {('7', 'donor4')},
            # #     {('9', 'donor0')},
            # #     {('9', 'donor5')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('7', 'donor7')}, ['N273', 'N274', 'N275', 'N276', 'N250', 'N251', 'N252', 'N255', 'N101']), 
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('5', 'donor3'), 
            #     ('6', 'donor6'), 
            #     ('2', 'donor9'), 
            #     ('1', 'donor8'), 
            #     ('8', 'donor1'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('3', 'donor7'), # reject pairing to any expected donor.
            #     ('5', 'donor7'), # reject pairing to any expected donor.
            #     ('unassigned', 'donor7'), # reject pairing to any expected donor.
            #     ('0', 'unassigned'), 
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('7', 'donor7', 'N250'), # looking at the VAF scatters, looks better for N250, and N250 is more expected than N101, so i think it is good enough.
            ],
            # 'rematch_anyway': True,
        },
        'demux_01_02_21_2': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('2', 'donor2'), ('4', 'donor2'), ('unassigned', 'donor2'), ('2', 'donor3'), ('4', 'donor3'), ('unassigned', 'donor3')},
            # ],

            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('2', 'donor2'), ('4', 'donor2'), ('unassigned', 'donor2'), ('2', 'donor3'), ('4', 'donor3'), ('unassigned', 'donor3')}, ['N273', 'N130','N133','N135','N232']),
            #     # ({('6', 'donor1'), ('unassigned', 'donor1'), ('doublet', 'donor1')}, ['N273', 'N274', 'N275', 'N276', 'N250', 'N251', 'N252', 'N255']), # just to see whether we get to the number of SNPs MIP-vireo has. didn't get to that number of SNPs, I guess due to soup's default MIN_ALT and MIN_REF arguments.
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('6', 'donor4'), 
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('2', 'unassigned'), # matches both N133 and N135 nicely.
                ('4', 'unassigned'), # matches both N130 and N135 nicely.
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('3', 'donor7', 'N251'), # what.
            #     ('5', 'donor7', 'N251'), # what.
            #     ('unassigned', 'donor7', 'N251'), # what.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_01_02_21_2_illumina': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('2', 'donor2'), ('4', 'donor2'), ('unassigned', 'donor2'), ('2', 'donor3'), ('4', 'donor3'), ('unassigned', 'donor3')},
            # ],

            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('4', 'donor5'), 
                ('5', 'donor1'), 
                ('0', 'donor4'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('2', 'unassigned'), # matches both N133 and N135 nicely.
            #     ('4', 'unassigned'), # matches both N130 and N135 nicely.
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('3', 'donor7', 'N251'), # what.
            #     ('5', 'donor7', 'N251'), # what.
            #     ('unassigned', 'donor7', 'N251'), # what.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_22_02_21_1': {
            

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('6', 'donor1')}, 
            # ],

            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('6', 'donor1'), ('unassigned', 'donor1'), ('doublet', 'donor1')}, ['N273', 'N274', 'N275', 'N276', 'N250', 'N251', 'N252', 'N255']), # just to see whether we get to the number of SNPs MIP-vireo has. didn't get to that number of SNPs, I guess due to soup's default MIN_ALT and MIN_REF arguments.
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('6', 'donor1'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('0', 'unassigned'), 
            #     ('2', 'unassigned'), 
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('3', 'donor7', 'N251'), # what.
            #     ('5', 'donor7', 'N251'), # what.
            #     ('unassigned', 'donor7', 'N251'), # what.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_22_02_21_2': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('6', 'donor1'), ('unassigned', 'donor1'), ('doublet', 'donor1')}, ['N273', 'N274', 'N275', 'N276', 'N250', 'N251', 'N252', 'N255']), # just to see whether we get to the number of SNPs MIP-vireo has. didn't get to that number of SNPs, I guess due to soup's default MIN_ALT and MIN_REF arguments.
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('6', 'donor4'), 
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('0', 'unassigned'), 
                ('2', 'unassigned'), 
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('3', 'donor7', 'N251'), # what.
            #     ('5', 'donor7', 'N251'), # what.
            #     ('unassigned', 'donor7', 'N251'), # what.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_22_02_21_2_illumina': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('0', 'donor6')},
            #     # {('unassigned', 'donor6')},
            #     {('0', 'donor6'), ('unassigned', 'donor6')},
            #     # {('0', 'unassigned')},
            #     # {('2', 'unassigned')},
            #     # {('4', 'unassigned')},
            #     # {('6', 'unassigned')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor5'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('0', 'unassigned'), 
            #     ('2', 'unassigned'), 
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('3', 'donor7', 'N251'), # what.
            #     ('5', 'donor7', 'N251'), # what.
            #     ('unassigned', 'donor7', 'N251'), # what.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_11_04_21_2': {
            # 1,unassigned doesnt look like N177, nor N184.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_11_04_21_2/vaf_df_csvs/4_donor6_vaf_df.csv', # N177
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_11_04_21_2/vaf_df_csvs/1_unassigned_vaf_df.csv', # N177? N184?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_11_04_21_2/vaf_df_csvs/0_donor2_vaf_df.csv', # N184
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('3', 'unassigned')}, ['N177', 'N178', 'N184', 'N200', 'N210']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor2'), 
                ('5', 'donor4'), 
                ('6', 'donor3'), 
                ('2', 'donor1'), 
                ('4', 'donor6'), 
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('1', 'unassigned'), # matches both N177 and N184 nicely.
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('3', 'unassigned', 'N184'), # NOTE: removed only on 230831. ugh.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_11_04_21_2_illumina': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],

            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('3', 'unassigned')}, ['N177', 'N178', 'N184', 'N200', 'N210']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('6', 'donor4'), 
                ('0', 'donor2'), 
                ('3', 'donor1'), 
                ('5', 'donor0'), 
                ('5', 'donor3'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('1', 'unassigned'), # matches both N177 and N184 nicely.
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('3', 'unassigned', 'N184'),
            # ],
            # 'rematch_anyway': True,
        },
        'N328_15_12_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('1', 'donor0'), 
                ('2', 'donor0'), 
                ('unassigned', 'donor0'), 
                ('0', 'donor1'), 
                ('1', 'donor1'), 
                ('2', 'donor1'), 
                ('unassigned', 'donor1'), 
                ('0', 'donor2'), 
                ('1', 'donor2'), 
                ('2', 'donor2'), 
                ('unassigned', 'donor2'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('unassigned', 'donor7'), # reject pairing to any expected donor.
            #     ('0', 'unassigned'), 
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N328'),
                ('1', 'unassigned', 'N328'),
                ('2', 'unassigned', 'N328'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_30_11_20_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('0', 'donor0'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('1', 'donor0'),
            #     ('2', 'donor1'),
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('2', 'unassigned', 'N328'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_01_01_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('5', 'donor0')}, ['N1', 'N196', 'N240', 'N323', 'N325']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('3', 'donor6'), 
            #     ('6', 'donor4'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('2', 'unassigned', 'N328'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_15_01_23_1': { # ugh. souporcell failed with an error.
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor2'), 
                ('0', 'donor3'), 
                ('0', 'donor4'), 
                ('0', 'donor5'), 
                ('0', 'donor6'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('2', 'unassigned', 'N328'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_05_02_23_1': {
            # wrote at some point, with regard to ('1', 'unassigned'), ('unassigned', 'donor6'), but later thought it might be dangerous to exclude these cells. maybe they are of a different type, thus expressing different genes and not clustering together? so in the end decided i should use them.
            # # both are probably N198, but vaf_dfs showed them to be somewhat noisy, and we already have enough cells for N198, so better be safe and discard them.

            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_05_02_23_1/vaf_df_csvs/6_donor6_vaf_df.csv', # N198
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_05_02_23_1/vaf_df_csvs/6_unassigned_vaf_df.csv', # N198? seems like yes
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_05_02_23_1/vaf_df_csvs/1_unassigned_vaf_df.csv', # N198? not good enough.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_05_02_23_1/vaf_df_csvs/unassigned_donor6_vaf_df.csv', # N198? ok? but a bit noisy...

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('6', 'donor6')}, ['N276','N198','N342','N343','N344','N31','N40','N53','N55','N56','N57','N71','N72']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('6', 'donor6'),
                # ('unassigned', 'donor6'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     # ('6', 'donor6'), # considered rejecting, because we already have here ~8k cells for N198.
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('1', 'unassigned', 'N198'),
            #     ('6', 'unassigned', 'N198'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_27_02_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_27_02_23_1/vaf_df_csvs/4_unassigned_vaf_df.csv', # N350? seems like yes.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_08_02_23_1/vaf_df_csvs/3_donor1_vaf_df.csv', # N350
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_08_02_23_1/vaf_df_csvs/4_donor7_vaf_df.csv', # N353
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('5', 'donor4')}, ['N354', 'N355', 'N356', 'N357', 'N260', 'N359', 'N360','N307','N362','N363']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('5', 'donor6'), 
                ('5', 'donor4'), 
                ('4', 'donor7'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('7', 'donor3'), # N225?? though this really doesnt make sense because N225 did not dominate their experiment... vaf_dfs didnt match.
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('4', 'unassigned', 'N350'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_27_02_23_2': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('5', 'donor4')}, ['N354', 'N355', 'N356', 'N357', 'N260', 'N359', 'N360','N307','N362','N363']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('0', 'donor0'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('1', 'donor1'),
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_28_02_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('5', 'donor4')}, ['N354', 'N355', 'N356', 'N357', 'N260', 'N359', 'N360','N307','N362','N363']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('0', 'donor0'), 
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('7', 'donor0'), # only 5 cells
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('0', 'donor0', 'N260'), # what (?).
            # ],
            # 'rematch_anyway': True,
        },
        'demux_26_03_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('0', 'donor0'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('0', 'donor0', 'N260'),
            # ],
            # 'rematch_anyway': True,
        },
        'T6_25_05_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('0', 'donor0'), 
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('0', 'donor0'), # just to avoid raising an exception. no MIP data so can't tell whether these cells belong to T6 or not (but seems most likely that they do belong to her).
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('0', 'donor0', 'N260'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_29_01_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_29_01_23_1/vaf_df_csvs/1_donor6_vaf_df.csv', # N335
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_29_01_23_1/vaf_df_csvs/1_unassigned_vaf_df.csv', # N335? i guess yes, though it is a bit noisy
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_29_01_23_1/vaf_df_csvs/0_donor2_vaf_df.csv', # N336
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_29_01_23_1/vaf_df_csvs/0_unassigned_vaf_df.csv', # N336? i guess yes, though it is a bit noisy
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_29_01_23_1/vaf_df_csvs/6_donor4_vaf_df.csv', # N337
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_29_01_23_1/vaf_df_csvs/6_unassigned_vaf_df.csv', # N337? i guess yes, though it is a bit noisy
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('5', 'donor4')}, ['N354', 'N355', 'N356', 'N357', 'N260', 'N359', 'N360','N307','N362','N363']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor2'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('1', 'unassigned', 'N335'),
            #     # ('0', 'unassigned', 'N336'), 
            #     ('6', 'unassigned', 'N337'),
            # ],
            # 'rematch_anyway': True,
        },
        'demux_08_02_23_1': {
            
            
            
            


            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('5', 'donor4')}, ['N354', 'N355', 'N356', 'N357', 'N260', 'N359', 'N360','N307','N362','N363']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('3', 'donor1'),
            ],
            
            # ('N350', 'donor0', 'convincing VAF scatter'),
            # ('N328', 'donor1', 'convincing VAF scatter'),
            # ('N351', 'donor2', 'convincing VAF scatter'),
            # ('N352', 'donor3', 'convincing VAF scatter'),
            # # donor4 is only ~800 cells, so i guess we give it a benefit of the doubt and don't call it a contamination.
            # ('N353', 'donor5', 'convincing VAF scatter'),
            
            'manually_rejected_soup_vireo_pairs': [
                ('0', 'donor5'), 
                ('3', 'donor3'), 
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N276'),
            ],
            # 'rematch_anyway': True,
        },
        'N279_26_10_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N279', 'N12']),
            # #     ({('2', 'unassigned')}, ['N279', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('2', 'donor0'), 
                ('2', 'donor1'), 
                ('2', 'donor2'), 
            ],
            # # 'manually_rejected_soup_vireo_pairs': [
            # #     ('4', 'unassigned'), 
            # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('2', 'unassigned', 'N279'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_06_11_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('unassigned', 'donor0')}, ['N277', 'N280', 'N281']),
            #     ({('0', 'donor2')}, ['N277', 'N280', 'N281']),
            #     ({('2', 'unassigned')}, ['N277', 'N280', 'N281']),
            #     ({('1', 'unassigned')}, ['N277', 'N280', 'N281']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('0', 'donor2'), # log ratio is a bit too low, but good enough, i think. especially if you look at vireo-mip for donor2, while ('0', 'donor2') is 2851 out of the ~3750 cells making donor2. unfortunately, we don't have N277 in another experiment so we can't compare vaf_dfs.
            #     ('2', 'donor0'), 
            #     ('unassigned', 'donor0'), 
            #     ('1', 'donor1'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('3', 'unassigned'), # fits kind of ok N19, but it doesn't sound plausible. doesn't fit anyone else well.
            #     ('4', 'unassigned'), # ugh. very similar log ratio for N280 and N281.
            # ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('2', 'unassigned', 'N280'),
            #     ('1', 'unassigned', 'N281'),
            # ],
            # 'rematch_anyway': True,
        },
        'N251_04_12_22_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_04_12_22_1/vaf_df_csvs/1_donor0_vaf_df.csv', # N234? N12? not N234. also, not good enough for N12.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_04_12_22_1/vaf_df_csvs/1_unassigned_vaf_df.csv', # N234? N12? not N234. anyway, looks like it is the same one as 1_donor0. doesnt look like N12.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_04_12_22_1/vaf_df_csvs/6_donor6_vaf_df.csv', # N12
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_14_02_22_1/vaf_df_csvs/2_donor3_vaf_df.csv', # N234

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('1', 'unassigned'), ('1', 'donor0')}, # ah?? N234??
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('unassigned', 'donor5'), 
                ('0', 'donor5'), 
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('1', 'donor0'), # 230711: these are 9.4k cells (!) is it N234? N12? a mix? none?
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('1', 'unassigned', 'N234'), # 230709: is it N234?? anyway, i guess we discard all accidents, so it doesnt matter much.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_13_11_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_06_06_22_1_ultima': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_13_11_22_2': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('unassigned', 'unassigned')}, # 230612: didnt fit any donor well
            #     # {('4', 'donor3'), ('0', 'donor3'), ('1', 'donor3'), ('2', 'donor3'), ('unassigned', 'donor3')},
            #     # {('0', 'donor5'), ('1', 'donor5'), ('2', 'donor5'), ('3', 'donor5'), ('4', 'donor5'), ('unassigned', 'donor5')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('5', 'unassigned')}, ['N282', 'N283', 'N284', 'N285', 'N250', 'N251', 'N252', 'N255', 'N12']),
            #     ({('4', 'unassigned')}, ['N282', 'N283', 'N284', 'N285', 'N250', 'N251', 'N252', 'N255', 'N12']),
            #     ({('3', 'unassigned')}, ['N282', 'N283', 'N284', 'N285', 'N250', 'N251', 'N252', 'N255', 'N12']),
            #     ({('0', 'donor4')}, ['N282', 'N283', 'N284', 'N285', 'N250', 'N251', 'N252', 'N255', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('5', 'donor7'),
                ('1', 'donor5'),
                ('2', 'donor2'),
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('2', 'unassigned'), # when compared to another soup-vireo donor that was clearly N283, the agreement did not look good enough.
            #     ('0', 'unassigned'), # when compared to another soup-vireo donor that was clearly N284, the agreement did not look good enough.
            # ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor1',
            #     'donor3', 
            # },
            'manually_added_soup_vireo_mip_triplets': [
                ('4', 'unassigned', 'N271'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_14_11_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('0', 'donor5'), ('1', 'donor5'), ('2', 'donor5'), ('3', 'donor5'), ('4', 'donor5'), ('unassigned', 'donor5')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('0', 'unassigned')}, ['N286', 'N287', 'N288', 'N289', 'N290']),
            #     ({('1', 'unassigned')}, ['N286', 'N287', 'N288', 'N289', 'N290']),
            #     ({('3', 'unassigned')}, ['N286', 'N287', 'N288', 'N289', 'N290']),
            #     ({('2', 'unassigned')}, ['N286', 'N287', 'N288', 'N289', 'N290']),
            #     ({('4', 'unassigned')}, ['N286', 'N287', 'N288', 'N289', 'N290']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('2', 'donor0'),
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('3', 'unassigned'),
                ('5', 'unassigned'),
            ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor3', 
            #     'donor5',
            # },
            # 'manually_added_soup_vireo_mip_triplets': [
            #     # ('3', 'unassigned', 'N286'), 
            # ],
            # 'rematch_anyway': True,
        },
        'demux_23_11_22_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_23_11_22_1/vaf_df_csvs/5_donor6_vaf_df.csv', # N295
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_23_11_22_1/vaf_df_csvs/5_unassigned_vaf_df.csv', # N295? not enough data

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('0', 'donor5'), ('1', 'donor5'), ('2', 'donor5'), ('3', 'donor5'), ('4', 'donor5'), ('unassigned', 'donor5')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('5', 'unassigned')}, ['N291', 'N292', 'N293', 'N294', 'N295']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('3', 'donor4'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('6', 'unassigned'), 
            # ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor3', 
            #     'donor5',
            # },
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('5', 'unassigned', 'N295'), # not enough data to feel confident (when compared to 5,donor6, which is clearly N295). anyway it isn't many cells
            # ],
            # 'rematch_anyway': True,
        },
        'NS20_29_09_21_1': {
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/NS20_29_09_21_1/vaf_df_csvs/0_donor0_vaf_df.csv' # NS20?
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/NS20_29_09_21_1/vaf_df_csvs/0_donor1_vaf_df.csv' # NS20?
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/NS20_29_09_21_1/vaf_df_csvs/0_donor2_vaf_df.csv' # NS20?
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/NS20_29_09_21_1/vaf_df_csvs/0_unassigned_vaf_df.csv' # NS20?

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('0', 'donor5'), ('1', 'donor5'), ('2', 'donor5'), ('3', 'donor5'), ('4', 'donor5'), ('unassigned', 'donor5')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('3', 'donor4'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('6', 'unassigned'), 
            # ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor3', 
            #     'donor5',
            # },
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'donor0', 'NS20'), 
                ('0', 'donor1', 'NS20'), 
                ('0', 'donor2', 'NS20'), 
                ('0', 'unassigned', 'NS20'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_04_12_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('6', 'donor6')},
            #     # {('0', 'donor5'), ('1', 'donor5'), ('2', 'donor5'), ('3', 'donor5'), ('4', 'donor5'), ('unassigned', 'donor5')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('6', 'donor6'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('6', 'unassigned'), 
            # ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor3', 
            #     'donor5',
            # },
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('6', 'donor6', 'N12'), 
            # ],
            # 'rematch_anyway': True,
        },
        'demux_08_12_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # #     {('10', 'donor6')},
            # #     # {('10', 'donor14')},
            # #     # {('10', 'donor8')},
            #     {('10', 'donor6')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('10', 'donor6')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            #     # ({('14', 'donor5')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
                
            #     # ({('10', 'donor0')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            #     # ({('10', 'donor9')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            #     ({('10', 'unassigned')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('12', 'donor8'), 
                ('9', 'donor1'), 
                ('11', 'donor7'), 
                ('14', 'donor2'), 
                ('8', 'donor3'), 
                ('10', 'donor6'),
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('10', 'donor0'), # i guess it is N282, but not convincing enough. also, match with N12 is too good, but N12 is weird...
                ('10', 'donor9'), # i guess it is N282, but not convincing enough. also, match with N12 is too good, but N12 is weird...
                ('10', 'unassigned'), # i guess it is N282, but not convincing enough. also, match with N12 is too good, but N12 is weird...
                ('14', 'donor5'), # 230710: the log ratio does not look good for any donor currently
            ],
            'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
                'donor5',
            },
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('6', 'donor0', 'N285'), 
            #     ('6', 'donor0', 'N285'), 
            # ],
            # 'rematch_anyway': True,
        },
        'demux_13_12_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('10', 'donor6')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('10', 'donor6')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            #     # ({('14', 'donor5')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
                
            #     # ({('10', 'donor0')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            #     # ({('10', 'donor9')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            #     ({('10', 'unassigned')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('12', 'donor8'), 
            #     ('9', 'donor1'), 
            #     ('11', 'donor7'), 
            #     ('14', 'donor2'), 
            #     ('14', 'donor5'), # the log ratio is not too bad, but for 1.5k cells it's a bit low. anyway, when looking at vireo-mip matching, donor5 looks pretty good, so i think it is good enough.
            #     ('10', 'donor6'), 
            #     ('8', 'donor3'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('10', 'donor0'), # i guess it is N282, but not convincing enough. also, match with N12 is too good, but N12 is weird...
            #     ('10', 'donor9'), # i guess it is N282, but not convincing enough. also, match with N12 is too good, but N12 is weird...
            #     ('10', 'unassigned'), # i guess it is N282, but not convincing enough. also, match with N12 is too good, but N12 is weird...
            # ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor3', 
            #     'donor5',
            # },
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('6', 'donor0', 'N285'), 
            #     ('6', 'donor0', 'N285'), 
            # ],
            # 'rematch_anyway': True,
        },
        'demux_23_11_20_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('10', 'donor6')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('10', 'donor6')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('12', 'donor8'), 
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('2', 'donor5'), 
            #     ('4', 'donor0'), 
            #     ('0', 'donor1'), 
            #     ('5', 'donor2'), 
            #     ('unassigned', 'donor5'), 
            # ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor3', 
            # 'manually_added_soup_vireo_mip_triplets': [
            # },
            #     ('6', 'donor0', 'N285'), 
            # ],
            # 'rematch_anyway': True,
        },
        'demux_07_03_21_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('10', 'donor6')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('10', 'donor6')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('5', 'unassigned'), 
                ('4', 'unassigned'), 
            ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor5',
            # },
            # 'manually_added_soup_vireo_mip_triplets': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_07_03_21_1_illumina': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('10', 'donor6')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     # ({('10', 'donor6')}, ['N269', 'N270', 'N271', 'N272', 'N282', 'N283', 'N284', 'N285', 'N312', 'N313', 'N314', 'N315', 'N316', 'N6', 'N12', 'N211']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('2', 'donor0'),
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('5', 'unassigned'), 
            # ],
            # 'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
            #     'donor5',
            # },
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('4', 'unassigned', 'N169'), 
            # ],
            # 'rematch_anyway': True,
        },

        'N48_bm_12_03_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            #     {('0', 'doublet')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('0', 'unassigned'), # adding the triplet later
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N48'),
            ],
            # 'rematch_anyway': True,
        },
        'N276_bm_05_02_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor4'), # pretty much the same as the vireo donor
                ('2', 'donor3'), # pretty much the same as the vireo donor
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('6', 'donor5'),
                ('unassigned', 'donor0'),
                ('1', 'unassigned'),
                ('2', 'donor1'),
            ],
            'vireo_donor_names_with_good_mip_match_but_without_good_soup_match': {
                'donor7',
            },
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('0', 'unassigned', 'N48'),
            # ],
            # 'rematch_anyway': True,
        },
        'N365_bm_12_03_23_1': {
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N365_bm_12_03_23_1/vaf_df_csvs/0_donor2_vaf_df.csv' # N365?
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N365_bm_12_03_23_1/vaf_df_csvs/0_donor1_vaf_df.csv' # N365?
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N365_bm_12_03_23_1/vaf_df_csvs/0_donor0_vaf_df.csv' # N365?
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N365_bm_12_03_23_1/vaf_df_csvs/0_unassigned_vaf_df.csv' # N365? seems like all donors are the same donor, so i guess it is indeed N365.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            #     {('0', 'doublet')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            # # 'manually_rejected_soup_vireo_pairs': [
            # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N365'),
            ],
            # 'rematch_anyway': True,
        },
        'N242_bm_12_03_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            # # 'manually_rejected_soup_vireo_pairs': [
            # #     ('0', 'unassigned'), # adding the triplet later
            # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N242'),
            ],
            # 'rematch_anyway': True,
        },
        'N279_bm_26_10_22_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N279_bm_26_10_22_1/vaf_df_csvs/2_donor0_vaf_df.csv', # N279
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N279_bm_26_10_22_1/vaf_df_csvs/2_unassigned_vaf_df.csv', # N279? yes.
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            #     {('1', 'doublet')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor0'), 
                ('1', 'donor1'), 
                ('1', 'donor2'), 
            ],
            # # # 'manually_rejected_soup_vireo_pairs': [
            # # #     ('0', 'unassigned'), # adding the triplet later
            # # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('1', 'unassigned', 'N279'),
            ],
            # 'rematch_anyway': True,
        },
        'N188_bm_15_01_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            # # 'manually_rejected_soup_vireo_pairs': [
            # #     ('0', 'unassigned'), # adding the triplet later
            # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N188'),
            ],
            # 'rematch_anyway': True,
        },
        'N257_bm_06_10_22_1_illumina': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N257'),
            ],
            # 'rematch_anyway': True,
        },
        'N279_26_10_22_1_illumina': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N279'),
            ],
            # 'rematch_anyway': True,
        },
        'N279_bm_26_10_22_1_illumina': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N279'),
            ],
            # 'rematch_anyway': True,
        },
        'N12_bm_04_12_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            #     {('0', 'doublet')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            # # # 'manually_rejected_soup_vireo_pairs': [
            # # #     ('0', 'unassigned'), # adding the triplet later
            # # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N251'),
            ],
            # 'rematch_anyway': True,
        },
        'N196_bm_01_01_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            #     {('0', 'doublet')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'), 
                ('0', 'donor1'), 
                ('0', 'donor2'), 
            ],
            # # # # 'manually_rejected_soup_vireo_pairs': [
            # # # #     ('0', 'unassigned'), # adding the triplet later
            # # # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'N196'),
            ],
            # 'rematch_anyway': True,
        },
        'N280_bm_06_11_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('all', 'all')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('3', 'donor7'), 
            #     ('6', 'donor7'), 
            #     ('unassigned', 'donor7'), 
            #     ('8', 'donor2'), 
            #     ('unassigned', 'donor2'), 
            #     ('2', 'donor3'), 
            # ],
            'manually_rejected_soup_vireo_pairs': [
                # ('0', 'unassigned'),
                # ('4', 'unassigned'),
                ('unassigned', 'donor1'), # maybe ok, but no reason to risk it, as we already have ~10x more cells of N280.
                ('unassigned', 'donor0'), # maybe ok, but no reason to risk it, as we already have ~10x more cells of N279.
                ('5', 'unassigned'), # maybe ok, but no reason to risk it, as we already have ~10x more cells of N279.
                
                # both of these potential N279 soup-vireo clusters did not look good enough when comparing to another soup-vireo donor which is N279 for sure. most worrying, they did not fully agree with each other.
                # ('3', 'unassigned'),
                # ('6', 'unassigned'),

                # ('1', 'unassigned'), # did not fully agree with another soup-vireo donor which is N220 for sure, yet there were only a few SNPs.
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('9', 'donor9', 'N262'),
                ('9', 'unassigned', 'N262'),
                ('2', 'unassigned', 'N245'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_02_03_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'N257_bm_06_10_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('2', 'donor0'), 
                ('2', 'donor1'), 
                ('2', 'donor2'), 
            ],
            # # # 'manually_rejected_soup_vireo_pairs': [
            # # #     ('0', 'unassigned'), # adding the triplet later
            # # # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('2', 'unassigned', 'N257'),
            ],
            # 'rematch_anyway': True,
        },
        'N334_bm_15_01_23_1': {
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment_before_230707/demux_03_07_22_2_ultima/vaf_df_csvs/2_donor0_vaf_df.csv' # N272
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment_before_230707/demux_03_07_22_2_ultima/vaf_df_csvs/0_donor5_vaf_df.csv' # N269
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment_before_230707/demux_03_07_22_1_ultima/vaf_df_csvs/2_donor0_vaf_df.csv' # N269
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N334_bm_15_01_23_1/vaf_df_csvs/2_donor1_vaf_df.csv' # N269? doesn't look like it.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            #     {('all', 'all')},
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor0'), 
                ('1', 'donor2'), 
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('2', 'donor1'), # looks a bit noisy for almost 1.5k cells. ugh.
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('1', 'unassigned', 'N334'),
            ],
            # 'rematch_anyway': True,
        },
        'N251_bm_04_12_22_1': {
            # looks like 1,donor3 and 1,donor6 are the same donor. yet when comparing to /dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_04_12_22_1/vaf_df_csvs/0_donor5_vaf_df.csv, which is surely N306, the vireo-soup donors here look relatively "dirty". but i guess they must have came from that very same experiment. so what's going on? they can't have different ambient noise. and it seems unlikely that each of them is a mix of the same two donors (one of which is N306). could it be that there is some barcode-clashing which introduces the noise that we see here??
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_donor6_vaf_df.csv', # N306? not good enough..
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_donor3_vaf_df.csv', # N306? not good enough..
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_unassigned_vaf_df.csv', # N306? not good enough..
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_04_12_22_1/vaf_df_csvs/0_donor5_vaf_df.csv', # N306
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_04_12_22_1/vaf_df_csvs/1_donor7_vaf_df.csv', # N307
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/7_unassigned_vaf_df.csv', # N307? not good enough..

            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/8_donor2_vaf_df.csv', # N12
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/unassigned_donor2_vaf_df.csv', # N12? maybe a bit dirty, and we have so many 

            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N334_bm_15_01_23_1/vaf_df_csvs/1_donor2_vaf_df.csv', # N334
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/9_donor1_vaf_df.csv', # N334? very good.

            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N12_bm_04_12_22_1/vaf_df_csvs/all_all_vaf_df.csv', # N251
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/5_unassigned_vaf_df.csv', # N251? seems like yes.

            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_unassigned_vaf_df.csv' # N306? N234? seems like N306, but somewhat borderline. but looks the same as 
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_14_02_22_1/vaf_df_csvs/0_donor0_vaf_df.csv' # N231
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_14_02_22_1/vaf_df_csvs/2_donor3_vaf_df.csv' # N234
            # vaf1_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_04_12_22_1/vaf_df_csvs/0_donor5_vaf_df.csv' # N306
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_donor6_vaf_df.csv' # N306? N234? not enough data for both.
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_unassigned_vaf_df.csv' # N306? N234? seems like N306, but somewhat borderline. but looks the same as 1,donor3, so i guess it is really N306.
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_donor3_vaf_df.csv' # N306? N234? seems like N306, but borderline
            # vaf2_df_csv_file_path = '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/N251_bm_04_12_22_1/vaf_df_csvs/1_donor3__1_unassigned_vaf_df.csv' # N306? N234? seems like N306, especially with min_total_ref_and_alt_depth=30.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     # {('all', 'all')},
            #     {('1', 'unassigned'), ('1', 'donor3')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({('all', 'all')}, ['N328', 'N12', 'N48']),
            # #     ({('0', 'unassigned')}, ['N328', 'N12']),
            # #     ({('1', 'unassigned')}, ['N328', 'N12']),
            # #     ({('2', 'unassigned')}, ['N328', 'N12']),
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     # ('8', 'donor2'), # already automatically approved.
            # ],
            'manually_rejected_soup_vireo_pairs': [
                # NOTE: N100 and N234 are the parents of N306
                ('1', 'donor6'), # probably N306, but seems like not enough data to be sure, and anyway i guess we have enough of his cells.
                ('1', 'donor3'), # probably N306, but when comparing to a vireo-soup donor who is surely N306, it doesnt look good enough, so better to not risk it, as anyway i guess we have enough of his cells.
                ('1', 'unassigned'), # probably N306, but when comparing to a vireo-soup donor who is surely N306, it doesnt look good enough, so better to not risk it, as anyway i guess we have enough of his cells.
                ('7', 'unassigned'), # probably N307, but looks too noisy when comparing to a vireo-soup donor who is surely N307.
                ('unassigned', 'donor2'), # maybe a bit dirty, and we have so many cells for N12, so no need to have these.
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('5', 'unassigned', 'N251'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_14_02_22_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_22_06_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_22_06_23_1/vaf_df_csvs/2_donor4_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_22_06_23_1/vaf_df_csvs/2_unassigned_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_22_06_23_1/vaf_df_csvs/4_unassigned_vaf_df.csv', # different from 2,donor4 and from 2,unassigned
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_22_06_23_1/vaf_df_csvs/5_donor1_vaf_df.csv', # different from 2,donor4 and from 4,unassigned and from 2,unassigned

            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_22_06_23_1/vaf_df_csvs/3_donor3_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_22_06_23_1/vaf_df_csvs/doublet_donor3_vaf_df.csv', # similar to 3,donor3, but not enough to be worth the danger of using something souporcell identified as doublet.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('doublet', 'donor3')},
            #     {('doublet', 'donor1')},
            #     {('4', 'doublet')},
            #     {('2', 'doublet')},
            # ],
            # 'soup_vireo_pair_sets_and_donor_ids_to_plot_VAF_scatter': [
            #     ({
            #         # ('unassigned', 'donor4'), ('5', 'donor4'), ('4', 'donor4'), 
            #         ('doublet', 'donor4'), 
            #     }, ['N365', 'N373', 'N374', 'N375', 'N376']),
            #     # ({
            #     #     ('unassigned', 'donor4'), ('5', 'donor4'), ('4', 'donor4'), 
            #     # }, ['N365', 'N373', 'N374', 'N375', 'N376']),
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('5', 'donor1'),
                
                ('doublet', 'donor3'), # matches the MIP genotyping data of N365 better than 3,donor3 matches it...
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('3', 'donor3'), # matches the MIP genotyping data of N365 worse than doublet,donor3 matches it... and has more monocytes, while this is (IIRC) the experiment with highest enrichment for mebemp monocyte doublets.
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     # ('2', 'unassigned', 'N375'), # probably N375, but lower UMIs per cell relative to 2,donor4, which is 21.6k cells, so better to give up on these, i think.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_22_06_23_2': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_15_08_23_1': {
            # ugh all three donors are similar to each other. but still better give up on the unassigned, i think, because of low UMI count... ugh.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_15_08_23_1/vaf_df_csvs/1_donor0_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_15_08_23_1/vaf_df_csvs/1_donor2_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_15_08_23_1/vaf_df_csvs/1_unassigned_vaf_df.csv',
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor0'),
                ('1', 'donor2'),
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            #     ('1', 'unassigned', 'N280'), # UMI count dist looks bad, so give up on it. ugh.
            # ],
            # 'rematch_anyway': True,
        },
        'demux_g_24_07_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/G8_13_07_23_1/vaf_df_csvs/0_unassigned_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_24_07_23_1/vaf_df_csvs/0_unassigned_vaf_df.csv', # seems like not G8.
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_24_07_23_1/vaf_df_csvs/0_donor0_vaf_df.csv', # seems like not G8, but not enough data. seems identical to 0,unassigned
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_24_07_23_1/vaf_df_csvs/0_donor1_vaf_df.csv', # seems like not G8, but not enough data. seems identical to 0,unassigned
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_24_07_23_1/vaf_df_csvs/0_donor2_vaf_df.csv', # seems like not G8, but not enough data. seems identical to 0,unassigned
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_24_07_23_1/vaf_df_csvs/0_donor4_vaf_df.csv', # seems like not G8, but not enough data. seems identical to 0,unassigned
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_24_07_23_1/vaf_df_csvs/0_donor5_vaf_df.csv', # seems like not G8, but not enough data. seems identical to 0,unassigned


            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('0', 'donor3')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'),
                ('0', 'donor1'),
                ('0', 'donor2'),
                ('0', 'donor3'),
                ('0', 'donor4'),
                ('0', 'donor5'),
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'unassigned', 'G9'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_12_09_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_04_12_22_1/vaf_df_csvs/1_donor7_vaf_df.csv', # N307
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_12_09_23_1/vaf_df_csvs/2_donor0_vaf_df.csv', # N307? i guess yes? but not enough data.

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('2', 'donor0')},
            # ],
            # 'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
            #     ('1', 'donor0'),
            #     ('1', 'donor2'),
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('2', 'donor0'), # reject and then accept later. TODO: change this if needed when we have genotype data.
                ('0', 'donor0'), # reject and then accept later. TODO: change this if needed when we have genotype data.
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_g_28_09_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('4', 'donor0'), # might be one of the donors, but too few cells anyway
            # ],
            # 'rematch_anyway': True,
        },
        'demux_g_07_12_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor4'),
                ('0', 'donor1'),
                ('0', 'donor3'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_18_12_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_18_12_23_1/vaf_df_csvs/0_donor0_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_18_12_23_1/vaf_df_csvs/0_donor1_vaf_df.csv', # N204?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_09_01_22_1/vaf_df_csvs/5_donor0_vaf_df.csv', # N204
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_09_01_22_1/vaf_df_csvs/0_donor3_vaf_df.csv', # N206
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('0', 'donor0'),
                ('0', 'donor1'),
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_n_18_12_23_2': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_18_12_23_2/vaf_df_csvs/2_donor1_vaf_df.csv', # N413?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_18_12_23_2/vaf_df_csvs/2_unassigned_vaf_df.csv', # N413?

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('2', 'unassigned'), # we already have 55k which belong (very probably) to the same donor, presumably N413. so i don't see why to risk it.
                ('0', 'donor0'), # reject and then accept later.
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('0', 'donor0', 'N281'),
                ('2', 'donor1', 'N413'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_19_09_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('2', 'donor0'),
                ('2', 'donor1'),
                ('2', 'donor2'),
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            'manually_added_soup_vireo_mip_triplets': [
                ('2', 'unassigned', 'N200'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_30_11_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_n_19_09_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_19_09_23_1/vaf_df_csvs/1_donor0_vaf_df.csv', # N403?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_19_09_23_1/vaf_df_csvs/1_donor1_vaf_df.csv', # N403?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_19_09_23_1/vaf_df_csvs/1_unassigned_vaf_df.csv', # N403?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_19_09_23_1/vaf_df_csvs/all_all_vaf_df.csv', # N200
            
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor0'),
                ('1', 'donor1'),
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_n_04_12_23_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_n_04_12_23_2': {
            
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_04_12_23_2/vaf_df_csvs/1_donor2_vaf_df.csv', # N367
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_04_12_23_2/vaf_df_csvs/1_donor3_vaf_df.csv', # N367
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_04_12_23_2/vaf_df_csvs/1_unassigned_vaf_df.csv', # N367?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_04_12_23_1/vaf_df_csvs/0_donor1_vaf_df.csv', # N367
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_04_12_23_1/vaf_df_csvs/3_donor0_vaf_df.csv', # N191
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor0'),
                ('1', 'donor2'),
                ('1', 'donor3'),
            ],
            # 'manually_added_soup_vireo_mip_triplets': [
            # NOTE: i really don't like these many N367 cells with lower umis per cell. what's going on? some kind of clog???
            #     ('1', 'unassigned', 'N367'), # ugh. much lower umis per cell. i guess better not risk it. we already have many cells of him in demux_n_04_12_23_1...
            # ],
            # 'manually_rejected_soup_vireo_pairs': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_g_26_09_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_26_09_23_1/vaf_df_csvs/1_donor0_vaf_df.csv', # G15?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_26_09_23_1/vaf_df_csvs/1_donor1_vaf_df.csv', # G15?
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_26_09_23_1/vaf_df_csvs/1_unassigned_vaf_df.csv', # G15?

            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('all', 'all')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor1'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_01_01_24_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor2'), # far from great, but i guess that's good enough.
            ],
            # 'manually_rejected_soup_vireo_pairs': [
            #     ('1', 'donor2'), 
            # ],
            # 'rematch_anyway': True,
        },
        'demux_n_12_02_24_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_rejected_soup_vireo_pairs': [
                ('1', 'donor0'), # reject and then accept later.
            ],
            'manually_added_soup_vireo_mip_triplets': [
                ('5', 'donor0', 'N339'),
                ('1', 'donor0', 'N418'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_g_28_12_23_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_28_12_23_1/vaf_df_csvs/6_donor0_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_28_12_23_1/vaf_df_csvs/6_donor6_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_28_12_23_1/vaf_df_csvs/2_donor4_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_g_28_12_23_1/vaf_df_csvs/6_unassigned_vaf_df.csv',
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('0', 'donor2')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('6', 'donor0'),
                ('6', 'donor6'),
            ],
            'manually_rejected_soup_vireo_pairs': [
                ('2', 'donor4'), # reject and then accept later. TODO: update after we get MIP genotyping data
            ],
            # 'rematch_anyway': True,
        },
        'demux_n_bm_12_02_24_1': {
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_12_02_24_1/vaf_df_csvs/0_donor2_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_12_02_24_1/vaf_df_csvs/0_unassigned_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_12_02_24_1/vaf_df_csvs/2_donor0_vaf_df.csv',
            # '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_n_bm_12_02_24_1/vaf_df_csvs/2_unassigned_vaf_df.csv',
            
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            # 'rematch_anyway': True,
        },
        'demux_04_01_21_2': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            #     {('doublet', 'donor3')},
            #     {('3', 'donor3')},
            #     {('doublet', 'donor3'), ('3', 'donor3')},
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                # far from great, but i guess that's good enough, and we really want N220...
                ('3', 'donor3'), 
                ('doublet', 'donor3'),
            ],
            # 'rematch_anyway': True,
        },
        'demux_01_03_21_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('6', 'donor1'), 
                ('6', 'donor4'), 
            ],
            # 'rematch_anyway': True,
        },
        'demux_04_01_21_1': {
            # 'soup_vireo_pair_sets_for_log_ratio_per_donor_plot': [
            # ],
            'manually_confirmed_soup_vireo_pairs': [ # i.e., assign according to vireo-mip matching
                ('1', 'donor1'), 
                ('1', 'donor6'), 
            ],
            # 'rematch_anyway': True,
        },
        
    },
    
    'exp_name_to_empty_droplet_umi_count_thresholds': {
        'N279_26_10_22_1_illumina': (2**4.5, 2**6),
        'N279_bm_26_10_22_1': (2**7.5, 2**8.5),
        'N334_bm_15_01_23_1': (2**5.5, 2**7.3),
        'N276_bm_05_02_23_1': (2**4, 2**6),
        'N365_bm_12_03_23_1': (2**7.4, None),
        'N188_bm_15_01_23_1': (2**6.2, 2**7.3),
        'demux_03_07_22_2_ultima': (2**4.2, 2**5.9),
        'demux_01_02_21_2': (2**5.5, 2**6.5),
        'demux_11_04_21_2': (2**4.8, 2**6),
        'demux_14_02_22_1': (2**4.6, 2**5.7),
        'demux_14_02_22_2': (2**4.6, 2**5.7),
        'demux_09_06_22_1': (2**4, 2**4.8),
        'demux_09_06_22_2': (2**4, 2**4.8),
    },
    'exp_name_to_soup_vireo_mip_triplets_to_discard_in_retrospect': {
    },

    # even if the SNP was not used for clustering, it might be used when comparing to MIP data. even if no soup/vireo donor has the alt allele, it could still help when comparing to all MIP donors. min_ref=10 to not use SNPs that pretty much don't exist in the data.
    'soup_min_ref_arg': 10,
    'soup_min_alt_arg': 0,

    'max_num_of_vireo_donor_cells_to_discard_silently': 1000,
    'exp_name_to_soup_and_vireo_donor_names_to_discard_despite_large_num_of_cells': {
        'demux_n_18_12_23_2': { 
            ('2', 'unassigned'), # 240311: that's probably N413, but we already have ~56k cells of N413, so i think there is no reason to push our luck.
        },
        'demux_n_19_09_23_1': { 
            ('1', 'unassigned'), # probably N403, but we already have ~8k cells of N403, so better to give up on these, i think.
        },
        'demux_g_03_07_23_1': { 
            # TODO: hopefully can remove this when we have more MIP genotype data.
            ('1', 'donor1', 'g_03_07_23_1_a'),
            ('0', 'donor2', 'g_03_07_23_1_b'),
            ('4', 'donor3', 'g_03_07_23_1_c'),
        },
        'demux_g_12_09_23_1': { 
            # TODO: hopefully can remove this when we have more MIP genotype data.
            # ('3', 'donor2', 'g_12_09_23_1_N320'),
            ('2', 'donor0', 'g_12_09_23_1_b'), # might be N307, but not enough data to tell. only 295 cells
            ('0', 'donor0', 'g_12_09_23_1_c'), # might be N307, but not enough data to tell. only 129 cells
        },
        'demux_g_28_09_23_1': { 
            # TODO: hopefully can remove this when we have more MIP genotype data.
            ('2', 'donor4', 'g_28_09_23_1_c'),
        },
        'demux_g_26_09_23_1': { 
            ('1', 'donor0'), # i guess this is G15, but this is the experiment with mice, and we already have 1625 cells for 1,donor1, so no reason to risk it.
            ('1', 'unassigned'), # i guess this is G15, but this is the experiment with mice, and we already have 1625 cells for 1,donor1, so no reason to risk it.
        },
        'demux_g_28_12_23_1': { 
            # TODO: hopefully can remove this when we have more MIP genotype data.
            ('2', 'donor4', 'g_28_12_23_1_a'),
            ('0', 'donor2', 'g_28_12_23_1_b'), # only 338 cells, and maybe it is also G31, but can't know for sure, so should wait for more genotype data maybe.
            ('6', 'unassigned'), # 240522: that's the same donor as 6,donor6 and 6,donor0, and we already have ~8k cells of that donor, so i think there is no reason to push our luck.
        },
        'demux_g_11_01_24_1': {
            # TODO: hopefully can remove this when we have more MIP genotype data.
            ('2', 'donor3', 'g_11_01_24_1_b'),
        },
        'demux_g_22_02_24_1': {
            # TODO: hopefully can remove this when we have more MIP genotype data.
            ('4', 'donor4', 'g_22_02_24_1_c'),
            ('6', 'donor3', 'g_22_02_24_1_e'),
        },
        'demux_n_bm_12_02_24_1': {
            ('0', 'unassigned'), # 240311: doesn't fit any of current donor MIP genotype much better than others.
            ('2', 'unassigned'), # 240311: doesn't fit any of current donor MIP genotype much better than others.
        },
        'N334_bm_15_01_23_1': {
            ('2', 'donor1'), # a contamination by ultima of a donor not in our cohort??? NOTE: empty droplets have RPS4Y1...
        },

        'demux_22_06_23_1': { 
            ('4', 'unassigned'), # 230918: doesn't fit any of current donor MIP genotype much better than others. 
            ('2', 'unassigned'), # probably N375, but lower UMIs per cell relative to 2,donor4, which is 21.6k cells, so better to give up on these, i think.
        },
        'demux_n_20_07_23_1': { 
            ('0', 'unassigned'), # 230918: fits N383 okay-ish, but we already have 47k cells of N383 in 0,donor1, so i think there is no reason to push our luck.
            ('7', 'unassigned'), # 230918: doesn't fit any of current donor MIP genotype much better than others. 
            ('unassigned', 'donor1'), # 230918: fits N383 okay-ish, but we already have 47k cells of N383 in 0,donor1, so i think there is no reason to push our luck.
        },
        'demux_n_20_07_23_2': { 
            ('4', 'unassigned'), # 230918: fits N383 okay-ish, but we already have 47k cells of N383 in 4,donor0, so i think there is no reason to push our luck.
            ('6', 'unassigned'), # 230918: doesn't fit any of current donor MIP genotype much better than others. 
            ('unassigned', 'donor0'), # 230918: fits N383 okay-ish, but we already have 47k cells of N383 in 4,donor0, so i think there is no reason to push our luck.
        },

        'demux_13_06_22_1_ultima': { 
            ('1', 'donor1'), # seems like this is (very probably) NS20. we already have enough cells of NS20, and we currently discard his cells anyway, because they were not enriched for CD34 (IIRC). TODO: make sure when we have NS20 MIP genotyping data.
            ('0', 'unassigned'), # 231111: doesn't fit any of current donor MIP genotype much better than others. 
        },
        'demux_28_02_22_2': {
            ('1', 'unassigned'), # seems like souporcell and vireo just failed to cluster N219 and N236 cells (i guess we had very little cells of N237). yet, it is quite weird that '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_2/vaf_df_csvs/1_unassigned_vaf_df.csv' and '/dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_28_02_22_1/vaf_df_csvs/0_unassigned_vaf_df.csv' are so similar. maybe they are really N219 but with very high ambient noise?? POTENTIAL TODO: compare to vaf_df of N219 in his previous experiment - in the healthy atlas.
        },
        'demux_28_02_22_1': {
            # only 589 cells, but here because i suspect it is also N219.
            ('0', 'unassigned'), # dito demux_28_02_22_2, ('1', 'unassigned')
        },
        'demux_01_02_21_2': { 
            # POTENTIAL TODO: compare scRNA-seq vaf_dfs, but i don't think it would do much help, as it isn't like we are missing data - it is just ambiguous, and i guess would be ambiguous also when comparing vaf_dfs.
            # each doesn't fit any donor extremely well, but fit two donors quite well, and can't tell which is right. i guess each contains cells of at least two donors? anyway, still got enough cells for each donor in the experiment.
            ('2', 'unassigned'),
            ('4', 'unassigned'),
        },
        'demux_22_02_21_2': { 
            # POTENTIAL TODO: compare scRNA-seq vaf_dfs, but i don't think it would do much help, as it isn't like we are missing data - it is just ambiguous, and i guess would be ambiguous also when comparing vaf_dfs.
            # each doesn't fit any donor very well. i guess each contains cells of at least two donors? anyway, still got enough cells for each donor in the experiment.
            ('0', 'unassigned'),
            ('2', 'unassigned'),
        },
        'demux_07_03_21_1': { 
            # POTENTIAL TODO: compare scRNA-seq vaf_dfs, but i don't think it would do much help, as it isn't like we are missing data - it is just ambiguous, and i guess would be ambiguous also when comparing vaf_dfs.
            ('5', 'unassigned'), # fits two donors very well, and can't tell which is right. i guess each contains cells of at least two donors? anyway, still got enough cells for each donor in the experiment.
            ('4', 'unassigned'), # fits one donor very well, yet also two other donors from the same experiment, which is suspicious. and we already have 9k of the same donor (N169) in this experiment, so better not risk it. NOTE: this donor has DNMT3A with very high VAF. a coincidence?
        },
        'demux_11_04_21_2': { 
            # i guess it contains cells of at least two donors? anyway, still got enough cells for each donor in the experiment.
            ('1', 'unassigned'),
            ('3', 'unassigned'),
        },
        'demux_27_02_23_1': {
            ('5', 'unassigned'), # 230724: doesn't fit any of current donor MIP genotype well.
        },
        'demux_06_11_22_1': { 
            ('3', 'unassigned'), # 230722: doesn't fit any of current donor MIP genotype well.
        },
        'N280_bm_06_11_22_1': { 
            # only donor0, donor1, and donor8 look good. everything else is a mess.
            ('8', 'unassigned'), # 230724: doesn't fit any of current donor MIP genotype well.
            ('6', 'unassigned'), # 230724: doesn't fit any of current donor MIP genotype well.
            ('8', 'donor6'), # 230724: doesn't fit any of current donor MIP genotype well.
            ('8', 'donor9'), # 230724: doesn't fit any of current donor MIP genotype well.
            ('unassigned', 'donor1'), # maybe ok, but no reason to risk it, as we already have ~10x more cells of N280.
            ('unassigned', 'donor0'), # maybe ok, but no reason to risk it, as we already have ~10x more cells of N279.
            ('5', 'unassigned'), # maybe ok, but no reason to risk it, as we already have ~10x more cells of N279.
        },
        'demux_13_11_22_2': { 
            ('0', 'unassigned'), # 230723: doesn't fit any of current donor MIP genotype much better than others.
            ('2', 'unassigned'), # 230723: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_23_11_22_1': { 
            ('5', 'unassigned'), # 230723: doesn't fit any of current donor MIP genotype much better than others.
            ('3', 'unassigned'), # 230723: doesn't fit any of current donor MIP genotype much better than others.
        },

        'N251_04_12_22_1': { 
            # both of these seemed somewhat similar to N12, but not good enough when comparing vaf_dfs. weirdly, the log ratio indicated similarity to N234, which doesn't make much sense, except that N306 is the son of N234 (and N100), and N306 is most dominant in demux_04_12_22_1... anyway, not important enough to figure out what happened here, i guess.
            ('1', 'donor0'),
            ('1', 'unassigned'),
        },

        'N251_bm_04_12_22_1': { 
            ('1', 'donor3'), # probably N306, but when comparing to a vireo-soup donor who is surely N306, it doesnt look good enough, so better to not risk it, as anyway i guess we have enough of his cells.
            ('1', 'unassigned'), # probably N306, but when comparing to a vireo-soup donor who is surely N306, it doesnt look good enough, so better to not risk it, as anyway i guess we have enough of his cells.
            ('7', 'unassigned'), # probably N307, but looks too noisy when comparing to a vireo-soup donor who is surely N307.
            ('unassigned', 'donor2'), # maybe a bit dirty, and we have so many cells for N12, so no need to have these.
            
            ('3', 'unassigned'), # 230710: doesn't fit any of current donor MIP genotype well.
        },

        'demux_03_07_22_1_ultima': { 
            ('1', 'unassigned'), # vaf_dfs show it is definitely not N269. it seems to fit N57, but this doesn't seem to make any sense, and it isn't the first time N57 seems to fit where it doesn't make sense. i guess N57 just has a lot of MIP genotyping data?
        },
        'demux_03_07_22_2_ultima': { 
            ('4', 'unassigned'), # vaf_dfs show it is definitely not N269. it seems to fit N30, but this doesn't seem to make any sense, and it isn't the first time N30 seems to fit where it doesn't make sense. i guess N30 just has a lot of MIP genotyping data? weirdly, it seems to match /dummy/dummy/dummy/raid/PBMC_CD34/nili_data/donor_assignment/demux_03_07_22_1_ultima/vaf_df_csvs/1_unassigned_vaf_df.csv.
        },
        'demux_08_12_22_1': { 
            ('10', 'unassigned'), # 230710: doesn't fit any of current donor MIP genotype much better than others.
            ('14', 'donor5'), # 230710: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_10_07_22_1_ultima': { 
            ('7', 'unassigned'), # 230710: doesn't fit any of current donor MIP genotype much better than others.
        },

        'demux_13_12_22_1': { 
            ('4', 'unassigned'), # 230710: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_29_01_23_1': { 
            ('2', 'unassigned'), # 230710: doesn't fit any of current donor MIP genotype much better than others.
            ('0', 'unassigned'), # 230710: doesn't fit any of current donor MIP genotype much better than others.
            
            # 231005: we have enough assigned soup-vireo cells for N337 and N335, so no reason to risk adding small soup vireo donors unassigned by soup/vireo.
            ('1', 'unassigned'),
            ('6', 'unassigned'),
        },
        'demux_19_12_21_1': { 
            ('5', 'unassigned'), # 230909: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_14_11_22_1': { 
            ('3', 'unassigned'), # 230910: doesn't fit any of current donor MIP genotype much better than others.
            ('5', 'unassigned'), # 230910: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_01_02_21_2_illumina': { 
            ('0', 'unassigned'), # 230713: doesn't fit any of current donor MIP genotype much better than others.
            ('4', 'unassigned'), # 230713: doesn't fit any of current donor MIP genotype much better than others.
            ('5', 'unassigned'), # 230713: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_07_03_21_1_illumina': { 
            ('0', 'unassigned'), # 230713: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_11_04_21_2_illumina': { 
            ('5', 'unassigned'), # 230713: doesn't fit any of current donor MIP genotype much better than others.
        },
        'demux_08_02_23_1': { 
            ('6', 'unassigned'), # 230723: doesn't fit any of current donor MIP genotype much better than others. as these are 5842 cells, i guess this is a mix of multiple donors.
        },
        'N276_bm_05_02_23_1': { 
            ('6', 'donor5'), # 230913: i guess this is a contamination by a donor whose genotype data is missing.
            ('unassigned', 'donor0'), # 230913: doesn't fit any of current donor MIP genotype much better than others.
            ('1', 'unassigned'), # 230913: doesn't fit any of current donor MIP genotype much better than others.
            ('2', 'donor1'), # 230913: i guess this is a contamination by a donor whose genotype data is missing.
        },
        'demux_05_02_23_1': { 
            # 231005: we have a soup vireo donor with 36k cells, so no reason to risk adding small soup vireo donors unassigned by soup/vireo. also, in the unified metacell model, cell state compositions were slightly different between some soup vireo donors of this exp.
            ('1', 'unassigned'),
            ('6', 'unassigned'),
            ('unassigned', 'donor6'),
        },
        'demux_02_01_22_1': { 
            # 231005: we have two soup vireo donor with 30k cells each, so no reason to risk adding small soup vireo donors unassigned by soup/vireo. also, in the unified metacell model, cell state compositions were slightly different between some soup vireo donors of this exp.
            ('unassigned', 'donor2'),
        },
        'demux_02_05_22_1': { 
            # 231005: we have enough assigned soup-vireo cells for N245, so no reason to risk adding small soup vireo donors unassigned by soup/vireo.
            ('unassigned', 'donor1'),
        },
        'demux_02_05_22_2': { 
            # 231005: we have enough assigned soup-vireo cells for N245, so no reason to risk adding small soup vireo donors unassigned by soup/vireo.
            ('unassigned', 'donor1'),
        },
        'demux_13_06_22_1': { 
            # 231005: we have enough assigned soup-vireo cells for N262, so no reason to risk adding small soup vireo donors unassigned by soup/vireo.
            ('unassigned', 'donor0'),
        },
        'demux_n_bm_15_08_23_1': { 
            # 231005: we have enough assigned soup-vireo cells for N262, so no reason to risk adding small soup vireo donors unassigned by soup/vireo.
            ('1', 'unassigned'),
        },
        'demux_n_04_12_23_2': { 
            # NOTE: looking at stacked composition plot, as well as s_phase_sig_across_hsc_mpp_mebemp_l, there seems to be a difference between N367 cells in demux_n_04_12_23_1 and demux_n_04_12_23_2. Could it actually be different clones?? could it be that demux_n_04_12_23_2 is a more abnormal clone???
            ('1', 'unassigned'), # maybe ok, but no reason to risk it, as we already have ~19k cells of N367 in demux_n_04_12_23_1
        },
        'demux_g_28_03_24_1': { 
            # 240522: we have enough assigned soup-vireo cells for each of the two donors, so no reason to push our luck, i think.
            ('2', 'unassigned'),
        },
    },

    'exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs': {
        ### start single donor exps
        'G8_13_07_23_1': { # TODO: replace all of this with "all_all in plot looked good", if indeed right.
            'donor_id': 'G8', # TODO: verify that it is ok when we have MIP genotyping.
        },
        'demux_n_24_07_23_1': {
            'donor_id': 'N387', # all_all in plot looked good.
            
            ### adding other attributes because the combination of donor_id and bleeding_date are not in clinical_df.
            # 'donor_sex': '',
            # 'donor_age': -3,
            # 'diagnosis': '',
        },
        'demux_n_bm_24_07_23_1': {
            'donor_id': 'N387', # all_all in plot looked good.
            
            ### adding other attributes because the combination of donor_id and bleeding_date are not in clinical_df.
            # 'donor_sex': '',
            # 'donor_age': -3,
            # 'diagnosis': '',
        },
        'N48_bm_12_03_23_1': {
            'donor_id': 'N48', # all_all in plot looked good.
        },
        'N242_bm_12_03_23_1': {
            'donor_id': 'N242', # all_all in plot looked good.
        },
        'N188_bm_15_01_23_1': {
            'donor_id': 'N188', # all_all in plot looked good.
        },
        'N12_bm_04_12_22_1': {
            'donor_id': 'N251', # all_all in plot looked good. also, 2,doublet looked good. this is not a mistake, N12_bm_04_12_22_1 and N251_bm_04_12_22_1 labels were swapped, it seems.
        },
        'N196_bm_01_01_23_1': {
            'donor_id': 'N196', # all_all in plot looked good.
        },
        'N257_bm_06_10_22_1_illumina': {
            'donor_id': 'N257', # all_all in plot looked good.
        },
        'N279_26_10_22_1_illumina': {
            'donor_id': 'N279', # all_all in plot looked good.
        },
        'N279_bm_26_10_22_1_illumina': {
            'donor_id': 'N279', # all_all in plot looked good.
        },
        'N365_bm_12_03_23_1': {
            'donor_id': 'N365', # all_all in plot looked good.
        },
        # NOTE: a (wild?) guess: souporcell fails in the doublet identification step if it is run with num_of_donors>1 but estimates all of the droplets to belong to the same donor.
        'N328_15_12_22_1': {
            'donor_id': 'N328', # all_all in plot looked good.
            
            ### no need to add other attributes because the combination of donor_id and bleeding_date are in clinical_df.
            # 'donor_sex': 'female',
            # 'donor_age': 2.5, 
            # 'diagnosis': 'JMML', 
        },
        'N257_bm_06_10_22_1': {
            'donor_id': 'N257', # all_all in plot looked good.
            
            ### adding other attributes because the combination of donor_id and bleeding_date are not in clinical_df.
            'donor_sex': 'female',
            # 'diagnosis': 'AML',
            'donor_age': 63,
        },
        'N279_26_10_22_1': {
            'donor_id': 'N279', # all_all in plot looked good.
            # 'donor_sex': 'female',
            # 'diagnosis': 'MDS',
            # 'donor_age': 82,
        },
        'N279_bm_26_10_22_1': {
            'donor_id': 'N279', # all_all in plot looked good.
        },
        'demux_n_bm_19_09_23_1': {
            'donor_id': 'N200', # all_all in plot looked good.
        },
        ### end single donor exps



        ### start special cases (e.g. samples that waited and their controls (which didn't wait))
        ('demux_22_06_23_2', 'N264'): {
            'donor_sex': 'female',
            'diagnosis': 'normal_24h',
            'donor_age': 32,
        },
        ('demux_22_06_23_2', 'N275'): {
            'donor_sex': 'male',
            'diagnosis': 'normal_36h',
            'donor_age': 41,
        },
        ('demux_22_06_23_2', 'N285'): {
            'donor_sex': 'female',
            'diagnosis': 'normal_48h',
            'donor_age': 41,
        },
        ('demux_n_15_08_23_2', 'N264'): {
            'donor_sex': 'female',
            'diagnosis': 'normal',
            'donor_age': 32,
        },
        ('demux_n_15_08_23_2', 'N1'): {
            'donor_sex': 'male',
            'diagnosis': 'normal',
            'donor_age': 48,
        },
        ('demux_n_15_08_23_3', 'N264'): {
            'donor_sex': 'female',
            'diagnosis': 'normal_24h',
            'donor_age': 32,
            'bleeding_date': '14.08.23', # NOTE: i guess it is likely this will make my code fail. maybe silently. ugh.
        },
        ('demux_n_15_08_23_3', 'N1'): {
            'donor_sex': 'male',
            'diagnosis': 'normal_24h',
            'donor_age': 47,
            'bleeding_date': '14.08.23', # NOTE: i guess it is likely this will make my code fail. maybe silently. ugh.
        },
        ('demux_g_24_07_23_1', 'G8'): {
            # 'donor_sex': 'male',
            # 'diagnosis': 'pmf',
            # 'donor_age': 73,
            'bleeding_date': '13.07.23', # NOTE: i guess it is likely this will make my code fail. maybe silently. ugh.
        },
        ('demux_g_12_09_23_1', 'N320'): {
            'donor_sex': 'female',
            'diagnosis': 'other_exp_condition_excluded',
            'donor_age': 38,
        },
        ('demux_n_18_12_23_2', 'N281'): {
            'donor_age': 68,
            'donor_sex': 'female',
            'diagnosis': 'other_exp_condition_excluded',
        },
        ('N257_bm_06_10_22_1_illumina', 'N257'): {
            'donor_sex': 'female',
            'donor_age': 63,
        },
        'T6_25_05_22_1': {
            'donor_id': 'T6', # NOTE: no way to verify, i think.

            ### adding other attributes because the combination of donor_id and bleeding_date are not in clinical_df.
            'donor_sex': 'female',
            'donor_age': 2.0, # NOTE: not sure this is right, but not relevant to us anyway, i think.
        },
        'NS20_29_09_21_1': {
            'donor_id': 'NS20', # TODO: verify that it is ok when we have MIP genotyping for NS20, though it seems highly likely that it is NS20.
            
            ### adding other attributes because the combination of donor_id and bleeding_date are not in clinical_df.
            'donor_sex': 'male',
            'donor_age': -3, # TODO: fix.
        },
        ### start special cases (e.g. samples that waited and their controls (which didn't wait))



        ### start missing from table but i guess would stay missing.
        ('demux_17_08_20_1', 'N35'): { # TODO: remove this. should be in the clinical_df
            'donor_sex': 'female',
            'diagnosis': 'normal', 
            'donor_age': 68,
        },
        ('demux_17_08_20_1', 'N37'): { # TODO: remove this. should be in the clinical_df
            'donor_sex': 'female',
            'diagnosis': 'normal',
            'donor_age': 79,
        },
        ### end missing from table but i guess would stay missing.


        ('demux_g_12_09_23_1', 'G12'): { # TODO: remove this. should be in the clinical_df
            'donor_age': 71,
        },
        ('demux_17_08_20_1', 'N199'): { # TODO: remove this. should be in the clinical_df
            'donor_sex': 'male',
            'diagnosis': 'PV', 
            'donor_age': 66,
        },
        ('demux_01_02_21_2', 'N232'): { # TODO: remove this. should be in the clinical_df
            'donor_sex': 'male',
            'diagnosis': 'cytopenia',
            'donor_age': 68,
        },
        ('demux_01_02_21_2_illumina', 'N232'): { # TODO: remove this. should be in the clinical_df
            'donor_sex': 'male',
            'diagnosis': 'cytopenia',
            'donor_age': 68,
        },
        ('demux_01_02_21_2', 'N135'): { # TODO: remove this. should be in the clinical_df
            'diagnosis': 'leukocytosis',
        },
        ('demux_01_02_21_2_illumina', 'N135'): { # TODO: remove this. should be in the clinical_df
            'diagnosis': 'leukocytosis',
        },
        
    },

    'generate_exp_metadata_csv': {
        'ultima_summary_file_infos': [
            {
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/data/users/aviezerl/proj/blood_aging/data/ultima/20220628',
                'fastq_summary_file_name': 'weizmann_library_summary_20220628.csv',
                'approx_date_of_fastq_files': '30_6_22',
            },
            {
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/data/users/aviezerl/proj/blood_aging/data/ultima/reads_24_libraries',
                'fastq_summary_file_name': 'weizmann_library_summary_20220620.csv',
                'approx_date_of_fastq_files': '20_6_22',
            },
        ],
        'individual_experiment_info': [
            {
                'exp_name': 'demux_02_05_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '020522_1',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/first/ni/nilisf/220704_A00929_0697_AHHFW5DRX2/HHFW5DRX2/outs/fastq_path/020522_1',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_02_05_22_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '020522_2',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/first/ni/nilisf/220704_A00929_0697_AHHFW5DRX2/HHFW5DRX2/outs/fastq_path/020522_2',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_06_06_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '060622_1',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/second/ni/nilisf/220704_A00929_0698_BH73T2DRX2/H73T2DRX2/outs/fastq_path/060622_1',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_06_06_22_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '060622_2',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/first/ni/nilisf/220704_A00929_0697_AHHFW5DRX2/HHFW5DRX2/outs/fastq_path/060622_2',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_09_06_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '090622_1',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/first/ni/nilisf/220704_A00929_0697_AHHFW5DRX2/HHFW5DRX2/outs/fastq_path/090622_1',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_09_06_22_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '090622_2',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/first/ni/nilisf/220704_A00929_0697_AHHFW5DRX2/HHFW5DRX2/outs/fastq_path/090622_2',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_12_06_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '120622',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/second/ni/nilisf/220704_A00929_0698_BH73T2DRX2/H73T2DRX2/outs/fastq_path/120622',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_13_06_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '130622',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/first/ni/nilisf/220704_A00929_0697_AHHFW5DRX2/HHFW5DRX2/outs/fastq_path/130622',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_19_06_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '190622_1',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/second/ni/nilisf/220704_A00929_0698_BH73T2DRX2/H73T2DRX2/outs/fastq_path/190622_1',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_19_06_22_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '190622_2',
                'approx_date_of_fastq_files': '5_7_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/5_7_22/first/ni/nilisf/220704_A00929_0697_AHHFW5DRX2/HHFW5DRX2/outs/fastq_path/190622_2',
                'is_ultima': False,
            },
            {
                'exp_name': 'demux_03_07_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '030722_1',
                'approx_date_of_fastq_files': '9_11_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/9_11_22/ni/nilisf/221108_A00929_0790_AHJ2FCDRX2/030722_1',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_03_07_22_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '030722_2',
                'approx_date_of_fastq_files': '9_11_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/9_11_22/ni/nilisf/221108_A00929_0790_AHJ2FCDRX2/030722_2',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_10_07_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '100722_1',
                'approx_date_of_fastq_files': '9_11_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/9_11_22/ni/nilisf/221108_A00929_0790_AHJ2FCDRX2/100722_1',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_10_07_22_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '100722_2',
                'approx_date_of_fastq_files': '9_11_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/9_11_22/ni/nilisf/221108_A00929_0790_AHJ2FCDRX2/100722_2',
                'is_ultima': False, 
            },
            
            {
                'exp_name': 'N257_bm_06_10_22_1_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '257_BM',
                'approx_date_of_fastq_files': '9_11_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/9_11_22/ni/nilisf/221108_A00929_0790_AHJ2FCDRX2/257_BM',
                'is_ultima': False, 
                'ids_in_exp': 'N257',
            },
            
            {
                'exp_name': 'N279_26_10_22_1_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '279_PB',
                'approx_date_of_fastq_files': '9_11_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/9_11_22/ni/nilisf/221108_A00929_0790_AHJ2FCDRX2/279_PB',
                'is_ultima': False, 
                'ids_in_exp': 'N279',
            },
            
            {
                'exp_name': 'N279_bm_26_10_22_1_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '279_BM',
                'approx_date_of_fastq_files': '9_11_22',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/9_11_22/ni/nilisf/221108_A00929_0790_AHJ2FCDRX2/279_BM',
                'is_ultima': False, 
                'ids_in_exp': 'N279',
            },
            
            {
                # 'exp_name': 'demux_tal_25_05_22_1', # NOTE: changed to T6_25_05_22_1 in retrospect. hopefully won't cause any problems.
                'exp_name': 'T6_25_05_22_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Tal_3',
                'approx_date_of_fastq_files': '03_07_22',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/PBMC_CD34/raw_data/tal_data/T6_sequencing_on_03_07_22/10x_on_25_05_22',
                'is_ultima': False, 
                'ids_in_exp': 'T6',
            },
            {
                'exp_name': 'demux_23_11_20_1_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'cDNA_2411',
                'approx_date_of_fastq_files': '04_12_20',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/4_12_20/ni/nilisf/201203_A00929_0215_AHVHJJDRXX/HVHJJDRXX/outs/fastq_path/cDNA_2411',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_30_11_20_1_illumina', # NOTE: 240523: my mistake, i think. seems like it really should have been demux_30_11_20_2_illumina, but it's too late now and we will keep it this way.
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'cDNA_3011p2',
                'approx_date_of_fastq_files': '04_12_20',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/4_12_20/ni/nilisf/201203_A00929_0215_AHVHJJDRXX/HVHJJDRXX/outs/fastq_path/cDNA_3011p2',
                'is_ultima': False, 
            },
            
            {
                'exp_name': 'demux_01_02_21_2_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '010221_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/010221_p2',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_22_02_21_2_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '220221_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/220221_p2',
                'ids_in_exp': 'N152,N153,N155,N156,N233',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_07_03_21_1_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '070321_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/070321_p1',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_11_04_21_2_illumina',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '110421_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/110421_p2',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_11_04_21_1', # was here at first only because i compared it to demux_11_04_21_1_elembio
                'cellranger_with_introns': False,
                'fastq_file_prefix': '110421_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/110421_p1',
                'ids_in_exp': 'N83,N182,N183,N244',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_04_01_21_2', # here because of N220 (AKA N116)
                'cellranger_with_introns': False,
                'fastq_file_prefix': '040121_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/040121_p2',
                'ids_in_exp': 'N108,N110,N220,N111',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_11_04_21_1_elembio',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Weizmann2',
                'approx_date_of_fastq_files': '11_12_23',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20231129_shlush_elembio/20231120_AV224504_EBSL-0444-OBPA/Weizmann2',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_07_03_21_2', # was here at first only because i compared it to demux_07_03_21_2_elembio
                'cellranger_with_introns': False,
                'fastq_file_prefix': '070321_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/070321_p2',
                'ids_in_exp': 'N174,N173,N172,N170,N231',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_07_03_21_2_elembio',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Weizmann1',
                'approx_date_of_fastq_files': '11_12_23',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20231129_shlush_elembio/20231120_AV224504_EBSL-0444-OBPA/Weizmann1',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_22_06_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '220623_p1',
                'approx_date_of_fastq_files': '5_7_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230622_nili/outs/fastq_path',
                'is_ultima': False, 
                'ids_in_exp': 'N373,N375,N376,N365', # need this here because different donors in the two experiments on that day.
            },
            {
                'exp_name': 'demux_22_06_23_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '220623_p2',
                'approx_date_of_fastq_files': '5_7_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230622_nili/outs/fastq_path',
                'is_ultima': False, 
                'ids_in_exp': 'N377,N285,N264,N275', # need this here because different donors in the two experiments on that day.
            },
            
            {
                'exp_name': 'demux_17_08_20_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili-17_8_20',
                'approx_date_of_fastq_files': '13_8_20',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/13_8_20/bcl_output/outs/fastq_path',
                'is_ultima': False, 
                'ids_in_exp': 'N35,N199,N37,N235', # N36 is N199. N38 is N235.
            },
            
            {
                'exp_name': 'demux_g_26_06_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Gal_mpn_1',
                'approx_date_of_fastq_files': '27_8_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230816_shlush_admera/unified_fastq_folder',
                'is_ultima': False, 
                'ids_in_exp': 'G1,G2,G3,G4',
            },
            {
                'exp_name': 'demux_g_03_07_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Gal_mpn_2',
                'approx_date_of_fastq_files': '27_8_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230816_shlush_admera/unified_fastq_folder',
                'is_ultima': False, 
                'ids_in_exp': 'G5,G6,G7',
            },
            {
                'exp_name': 'G8_13_07_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Gal_mpn_3',
                'approx_date_of_fastq_files': '27_8_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230816_shlush_admera/unified_fastq_folder',
                'is_ultima': False, 
                'ids_in_exp': 'G8',
            },
            
            {
                'exp_name': 'demux_n_20_07_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_20_07_23_1',
                'approx_date_of_fastq_files': '27_8_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230816_shlush_admera/unified_fastq_folder',
                'is_ultima': False, 
                'ids_in_exp': 'N381,N382,N383,N384,N385,N386', # why not be on the safe side
            },
            {
                'exp_name': 'demux_n_20_07_23_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_20_07_23_2',
                'approx_date_of_fastq_files': '27_8_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230816_shlush_admera/unified_fastq_folder',
                'is_ultima': False, 
                'ids_in_exp': 'N381,N382,N383,N384,N385,N386', # why not be on the safe side
            },
            {
                'exp_name': 'demux_n_bm_24_07_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_24_07_23_BM',
                'approx_date_of_fastq_files': '27_8_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230816_shlush_admera/unified_fastq_folder',
                'is_ultima': False, 
                'ids_in_exp': 'N387', 
            },
            {
                'exp_name': 'demux_n_24_07_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_24_07_23_PB',
                'approx_date_of_fastq_files': '27_8_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20230816_shlush_admera/unified_fastq_folder',
                'is_ultima': False, 
                'ids_in_exp': 'N387', 
            },
            
            
            {
                'exp_name': 'demux_21_02_21_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '210221_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_nimrod_blood_aging_dir/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/210221_p1',
                'is_ultima': False, 
                'ids_in_exp': 'N144,N148,N149,N151',
            },
            {
                'exp_name': 'demux_22_02_21_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '220221_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_nimrod_blood_aging_dir/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/220221_p1',
                'is_ultima': False, 
                'ids_in_exp': 'N48,N154,N158,N160,N215',
            },


            {
                'exp_name': 'demux_n_13_08_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_13_08_23',
                'approx_date_of_fastq_files': '8_11_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20231010_shlush_admera',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_bm_15_08_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_15_08_23_1',
                'approx_date_of_fastq_files': '8_11_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20231010_shlush_admera',
                'is_ultima': False, 
                'ids_in_exp': 'N392,N280', # N393 == N280
            },
            {
                'exp_name': 'demux_n_15_08_23_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_15_08_23_2',
                'approx_date_of_fastq_files': '8_11_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20231010_shlush_admera',
                'is_ultima': False, 
                'ids_in_exp': 'N392,N280,N264,N1', # N393 == N280
            },
            {
                'exp_name': 'demux_n_15_08_23_3',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_15_08_23_3',
                'approx_date_of_fastq_files': '8_11_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20231010_shlush_admera',
                'is_ultima': False, 
                'ids_in_exp': 'N264,N1,N252',
            },
            {
                'exp_name': 'demux_g_24_07_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'gal_mpn_pool_4',
                'approx_date_of_fastq_files': '8_11_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20231010_shlush_admera',
                'is_ultima': False, 
                'ids_in_exp': 'G8,G9,G10,G11', # NOTE: G8 is special (so must specify him here explicitly) because these are frozen cells whose bleeding date is actually 230713. i guess this would cause my code to fail somewhere...
            },
            {
                'exp_name': 'demux_g_12_09_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'gal_mpn_pool_5',
                'approx_date_of_fastq_files': '8_11_23',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20231010_shlush_admera',
                'ids_in_exp': 'G12,G13,G14,N320', 
                'is_ultima': False, 
            },
            
            {
                'exp_name': 'demux_g_28_09_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Gal_MPN_6',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_g_07_12_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Gal_MPN_7',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                # 'ids_in_exp': 'G21,G22,G23',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_18_12_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_18_12_23_1',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                'ids_in_exp': 'N204,N413',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_18_12_23_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_18_12_23_2',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                'ids_in_exp': 'N413,N281',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_bm_19_09_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_19_9_23_BM',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                'ids_in_exp': 'N200',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_19_09_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_19_9_23_PB',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_30_11_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_30_11_23_1',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_04_12_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_4_12_23_1',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_04_12_23_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_4_12_23_2',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_05_12_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_5_12_23',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                'ids_in_exp': 'N382,N410,N192',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_g_26_09_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Orly_260923',
                'approx_date_of_fastq_files': '23_2_24',
                'fastq_dir_path': '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tg_seq_runs/20240131_shlush_admera/unified_fastq_folder',
                'ids_in_exp': 'G15', # TODO: update? not really needed if we force num of donors...
                'is_ultima': False, 
            },
            
            
            {
                'exp_name': 'demux_n_01_01_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_010124',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240303_shlush_lab_nili_furer/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_bm_12_02_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_120224_BM',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240303_shlush_lab_nili_furer/outs/fastq_path',
                'ids_in_exp': 'N418,N419',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_12_02_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_120224_PB',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240303_shlush_lab_nili_furer/outs/fastq_path',
                'ids_in_exp': 'N418,N419,N339,N421',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_n_14_02_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'Nili_140224',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240303_shlush_lab_nili_furer/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            
            {
                'exp_name': 'demux_g_28_12_23_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'mpn_9',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240229_shlush_lab_gal_dadi/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_g_08_01_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'mpn_10',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240229_shlush_lab_gal_dadi/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_g_11_01_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'mpn_11',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240229_shlush_lab_gal_dadi/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_g_01_02_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'mpn_12',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240229_shlush_lab_gal_dadi/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_g_08_02_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'mpn_13',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240229_shlush_lab_gal_dadi/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_g_22_02_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'mpn_14',
                'approx_date_of_fastq_files': '4_3_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240229_shlush_lab_gal_dadi/outs/fastq_path',
                # 'ids_in_exp': '',
                'is_ultima': False, 
            },

            
            # generated ids_in_exp by:
            # set(n_c_ad.obs['exp_name'].unique()) - set(c_ad.obs['exp_name'].unique())
            # donor_ids = sorted(n_c_ad.obs.loc[n_c_ad.obs['exp_name'] == 'xxxxx', 'donor_id'].unique())
            # donor_ids = sorted(get_sc_rna_seq_preprocessing_params()['other_donor_id_to_donor_id_i_use'].get(x, x) for x in donor_ids)
            # ','.join(donor_ids)
            {
                'exp_name': 'demux_01_02_21_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '010221_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/010221_p1',
                'ids_in_exp': 'N113,N132,N134,N223',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_01_03_21_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '010321_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/010321_p1',
                'ids_in_exp': 'N161,N162,N163,N225,N226',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_04_01_21_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '040121_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/040121_p1',
                'ids_in_exp': 'N107,N112,N114,N115,N238',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_07_12_20_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '071220_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/071220_p1',
                'ids_in_exp': 'N269,N84,N86,N89',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_07_12_20_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '071220_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/071220_p2',
                'ids_in_exp': 'N213,N59,N85,N87',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_11_01_21_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '110121_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/110121_p1',
                'ids_in_exp': 'N117,N119,N121,N18',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_11_01_21_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '110121_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/110121_p2',
                'ids_in_exp': 'N118,N120,N123,N254',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_21_01_21_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '210121_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/210121_p1',
                'ids_in_exp': 'N126,N127,N16,N17',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_21_01_21_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '210121_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/210121_p2',
                'ids_in_exp': 'N124,N125,N128,N129',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_21_02_21_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '210221_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/210221_p2',
                'ids_in_exp': 'N143,N145,N146,N219,N273',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_21_12_20_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '211220_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/211220_p1',
                'ids_in_exp': 'N268,N90,N93,N94,N92', 
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_21_12_20_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '211220_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/211220_p2',
                'ids_in_exp': 'N216,N92,N95,N96',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_28_12_20_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '281220_p1',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/281220_p1',
                'ids_in_exp': 'N100,N102,N103,N234',
                'is_ultima': False, 
            },
            {
                'exp_name': 'demux_28_12_20_2',
                'cellranger_with_introns': False,
                'fastq_file_prefix': '281220_p2',
                'approx_date_of_fastq_files': '21_4_21',
                'fastq_dir_path': '/home/nimrodra/proj/blood_aging/our_data/21_4_21/ni/nilisf/210420_A00929_0314_AHYJ3NDMXX/HYJ3NDMXX/outs/fastq_path/281220_p2',
                'ids_in_exp': 'N101,N104,N105,N106',
                'is_ultima': False, 
            },
            
            {
                'exp_name': 'demux_g_28_03_24_1',
                'cellranger_with_introns': False,
                'fastq_file_prefix': 'mpn15',
                'approx_date_of_fastq_files': '21_5_24',
                'fastq_dir_path': '/dummy/dummy/dummy/raid/symlink_to_tgdb/tg_seq_runs/20240521_shlush_lab_gal_dadi/fastq_path/',
                'ids_in_exp': 'G49,G50',
                'is_ultima': False, 
            },
        ],


        'metadata_key_vals_to_add': {
            'transcriptome_for_cellranger_dir_path': '/dummy/dummy/dummy/raid/human_genome_files/refdata-cellranger-GRCh38-3.0.0',

            # these are SNPs that at the time (when nimrodra downloaded the file), had an minor allele frequency > 0.05 (i.e., 5e-2) - i guess in the thousand genomes project, phase 3.
            # looks like it was downloaded from https://sourceforge.net/projects/cellsnp/files/SNPlist/, which says about this file: "7.4M SNPs with minor allele frequency (MAF) > 0.05" 
            'souporcell_and_cellsnp_lite_candidate_snps_path': '/net/dummy/dummy/dummy/export/tgdata/users/orenmil/human_genome_files/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf', # NOTE: 240715: downloaded presumably the same file from https://sourceforge.net/projects/cellsnp/, but it was slightly different: ~3k more lines. So i guess it is being updated once in a while?
            'cellsnp_lite_minCOUNT': 20,
            
            # hg38
            'souporcell_ref_file_path': '/dummy/dummy/dummy/raid/human_genome_files/genome.fa',
        },

        'exp_name_to_forced_ids_in_exp': {
            # this isn't really useful for experiments in individual_experiment_info - just set ids_in_exp there...
            # NOTE: this comes before presumably_swapped_label_exp_name_pairs and adding ids according to identified accidents (in generate_sc_rna_seq_exp_metadata_csv.ipynb).
            'demux_04_12_22_1': 'N12,N306,N307,N308,N309',
            'demux_22_02_21_2': 'N152,N153,N155,N156,N233',
        },
        'presumably_swapped_label_exp_name_pairs': {
            # NOTE: this comes before exp_name_to_forced_key_vals (in generate_sc_rna_seq_exp_metadata_csv.ipynb).
            frozenset({'demux_27_02_23_1', 'demux_28_02_23_1'}),
            frozenset({'N12_bm_04_12_22_1', 'N251_bm_04_12_22_1'}),
        },
        'exp_name_to_forced_key_vals': {
            # NOTE: this comes last in generate_sc_rna_seq_exp_metadata_csv.ipynb, so overwrites everything else.
            
        },
    },
    'test_set_exp_names': [

    ],

    'ultima_sequenced_after_illumina_experiment_names': [
        'demux_02_05_22_1',
        'demux_02_05_22_2',
        'demux_03_07_22_1',
        'demux_03_07_22_2',
        'demux_06_06_22_1',
        'demux_06_06_22_2',
        'demux_09_06_22_1',
        'demux_09_06_22_2',
        'demux_10_07_22_1',
        'demux_10_07_22_2',
        'demux_12_06_22_1',
        'demux_13_06_22_1',
        'demux_19_06_22_1',
        'demux_19_06_22_2',
    ],

    'technical_repeat_experiment_names': [
        'demux_23_11_20_1', 'demux_23_11_20_1_illumina', 
        'demux_30_11_20_1', # NOTE: 240523: my mistake, i think. seems like it really should have been demux_30_11_20_2, but it's too late now and we will keep it this way.
        'demux_30_11_20_1_illumina', # NOTE: 240523: my mistake, i think. seems like it really should have been demux_30_11_20_2_illumina, but it's too late now and we will keep it this way.
        'demux_01_02_21_2', 'demux_01_02_21_2_illumina',
        'demux_22_02_21_2', 'demux_22_02_21_2_illumina', 'demux_07_03_21_1', 'demux_07_03_21_1_illumina', 'demux_11_04_21_2', 'demux_11_04_21_2_illumina',
    ],

    'donors_with_any_sample_that_waited_ids': [
        'N1', 'N264', 'N275', 'N285',
    ],
    'exps_with_donors_with_any_sample_that_waited_names': [
        'demux_01_01_23_1',
        'demux_n_15_08_23_2',
        'demux_n_15_08_23_3',
        'demux_22_06_23_2',
        'demux_10_07_22_1',
        'demux_10_07_22_2',
        'demux_10_07_22_1_ultima',
        'demux_10_07_22_2_ultima',
        'demux_13_11_22_1',
        'demux_13_11_22_2',
    ],

    'very_high_expr_range_gene_names': [
        'HBB', 'HBA2', 'HBA1', 'S100A9', 'JCHAIN', 'S100A8', 'LYZ', 'DLK1',
        'GNLY', 'PSTPIP2', 'CLC', 'NKG7', 'IGHA1', 'MPO', 'S100A12',
        'CST3', 'IGKC', 'CCL4', 'VCAN', 'LGALS1', 'MS4A1', 'DNTT', 'MGP',
        'AIF1', 'MNDA', 'CTSS', 'HLA-DRA', 'IL7R', 'TYROBP', 'S100A6',
        'LMO4', 'SOX4', 'LST1', 'FCER1G', 'MS4A2', 'TPSB2', 'SPINK2',
        'HBD', 'GZMA', 'FCN1', 'S100A11', 'FAM30A', 'CCL1', 'CCL5', 'MZB1',
        'LTB', 'S100A10', 'CAV1', 'PRTN3', 'FCGR3A', 'IGFBP5', 'KLRD1',
        'ID2', 'IL32', 'S100A4', 'IGHM', 'HPGD', 'PTPRC', 'VIM', 'SRGN',
        'ELANE', 'FGL2', 'GZMB', 'ANXA1', 'KLRB1', 'TPSAB1', 'HLA-DPA1',
        'HDC', 'FCRL5', 'AGR2', 'CYBB', 'MS4A7', 'CCL3', 'SATB1',
        'HLA-DRB1', 'RCBTB2', 'SPARCL1', 'SLC40A1', 'TCF4', 'CD2',
        'HLA-DPB1', 'CXCL8', 'CD79A', 'CD74', 'RALGPS2', 'IFI27', 'PRF1',
        'CST7', 'GZMK', 'CDK6', 'IGLL1', 'HSP90B1', 'FGFBP2', 'ACY3',
        'SEC11C', 'PCDH9', 'COTL1', 'C1QA', 'GZMH', 'BIRC3', 'MEF2C',
        'TYMP', 'C1QB', 'AZU1', 'MEG3', 'MS4A6A', 'TCL1A', 'CD52',
        'GIMAP7', 'MAF', 'TSC22D1', 'VPREB1', 'PLEK', 'HOPX', 'DUSP6',
        'KLF2', 'CD3D', 'CFD', 'EVI2B', 'CPA3', 'IFITM2', 'SLC25A37',
        'LINC01781', 'IFITM3', 'FCRL3', 'CCDC50', 'AHNAK', 'IRF8', 'ANXA2',
        'MME', 'KIAA0087', 'BANK1', 'SYNE2', 'AVP', 'HLA-DQB1', 'PSAP',
        'BASP1', 'NEAT1', 'CD3G', 'SAMHD1', 'CLU', 'CYP1B1', 'AC084033.3',
        'TPGS2', 'ZEB2', 'IQGAP2', 'SAT1', 'CD36', 'HEMGN', 'TNFRSF17',
        'ZFP36L1', 'APLP2', 'HLA-DQA1', 'GSTP1', 'NRIP1', 'IGSF6', 'ITGB2',
        'PLIN2', 'CLEC2D', 'IKZF3', 'MARCKS', 'SQLE', 'LGALS2', 'CSTA',
        'SLC2A5', 'CXCR4', 'AC020656.1', 'TM4SF1', 'PLAC8', 'RB1', 'HCST',
        'BCL2', 'XRN2', 'CEBPD', 'VMP1', 'FCER1A', 'CD247', 'ANGPT1',
        'SELL', 'MARCH1', 'RNASE2', 'GPR183', 'PLPP5', 'VPREB3', 'XIST',
        'CPVL', 'PPIB', 'BCL11A', 'ISG20', 'NAMPT', 'UGCG', 'CD14',
        'PRKCB', 'FYB1', 'GUCY1A1', 'ENO1', 'CYTL1', 'MPEG1', 'AGPAT5',
        'CASP3', 
        'PF4', 'PPBP',
    ],
    
}

all_comments = set()
for curr_comments in SC_RNA_SEQ_PREPROCESSING_PARAMS['matching_mip_to_vireo_manual_examination_comments'].values():
    all_comments |= curr_comments
SC_RNA_SEQ_PREPROCESSING_PARAMS['matching_mip_to_vireo_manual_examination_comments']['all'] = all_comments


def fix_exp_names_for_10x_batch_mixing_accidents(data_ad):
    assert 'exp_name_before_fixing_to_accident' not in data_ad.obs.columns, 'exp_name_before_fixing_to_accident already in obs.columns'
    exp_name_to_other_exp_date_part_presumably_accidentally_mixed_and_donor_ids = SC_RNA_SEQ_PREPROCESSING_PARAMS[
        'exp_name_to_other_exp_date_part_presumably_accidentally_mixed_and_donor_ids']
    data_ad.obs['exp_name'] = data_ad.obs['exp_name'].astype(str)
    data_ad.obs['exp_name_before_fixing_to_accident'] = data_ad.obs['exp_name'].astype(str)
    
    for exp_name, other_exp_date_part_presumably_accidentally_mixed_and_donor_ids in (
            exp_name_to_other_exp_date_part_presumably_accidentally_mixed_and_donor_ids.items()):
        for other_exp_date_part, accidental_donor_ids in other_exp_date_part_presumably_accidentally_mixed_and_donor_ids:
            cell_mask = (data_ad.obs['exp_name'] == exp_name) & data_ad.obs['donor_id'].isin(accidental_donor_ids)
            if cell_mask.any():
                accident_exp_name = f'{exp_name}__{other_exp_date_part}_accident'
                print(f'fixing exp name for {accidental_donor_ids} in {exp_name} ({cell_mask.sum()} cells) to {accident_exp_name}')
                data_ad.obs.loc[cell_mask, 'exp_name'] = accident_exp_name
                data_ad.obs.loc[cell_mask, 'exp_name_before_fixing_to_accident'] = exp_name

def get_exp_names_involved_in_10x_batch_mix_accidents(cells_ad):
    exp_names = set(cells_ad.obs['exp_name'].unique())
    return sorted(
        exp_names & (
            set(get_orig_exp_name_to_accidental_donor_ids()) |
            {x for x in exp_names if x.endswith('_accident')}
        )
    )

def get_exp_names_not_involved_in_10x_batch_mix_accidents(cells_ad):
    exp_names_involved_in_10x_batch_mix_accidents = get_exp_names_involved_in_10x_batch_mix_accidents(cells_ad)
    return sorted(x for x in cells_ad.obs['exp_name'].unique() if x not in exp_names_involved_in_10x_batch_mix_accidents)

def get_donor_id_and_bleeding_date_and_numbered_donor_id_df():
    df = pd.DataFrame(
        SC_RNA_SEQ_PREPROCESSING_PARAMS['list_of_donor_id_and_bleeding_date_and_numbered_donor_id'], columns=['donor_id', 'bleeding_date', 'numbered_donor_id'])
    assert not df.duplicated().any()
    return df

def get_df_with_numbered_donor_id(
        df, add_bm_according_to_exp_name=True, consider_different_exp_names_on_same_date=True, consider_different_soup_vireo_donors=False):
    if 'numbered_donor_id' in df.columns:
        print('NOTE: numbered_donor_id is already in df. overwriting it')
        df.drop(columns='numbered_donor_id', inplace=True)
    df = generic_utils.merge_preserving_df1_index_and_row_order(
        df, get_donor_id_and_bleeding_date_and_numbered_donor_id_df(), on=['donor_id', 'bleeding_date'], how='left')
    numbered_donor_id_mask = df['numbered_donor_id'].notna()
    df.loc[~numbered_donor_id_mask, 'numbered_donor_id'] = df.loc[~numbered_donor_id_mask, 'donor_id']
    df['numbered_donor_id'] = df['numbered_donor_id'].astype(str)
    
    if add_bm_according_to_exp_name:
        df.loc[df['exp_name'].str.contains('_bm_'), 'numbered_donor_id'] += '_bm'
        # print(df.loc[df['exp_name'].str.contains('_bm_'), 'numbered_donor_id'].unique())
    else:
        RuntimeError('if you can, please use add_bm_according_to_exp_name. i think it might prevent some future bugs/mistakes.')
    
    if consider_different_exp_names_on_same_date:
        # NOTE: as far as I can tell, presumably_swapped_label_exp_name_pairs shouldn't affect this, as long as we have only the two pairs we had on 230808.
        donor_exp_df = df[['numbered_donor_id', 'bleeding_date', 'exp_name']].drop_duplicates().astype(str)
        
        # numbered_donor_ids_with_multiple_exps = sorted((donor_exp_df['numbered_donor_id'].value_counts() > 1).loc[lambda x: x].index)
        list_of_old_and_new_numbered_donor_id_and_exp_name = []
        # for numbered_donor_id in numbered_donor_ids_with_multiple_exps:
        for numbered_donor_id in sorted(donor_exp_df['numbered_donor_id'].unique()): # changed to this on 240509. hopefully won't mess anything up
            curr_donor_exp_df = donor_exp_df[donor_exp_df['numbered_donor_id'] == numbered_donor_id]
            bleeding_date = curr_donor_exp_df['bleeding_date'].unique()
            assert len(bleeding_date) == 1, f'{numbered_donor_id}\n' + str(bleeding_date) + '\n' + str(curr_donor_exp_df)
            curr_exp_names = sorted(curr_donor_exp_df['exp_name'])
            list_of_old_and_new_numbered_donor_id_and_exp_name.extend(
                (numbered_donor_id, numbered_donor_id + '_' + get_exp_name_suffix(exp_name), exp_name)
                for exp_name in curr_exp_names
            )

        new_numbered_donor_id_df = pd.DataFrame(
            list_of_old_and_new_numbered_donor_id_and_exp_name, columns=['numbered_donor_id', 'new_numbered_donor_id', 'exp_name'])
        new_numbered_donor_id_df['new_numbered_donor_id'] = new_numbered_donor_id_df['new_numbered_donor_id'].str.replace(
            'illumina', 'ill').str.replace('ultima', 'ult')
        assert 'new_numbered_donor_id' not in df.columns
        df = generic_utils.merge_preserving_df1_index_and_row_order(df, new_numbered_donor_id_df, on=['numbered_donor_id', 'exp_name'], how='left')
        new_mask = df['new_numbered_donor_id'].notna()
        df.loc[new_mask, 'numbered_donor_id'] = df.loc[new_mask, 'new_numbered_donor_id']
        df.drop(columns='new_numbered_donor_id', inplace=True)
        
        hopefully_no_duplicate_df = df[['numbered_donor_id', 'bleeding_date', 'exp_name']].drop_duplicates(subset=['bleeding_date', 'exp_name'])
        problematic_numbered_donor_ids = sorted((hopefully_no_duplicate_df['numbered_donor_id'].value_counts() > 1).loc[lambda x: x].index)
        assert not problematic_numbered_donor_ids, str(
            hopefully_no_duplicate_df[hopefully_no_duplicate_df['numbered_donor_id'].isin(problematic_numbered_donor_ids)].sort_values('numbered_donor_id'))

    # print('3', df['numbered_donor_id'].unique())
    if consider_different_soup_vireo_donors:
        df['numbered_donor_id'] = df['numbered_donor_id'].astype(str) + '_' + df['soup_donor_name'].astype(str) + '_' + df['vireo_donor_name'].astype(str)
    # print('3', df['numbered_donor_id'].unique())
    donor_id_and_bleeding_date_df = df[['numbered_donor_id', 'bleeding_date']].drop_duplicates()
    unexpected_duplicate_mask = donor_id_and_bleeding_date_df.duplicated('numbered_donor_id', keep=False)
    if unexpected_duplicate_mask.any():
        print('donor_id_and_bleeding_date_df[unexpected_duplicate_mask]')
        print(donor_id_and_bleeding_date_df[unexpected_duplicate_mask])
        raise RuntimeError('should probably update list_of_donor_id_and_bleeding_date_and_numbered_donor_id accordingly')

    df['numbered_donor_id_without_tech_rep_suffix'] = get_numbered_donor_id_series_without_tech_rep_suffix(df['numbered_donor_id'])

    return df

def add_df_with_all_donor_bleeding_dates_and_with_donor_id_replaced_by_numbered_donor_id(df, exp_donor_df=None):
    assert 'bleeding_date' not in df.columns
    
    if exp_donor_df is None:
        df = df.merge(exp_donor_df, how='left')
    
    assert df['donor_id'].is_unique
    df = df.merge(get_donor_id_and_bleeding_date_and_numbered_donor_id_df()[['donor_id', 'bleeding_date']], how='left')
    df = generic_utils.merge_preserving_df1_index_and_row_order(
        df, get_nili_experiment_date_and_exp_name_df().rename(columns={'experiment date': 'bleeding_date'}))
    
    return get_df_with_numbered_donor_id(df)

def verify_donor_id_couples_in_other_donor_id_to_donor_id_i_use(df, donor_id_column_name, previous_donor_id_column_name):
    other_donor_id_to_donor_id_i_use = SC_RNA_SEQ_PREPROCESSING_PARAMS['other_donor_id_to_donor_id_i_use']
    other_donor_id_to_donor_id_i_use_as_couple_sets = {frozenset(x) for x in other_donor_id_to_donor_id_i_use.items()}
    

    previous_num_mask = ~df[previous_donor_id_column_name].isna()

    df = df[previous_num_mask].copy()
    df[previous_donor_id_column_name].replace({'N200, N180': 'N200'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N180, N200, N372': 'N200'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N367, NS21 (both ARCH only)': 'N367'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N339 (ARCH only)': 'N339'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N251, N311': 'N251'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N372, N200, N180': 'N200'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N10, N264, N379, N396': 'N264'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N10, N264, N379': 'N264'}, inplace=True)
    df[previous_donor_id_column_name].replace({'N1, N324': 'N1'}, inplace=True)

    curr_and_previous_donor_ids_as_couple_sets = {
        frozenset(x) for x in df.loc[previous_num_mask, [donor_id_column_name, previous_donor_id_column_name]].to_records(index=False).tolist()}
    # print(curr_and_previous_donor_ids_as_couple_sets)

    my_manual_donor_id_couple_sets_not_in_nilis = other_donor_id_to_donor_id_i_use_as_couple_sets - curr_and_previous_donor_ids_as_couple_sets
    nilis_donor_id_couple_sets_not_in_my_manual = curr_and_previous_donor_ids_as_couple_sets - other_donor_id_to_donor_id_i_use_as_couple_sets
    # print(nilis_donor_id_couple_sets_not_in_my_manual)

    for curr_and_previous_donor_ids_set in set(nilis_donor_id_couple_sets_not_in_my_manual):
        donor1_id, donor2_id = sorted(curr_and_previous_donor_ids_set)
        
        donor1_id_and_other_ids = set()
        donor2_id_and_other_ids = set()
        for x in other_donor_id_to_donor_id_i_use_as_couple_sets:
            if donor1_id in x:
                donor1_id_and_other_ids |= x
            if donor2_id in x:
                donor2_id_and_other_ids |= x
        
        if donor1_id_and_other_ids & donor2_id_and_other_ids:
            nilis_donor_id_couple_sets_not_in_my_manual.remove(curr_and_previous_donor_ids_set)
    

    my_manual_donor_id_couple_sets_not_in_nilis -= {
        # here because it won't enter the MDS open file anyway, i think.
        frozenset({'N285', 'N378'}),
        frozenset({'N380', 'N275'}),
        frozenset({'N379', 'N264'}),
        frozenset({'N396', 'N264'}),
        frozenset({'N398', 'N264'}),
        frozenset({'T34', 'N264'}),
        frozenset({'T35', 'N320'}),
        frozenset({'G21', 'N259'}),
        frozenset({'N163', 'N221'}),
        frozenset({'N162', 'N222'}),
        frozenset({'G13', 'N400'}),

        # TODO: remove the following. here only because haven't uploaded the mds open file yet.
        # nothing here currently.
    } 


    for curr_and_previous_donor_ids_set in set(my_manual_donor_id_couple_sets_not_in_nilis):
        donor1_id, donor2_id = sorted(curr_and_previous_donor_ids_set)
        
        donor1_id_and_other_ids = set(df.loc[
            (df[donor_id_column_name] == donor1_id) |
            (df[previous_donor_id_column_name] == donor1_id),
            [donor_id_column_name, previous_donor_id_column_name]
        ].stack())
        donor2_id_and_other_ids = set(df.loc[
            (df[donor_id_column_name] == donor2_id) |
            (df[previous_donor_id_column_name] == donor2_id),
            [donor_id_column_name, previous_donor_id_column_name]
        ].stack())
        # print(donor1_id, donor2_id, donor1_id_and_other_ids, donor2_id_and_other_ids)
        if donor1_id_and_other_ids & donor2_id_and_other_ids:
            my_manual_donor_id_couple_sets_not_in_nilis.remove(curr_and_previous_donor_ids_set)

    if nilis_donor_id_couple_sets_not_in_my_manual:
        print(f'nilis_donor_id_couple_sets_not_in_my_manual: {nilis_donor_id_couple_sets_not_in_my_manual}')
        raise RuntimeError(f'unexpected nilis_donor_id_couple_sets_not_in_my_manual ({nilis_donor_id_couple_sets_not_in_my_manual})')

    if my_manual_donor_id_couple_sets_not_in_nilis:
        print(f'my_manual_donor_id_couple_sets_not_in_nilis: {my_manual_donor_id_couple_sets_not_in_nilis}')
        raise RuntimeError(f'unexpected my_manual_donor_id_couple_sets_not_in_nilis ({my_manual_donor_id_couple_sets_not_in_nilis})')

def get_exp_name_to_swapped_label_exp_name():
    exp_name_to_swapped_label_exp_name = {}
    for exp_name_pair in SC_RNA_SEQ_PREPROCESSING_PARAMS['generate_exp_metadata_csv']['presumably_swapped_label_exp_name_pairs']:
        exp1_name, exp2_name = sorted(exp_name_pair)
        assert exp1_name not in exp_name_to_swapped_label_exp_name
        assert exp2_name not in exp_name_to_swapped_label_exp_name
        exp_name_to_swapped_label_exp_name[exp1_name] = exp2_name
        exp_name_to_swapped_label_exp_name[exp2_name] = exp1_name
    return exp_name_to_swapped_label_exp_name

def get_exp_name_suffix(exp_name):
    # NOTE: hopefully, the g or n after demux_ would not mess things up...
    match_obj = re.match(r'^(nimrod_)?demux_([ng]_)?' + DATE_PART_IN_EXP_NAME_REGEX, exp_name)
    if not match_obj:
        match_obj = re.match('^' + DONOR_ID_REGEX + '_' + DATE_PART_IN_EXP_NAME_REGEX, exp_name)
    if not match_obj:
        match_obj = re.match('^' + DONOR_ID_REGEX + '_bm_' + DATE_PART_IN_EXP_NAME_REGEX, exp_name)
    if not match_obj:
        raise RuntimeError(f'exp_name {exp_name} does not match any of the expected formats')
    return exp_name[len(match_obj.group()):]

def get_exp_name_to_exp_date(exp_names):
    exp_name_to_exp_date = {}
    exp_name_to_swapped_label_exp_name = get_exp_name_to_swapped_label_exp_name()
    for exp_name in exp_names:
        if exp_name in exp_name_to_swapped_label_exp_name:
            exp_name_without_suffix = exp_name_to_swapped_label_exp_name[exp_name]
        else:
            exp_name_without_suffix = exp_name
        if '_ult_batch_' in exp_name_without_suffix:
            exp_name_without_suffix, _, ultima_batch_name = exp_name_without_suffix.partition('_ult_batch_')
            assert re.fullmatch(r'[0-9]{6}', ultima_batch_name)

        for suffix_to_remove in ('_accident', '_ultima', '_illumina', '_elembio'):
            if exp_name_without_suffix.endswith(suffix_to_remove):
                exp_name_without_suffix = exp_name_without_suffix.partition(suffix_to_remove)[0]
                if suffix_to_remove == '_accident':
                    assert '__' in exp_name_without_suffix
                    exp_name_without_suffix = '__' + exp_name_without_suffix.partition('__')[2] # add '__' is an ugly hack to make thing look normal below...


        if re.fullmatch(r'_[123]', exp_name_without_suffix[-2:]):
            exp_name_without_suffix = exp_name_without_suffix[:-2]
        date_part_template = '_dd_mm_yy'
        date_part_len = len(date_part_template)
        # print(exp_name)
        # print(exp_name_without_suffix)
        if exp_name_without_suffix == '__unknown':
            bleeding_date = 'unknown'
        else:
            assert len(exp_name_without_suffix) > date_part_len
            exp_name_date_part = exp_name_without_suffix[-date_part_len:]
        
            # print(exp_name, exp_name_date_part)
            assert re.fullmatch(r'_[0-9]{2}_[0-9]{2}_[0-9]{2}', exp_name_date_part), exp_name_date_part
            bleeding_date = exp_name_date_part[1:].replace('_', '.')
        exp_name_to_exp_date[exp_name] = bleeding_date
    return exp_name_to_exp_date

def add_exp_date_column_to_df(df):
    exp_names = list(df['exp_name'].unique())
    df['exp_date'] = df['exp_name'].map(get_exp_name_to_exp_date(exp_names))
    assert not df['exp_date'].isna().any()
    df['exp_date'] = df['exp_date'].astype(str)


def get_orig_exp_name_to_accidental_donor_ids():
    orig_exp_name_to_accidental_donor_ids = {}
    for orig_exp_name, other_exp_date_part_presumably_accidentally_mixed_and_donor_ids in SC_RNA_SEQ_PREPROCESSING_PARAMS[
            'exp_name_to_other_exp_date_part_presumably_accidentally_mixed_and_donor_ids'].items():
        accidental_donor_ids = set()
        for _, donor_ids in other_exp_date_part_presumably_accidentally_mixed_and_donor_ids:
            accidental_donor_ids |= set(donor_ids)
        orig_exp_name_to_accidental_donor_ids[orig_exp_name] = sorted(accidental_donor_ids)
    return orig_exp_name_to_accidental_donor_ids

def get_orig_exp_name_to_accidental_exp_date_parts():
    orig_exp_name_to_accidental_exp_date_parts = {}
    for orig_exp_name, other_exp_date_part_presumably_accidentally_mixed_and_donor_ids in SC_RNA_SEQ_PREPROCESSING_PARAMS[
            'exp_name_to_other_exp_date_part_presumably_accidentally_mixed_and_donor_ids'].items():
        accidental_exp_date_parts = set()
        for accidental_exp_date_part, _ in other_exp_date_part_presumably_accidentally_mixed_and_donor_ids:
            accidental_exp_date_parts.add(accidental_exp_date_part)
        orig_exp_name_to_accidental_exp_date_parts[orig_exp_name] = sorted(accidental_exp_date_parts)
    return orig_exp_name_to_accidental_exp_date_parts

def add_exp_and_bleeding_date_according_to_exp_name_and_donor_id(df):
    assert 'exp_name' in df.columns
    assert 'donor_id' in df.columns
    
    df['exp_date'] = df['exp_name'].map(get_exp_name_to_exp_date(df['exp_name'].unique()))
    assert not df['exp_date'].isna().any()
    df['bleeding_date'] = df['exp_date']

    for exp_name_or_exp_name_and_donor_id, manually_added_cell_attrs in SC_RNA_SEQ_PREPROCESSING_PARAMS[
            'exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs'].items():
        if 'bleeding_date' not in manually_added_cell_attrs:
            continue
        if isinstance(exp_name_or_exp_name_and_donor_id, str):
            exp_name = exp_name_or_exp_name_and_donor_id
            mask = df['exp_name'] == exp_name
        else:
            assert isinstance(exp_name_or_exp_name_and_donor_id, tuple)
            assert len(exp_name_or_exp_name_and_donor_id) == 2
        
            exp_name, donor_id = exp_name_or_exp_name_and_donor_id

            mask = (df['exp_name'] == exp_name) & (df['donor_id'] == donor_id)
        
        df.loc[mask, 'bleeding_date'] = manually_added_cell_attrs['bleeding_date']
    

def get_numbered_donor_id_series_without_tech_rep_suffix(numbered_donor_id_series):
    return numbered_donor_id_series.str.partition('__')[0]


def get_tech_rep_numbered_donor_ids_to_discard_by_cell_count(c_ad, c_mask=None):
    if c_mask is None:
        c_mask = np.full(c_ad.n_obs, True)

    c_ad_obs = c_ad.obs.loc[c_mask, ['donor_id', 'bleeding_date', 'diagnosis', 'exp_name']]
    c_ad_obs = get_df_with_numbered_donor_id(c_ad_obs, add_bm_according_to_exp_name=True, consider_different_exp_names_on_same_date=True)
    numbered_donor_id_to_cell_count = c_ad_obs['numbered_donor_id'].astype(str).value_counts().to_dict()
    
    df = c_ad_obs.drop_duplicates().copy()
    # print(df)
    df['exp_name'] = df['exp_name'].astype(str)
    common_tech_rep_cols = ['donor_id', 'bleeding_date', 'diagnosis']
    df[common_tech_rep_cols] = df[common_tech_rep_cols].astype(str)
    multiple_exp_donor_id_bleeding_date_diagnosis_df = df.groupby(common_tech_rep_cols)['exp_name'].nunique().loc[
        lambda x: x > 1].index.to_frame().reset_index(drop=True)
    # print(multiple_exp_donor_id_bleeding_date_diagnosis_df)
    df = generic_utils.get_df1_masked_by_inner_join_with_df2_preserving_df1_row_order(df, multiple_exp_donor_id_bleeding_date_diagnosis_df)
    numbered_donor_ids_to_discard_candidates = list(df['numbered_donor_id'])
    assert len(numbered_donor_ids_to_discard_candidates) == len(set(numbered_donor_ids_to_discard_candidates))
    
    df['cell_count'] = df['numbered_donor_id'].map(numbered_donor_id_to_cell_count)
    assert df['cell_count'].notna().all()
    
    best_tech_rep_numbered_donor_ids = set(df.sort_values('cell_count', ascending=False).drop_duplicates(subset=common_tech_rep_cols, keep='first')['numbered_donor_id'])
    
    return sorted(set(numbered_donor_ids_to_discard_candidates) - best_tech_rep_numbered_donor_ids)

def get_ultima_exp_name_to_same_10x_lib_non_ultima_exp_names(
        exp_names, extra_ultima_exp_name_to_same_10x_lib_illumina_exp_names={}, 
        each_value_is_a_single_exp=False,
):
    ultima_exp_name_to_same_10x_lib_illumina_exp_names = collections.defaultdict(set)
    for exp_name1, exp_name2 in itertools.product(exp_names, repeat=2):
        if (exp_name1 == exp_name2 + '_ultima') or (exp_name2 == exp_name1 + '_illumina'):
            ultima_exp_name_to_same_10x_lib_illumina_exp_names[exp_name1].add(exp_name2)

    for exp_name1, exp_name2 in extra_ultima_exp_name_to_same_10x_lib_illumina_exp_names.items():
        ultima_exp_name_to_same_10x_lib_illumina_exp_names[exp_name1] |= exp_name2

    ultima_exp_name_to_same_10x_lib_illumina_exp_names = dict(ultima_exp_name_to_same_10x_lib_illumina_exp_names) # I don't want a defaultdict moving around.

    if each_value_is_a_single_exp:
        assert all(len(x) == 1 for x in ultima_exp_name_to_same_10x_lib_illumina_exp_names.values())
        ultima_exp_name_to_same_10x_lib_illumina_exp_names = {k: next(iter(v)) for k, v in ultima_exp_name_to_same_10x_lib_illumina_exp_names.items()}


    return ultima_exp_name_to_same_10x_lib_illumina_exp_names

def get_all_ult_ill_tech_rep_exp_names(exp_names, extra_ultima_exp_name_to_same_10x_lib_illumina_exp_names={}):
    ult_to_ill = get_ultima_exp_name_to_same_10x_lib_non_ultima_exp_names(exp_names, extra_ultima_exp_name_to_same_10x_lib_illumina_exp_names)
    ill_exp_names = set()
    for x in ult_to_ill.values():
        assert len(x) == 1
        ill_exp_names |= x
    return sorted(set(ult_to_ill) | ill_exp_names)

def get_mask_of_cells_allowing_pooling_all_pb_cells_per_donor(c_ad_obs, always_choose_non_ultima=True):
    mask = (
        (~c_ad_obs['numbered_donor_id'].isin(SC_RNA_SEQ_PREPROCESSING_PARAMS['delayed_sample_numbered_donor_ids'])) 
        & (~c_ad_obs['numbered_donor_id'].str.contains('_bm'))
    )
    exp_name_to_median_num_of_non_excluded_umis = c_ad_obs[mask].groupby('exp_name')['num_of_non_excluded_umis'].median().to_dict()
    for exp_name, non_ult_seq_tech_rep_exp_names in get_ultima_exp_name_to_same_10x_lib_non_ultima_exp_names(sorted(c_ad_obs.loc[mask, 'exp_name'].unique())).items():
        all_curr_exp_names = [exp_name, *non_ult_seq_tech_rep_exp_names]
        candidate_exp_names = non_ult_seq_tech_rep_exp_names if always_choose_non_ultima else all_curr_exp_names
        chosen_exp_name = max(candidate_exp_names, key=lambda x: exp_name_to_median_num_of_non_excluded_umis[x])
        not_chosen_exp_names = [x for x in all_curr_exp_names if x != chosen_exp_name]
        print(exp_name, chosen_exp_name, mask.sum(), candidate_exp_names, not_chosen_exp_names)
        mask &= ~(c_ad_obs['exp_name'].isin(not_chosen_exp_names))
        print(exp_name, chosen_exp_name, mask.sum())
    return mask

def get_exp_name_to_cell_ranger_out_path(metadata_df, exp_names):
    assert metadata_df['exp_name'].is_unique
    exp_name_to_raw_feature_mat_path = {
        exp_name: metadata_df.loc[metadata_df['exp_name'] == exp_name, 'tenx_path'].iloc[0]
        for exp_name in exp_names
    }
    assert len(set(exp_name_to_raw_feature_mat_path.values())) == len(exp_names)

    return exp_name_to_raw_feature_mat_path

def get_exp_name_to_manual_biased_composition_type():
    exp_name_to_manual_biased_composition_type = {}
    for bias_type, exp_name_to_stuff in SC_RNA_SEQ_PREPROCESSING_PARAMS['type_to_exp_name_to_forced_num_of_filtered_barcodes_and_min_umi_count'].items():
        for exp_name in exp_name_to_stuff:
            exp_name_to_manual_biased_composition_type[exp_name] = bias_type

    return exp_name_to_manual_biased_composition_type

def get_exp_name_to_min_umi_count():
    exp_name_to_min_umi_count = {}

    all_types = (
        SC_RNA_SEQ_PREPROCESSING_PARAMS['umi_count_per_barcode_exp_types_that_can_be_used_when_analyzing_composition']
        | SC_RNA_SEQ_PREPROCESSING_PARAMS['umi_count_per_barcode_exp_types_to_ignore_when_analyzing_composition']
    )
    missing_from_all_min_fn_and_min_fp_types = (
        set(SC_RNA_SEQ_PREPROCESSING_PARAMS['type_to_exp_name_to_forced_num_of_filtered_barcodes_and_min_umi_count']) -
        all_types
    )
    assert not missing_from_all_min_fn_and_min_fp_types, str(missing_from_all_min_fn_and_min_fp_types)

    for problematic_type, type_dict in SC_RNA_SEQ_PREPROCESSING_PARAMS['type_to_exp_name_to_forced_num_of_filtered_barcodes_and_min_umi_count'].items():
        if problematic_type == 'give_up_on':
            print('sc_rna_seq_preprocessing_params.py: get_exp_name_to_min_umi_count: skipping give_up_on exps. make sure they will be excluded later')
            continue
        # print(problematic_type)
        for exp_name, (_, min_umi_count) in type_dict.items():
            exp_name_to_min_umi_count[exp_name] = min_umi_count

    return exp_name_to_min_umi_count

def get_biased_composition_due_to_discarded_low_umi_count_barcodes_exp_names():
    assert not (
        SC_RNA_SEQ_PREPROCESSING_PARAMS['umi_count_per_barcode_exp_types_to_ignore_when_analyzing_composition'] &
        SC_RNA_SEQ_PREPROCESSING_PARAMS['umi_count_per_barcode_exp_types_that_can_be_used_when_analyzing_composition']
    )
    assert (
        SC_RNA_SEQ_PREPROCESSING_PARAMS['umi_count_per_barcode_exp_types_to_ignore_when_analyzing_composition'] |
        SC_RNA_SEQ_PREPROCESSING_PARAMS['umi_count_per_barcode_exp_types_that_can_be_used_when_analyzing_composition']
    ) == set(SC_RNA_SEQ_PREPROCESSING_PARAMS['type_to_exp_name_to_forced_num_of_filtered_barcodes_and_min_umi_count'])
    biased_composition_exp_names = set()
    for umi_count_type in SC_RNA_SEQ_PREPROCESSING_PARAMS['umi_count_per_barcode_exp_types_to_ignore_when_analyzing_composition']:
        biased_composition_exp_names |= set(SC_RNA_SEQ_PREPROCESSING_PARAMS[
            'type_to_exp_name_to_forced_num_of_filtered_barcodes_and_min_umi_count'][umi_count_type])
    return sorted(biased_composition_exp_names)

# NOTE: seems like we don't need this...            
# def get_donor_id_and_exp_date_to_forced_attr_dict():
#     donor_id_and_bleeding_date_to_forced_attr_dict = {
#         k: {'diagnosis': v}
#         for k,v in SC_RNA_SEQ_PREPROCESSING_PARAMS['clinical_table_donor_id_and_exp_date_to_forced_diagnosis'].items()
#     }
    
#     for exp_name_or_exp_name_and_donor_id, manually_added_cell_attrs in SC_RNA_SEQ_PREPROCESSING_PARAMS[
#             'exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs'].items():
#         # if (not isinstance(exp_name_or_exp_name_and_donor_id, tuple)) or ('diagnosis' not in manually_added_cell_attrs): # 231109: why did i add the diagnosis part here??
#         if not isinstance(exp_name_or_exp_name_and_donor_id, tuple):
#             assert isinstance(exp_name_or_exp_name_and_donor_id, str)
#             if 'bleeding_date' in manually_added_cell_attrs:
#                 raise RuntimeError('please manually add bleeding_date per donor, not for the whole experiment.')
#             continue
#         assert len(exp_name_or_exp_name_and_donor_id) == 2
        
#         exp_name, donor_id = exp_name_or_exp_name_and_donor_id
#         if 'bleeding_date' in manually_added_cell_attrs:
#             bleeding_date = manually_added_cell_attrs['bleeding_date']
#         else:
#             bleeding_date = get_exp_name_to_bleeding_date([exp_name])[exp_name]
#         donor_id_and_bleeding_date = (donor_id, bleeding_date)
#         if donor_id_and_bleeding_date not in donor_id_and_bleeding_date_to_forced_attr_dict:
#             donor_id_and_bleeding_date_to_forced_attr_dict[donor_id_and_bleeding_date] = {}
        
#         for attr, attr_val in manually_added_cell_attrs.items():
#             if attr == 'bleeding_date':
#                 continue
#             if attr in donor_id_and_bleeding_date_to_forced_attr_dict[donor_id_and_bleeding_date]:
#                 assert donor_id_and_bleeding_date_to_forced_attr_dict[donor_id_and_bleeding_date][attr] == attr_val, f'{donor_id_and_bleeding_date}, {attr}, {donor_id_and_bleeding_date_to_forced_attr_dict[donor_id_and_bleeding_date][attr]}, {attr_val}'
#             else:
#                 donor_id_and_bleeding_date_to_forced_attr_dict[donor_id_and_bleeding_date][attr] = attr_val
#     return donor_id_and_bleeding_date_to_forced_attr_dict

def get_single_donor_exp_names():
    exp_names = set()
    for exp_name_or_exp_name_and_donor_id, manually_added_cell_attrs in (
        SC_RNA_SEQ_PREPROCESSING_PARAMS['exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs'].items()
    ):
        if not isinstance(exp_name_or_exp_name_and_donor_id, str):
            continue
        if 'donor_id' in manually_added_cell_attrs:
            exp_names.add(exp_name_or_exp_name_and_donor_id)

    return sorted(exp_names)

def overwrite_according_to_exp_name_and_donor_id(
    df, 
    exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs=SC_RNA_SEQ_PREPROCESSING_PARAMS[
        'exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs'], 
    column_names_to_set=None,
):
    df = df.copy()
    if column_names_to_set != {'donor_id'}:
        # because this might affect later overwriting.
        overwrite_according_to_exp_name_and_donor_id(
            df, 
            exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs=exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs, 
            column_names_to_set={'donor_id'},
        )
    for exp_name_or_exp_name_and_donor_id, manually_added_cell_attrs in exp_name_or_exp_name_and_donor_id_to_manually_added_cell_attrs.items():
        if isinstance(exp_name_or_exp_name_and_donor_id, tuple):
            assert len(exp_name_or_exp_name_and_donor_id) == 2
            exp_name, donor_id = exp_name_or_exp_name_and_donor_id
            cell_mask = (df['exp_name'] == exp_name) & (df['donor_id'] == donor_id)
        else:
            assert isinstance(exp_name_or_exp_name_and_donor_id, str)
            exp_name = exp_name_or_exp_name_and_donor_id
            cell_mask = df['exp_name'] == exp_name
            donor_id = None
        if not cell_mask.any():
            continue
        for column_name, val in manually_added_cell_attrs.items():
            if (column_names_to_set is not None) and (column_name not in column_names_to_set):
                continue
            if df[column_name].dtype.name == 'category':
                assert not df[column_name].isna().any()
                print(f'NOTE: converting {column_name} to str')
                df[column_name] = df[column_name].astype(str)
            df.loc[cell_mask, column_name] = val
            if column_name == 'donor_id':
                assert donor_id is None
                # print(f'for {exp_name} cells, forcing donor_id={val}, so accordingly setting soup_donor_name and vireo_donor_name, soup_or_vireo_doublet, soup_or_vireo_doublet_or_unassigned_but_not_manually_confirmed and mip_to_10x_match_confidence.')
                # NOTE: i don't see why to overwrite soup_donor_name and vireo_donor_name, so i'm commenting it out.
                # if 'soup_donor_name' in df.columns:
                #     df.loc[cell_mask, 'soup_donor_name'] = 'dummy'
                # if 'vireo_donor_name' in df.columns:
                #     df.loc[cell_mask, 'vireo_donor_name'] = 'dummy'
                if 'genotype_doublet' in df.columns:
                    print(f'for {exp_name} cells, forcing donor_id={val}, so accordingly setting genotype_doublet=False')
                    df.loc[cell_mask, 'genotype_doublet'] = False
                if 'mip_to_10x_match_confidence' in df.columns:
                    print(f'for {exp_name} cells, forcing donor_id={val}, so accordingly setting mip_to_10x_match_confidence=forced_donor_id')
                    df.loc[cell_mask, 'mip_to_10x_match_confidence'] = 'forced_donor_id'
    
    return df


def get_barcode_set_from_barcode_file(barcode_file_path):
    import csv
    import gzip
    if barcode_file_path.endswith('.tsv.gz'):
        return {row[0] for row in csv.reader(gzip.open(barcode_file_path, mode="rt"), delimiter="\t")}
    return set(generic_utils.read_text_file(barcode_file_path).split())

def get_c_ad_from_pipseq_filtered_matrix(filtered_matrix_dir_path):
    import csv
    import gzip
    import scipy.io
    import anndata as ad

    barcode_gz_file_path = os.path.join(filtered_matrix_dir_path, 'barcodes.tsv.gz')
    barcode_df = pd.DataFrame(csv.reader(gzip.open(barcode_gz_file_path, mode="rt"), delimiter="\t"), columns=['barcode'])
    barcode_df.index = barcode_df['barcode'].to_numpy()

    feature_gz_file_path = os.path.join(filtered_matrix_dir_path, 'features.tsv.gz')
    feature_df = pd.DataFrame(csv.reader(gzip.open(feature_gz_file_path, mode="rt"), delimiter="\t"), columns=['gene_id', 'gene', 'type'])
    assert (feature_df['type'] == 'Gene Expression').all()
    feature_df.drop(columns=['type'], inplace=True)
    feature_df.index = feature_df['gene'].to_numpy()

    mtx_gz_file_path = os.path.join(filtered_matrix_dir_path, 'matrix.mtx.gz')
    print('starting to run scipy.io.mmread()')
    mat = scipy.io.mmread(mtx_gz_file_path).T.tocsr()

    c_ad = ad.AnnData(X=mat, obs=barcode_df, var=feature_df)
    c_ad.var_names_make_unique(join='---')
    assert c_ad.X.has_canonical_format 

    return c_ad

def get_ordered_exp_names(exp_names):
    assert len(set(exp_names)) == len(exp_names)

    exp_name_to_date_tuple = {}
    for exp_name in exp_names:
        date_repr = re.compile('_(\d\d_\d\d_\d\d)').search(exp_name).group(1)
        # print(exp_name, date_repr)
        date_tuple = tuple(int(x) for x in date_repr.split('_')[::-1])
        exp_name_to_date_tuple[exp_name] = date_tuple
    return sorted(exp_name_to_date_tuple, key=lambda x: (*exp_name_to_date_tuple[x], x))
    # re.compile('\["Chemistry","([^\]]*)"\]').search(cell_ranger_summary_html_as_str).group(1)

def get_ordered_exp_names_from_df_with_bleeding_dates(df):
    exp_name_to_date_tuple = {}
    for _, row in df[['bleeding_date', 'exp_name']].drop_duplicates().iterrows():
        exp_name = row['exp_name']
        bleeding_date = row['bleeding_date']
        date_tuple = tuple(int(x) for x in bleeding_date.split('.')[::-1])
        assert exp_name not in exp_name_to_date_tuple, f'{exp_name} is already in exp_name_to_date_tuple ({exp_name_to_date_tuple})'
        exp_name_to_date_tuple[exp_name] = date_tuple
    return sorted(exp_name_to_date_tuple, key=lambda x: (*exp_name_to_date_tuple[x], x))


def get_misc_numbered_donor_info(c_ad):
    import matplotlib
    import matplotlib.pyplot as plt

    exp_names = get_ordered_exp_names(list(c_ad.obs['exp_name'].unique()))
    tech_rep_numbered_donor_ids_to_discard = get_tech_rep_numbered_donor_ids_to_discard_by_cell_count(c_ad)

    numbered_donor_df = c_ad.obs[[
        'numbered_donor_id', 'donor_id', 'bleeding_date', 'exp_name', 'is_bm', 'donor_age', 'donor_sex',
        'numbered_donor_id_without_tech_rep_suffix', 'is_ultima', 'diagnosis', 'diagnosis_class',
    ]].drop_duplicates().reset_index(drop=True)
    numbered_donor_df = generic_utils.merge_preserving_df1_index_and_row_order(
        numbered_donor_df,
        c_ad.obs['numbered_donor_id'].value_counts().reset_index(name='cell_count'),
    )
    numbered_donor_df['log_cell_count'] = np.log2(numbered_donor_df['cell_count'])
    
    numbered_donor_id_to_exp_name = generic_utils.get_dict_mapping_one_df_column_to_other(numbered_donor_df, 'numbered_donor_id', 'exp_name')
    numbered_donor_id_to_is_ultima = generic_utils.get_dict_mapping_one_df_column_to_other(numbered_donor_df, 'numbered_donor_id', 'is_ultima')
    numbered_donor_id_to_diagnosis = generic_utils.get_dict_mapping_one_df_column_to_other(numbered_donor_df, 'numbered_donor_id', 'diagnosis')
    numbered_donor_id_to_diagnosis_class = generic_utils.get_dict_mapping_one_df_column_to_other(numbered_donor_df, 'numbered_donor_id', 'diagnosis_class')
    numbered_donor_id_to_numbered_donor_id_without_tech_rep_suffix = generic_utils.get_dict_mapping_one_df_column_to_other(
        numbered_donor_df, 'numbered_donor_id', 'numbered_donor_id_without_tech_rep_suffix')
    
    numbered_donor_id_to_is_ultima_color = {k: 'black' if v else 'white' for k, v in numbered_donor_id_to_is_ultima.items()}

    add_exp_date_column_to_df(numbered_donor_df)
    numbered_donor_df['exp_days_since_2020'] = generic_utils.convert_two_dot_two_dot_two_date_col_to_days_since_ref_date(
        numbered_donor_df, 'exp_date', (2020, 1, 1))
    binary_cmap = plt.get_cmap('binary')
    numbered_donor_df['norm_exp_days_since_2020'] = numbered_donor_df['exp_days_since_2020']
    generic_utils.normalize_df_columns(
        numbered_donor_df, column_names=['norm_exp_days_since_2020'], 
        col_to_min_and_max={'norm_exp_days_since_2020': (0, 365*5)},
        # allow_auto_min_max=True,
        inplace=True,
    )
    numbered_donor_df['exp_days_since_2020_color'] = numbered_donor_df['norm_exp_days_since_2020'].apply(
        lambda x: matplotlib.colors.to_hex(binary_cmap(x)))
    numbered_donor_id_to_exp_date_color = generic_utils.get_dict_mapping_one_df_column_to_other(
        numbered_donor_df, 'numbered_donor_id', 'exp_days_since_2020_color')

    exp_name_to_color = {
        k: v
        for k,v in zip(exp_names, generic_utils.get_n_colors(len(exp_names)))
    }
    numbered_donor_id_to_exp_color = {k: exp_name_to_color[v] for k,v in numbered_donor_id_to_exp_name.items()}

    numbered_donor_df['donor_age_color'] = numbered_donor_df['donor_age'].apply(
        lambda x: matplotlib.colors.to_hex(binary_cmap(x)))
    # {'male': 'blue', 'female': 'red'}

    color_df = generic_utils.get_color_df(
        numbered_donor_df[['donor_age', 'donor_sex', 'log_cell_count']], allow_nans=True, 
        col_to_cmap_and_min_and_max={
            'donor_age': (binary_cmap, 30, 80),
            'log_cell_count': (binary_cmap, 8, 13),
        },
        col_to_val_to_color={
            'donor_sex': {'male': 'blue', 'female': 'red'},
        },
    )
    color_df['numbered_donor_id'] = numbered_donor_df['numbered_donor_id']
    numbered_donor_id_to_age_color = generic_utils.get_dict_mapping_one_df_column_to_other(
        color_df, 'numbered_donor_id', 'donor_age')
    numbered_donor_id_to_sex_color = generic_utils.get_dict_mapping_one_df_column_to_other(
        color_df, 'numbered_donor_id', 'donor_sex')
    numbered_donor_id_to_log_cell_count_color = generic_utils.get_dict_mapping_one_df_column_to_other(
        color_df, 'numbered_donor_id', 'log_cell_count')
    numbered_donor_df['tech_rep_to_discard'] = numbered_donor_df['numbered_donor_id'].isin(tech_rep_numbered_donor_ids_to_discard)
    numbered_donor_df.drop(columns='numbered_donor_id', inplace=True) # numbered_donor_id is dangerous because it might not be the same if you take different tech reps... 

    return dict(
        numbered_donor_df=numbered_donor_df,

        numbered_donor_id_to_exp_name=numbered_donor_id_to_exp_name,
        numbered_donor_id_to_age_color=numbered_donor_id_to_age_color,
        numbered_donor_id_to_sex_color=numbered_donor_id_to_sex_color,
        numbered_donor_id_to_log_cell_count_color=numbered_donor_id_to_log_cell_count_color,
        numbered_donor_id_to_is_ultima=numbered_donor_id_to_is_ultima,
        numbered_donor_id_to_is_ultima_color=numbered_donor_id_to_is_ultima_color,
        numbered_donor_id_to_diagnosis=numbered_donor_id_to_diagnosis,
        numbered_donor_id_to_diagnosis_class=numbered_donor_id_to_diagnosis_class,
        numbered_donor_id_to_numbered_donor_id_without_tech_rep_suffix=numbered_donor_id_to_numbered_donor_id_without_tech_rep_suffix,
        numbered_donor_id_to_exp_date_color=numbered_donor_id_to_exp_date_color,
        numbered_donor_id_to_exp_color=numbered_donor_id_to_exp_color,
    )

def combine_misc_numbered_donor_infos(misc_numbered_donor_infos):
    numbered_donor_df = pd.concat([x['numbered_donor_df'] for x in misc_numbered_donor_infos], ignore_index=True)
    
    assert not numbered_donor_df[['donor_id', 'exp_name']].duplicated().any(), 'misc_numbered_donor_infos must be mutual exclusive'
    # assert numbered_donor_df['numbered_donor_id'].is_unique, 

    combined_misc_numbered_donor_info = dict(
        numbered_donor_df=numbered_donor_df,
    )
    for dict_name in misc_numbered_donor_infos[0]:
        if dict_name == 'numbered_donor_df':
            continue
        combined_misc_numbered_donor_info[dict_name] = generic_utils.merge_mutual_exclusive_dicts([x[dict_name] for x in misc_numbered_donor_infos])
    return combined_misc_numbered_donor_info

def get_donor_id_to_latest_diagnosis_class(combined_numbered_donor_df):
    assert 'temp_bleeding_date_as_date' not in combined_numbered_donor_df.columns
    combined_numbered_donor_df['temp_bleeding_date_as_date'] = pd.to_datetime(combined_numbered_donor_df['bleeding_date'], format='%d.%m.%y', errors='raise')
    res_dict = generic_utils.get_dict_mapping_one_df_column_to_other(
        combined_numbered_donor_df.sort_values('temp_bleeding_date_as_date').drop_duplicates(subset='donor_id', keep='last'),
        'donor_id', 'diagnosis_class',
    )
    combined_numbered_donor_df.drop(columns='temp_bleeding_date_as_date', inplace=True)
    return res_dict

def add_latest_diagnosis_class_col(combined_numbered_donor_df):
    combined_numbered_donor_df['latest_diagnosis_class'] = combined_numbered_donor_df['donor_id'].map(get_donor_id_to_latest_diagnosis_class(combined_numbered_donor_df))

def write_misc_numbered_donor_info_pickle(c_ad, pickle_file_path):
    misc_numbered_donor_info = get_misc_numbered_donor_info(c_ad)
    with open(pickle_file_path, 'wb') as f:
        pickle.dump(misc_numbered_donor_info, f, protocol=5)


def get_earliest_and_latest_bio_rep_df(numbered_donor_df, allow_tech_reps_in_received_donor_df=False):
    if not allow_tech_reps_in_received_donor_df:
        assert numbered_donor_df['numbered_donor_id_without_tech_rep_suffix'].is_unique, f'set allow_tech_reps_in_received_donor_df=True if you wish...'
    # NOTE: if you wish to compare 2nd to 3rd bio reps, you could, for example, only keep the latest 2 of each donor, before calling this func.
    donors_with_bio_reps_ids = (numbered_donor_df.groupby('donor_id')['numbered_donor_id_without_tech_rep_suffix'].nunique() >= 2).loc[lambda x: x].index
    numbered_donor_df = numbered_donor_df[numbered_donor_df['donor_id'].isin(donors_with_bio_reps_ids)].copy()
    numbered_donor_df['bio_rep'] = numbered_donor_df['numbered_donor_id_without_tech_rep_suffix'].str.rpartition('_', expand=True)[2].astype(int)
    earliest_bio_rep_df = numbered_donor_df.sort_values('bio_rep', ascending=True).drop_duplicates(subset='donor_id', keep='first')
    latest_bio_rep_df = numbered_donor_df.sort_values('bio_rep', ascending=False).drop_duplicates(subset='donor_id', keep='first')
    bio_rep_df = generic_utils.merge_preserving_df1_index_and_row_order(earliest_bio_rep_df, latest_bio_rep_df, on='donor_id', suffixes=('_early', '_late'))
    return bio_rep_df

def get_best_or_random_numbered_donor_id_per_donor_or_donor_sample(
        numbered_donor_df, sort_by_cols=['biased_composition_due_to_discarded_low_umi_count_barcodes', 'cell_count'], sort_by_ascending=[True, False], numbered_donor_ids=None, single_nd_id_per_donor=False, random_state=0):
    assert numbered_donor_df['numbered_donor_id'].is_unique
    curr_df = numbered_donor_df[['numbered_donor_id', 'numbered_donor_id_without_tech_rep_suffix', 'donor_id'] + sort_by_cols]
    if numbered_donor_ids is not None:
        curr_df = curr_df[curr_df['numbered_donor_id'].isin(numbered_donor_ids)]
    
    # print(curr_df)
    if sort_by_cols is not None:
        assert len(sort_by_cols) == len(sort_by_ascending)
        curr_df = curr_df.sort_values(sort_by_cols, ascending=sort_by_ascending)
    else:
        curr_df = curr_df.sample(frac=1, random_state=random_state)

    subset_col = 'donor_id' if single_nd_id_per_donor else 'numbered_donor_id_without_tech_rep_suffix'
    return list(curr_df.drop_duplicates(subset=subset_col, keep='first')['numbered_donor_id'])

def get_illu_ult_tech_rep_exp_sets(illu_ult_numbered_donor_df):
    nd_ids_without_suffix_with_tech_reps = (illu_ult_numbered_donor_df[['numbered_donor_id_without_tech_rep_suffix', 'is_ultima']].drop_duplicates()[
        'numbered_donor_id_without_tech_rep_suffix'].value_counts() >= 2).loc[lambda x: x].index

    illu_ult_tech_rep_exp_sets = []
    for nd_id_without_suffix in nd_ids_without_suffix_with_tech_reps:
        curr_exps = set(illu_ult_numbered_donor_df.loc[illu_ult_numbered_donor_df['numbered_donor_id_without_tech_rep_suffix'] == nd_id_without_suffix, 'exp_name'])
        intersecting_set_i = None
        for i, exp_set in enumerate(illu_ult_tech_rep_exp_sets):
            if curr_exps & exp_set:
                intersecting_set_i = i
                break
        if intersecting_set_i is not None:
            illu_ult_tech_rep_exp_sets[intersecting_set_i] |= curr_exps
        else:
            illu_ult_tech_rep_exp_sets.append(curr_exps)
    # print(f'illu_ult_tech_rep_exp_sets: {illu_ult_tech_rep_exp_sets}')
    for exp_set in illu_ult_tech_rep_exp_sets:
        # print(exp_set, '\n', illu_ult_numbered_donor_df[illu_ult_numbered_donor_df['exp_name'].isin(exp_set)].drop_duplicates(subset='donor_id')['diagnosis_class'].value_counts())
        assert exp_set in [
            {
                'demux_02_05_22_2',
                'demux_02_05_22_1',
                'demux_02_05_22_2_ultima',
            },
            {
                'demux_07_03_21_1_illumina', 
                'demux_07_03_21_1',
            },
            {
                'demux_01_02_21_2',
                'demux_01_02_21_2_illumina',
            },
            {
                'demux_22_02_21_2_illumina',
                'demux_22_02_21_2',
            },
            {
                'demux_11_04_21_2',
                'demux_11_04_21_2_illumina',
            },
            {
                'demux_06_06_22_2_ultima',
                'demux_06_06_22_1_ultima',
                'demux_06_06_22_1',
                'demux_06_06_22_2',
            },
            {
                'demux_12_06_22_1',
                'demux_12_06_22_1_ultima',
            },
            {
                'demux_03_07_22_2_ultima',
                'demux_03_07_22_2',
            },
            {
                'demux_09_06_22_1',
                'demux_09_06_22_2_ultima',
                'demux_09_06_22_2',
            },
        ]
    return illu_ult_tech_rep_exp_sets

def get_illu_and_ult_blacklist_nd_ids(
        illu_ult_numbered_donor_df, target_training_set_frac=0.7, train_platform='ult', illu_ult_tech_rep_exp_sets_to_assign_to_ult_frac=1, 
        blacklist_anyway_numbered_donor_ids=None,
        random_state=0,
        verbose=True,
):
    # logic:
    # 1. discard blacklist_anyway_numbered_donor_ids (verify that the same donors with other nd ids are not left after this exclusion) 
    # 2. handle same day (and exp condition) donor samples with both illu and ult 10x libs. make sure that each (donor_id, date) after exclusion is only in ult or only in illu
    # 3. identify which donors still exist in both illu and ult, and assign according to wanted target_training_set_frac.
    
    illu_ult_tech_rep_exp_sets = get_illu_ult_tech_rep_exp_sets(illu_ult_numbered_donor_df)
    donor_df = illu_ult_numbered_donor_df.copy()
    
    ult_blacklist_nd_ids = []
    illu_blacklist_nd_ids = []

    if blacklist_anyway_numbered_donor_ids is not None:
        ult_blacklist_nd_ids.extend(
            donor_df.loc[donor_df['numbered_donor_id'].isin(blacklist_anyway_numbered_donor_ids) & donor_df['is_ultima'], 'numbered_donor_id'])
        illu_blacklist_nd_ids.extend(
            donor_df.loc[donor_df['numbered_donor_id'].isin(blacklist_anyway_numbered_donor_ids) & ~donor_df['is_ultima'], 'numbered_donor_id'])
        
        blacklist_anyway_donor_ids = set(donor_df.loc[donor_df['numbered_donor_id'].isin(blacklist_anyway_numbered_donor_ids), 'donor_id'])
        donor_df = donor_df[~donor_df['numbered_donor_id'].isin(blacklist_anyway_numbered_donor_ids)]
        donor_ids_partially_blacklisted = blacklist_anyway_donor_ids & set(donor_df['donor_id'])
        donor_ids_fully_blacklisted = blacklist_anyway_donor_ids - set(donor_df['donor_id'])
        print(f'NOTE: donor_ids_partially_blacklisted: {donor_ids_partially_blacklisted}')
        print(f'NOTE: donor_ids_fully_blacklisted: {donor_ids_fully_blacklisted}')
        # assert not donor_ids_partially_blacklisted, f'donor_ids_partially_blacklisted: {donor_ids_partially_blacklisted}'

    num_of_common_tech_rep_libs = len(illu_ult_tech_rep_exp_sets)
    num_of_common_tech_rep_libs_to_assign_to_ult = int(illu_ult_tech_rep_exp_sets_to_assign_to_ult_frac * num_of_common_tech_rep_libs)
    assert 0 <= num_of_common_tech_rep_libs_to_assign_to_ult <= num_of_common_tech_rep_libs
    random.seed(random_state - 1)
    illu_ult_tech_rep_exp_sets_to_assign_to_ult = random.sample(illu_ult_tech_rep_exp_sets, num_of_common_tech_rep_libs_to_assign_to_ult)
    illu_ult_tech_rep_exp_sets_to_assign_to_illu = [x for x in illu_ult_tech_rep_exp_sets if x not in illu_ult_tech_rep_exp_sets_to_assign_to_ult]
    assert len(illu_ult_tech_rep_exp_sets_to_assign_to_ult) + len(illu_ult_tech_rep_exp_sets_to_assign_to_illu) == num_of_common_tech_rep_libs
    
    illu_ult_tech_rep_exp_names_to_assign_to_ult = donor_df.loc[donor_df['is_ultima'] & donor_df['exp_name'].isin(
        sorted(set.union(set(), *illu_ult_tech_rep_exp_sets_to_assign_to_ult))), 'exp_name'].drop_duplicates()
    illu_ult_tech_rep_exp_names_to_not_assign_to_ult = donor_df.loc[donor_df['is_ultima'] & donor_df['exp_name'].isin(
        sorted(set.union(set(), *illu_ult_tech_rep_exp_sets_to_assign_to_illu))), 'exp_name'].drop_duplicates()
    illu_ult_tech_rep_exp_names_to_assign_to_illu = donor_df.loc[~donor_df['is_ultima'] & donor_df['exp_name'].isin(
        sorted(set.union(set(), *illu_ult_tech_rep_exp_sets_to_assign_to_illu))), 'exp_name'].drop_duplicates()
    illu_ult_tech_rep_exp_names_to_not_assign_to_illu = donor_df.loc[~donor_df['is_ultima'] & donor_df['exp_name'].isin(
        sorted(set.union(set(), *illu_ult_tech_rep_exp_sets_to_assign_to_ult))), 'exp_name'].drop_duplicates()
    assert not (set(illu_ult_tech_rep_exp_names_to_assign_to_ult) & set(illu_ult_tech_rep_exp_names_to_assign_to_illu))

    illu_ult_tech_rep_donor_ids_to_assign_to_ult = donor_df.loc[donor_df['exp_name'].isin(illu_ult_tech_rep_exp_names_to_assign_to_ult), 'donor_id'].drop_duplicates()
    illu_ult_tech_rep_donor_ids_to_assign_to_illu = donor_df.loc[donor_df['exp_name'].isin(illu_ult_tech_rep_exp_names_to_assign_to_illu), 'donor_id'].drop_duplicates()
    illu_ult_tech_rep_donor_ids_to_assign_to_ult_and_illu = set(illu_ult_tech_rep_donor_ids_to_assign_to_ult) & set(illu_ult_tech_rep_donor_ids_to_assign_to_illu)
    assert not illu_ult_tech_rep_donor_ids_to_assign_to_ult_and_illu, 'ugh. this might happen if a donor appears in more than one illu_ult_tech_rep_exp_set. ({illu_ult_tech_rep_donor_ids_to_assign_to_ult_and_illu})'

    illu_blacklist_nd_ids.extend(donor_df.loc[donor_df['donor_id'].isin(illu_ult_tech_rep_donor_ids_to_assign_to_ult) & ~donor_df['is_ultima'], 'numbered_donor_id'].to_list())
    ult_blacklist_nd_ids.extend(donor_df.loc[donor_df['donor_id'].isin(illu_ult_tech_rep_donor_ids_to_assign_to_illu) & donor_df['is_ultima'], 'numbered_donor_id'].to_list())
    
    illu_blacklist_nd_ids.extend(donor_df.loc[donor_df['exp_name'].isin(illu_ult_tech_rep_exp_names_to_not_assign_to_illu) & ~donor_df['is_ultima'], 'numbered_donor_id'].to_list())
    ult_blacklist_nd_ids.extend(donor_df.loc[donor_df['exp_name'].isin(illu_ult_tech_rep_exp_names_to_not_assign_to_ult) & donor_df['is_ultima'], 'numbered_donor_id'].to_list())

    donor_df = donor_df.loc[~donor_df['numbered_donor_id'].isin(illu_blacklist_nd_ids + ult_blacklist_nd_ids)]

    # print(f'ult_exp_names_with_illu_tech_rep: {ult_exp_names_with_illu_tech_rep}')
    # print(f'illu_exp_names_with_ult_tech_rep: {illu_exp_names_with_ult_tech_rep}')
    # raise

    donor_ids_sampled_in_illu_and_ult = sorted((donor_df[['donor_id', 'is_ultima']].drop_duplicates()['donor_id'].value_counts() >= 2).loc[lambda x: x].index)
    # test_platform = 'ult' if train_platform == 'illu' else 'illu'
    add_latest_diagnosis_class_col(donor_df)
    donor_id_to_latest_diagnosis_class = get_donor_id_to_latest_diagnosis_class(donor_df)

    platform_to_exclusive_donor_ids = {
        'ult': sorted(set(donor_df.loc[donor_df['is_ultima'], 'donor_id']) - set(donor_ids_sampled_in_illu_and_ult)),
        'illu': sorted(set(donor_df.loc[~donor_df['is_ultima'], 'donor_id']) - set(donor_ids_sampled_in_illu_and_ult)),
    }

    diagnosis_class_counts_df = donor_df.drop_duplicates(subset='donor_id')['latest_diagnosis_class'].value_counts().reset_index(name='donor_count')
    for name, curr_donor_ids in (
        ('common_illu_ult', donor_ids_sampled_in_illu_and_ult),
        ('exclusive_illu', platform_to_exclusive_donor_ids['illu']),
        ('exclusive_ult', platform_to_exclusive_donor_ids['ult']),
    ):
        diagnosis_class_counts_df = generic_utils.merge_preserving_df1_index_and_row_order(
            diagnosis_class_counts_df, 
            pd.Series([donor_id_to_latest_diagnosis_class[x] for x in curr_donor_ids], name='latest_diagnosis_class').value_counts().reset_index(name=f'{name}_donor_count'),
            how='left',
        )

    diagnosis_class_counts_df.fillna(0, inplace=True)
    print(f'diagnosis_class_counts_df (ignoring any that were excluded up to this point - due to blacklist anyway and due to illu_ult_tech_rep_exp_sets and donors that were in illu_ult_tech_rep_exp_sets and were assigned to the other (illu or ult)):\n{diagnosis_class_counts_df}')

    assert (
        diagnosis_class_counts_df['donor_count'] == 
        diagnosis_class_counts_df[['common_illu_ult_donor_count','exclusive_illu_donor_count','exclusive_ult_donor_count']].sum(axis=1)
    ).all()

    diagnosis_class_counts_df[f'target_train_donor_count'] = (diagnosis_class_counts_df['donor_count'] * target_training_set_frac).astype(int)
    diagnosis_class_counts_df[f'wanted_donor_count_to_add_to_train'] = (
        diagnosis_class_counts_df['target_train_donor_count'] - diagnosis_class_counts_df[f'exclusive_{train_platform}_donor_count'])
    diagnosis_class_counts_df[f'common_donor_count_to_assign_to_train'] = diagnosis_class_counts_df[
        ['wanted_donor_count_to_add_to_train', 'common_illu_ult_donor_count']].min(axis=1).astype(int)
    # return diagnosis_class_counts_df

    common_donor_ids_to_add_to_train = []
    common_donor_ids_to_add_to_test = []
    for i, row in diagnosis_class_counts_df.iterrows():
        shuffled_common_donor_ids = list(donor_df.loc[
            (donor_df['latest_diagnosis_class'] == row['latest_diagnosis_class']) & 
            donor_df['donor_id'].isin(donor_ids_sampled_in_illu_and_ult), 'donor_id'
        ].drop_duplicates().sample(frac=1, random_state=(random_state + i)))
        if row['common_donor_count_to_assign_to_train'] > 0:
            common_donor_ids_to_add_to_train.extend(shuffled_common_donor_ids[:row['common_donor_count_to_assign_to_train']])
            common_donor_ids_to_add_to_test.extend(shuffled_common_donor_ids[row['common_donor_count_to_assign_to_train']:])
        else:
            common_donor_ids_to_add_to_test.extend(shuffled_common_donor_ids)

    common_nd_ids_to_add_to_test = donor_df.loc[donor_df['donor_id'].isin(common_donor_ids_to_add_to_test), 'numbered_donor_id']
    common_nd_ids_to_add_to_train = donor_df.loc[donor_df['donor_id'].isin(common_donor_ids_to_add_to_train), 'numbered_donor_id']

    if train_platform == 'illu':
        illu_blacklist_nd_ids.extend(donor_df.loc[donor_df['numbered_donor_id'].isin(common_nd_ids_to_add_to_test) & ~donor_df['is_ultima'], 'numbered_donor_id'])
        ult_blacklist_nd_ids.extend(donor_df.loc[donor_df['numbered_donor_id'].isin(common_nd_ids_to_add_to_train) & donor_df['is_ultima'], 'numbered_donor_id'])
    else:
        illu_blacklist_nd_ids.extend(donor_df.loc[donor_df['numbered_donor_id'].isin(common_nd_ids_to_add_to_train) & ~donor_df['is_ultima'], 'numbered_donor_id'])
        ult_blacklist_nd_ids.extend(donor_df.loc[donor_df['numbered_donor_id'].isin(common_nd_ids_to_add_to_test) & donor_df['is_ultima'], 'numbered_donor_id'])
    
    illu_blacklist_nd_ids = sorted(set(illu_blacklist_nd_ids))
    ult_blacklist_nd_ids = sorted(set(ult_blacklist_nd_ids))

    unexpected_common_blacklist_nd_ids = (
        (set(illu_blacklist_nd_ids) & set(ult_blacklist_nd_ids)) - 
        (set(blacklist_anyway_numbered_donor_ids) if (blacklist_anyway_numbered_donor_ids is not None) else set())
    )
    assert not unexpected_common_blacklist_nd_ids, str(f'unexpected_common_blacklist_nd_ids: {unexpected_common_blacklist_nd_ids}')

    orig_donor_df = illu_ult_numbered_donor_df.copy()
    add_latest_diagnosis_class_col(orig_donor_df)
    final_illu_mask = ~orig_donor_df['is_ultima'] & ~orig_donor_df['numbered_donor_id'].isin(illu_blacklist_nd_ids)
    final_ult_mask = orig_donor_df['is_ultima'] & ~orig_donor_df['numbered_donor_id'].isin(ult_blacklist_nd_ids)

    final_illu_donor_ids = set(orig_donor_df.loc[final_illu_mask, 'donor_id'])
    final_ult_donor_ids = set(orig_donor_df.loc[final_ult_mask, 'donor_id'])
    final_common_donor_ids = final_illu_donor_ids & final_ult_donor_ids
    assert not final_common_donor_ids, set(final_common_donor_ids)

    illu_donor_diagnosis_class_counts = orig_donor_df.loc[final_illu_mask].drop_duplicates(subset='donor_id')['latest_diagnosis_class'].value_counts()
    ult_donor_diagnosis_class_counts = orig_donor_df.loc[final_ult_mask].drop_duplicates(
        subset='donor_id')['latest_diagnosis_class'].value_counts()
    print(f'illu_donor_diagnosis_class_counts:\n{illu_donor_diagnosis_class_counts}')
    print(f'ult_donor_diagnosis_class_counts:\n{ult_donor_diagnosis_class_counts}')

    assert not illu_ult_numbered_donor_df.loc[illu_ult_numbered_donor_df['numbered_donor_id'].isin(illu_blacklist_nd_ids), 'is_ultima'].any()
    assert illu_ult_numbered_donor_df.loc[illu_ult_numbered_donor_df['numbered_donor_id'].isin(ult_blacklist_nd_ids), 'is_ultima'].all()

    assert set(ult_blacklist_nd_ids) <= set(illu_ult_numbered_donor_df['numbered_donor_id'])
    assert set(illu_blacklist_nd_ids) <= set(illu_ult_numbered_donor_df['numbered_donor_id'])

    return illu_blacklist_nd_ids, ult_blacklist_nd_ids

def add_approx_birth_date_by_bleeding_date_and_age(df):
    bleeding_date_as_date = pd.to_datetime(df['bleeding_date'], format='%d.%m.%y', errors='raise')
    df['approx_birth_date'] = bleeding_date_as_date - pd.to_timedelta(df['donor_age'].astype(float) * 365, unit='days', errors='raise')
def add_age_by_approx_birth_date_and_bleeding_date(df):
    bleeding_date_as_date = pd.to_datetime(df['bleeding_date'], format='%d.%m.%y', errors='raise')
    df['donor_age'] = ((bleeding_date_as_date - df['approx_birth_date']).dt.days / 365).astype(int)



print('sc_rna_seq_preprocessing_params was loaded/reloaded')
