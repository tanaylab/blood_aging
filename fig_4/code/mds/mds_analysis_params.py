import numpy as np
import pandas as pd
import os.path
import pathlib
import re
from generic import generic_utils
# from generic import gtf_interface_and_utils
from generic import hg38_utils
from generic import mc_utils
from mds import pb_cd34_c_score_threshs
from mds import pb_cd34_mc_score_threshs
from mds import k1k_pb_mc_score_threshs
from mds import lateral_and_noisy_genes
from mds import mds_in_out_dir_paths
import collections
import pickle
import itertools

# https://stackoverflow.com/questions/16424493/pandas-setting-no-of-max-rows/16433953#16433953
# pd.describe_option('display')
pd.set_option('display.large_repr', 'truncate')
pd.set_option('display.max_rows', 40)
pd.set_option('display.min_rows', 40)


DEFAULT_CELL_STATE_COLORS_CSV_FILE_PATH = os.path.join(mds_in_out_dir_paths.MISC_INPUT_DIR_PATH, 'cell_type_colors.csv')

PLATELET_POOLED_DOWNSAMPLED_UMIS_DF_CSV_FILE_PATH = '/dummy/dummy/dummy/raid/mds/all_data/intermediate_output/230227_platelet_pooled_downsampled_umis_df.csv'


SPECIFIC_OUTLIER_EXP_DONORS_NAME_PAIRS = [
    ('demux_07_11_21_1', 'N186'),
    ('demux_09_01_22_1', 'N205'),
    ('demux_28_11_21_1', 'N191'),
    
    # N257 dominates some CLP-E metacells (and is AML)
    # ('N257_bm_06_10_22_1', 'N257'),
    # ('N257_bm_06_10_22_1_illumina', 'N257'),
    ('demux_06_06_22_2_ultima', 'N257'),
    ('demux_06_06_22_2', 'N257'),
    ('demux_06_06_22_1_ultima', 'N257'),
    ('demux_06_06_22_1', 'N257'),

    ('G8_13_07_23_1', 'G8'), # dominating some MEBEMP-L metacells
    ('demux_22_02_21_1', 'N154'), # dominating some B metacells, and healthy, so we can exclude her, i guess

    # ('', 'N235'), # dominating some GMP-L metacells (and also MGST1 in GMP-E). could just reduce the number of cells only for the 28_02_22 exps (say, each to 1.5k). i guess this would be good enough. but first MCNoise.
    # ('demux_21_02_21_1', 'N151'), # dominating some BEMPs # but maybe because this sample is illumina and old sample? check again after MCNoise. also, could just reduce the number of cells (from ~6.7k to 3k or even 2k). i guess this would be good enough. but first MCNoise.


    # i guess don't exclude
    # ('', 'N376'), # somewhat dominating some DCs and Monocytes # so i guess leave in the model...
    # ('demux_g_03_07_23_1', 'g_03_07_23_1_a'), # dominating some Monocytes. though maybe we don't really care about Monocytes so ok to leave in the model??
    # ('demux_g_03_07_23_1', 'g_03_07_23_1_b'), # dominating some Monocytes. though maybe we don't really care about Monocytes so ok to leave in the model?? maybe the same donor as g_03_07_23_1_a????
]

SPECIFIC_EXPS_TO_EXCLUDE_NAMES = [
    # previously had here more stuff. they moved to SPECIFIC_EXPS_TO_EXCLUDE_NAMES_EXCEPT_FOR_NIMROD in sc_rna_seq_preprocessing_params.py
    
]

SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS = [
    ('demux_10_01_22_1', 'N208'), # high interferon
    ('demux_g_12_09_23_1', 'N320'), 
]



# NOTE: https://www.ensembl.info/2021/03/15/retirement-of-clone-based-gene-names/: "Human, mouse, rat and zebrafish genes without an official gene symbol, have been conventionally named after BAC clones in Ensembl. We plan to remove these names in Ensembl release 104 and resort to using the Ensembl gene stable IDs instead, in line with practices adopted for all other vertebrate species". at some point I considered to mark all of these as lateral, but decided against it.
AC_AND_AL_CLONE_BASED_GENE_NAMES_FILE_PATH = '/dummy/dummy/dummy/tanay_group/mds/AC_and_AL_clone_based_gene_names.txt'

DONORS_WITH_HIGH_FACS_BLAST_IDS_AND_ENOUGH_CLP_E_CELLS = [
    'N257',
    'N205',
    'N191',
    'N280',
]
DONORS_WITH_HIGH_FACS_BLAST_IDS = [
    *DONORS_WITH_HIGH_FACS_BLAST_IDS_AND_ENOUGH_CLP_E_CELLS,
    'N265',
]

# NON_STEM_CELL_STATE_NAMES = [
#     # 'Monocytes-FCGR3A-high',
#     'cMonocyte',
#     'B',
#     'B-high-SUB1',
#     'NKT',
#     'DC',
#     # 'DC-IRF8-low?',
#     'Endothel',
#     'HBB-MPO-high',
#     'post-GMP',
#     'post-EP',
#     'post-NKTDP-L',
#     'post-BEMP',
#     'Dissimilar-CD34-low',
#     'Dissimilar-CD34-low-JUN-high',


# ]

MYELOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB = [
    'BEMP',
    'EP',
    'ERYP',
    'high_MPP_sig_ERYP', # ?
    'MEBEMP-L',
    'MEBEMP-M', # TODO: remove
    'MEBEMP-E',
    'GMP-L',
    'GMP-E',
]

NON_LYMPHOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB = [
    *MYELOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB,
    'MPP',
    'HSC',
    'HSC_MPP',
]
PB_NON_LYMPHOID_HSPC_STATE_NAMES = [
    *NON_LYMPHOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB,
    'MKP',
    'high_MPP_sig_MKP',
    'CFD_tryptase_GMP-L', # only in 1 patient? (G12)
]
LYMPHOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB = [
    'CLP',
    'CLP-E', # TODO: remove
    'CLP-M', # TODO: remove
    'CLP-L', # TODO: remove
    'CLP_NKTDP_intermediate?', # TODO: remove
    'NKTDP',
    'pro-B?',
    'pre-B?',
]
PB_LYMPHOID_HSPC_STATE_NAMES = [
    *LYMPHOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB,
]

HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB = [
    *NON_LYMPHOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB,
    *LYMPHOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB,
]
PB_HSPC_STATE_NAMES = [
    *PB_NON_LYMPHOID_HSPC_STATE_NAMES,
    *PB_LYMPHOID_HSPC_STATE_NAMES,
    
    # 'DC-prog?',
    # 'cDC2-prog?', # NOTE: non-zero CD34, but seems batchy, so better not have it here, i think
    # 'high_cMonocyte_sig_cDC2-prog?', # NOTE: non-zero CD34, but seems batchy, so better not have it here, i think

    # # not sure that's ok that i commented them out...
    # 'N186-MDS',
    # 'N205-MDS',
    # 'N191-AML',
    # 'N257-AML?',
]

MATURE_B_CELL_STATE_NAMES = [
    'B',
    'aged-B?',
    'naive-B',
    'memory-B',
    'CCDC88A-high-B',
]
B_CELL_STATE_NAMES = [
    'pro-B?',
    'pre-B?',
    *MATURE_B_CELL_STATE_NAMES,
]

ORDERED_CELL_STATE_NAMES = [
    'Unknown',
    'post-BEMP',
    'post-EP',
    'post-GMP',
    'MKP', 
    'high_MPP_sig_MKP', 
    'BEMP',
    'HBB-MPO-high',
    'EP',
    'high_MKP_sig_ERYP', 
    'ERYP',
    'high_MPP_sig_ERYP',
    'MEBEMP-L',
    'MEBEMP-M',
    'MEBEMP-E',
    'CFD_tryptase_GMP-L',
    'GMP-L',
    'GMP-E',
    'MPP',
    'HSC',
    'HSC_MPP',
    'CLP-E',
    'CLP-M',
    'CLP',
    'CLP-L',
    'CLP-IL7R',
    'CLP_NKTDP_intermediate?',
    'NKTDP-E',
    'NKTDP-L',
    'NKTDP',
    'post-NKTDP-L',
    'pro-B?',
    'pre-B?',
    'B',
    'aged-B?',
    'naive-B',
    'memory-B',
    'plasmablast_igha_iglc',
    'plasmablast_ighm_ighg',
    'plasmablast',
    'CCDC88A-high-B',
    'B-high-SUB1',
    'NKT-prog?',
    'NKT',
    'NK',
    'T',
    # 'DC-IRF8-low?',
    'DC',
    'pDC',
    'cDC', # TODO: remove
    'cDC1',
    'cDC?', # TODO: remove
    'AS-DC',
    'mregDC', 
    'DC-prog?',
    'cDC2-prog?',
    'high_cMonocyte_sig_cDC2-prog?', 
    # 'Monocytes-FCGR3A-high',
    'Monocytes', # TODO: remove
    'Monocyte', # TODO: remove?
    'cMonocyte',
    'intermMonocyte',
    'ncMonocyte',
    'ncMonocytes', # TODO: remove
    'Endothel',
    'LAMP3_outlier',
    'Dissimilar-high-IFN',
    'Dissimilar-CD34-low',
    'Dissimilar-CD34-low-JUN-high',
    # 'N224-high-ribo',
    'N186-MDS',
    'N205-MDS',
    'N191-AML',
    'N257-AML?',
    'Dissimilar',
    'Outliers',
    'Mixture',
    'Doublet',
    'monocyte_B_doublet',
    'mono-neutrophil-Doublet?',
    'unknown-neutrophil-Doublet?',
    'state_unassigned',
    'unassigned',
    'nan',
]
assert set(PB_HSPC_STATE_NAMES) <= set(ORDERED_CELL_STATE_NAMES)

ORDERED_PB_HSPC_STATE_NAMES = [x for x in ORDERED_CELL_STATE_NAMES if x in PB_HSPC_STATE_NAMES]
ORDERED_PB_NON_HSPC_STATE_NAMES = [x for x in ORDERED_CELL_STATE_NAMES if x not in PB_HSPC_STATE_NAMES]

NEIGHBOR_STATE_PAIRS = [
    ('CLP', 'NKTDP'),
    ('CLP', 'pro-B?'),
    ('MEBEMP-L', 'BEMP'),
    ('MEBEMP-L', 'ERYP'),
    ('MEBEMP-L', 'MKP'),
    ('HSC_MPP', 'MEBEMP-L'),
    ('HSC_MPP', 'GMP-L'),
]
NEIGHBOR_STATE_LOG_RATIO_COLS = [f'log_ratio_c_{state2}_{state1}' for state1, state2 in NEIGHBOR_STATE_PAIRS]

MREG_DC_RULES = [
    # https://www.tandfonline.com/doi/full/10.1080/2162402X.2023.2294564: "population known as mregDCs (also called LAMP3+DCs).Citation32,Citation33 Remarkably, this population shows the highest expression of maturation markers, such as, LAMP3, CCR7, CD83, BIRC3 and MARCKSL1"
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9937888/ - specifically table 2.
    [('metadata', 'B', [False])],
    [('metadata', 'monocyte', [False])],
    [('metadata', 'lamp3_outlier_specific', pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH, np.inf)],
]

GENES_THAT_WERE_AT_LEAST_MINUS_3_IN_ANY_CELL = [
    'S100A9', 'S100A8', 
    'PPBP', 'TMSB4X', 
    'MALAT1',
    'IGHG4', 'IGHA1', 'IGHM', 'IGKC', 
    'HBA2', 'HBA1', 'HBB',
    'CCL4', 
    'FTL', 
    'IGLC2', 'IGLC3', 'IGLC7', 
    'DEFA3', # 240219 was extreme in <10 out of 1.3M cells (ill vs ult c_ad). 240219: in final_all_assigned_pb_except_problematics, each cell with extreme DEFA3 was an outlier (only 1, or three if lowering the threshold further)
    'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP6', 'MT-CO3', 'MT-ND3', 'MT-ND4', 'MT-CYB',
]
PRESUMABLY_BEMP_GENE_NAMES = [
    'LMO4', 'HPGD', 'CNKSR3', 'MS4A3', 'MS4A2', 'VWA5A', 'SLC18A2',
    'CLEC12A', 'ALOX5AP', 'GPR183', 'HDC', 'TPSB2', 'TPSAB1', 'ERN1',
    'FOXJ1', 'CLC', 'GRAP2', 'FAM83F', 'BACE2',
]

MAYBE_NOISY_6P_HLA_GENE_NAMES = lateral_and_noisy_genes.MAYBE_NOISY_6P_HLA_GENE_NAMES


NOISY_6P_HIST_GENE_NAMES = lateral_and_noisy_genes.NOISY_6P_HIST_GENE_NAMES


CORRELATED_WITH_VIM_AROUND_CLP_E_GENE_NAMES = [
    'VIM', 'ANXA1', 'TSPAN2', 'AHNAK', 'MYADM', 'LMNA', 
    ### 'S100A10',
       'DPYSL3', 'MTURN', 'ATP2B1', 'TUBA1A', 
    ###    'FHL1', 
       'PBXIP1', 'TAGLN2',
       'IDS', 'UCP2', 
    ###    'ANKRD28', 
    ###    'SEPT11', 
       'CEP126', 
    ###    'WDR49',
]


# FILTERED_RP_MOST_ANTI_CORR_WITH_DNTT_IN_HSC_TO_CLP_M = ['RPL12', 'RPL5', 'RPS6', 'RPL10A', 'RPS18', 'RPS3A', 'RPS8', 'RPL13', 'RPS4X', 'RPL37A', 'RPL3', 'RPL32', 'RPS23', 'RPS2', 'RPS14', 'RPS12', 'RPL7A', 'RPL15', 'RPS24', 'RPL23', 'RPL29', 'RPL36', 'RPS3', 'RPL18A', 'RPS25', 'RPL35', 'RPS5', 'RPLP0', 'RPL22', 'RPS19', 'RPL8', 'RPL36A', 'RPS21', 'RPL38', 'RPL31', 'RPL4', 'RPLP1']
FILTERED_RP_ANTI_CORR_WITH_DNTT_IN_HSC_TO_CLP_M = [
    'RPL9', 'RPL12', 'RPL34', 'RPL5', 'RPS6', 'RPL10A', 'RPS18', 'RPS3A', 'RPS8', 'RPL7', 'RPL13', 'RPL39', 'RPL37A', 'RPL3', 'RPL32', 'RPS23', 'RPL21', 'RPS2', 'RPL6', 'RPL26', 'RPS14', 'RPS12', 'RPL7A', 'RPL15', 'RPS24', 'RPL23', 'RPL29', 'RPS7', 'RPL36', 'RPL35A', 'RPL11', 'RPS3', 'RPL18A', 'RPS15A', 'RPS25', 'RPL13A', 'RPS29', 'RPL35', 'RPS5', 'RPL10', 'RPLP0', 'RPL22', 'RPS19', 'RPL37', 'RPL41', 'RPSA', 'RPL8', 'RPL36A', 'RPS9', 'RPS21', 'RPL38', 'RPL27A', 'RPL30', 'RPL28', 'RPL31', 'RPL24', 'RPL4', 'RPS15', 'RPS20', 'RPL18', 'RPL14', 'RPS28', 'RPS13', 'RPS27A', 'RPS16', 'RPLP1', 'RPL19', 'RPL27', 'RPS27', 'RPS26', 'RPS11', 'RPLP2', 'RPS10',
]

CORR_WITH_PCNA_IN_MEBEMP_L_GENES__1 = [
    # chosen by donor gene heatmap showing pooled expression of all CORR_WITH_PCNA_IN_MEBEMP_L_GENES across early (but not earliest) cells in MEBEMP_ERYP_traj.
    'HSP90B1', 'CENPV', 'CDCA7', 'MANF', 'CENPF', 'MKI67', 'TOP2A', 'MCM2', 'SVIP', 'HELLS', 'CBX5', 'MCM6', 'DHFR', 'ATAD5', 'BRCA1', 'XRCC2', 'WDR76', 'BOLA3', 'MAD2L1', 'DTL', 'CENPH', 'RPA3', 'NDUFA6', 'FEN1', 'HIRIP3', 'CLSPN', 'ORC6', 'RRM1', 'CHEK1', 'CENPK', 'MCM3', 'TFDP1', 'CKS1B', 'POLE4', 'NPW', 'CTNNAL1', 'DTYMK', 'SLC1A5', 'POP7', 'CENPM', 'UHRF1', 'DSCC1', 'CDC6', 'E2F1', 'UBE2C', 'MCM10', 'RAD51AP1', 'MND1', 'CDKN3', 'TK1', 'BIRC5', 'GGH', 'CCNB2', 'PRIM1', 'PTMA', 'ESCO2', 'CDK1', 'KIF15', 'NCAPG', 'MYBL2', 'TIMELESS', 'FAM111B', 'ASF1B', 'CDC45', 'CDCA5', 'CCNE1', 'NUSAP1', 'CCDC167', 'LIG1', 'SLC29A1', 'ATAD2', 'FANCI', 'ASPM', 'BRCA2',
]
CORR_WITH_PCNA_IN_MEBEMP_L_GENES111 = [
    'CDCA7', 'HELLS', 'CBX5', 'NASP', 'PRKDC', 
    ###'NCL', # pretty strong
    'MSH6', 'TMPO', 'SUPT16H', 'DNMT1', 'SVIP', 'MCM6', 'MCM2', 'TFDP1', 'SLC29A1', 'SLC1A5', 'ATAD2', 'BRCA1', 'SMC2', 'BRCA2', 'XRCC2', 'ATAD5', 
    ###'PTMA', # very strong
    'DTL', 'WDR76', 'FANCI', 'LIG1', 'PRIM1', 'DTYMK', 'RFC4', 'CCDC167', 'LMNB1', 'CENPF', 'MCM3', 'CENPK', 'NUSAP1', 'GMNN', 'CLSPN', 'ORC6', 'TOP2A', 'ASPM', 'MKI67', 'CENPH', 'DHFR', 'RRM1', 'HIRIP3',

    'TK1', 'CDKN3', 'BIRC5', 'CCNB2', 'UBE2C', 'DSCC1', 'MCM10', 'MND1', 'CENPM', 'PXMP2', 'UHRF1', 'TIMELESS', 'NCAPG', 'KIF15', 'MYBL2', 'CDCA5', 'CCNE1', 'ESCO2', 'ASF1B', 'FAM111B', 'E2F1', 'CTNNAL1', 'CDK1', 'CDC45', 'RAD51AP1', 'CDC6',
]
CORR_WITH_PCNA_IN_MEBEMP_L_GENES11 = [
    'CDCA7',
 'HELLS',
 'CBX5',
 'NASP',
 'PRKDC',
 'NCL',
 'MSH6',
 'TMPO',
 'SUPT16H',
 'DNMT1',
 'SVIP',
 'MCM6',
 'MCM2',
 'TFDP1',
 'SLC29A1',
 'SLC1A5',
 'ATAD2',
 'BRCA1',
 'SMC2',
 'BRCA2',
 'XRCC2',
 'ATAD5',
 'PTMA',
 'DTL',
 'WDR76',
 'FANCI',
 'LIG1',
 'PRIM1',
 'DTYMK',
 'RFC4',
 'CCDC167',
 'LMNB1',
 'CENPF',
 'MCM3',
 'CENPK',
 'NUSAP1',
 'GMNN',
 'CLSPN',
 'ORC6',
 'TOP2A',
 'ASPM',
 'MKI67',
 'CENPH',
 'DHFR',
 'RRM1',
 'HIRIP3',
 'TK1',
 'CDKN3',
 'BIRC5',
 'CCNB2',
 'UBE2C',
 'DSCC1',
 'MCM10',
 'MND1',
 'CENPM',
 'PXMP2',
 'UHRF1',
 'TIMELESS',
 'NCAPG',
 'KIF15',
 'MYBL2',
 'CDCA5',
 'CCNE1',
 'ESCO2',
 'ASF1B',
 'FAM111B',
 'E2F1',
 'CTNNAL1',
 'CDK1',
 'CDC45',
 'RAD51AP1',
 'CDC6',
 'DEK',
 'HNRNPAB',
 'NDUFA6',
 'POLE3',
 'PA2G4',
 'CENPX',
 'DCTPP1',
 'PGP',
 'TOMM40',
 'POLR3K',
 'SLBP',
 'POP7',
 'SNRNP25',
 'PSMG1',
 'CENPU',
 'MCM7',
 'MCM5',
 'PAICS',
 'CDT1',
 'LY6E',
 'MT2A',
 'NPW',
 'MANF',
 'PPIF',
 'SDF2L1',
 'MRTO4',
 'TMEM106C',
 'HSP90B1',
 'EBNA1BP2',
 'CENPV',
 'BOLA3',
 'RMI2',
 'POLE4',
 'MAGOHB',
 'UBE2T',
 'GGH',
 'ZWINT',
 'MAD2L1',
 'CENPW',
 'RPA3',
 'TMEM97',
 'TYMS',
 'MCM4',
 'CHEK1',
 'PCNA',
 'CKS1B',
 'UNG',
 'FEN1',
 'GINS2',
]
CORR_WITH_PCNA_IN_MEBEMP_L_GENES12 = [
    'RNASEH2C',
 'METTL26',
 'EIF4EBP1',
 'PHPT1',
 'NHP2',
 'UQCC2',
 'NDUFS6',
 'POLR2L',
 'NDUFS5',
 'COX5B',
 'ATP5MF',
 'NDUFB2',
 'HSP90AA1',
 'TMA7',
 'COX6C',
 'PGK1',
 'PRELID1',
 'PDIA3',
 'PSME2',
 'PSMA7',
 'ENO1',
 'HSPA5',
 'TUBB4B',
 'EIF5A',
 'ODC1',
 'CYC1',
 'CDK4',
 'PRMT1',
 'UQCRQ',
 'CYCS',
 'ATP1B3',
 'PGAM1',
 'SNRPB',
 'SNU13',
 'HSPA8',
 'LDHA',
 'SRM',
 'PHB',
 'NME4',
 'MRPL12',
 'TPI1',
 'CCT5',
 'YBX1',
 'MRPS34',
 'DNPH1',
 'PRDX4',
 'NUDT1',
 'VDAC1',
 'NME1',
 'SSBP1',
 'PPIA',
 'MRPL3',
 'GCSH',
 'NUDC',
 'PPM1G',
 'PARK7',
 'HSPD1',
]
CORR_WITH_PCNA_IN_MEBEMP_L_GENES1 = [
    'CDCA7',
 'HELLS',
 'CBX5',
 'NASP',
 'PRKDC',
 'NCL',
 'MSH6',
 'TMPO',
 'SUPT16H',
 'DNMT1',
 'SVIP',
 'MCM6',
 'MCM2',
 'TFDP1',
 'SLC29A1',
 'SLC1A5',
 'ATAD2',
 'BRCA1',
 'SMC2',
 'BRCA2',
 'XRCC2',
 'ATAD5',
 'PTMA',
 'DTL',
 'WDR76',
 'FANCI',
 'LIG1',
 'PRIM1',
 'DTYMK',
 'RFC4',
 'CCDC167',
 'LMNB1',
 'CENPF',
 'MCM3',
 'CENPK',
 'NUSAP1',
 'GMNN',
 'CLSPN',
 'ORC6',
 'TOP2A',
 'ASPM',
 'MKI67',
 'CENPH',
 'DHFR',
 'RRM1',
 'HIRIP3',
 'TK1',
 'CDKN3',
 'BIRC5',
 'CCNB2',
 'UBE2C',
 'DSCC1',
 'MCM10',
 'MND1',
 'CENPM',
 'PXMP2',
 'UHRF1',
 'TIMELESS',
 'NCAPG',
 'KIF15',
 'MYBL2',
 'CDCA5',
 'CCNE1',
 'ESCO2',
 'ASF1B',
 'FAM111B',
 'E2F1',
 'CTNNAL1',
 'CDK1',
 'CDC45',
 'RAD51AP1',
 'CDC6',
 'DEK',
 'HNRNPAB',
 'NDUFA6',
 'POLE3',
 'PA2G4',
 'CENPX',
 'DCTPP1',
 'PGP',
 'TOMM40',
 'POLR3K',
 'SLBP',
 'POP7',
 'SNRNP25',
 'PSMG1',
 'CENPU',
 'MCM7',
 'MCM5',
 'PAICS',
 'CDT1',
 'LY6E',
 'MT2A',
 'NPW',
 'MANF',
 'PPIF',
 'SDF2L1',
 'MRTO4',
 'TMEM106C',
 'HSP90B1',
 'EBNA1BP2',
 'CENPV',
 'BOLA3',
 'RMI2',
 'POLE4',
 'MAGOHB',
 'UBE2T',
 'GGH',
 'ZWINT',
 'MAD2L1',
 'CENPW',
 'RPA3',
 'TMEM97',
 'TYMS',
 'MCM4',
 'CHEK1',
 'PCNA',
 'CKS1B',
 'UNG',
 'FEN1',
 'GINS2',
 'RNASEH2C',
 'METTL26',
 'EIF4EBP1',
 'PHPT1',
 'NHP2',
 'UQCC2',
 'NDUFS6',
 'POLR2L',
 'NDUFS5',
 'COX5B',
 'ATP5MF',
 'NDUFB2',
 'HSP90AA1',
 'TMA7',
 'COX6C',
 'PGK1',
 'PRELID1',
 'PDIA3',
 'PSME2',
 'PSMA7',
 'ENO1',
 'HSPA5',
 'TUBB4B',
 'EIF5A',
 'ODC1',
 'CYC1',
 'CDK4',
 'PRMT1',
 'UQCRQ',
 'CYCS',
 'ATP1B3',
 'PGAM1',
 'SNRPB',
 'SNU13',
 'HSPA8',
 'LDHA',
 'SRM',
 'PHB',
 'NME4',
 'MRPL12',
 'TPI1',
 'CCT5',
 'YBX1',
 'MRPS34',
 'DNPH1',
 'PRDX4',
 'NUDT1',
 'VDAC1',
 'NME1',
 'SSBP1',
 'PPIA',
 'MRPL3',
 'GCSH',
 'NUDC',
 'PPM1G',
 'PARK7',
 'HSPD1',
]
CORR_WITH_PCNA_IN_MEBEMP_L_GENES2 = [
    'SEC61G',
 'RBX1',
 'LSM5',
 'SNRPD3',
 'PFDN2',
 'SERF2',
 'SNRPG',
 'NOP10',
 'BANF1',
 'SEC61B',
 'ERH',
 'SPCS2',
 'PTGES3',
 'CHCHD2',
 'CFL1',
 'HNRNPF',
 'SRSF9',
 'FABP5',
 'H2AFZ',
 'DUT',
 'ATP5MC3',
 'RAN',
 'PEBP1',
 'C1QBP',
 'PRDX1',
 'ACTB',
 'PFN1',
 'PPP1R14B',
 'PSMB5',
 'HPRT1',
 'GLRX5',
 'PSMB3',
 'PSMB2',
 'STOML2',
 'ALYREF',
 'PPIH',
 'MRPL4',
 'MRPS26',
 'MRPS6',
 'PCLAF',
 'HNRNPR',
 'SRSF7',
 'HNRNPD',
 'SRSF2',
 'NDUFAF8',
 'SNRPE',
 'SNRPF',
 'HMGN2',
 'HMGB2',
 'TUBA1B',
 'MINOS1',
 'SIVA1',
 'SNRPD1',
 'HSPE1',
 'ATP5MC1',
 'RANBP1',
 'PDIA6',
 'GTF3A',
 'PSMA3',
 'LSM4',
 'EEF1E1',
 'MRPL11',
 'POMP',
 'PDCD5',
 'CACYBP',
 'SNRPA1',
 'TXNDC17',
 'LSM3',
 'SMS',
 'COA4',
 'MRPL20',
 'UBE2N',
 'ANAPC11',
 'FKBP3',
 'SNRPC',
 'POLR2F',
 'CALM1',
 'STMN1',
 'GTF3C6',
 'NDUFAB1',
]
CORR_WITH_PCNA_IN_MEBEMP_L_GENES = [
    # 240408: corr > 0.75 with PCNA across MEBEMP-L MCs.
    'PCNA', 'TYMS', 'RANBP1', 'H2AFZ', 'CLSPN', 'GINS2', 'TUBA1B',
       'MCM4', 'RAN', 'SNRPD1', 'SIVA1', 'NME1', 'DUT', 'SNRPB', 'MCM3',
       'FEN1', 'HELLS', 'CHEK1', 'HSPE1', 'FABP5', 'HMGN2', 'ORC6',
       'RAD51AP1', 'LSM4', 'CENPM', 'MCM10', 'ATP5MC1', 'ZWINT',
       'PPP1R14B', 'ACTB', 'UNG', 'PTMA', 'SNRPG', 'GMNN', 'TK1', 'RPA3',
       'NDUFAF8', 'DNMT1', 'CDC6', 'CDT1', 'CDK4', 'ALYREF', 'BRCA1',
       'PCLAF', 'PSMA7', 'SRM', 'CDC45', 'DCTPP1', 'MRPS26', 'CKS1B',
       'NPW', 'CHCHD2', 'DTYMK', 'HIRIP3', 'PRELID1', 'PPIH', 'EIF4EBP1',
       'CENPK', 'MCM6', 'ODC1', 'VDAC1', 'FAM111B', 'PRKDC', 'ERH',
       'HSPD1', 'CACYBP', 'C1QBP', 'MAD2L1', 'ASPM', 'NCL', 'YBX1',
       'NASP', 'MCM7', 'SLBP', 'GTF3A', 'TOP2A', 'NHP2', 'PSME2', 'UBE2T',
       'TMEM106C', 'NUDT1', 'POLR2L', 'POMP', 'CENPW', 'LSM5', 'UQCRQ',
       'E2F1', 'GGH', 'PDIA6', 'ATP5MC3', 'CFL1', 'DHFR', 'RMI2',
       'SNRNP25', 'UHRF1', 'ASF1B', 'PSMG1', 'LDHA', 'MCM2', 'SRSF2',
       'LMNB1', 'SRSF7', 'NDUFAB1', 'EBNA1BP2', 'PRDX4', 'KIF15', 'CENPH',
       'NDUFB2', 'BOLA3', 'BANF1', 'PHPT1', 'MCM5', 'EIF5A', 'PAICS',
       'HNRNPD', 'CDCA5', 'PEBP1', 'MRPS34', 'RRM1', 'DEK', 'LY6E',
       'PRMT1', 'LIG1', 'CCT5', 'METTL26', 'ATAD2', 'SNRPF', 'CDK1',
       'XRCC2', 'CDCA7', 'CENPF', 'PSMB3', 'SNRPE', 'CYCS', 'HNRNPF',
       'PTGES3', 'NUDC', 'TOMM40', 'MRPL11', 'HSP90AA1', 'PDCD5', 'CENPU',
       'PA2G4', 'UBE2C', 'MRPL4', 'SNRPD3', 'SDF2L1', 'HSPA8', 'PRDX1',
       'MINOS1', 'ATP5MF', 'PGP', 'MRPL12', 'MKI67', 'SVIP', 'MYBL2',
       'TFDP1', 'HMGB2', 'WDR76', 'PPIF', 'CALM1', 'MRPL20', 'EEF1E1',
       'NDUFS5', 'FKBP3', 'SERF2', 'PXMP2', 'DSCC1', 'SMS', 'SMC2',
       'RFC4', 'ATP1B3', 'LSM3', 'SRSF9', 'CTNNAL1', 'SEC61B', 'PFN1',
       'MANF', 'BRCA2', 'PARK7', 'MAGOHB', 'SLC29A1', 'DTL', 'SSBP1',
       'MT2A', 'HPRT1', 'HSP90B1', 'NDUFS6', 'POLR3K', 'BIRC5', 'ESCO2',
       'PGAM1', 'DNPH1', 'TMPO', 'MRPL3', 'HSPA5', 'PFDN2', 'SLC1A5',
       'CCNE1', 'TIMELESS', 'ANAPC11', 'SNU13', 'POLR2F', 'TUBB4B',
       'STOML2', 'PPIA', 'MSH6', 'HNRNPAB', 'NUSAP1', 'SNRPC', 'CCNB2',
       'COA4', 'RBX1', 'NME4', 'UQCC2', 'SPCS2', 'CDKN3', 'MND1', 'CYC1',
       'PDIA3', 'GLRX5', 'PPM1G', 'GCSH', 'STMN1', 'CENPX', 'NDUFA6',
       'CENPV', 'TMEM97', 'SUPT16H', 'NOP10', 'NCAPG', 'POLE4',
       'RNASEH2C', 'COX5B', 'SNRPA1', 'ATAD5', 'HNRNPR', 'ENO1',
       'CCDC167', 'FANCI', 'PSMB5', 'PHB', 'GTF3C6', 'MRTO4', 'COX6C',
       'PSMB2', 'TMA7', 'PGK1', 'PSMA3', 'POLE3', 'UBE2N', 'SEC61G',
       'MRPS6', 'PRIM1', 'TPI1', 'POP7', 'CBX5', 'TXNDC17',
]

GENES_HIGHER_IN_CCDC88A_HIGH_CD34_ENRICHED_PB_B = [
    # NOTE: for our data, maybe better to just use CCDC88A alone. it is quite strong, and more specific than when using all of these genes.
    
    # filtered by marker heatmap:
    'CCDC88A',
    ###'RHOC',
    ###'TRAPPC2L',
    'RGS10',
    ###'DBI',
    'MYL9','LINC01480','CLNK','ADTRP','ABCA6','FMOD','APOD','PIGR',

    'FILIP1L', 'TSHZ2', 'WNT3', 'LEF1', 
    ###'AL591845.1', # doesn't appear in oneK1K
    ###'INPP5F',
    'SFTPB', 'KSR2', 'CHDH', 'RASGRF1', 
    ###'PTPN7', 
    'SPATS2L', 
    ###'FCMR',
    'FBXO27', 'DRAIC', 
    ###'HCST', 
    'ATP2A1', 
    ###'NOSIP', 
    'ROR1', 
    ###'CMTM3',
    'ADAM29', 
    'ABCA9', 
    'TRAC', 'IGFBP4', 
    ###'CYTH1',

    ###'TSTD1', 
    'CHST2', 
    ###'SRGN', 
    ###'SEPT7', 
    ###'YBX3', 
    'HACD1', 'MRC2',
    ###'CCDC141', 
    ###'FCRL2', 
    ###'CREB3L2', 
    ###'GLO1', 
    ###'AC010524.1', # doesn't appear in oneK1K
    'SYNE2',
    'ADAMTS6',
    'GPR34',
    ###'AC008105.3', # doesn't appear in oneK1K
    'ALDH5A1',
    'SFMBT1',
    # not: 'TRAPPC2L','CYTH1','CCDC141','GLO1','RHOC','SRGN','HCST','DBI','CMTM3','PTPN7','NOSIP','FCRL2','SEPT7','FCMR','INPP5F'
    # not: 'TSTD1','YBX3','CREB3L2','IGFBP4'

    # orig list
    # 'CCDC88A', 'FCER2', 'S100A6', 'TCF4', 'ADTRP', 'FCMR', 'CD27',
    #    'LINC01857', 'VOPP1', 'ABCA6', 'S100A4', 'CLNK', 'PIGR', 'APOD',
    #    'LINC01480', 'IGLC1', 'TRAC', 'PEBP1', 'HCST', 'RGS10', 'VIM',
    #    'ARHGAP24', 'SMARCB1', 'SRGN', 'MYL9', 'NOSIP', 'FMOD', 'DBI',
    #    'RHOC', 'YBX3', 'LILRA4', 'ID3', 'IGHG2', 'IGFBP4', 'FCGRT',
    #    'TRAPPC2L', 'SH3BP5', 'CTSH', 'GSTP1', 'IGFLR1', 'TSTD1',
    #    'APOBEC3G', 'IGLC2', 'GPX1', 'IGHA1', 'BCL2', 'MIR155HG', 'PIM1',
    #    'PMAIP1', 'CKAP4', 'ZBTB32', 'ITGB1', 'ATOX1', 'CTLA4', 'ACTG1',
    #    'ARPC1B', 'LITAF', 'AC074289.1', 'LSM2', 'SEPTIN9', 'ISCU',
    #    'CORO1B',
]
GENES_HIGHER_IN_ZEB2_HIGH_B_IN_ONEK1K = [
    # >2
    'CIB1', 
    ###'CRIP1', 
    'FCRL5', 
    ###'IFI30', 
    'EMP3', 
    'FGR', 
    ###'NEAT1', 
    
    'ZEB2', 
    'IGHG3', 
    ###'S100A11', 
    'HCK', 
    'HSPB1',

    # 1.5-2
    'MPP6', 
    ###'COTL1',
    'ITGB2', 
    ###'MS4A1', # maybe should uncomment...
    'LY6E', 
    'MACROD2', 
    'PSAP', 'GRN', 
    
    'DAPP1', 'RHOB', 
    
    
    ###'IGHG1', 
    ###'ITGB2-AS1',
    'TNFRSF1B', 
    ###'SAT1', 
    'MAP3K8', 
    ###'RP5-1171I10.5', 
    'GABARAPL2', 
    'SYT1',

    # manually added
    'ITGAX',
    'CD19',

]
GENES_HIGHER_IN_ZEB2_HIGH_B = [
    *GENES_HIGHER_IN_ZEB2_HIGH_B_IN_ONEK1K,
    'SCIMP', 

    'FCRL3', 
    'PCDH9', 
    ###'POU2F2', 
    'TENT5A', 
    ###'ITGB2-AS1', 
    'CD72',
    'LAPTM5', 'SLC11A1', 'FGD4', 'MS4A7',

    # higher when comparing cells with highest vs lowest sum of above genes
    'LRMP',
    'PLEK',
    'GDI2',
    'TBX21',
    'ADGRE5',
    'CYTH4',
    'PLEKHO1',
    'TSPAN3',
    'LITAF',
    'MTSS1',
    'DTX1',
    'PLXNC1',

    
    # # Strong over expression
    # 'FCRL5',
    # 'ZEB2',
    # 'DAPP1',
    # # Weaker overexpression
    # 'CD22', 'SCIMP', 'EMP3', 'MPP6', 'FGR', 'CD19',
    # # Strong down regulation
    # 'LTB',
]
GENES_LOWER_IN_ZEB2_HIGH_B_IN_ONEK1K = [
    ###'CXCR4', # looks meh in oneK1K
    'VPREB3', 'LTB', 
    
    'JCHAIN', 
    ###'TCL1A', # looks like a mistake
    'CD24', 'SELL', 
    'MARCKS',
]
GENES_LOWER_IN_ZEB2_HIGH_B = [
    *GENES_LOWER_IN_ZEB2_HIGH_B_IN_ONEK1K,

    ###'S100A4',
    ###'LINC01781',
    'ID3', 'METTL8', 
    ###'TCF4', 
    'QRSL1', 'TRAF5',
    ###'CXCR4',
]
GENES_HIGHER_IN_CLP_M_PRSS2_COMPARED_TO_OTHER_CLP_M = [
    # 'PRSS2',
    'LYN',
    'LIMS1', 
    'GCSAML',
    'MZB1',
    'TCTEX1D1',
    'SVIL',
    'LAPTM4B',
    'ADGRG1',
    # 'LINC01122', # not a great separator
    # 'CYTL1', # not a great separator, but might hint at a potential novel state.
]



HIGHER_IN_HSC_THAN_MEBEMP_L_AND_ALSO_CLP_M = [
    # >1 in both
    'AVP', 'SELL', 'MDK', 
    ###'CSF3R', 
    'CRHBP', 
    ###'C1QTNF4', 
    'BEX1',
       ###'AC011139.1', 
       'MEG3', 'AC026369.3', 
       ###'HSPB1', 
       'MLLT3', 'IFITM3',
       'MSRB3', 'SOCS2', 'PRDX1',
    # >0.6 in both (but not >1 in both)
    'PTMS', 'SKAP1', 'LAPTM4B', 'NPDC1',
       ###'VIM', 
       'CLEC9A', 'SELENOP', 'TM4SF1', 'NOG', 'PREX2', 
       ###'PHLDB2',
       'SOCS3', 
       ###'SCN9A', 
       ###'XBP1', 
       'NPR3', 
       ###'ANXA5',
       'NFKBIZ', 
       ###'ATP8B4',
       'BSPRY', 
       ###'MGST1', 
       'AL590226.1', 'LIMCH1', 'BST2', 'HLF', 'BEX5',
       'CDH7', 'KRT8', 
       ###'SORL1',
    # >0.5 in both (but not >0.6 in both)
    'CRYGD', 'PCBD1', 'CCDC175', 'BEX2',
       ###'PCDH9', 
       'TLE4', # also from manuscript fig1c
       ###'ANKRD28', 
       'UCHL1', 
       ###'TFEC', 
       ###'DLK1', 
#        'RPL13',
       'RABGAP1', 'CCDC42', 
#        'HINT1', 
       'CALN1', 'TCEAL2', 
       ###'EVA1B',
       'TMEM163', 
       ###'PPA1', 
       'CLDN15',
       
       
       'HOXB5', 'GATA3', # from manuscript fig1c
]











RIBO_GENES_SLIGHTLY_HIGHER_IN_HIGH_HLF_HSC_COMPARED_TO_MPP = [
    'RPL34', 'RPL39', 'RPL12', 'RPS18', 'RPS10', 'RPL13', 'RPL32',
       'RPL11', 'RPS3A', 'RPS12', 'RPL21', 'RPL26', 'RPL30', 'RPL9',
       'RPS4X', 'RPS9', 'RPL38', 'RPS19', 'RPL35A', 'RPS27A', 'RPS13',
       'RPL36', 
       ###'RPL39L', 
       'RPS24', 'RPL23', 'RPS6', 'RPS23',
]

MOST_ANTI_CORR_WITH_BEMP_GENES_ACROSS_BEMP_AND_MEBEMP_L = [
    # (0.8, 0.85) corr
    'YBX3', 'EIF4A2', 'ADA', 'EREG', 'SLC25A3', 'MPIG6B', 'PPA1',
       'SOD2', 'RIPOR3', 'FAM30A', 'LINC02573', 'CD164', 'MGLL', 'NPM1',
       'SEPT6', 'LYSMD2', 'MPP6', 'IGSF10', 'MICAL2', 'RGS10', 'KIFAP3',
       'CDCA7L', 'ST6GAL2', 'GUCY1B1', 'SH3BGRL', 'COMT', 'SNCA',
       'AKR7A2', 'RAP1B', 'ADGRG1', 'MYCT1', 'CAVIN1',
    
    # >0.85 corr
    'MED12L', 'TALDO1', 'FAM69B', 'STOM', 'PKIG', 'NAA38', 'GSTP1',
       'CYTL1', 'HLA-DPB1', 'C2orf88', 'HLA-DPA1', 'HTR1F', 'HINT1',
       'LSM7', 'NAP1L1', 'EGFL7', 'TPM1', 'CAT', 'BCL11A', 'B3GNT2',
       'ARHGAP22', 'PLAC8', 'CD34', 'MLLT3', 'IFI16', 'H2AFY', 'HLA-DRB1',
       'NPR3', 'CD74', 'PLXNC1', 'LAPTM4B', 'HLA-DRA', 'MEF2C', 'LMO2',
]

HIGHER_IN_LATEST_EP_THAN_MEBEMP_L_GENES = [
    # 240107: mc_utils.get_diff_exp_df(n_mc_ad, (n_mc_ad.obs['state'] == 'ERYP') & (mc_utils.get_genes_expr(n_mc_ad, 'BLVRB') > -11), n_mc_ad.obs['state'] == 'MEBEMP-L')
    # log_ratio > 2
    'HBB', 'APOC1', 'BLVRB', 'HBD', 'MYC', 'CA1', 'TMEM14C', 'MPC2',
       'REXO2', 'APOE', 'ANK1', 
    #    'S100A6', 
       'PVT1', 'FAM178B', 'NFIA',
       'CD36', 'TMEM14B',
]
MANUALLY_CHOSEN_ERY_SPECIFIC_GENES = [
    'HBB', 'HBD', 'AHSP',
    'HBA1', 'HBA2',
    'HBG1', 'HBG2',
]

APOPTOTIC_GENES_FROM_NIMROD = [
    # IIUC, this is specific to the experiments in which they mistakenly used a media which was toxic to the cells?
    'CDC42', 'PAQR7', 'FAM76A', 'AL353622.1', 'RBBP4', 'SYNC', 'AGO3', 'NFIA', 'FUBP1', 'PTBP2', 'CAPZA1', 'HIPK1', 'SEC22B', 'ARNT', 'C1orf56', 'CDC42SE1', 'KRTCAP2', 'ATP1B1', 'TOR1AIP2', 'RNF2', 'TROVE2', 'MDM4', 'NUCKS1', 'AL390728.6', 'IAH1', 'ADAM17', 'DNMT3A', 'PPP1CB', 'WDR43', 'YIPF4', 'CRIM1', 'SPTBN1', 'CCDC88A', 'BCL11A', 'MEIS1', 'MEIS1-AS2', 'GGCX', 'C2orf68', 'RGPD5', 'TLK1', 'CERKL', 'CNOT9', 'ANKRD28', 'EPM2AIP1', 'CTNNB1', 'RBM6', 'FOXP1', 'NFKBIZ', 'CD47', 'DPPA4', 'PARP15', 'MED12L', 'GPR171', 'IGSF10', 'MBNL1', 'RAP2B', 'AC106707.1', 'TRIM59', 'GOLIM4', 'TBL1XR1', 'DCUN1D1', 'PDE6B', 'CTBP1', 'TMEM33', 'PPP3CA', 'MARCH6', 'PURPL', 'NADK2', 'MAP3K1', 'OCLN', 'SERINC5', 'ARRDC3', 'SLC12A2', 'FNIP1', 'G3BP1', 'SFXN1', 'HNRNPH1', 'SERPINB9', 'RREB1', 'PHACTR1', 'C6orf62', 'HIST1H2BD', 'RUNX2', 'TSPYL1', 'RNF217', 'MYCT1', 'ARID1B', 'WTAP', 'FGFR1OP', 'AC092171.5', 'ACTB', 'AHR', 'TRA2A', 'HOXA3', 'AC004951.1', 'SUMF2', 'PHKG1', 'GIGYF1', 'THAP5', 'AC058791.1', 'CASP2', 'EIF2S3', 'ZFX', 'DDX3X', 'AC234772.3', 'PIM2', 'OGT', 'TSIX', 
    'JPX', 'FTX', 'RAP2C', 'PRRG3', 'ZDHHC2', 'CNOT7', 'DOCK5', 'NECAB1', 'TP53INP1', 'CBWD1', 'B4GALT1', 'CBWD3', 'ABCA1', 'SCAI', 'SET', 'TRIM5', 'WEE1', 'BEST1', 'MTA2', 
    'AP001157.1', 'SESN3', 'DDX6', 'FLI1', 'CELF2', 'CACNB2', 'PIP4K2A', 'MTPAP', 'SAR1A', 'PTEN', 'PDCD4', 'TIAL1', 'CCND2', 'ETNK1', 'CPNE8', 'ARID2', 'MAP3K12', 'RASSF3', 'MDM2', 'NUDT4', 'MED13L', 'KMT5A', 'ZNF26', 'PAN3', 'AL138963.3', 'INTS6', 'TPP2', 'ARGLU1', 'ANKRD10', 'ARF6', 'ATP5S', 'PSMA3-AS1', 'PCNX4', 'STON2', 'CHP1', 'MYO5C', 'TLE3', 'MPV17L', 'THUMPD1', 'TMEM159', 'SPN', 'C16orf54', 'AC026471.4', 'AMFR', 'EXOSC6', 'AP1G1', 'PLCG2', 'YWHAE', 'EIF5A', 'TP53', 'LRRC75A', 'SMCR5', 'WSB1', 'AC010761.1', 'DCAKD', 'TOB1', 'NME2', 'NOG', 'MSI2', 'VMP1', 'AC015802.6', 'AFMID', 'ZNF521', 'TCF4', 'SRSF6', 'INSR', 'ZNF121', 'ZNF562', 'CALR', 'AC092069.1', 'JAK3', 'ZNF431', 'ZC3H4', 'LENG8', 'LENG9', 'MIF', 'CCDC117', 'APOBEC3C', 'TNRC6B', 'MATR3', 'POLR2J3', 'LINC01355', 'EPS15', 'SLC30A7', 'ANP32E', 'RC3H1', 'HIST3H2A', 'TMSB10', 'GLS', 'EEF1B2', 'ARPC4', 'OSBPL10', 'ZNF662', 'ZBTB20-AS2', 'H1FX', 'ATP13A3', 'FAM184B', 'HNRNPDL', 'STARD4-AS1', 'PGGT1B', 'CDC42SE2', 'NPM1', 'RACK1', 'SERPINB1', 'ATXN1', 'EEF1A1', 'MMS22L', 'ZNF117', 'RASA4', 'PRPS2', 'TMSB4X', 'TRAPPC2', 'ANGPT2', 'YWHAZ', 'CBWD5', 'CARNMT1', 'AAED1', 'ZNF462', 'SNHG7', 'CPT1A', 'ZNF518A', 'ETV6', 'SLC48A1', 'HNRNPA1', 'AC084033.3', 'CPSF6', 'AC055713.1', 'UBC', 'NUP58', 'TPT1', 'PRPF39', 'SRSF5', 'ELMSAN1', 'CRIP1', 'TBC1D24', 'C16orf72', 'ZNF785', 'KDM6B', 'FLCN', 'SGK494', 'EIF1', 'H3F3B', 'ACOX1', 'ACTG1', 'PIAS2', 'SMAD2', 'ZFAS1', 'LYL1', 'ZNF506', 'HNRNPL', 'FTL', 'NDUFA3', 'CBX6', 
    # 'XIST', 
    # 'NEAT1', 
]

HIGHER_IN_MEBEMP_E_THAN_MONOCYTE_GENES = [
    # >3, and also monocytes expression <-13
    'CDK6', 'AC084033.3', 'PRSS57', 'STMN1', 'LDHB', 'ANKRD28', 'SOX4',
    'MSI2', 'MYB', 'SPINK2', 'CYTL1', 'TFDP2', 'SMIM24', 'SLC40A1',
    'NPR3', 'SNHG8', 'PCDH9', 'EGFL7', 'NUCB2', 'TESPA1', 'SNHG7',
    # 2.5-3, and also monocytes expression <-13
    'FAM117A', 'PRKACB', 'SEPT6', 'EBPL', 'LAPTM4B', 'GATA2', 'CD34',
    'FAM30A', 'C12orf57', 'CPA3', 'BEX3', 'ZNF385D', 'PBX1', 'MTURN',
    'AL157895.1',
]
BEST_HIGHER_IN_MEBEMP_L_THAN_MPP_GENES = [
    # 240107: mc_utils.get_diff_exp_df(n_mc_ad, n_mc_ad.obs['state'] == 'MEBEMP-L', n_mc_ad.obs['state'] == 'MPP')
    # log_ratio > 1.2
    'SLC40A1', 
    ###'HBD',
    'CNRIP1', 'KLF1', 'TFRC', 'CSF2RB', 
     ###'S100A4',
     ###'S100A6', 
       'GATA1', 'PDZD8', 'PLEK', 'CTNNBL1', 'CPA3', 'FCER1A',
       'MINPP1', 'ZEB2', 'CYTL1', 'CD84', 'ZMYND8', 'PRKAR2B', 'BMP2K',
       'ABCC4', 'GATA2', 'TLN1', 'CYTOR', 'CASP3', 'MYC', 'PDCD4', 'PKIG',
     ###'AL034397.3', 'AL157895.1', 
       'PSTPIP2',
]
HIGHER_IN_NKTDP_THAN_MPP_GENES = [
    # log_ratio > 3
    'ID2', 'ACY3', 'SAMHD1', 'LGALS1', 'JCHAIN', 'RUNX2', 'RUNX3',
       'CORO1A', 'LSP1', 'TTN', 'DDIT4', 'MEF2A', 'CYTH4', 
       ###'CALM2',
       'NFIL3', 'ITGB7', 'RAB11FIP1', 'NFKB1', 'TYROBP', 'SCPEP1',
       'GPR183', 'LCP1', 'TESC', 'CAPG', 'MEF2C', 'TNFRSF18', 'PLD4',
       'FYB1', 'BLNK', 'ATP10A', 'SLC38A1', 'MAF', 'NCF1', 'PIK3AP1',
]





BM_DC_PROG_SOMEWHAT_SPECIFIC_GENE_NAMES = [
    # ugh. might be a bit high in CLP-M??

    'LGALS1', 'IRF8', 'RNASE6', 'PLD4', 'SAMHD1', 'C12orf75', 'CCDC50',
       'CST3', 'PLAC8', 
       ###'CORO1A', # appears also in mds_analysis_params.HIGHER_IN_CLP_M_THAN_HSC 
       'PLP2', 
       'IGKC', 
       'HERPUD1', 'FCER1G',
       'LSP1',
]

BM_POST_GMP_MONOCYTE_SOMEWHAT_SPECIFIC_GENE_NAMES = [
    'LYZ', 'MPO', 'LGALS1', 'CSTA', 'AZU1', 'ELANE', 'CTSG', 'PRTN3',
       'SAMHD1', 'CFD', 'AC020656.1', 'SRGN', 'CEBPD', 'RNASE2', 'ANXA2',
       'IRF8', 'TYROBP',
]

BM_SELL_OUTLIER_GENE_NAMES = [
    'SELL', 
    # 'HSPA5',
    'SLC3A2', 'ICAM3', 'ANKRD28', 'APLP2', 'CTSD', 'SPNS3',
]

HCA_BM_HIGHER_IN_CLP_E_CANDIDATE_THAN_CLOSE_HSC_AND_ALSO_CLP_M = [
    # NOTE: 231226: took these genes and filtered them to get rid of problematic ones. looked at rp_s vs dntt_s plot for BM CD34+ palantir relevant cells and these genes were not higher where i expect to see CLP-E. and also there are outliers in palantir that are high in these genes... though some celebrities here (SF3B1, RUNX1, LMNA, FLT3). anyway, too few cells. could try to look at MPP and MLP in https://www.nature.com/articles/s41467-019-10291-0/figures/1... maybe CLP-E could be found there? though maybe best to just assume it is a transient state, and move on.
    'FLT3', 
    'CCNL1', 'CYTIP', 'IDH2', 'TARBP1', 'FAM172A', 'GPATCH8',
       'SF3B1', 
       'LMNA', 
       'CALM1', 'NDUFS7', 'RGL4', 'SEPT9', 
       ###'FAM103A1', # missing from palantir h5ad
       'NUPL2', 'THUMPD3', 'TRA2A', 'GNB1', 'SEPT1', 'GGNBP2', 'ORMDL1',
       'GON4L', 'VAV1', 'RCSD1', 'PAN3', 'SH3GLB1', 'RHBDD2', 'TM2D1',
       'C9orf78', 'SETDB1', 'ICA1', 
       ###'MYL12A', # a bit too high basal level
       'DGUOK', 
       'TCF25', 
       'DDX46',
       'KDM2A', 'PSIP1', 'SQSTM1', 'TM7SF2', 'LUZP1', 'ATF7IP2', 'DOCK2',
       'MANBA', 'BLCAP', 'RGS14', 'ILF3', 'CDK2AP1', 
       ###'HSPA5', # presumably lateral function
       ###'JUNB',
       'GNG11', 'CNIH1', 'FAM111A', 'KCNAB2', 'GOLGA2', 'TMEM187',
       'TIMM17B', 'LGALS3BP', 'TMEM91', 'VANGL1', 'LDLRAD4', 'CREB1',
       'IDE', 'VMA21', 'AP3S1', 'IFT57', 'TERF2IP', 'JUP', 'PLA2G12A',
       'TMEM256', 'SREK1', 
       'RUNX1', 
       'SMIM12', 'KRTCAP2', 'REST',
       ###'COX7A2L', # a bit too high basal level
       'RNF167', 'AARSD1', 
       ###'ATP5C1', # missing from palantir h5ad
       'SACS', 'HOOK3', 'NSMAF',
       'LEMD2', 'SLC30A9', 'WWOX', 'MTRF1L', 'NT5C', 'PTPN4', 'THAP1',
       'MFSD14A', 'HERC1', 
       ###'BRK1', # a bit too high basal level
       ###'MZT2B', # a bit too high basal level
       'ABCA2', 
       ###'C11orf84', # missing from palantir h5ad
       'GTF2I',
       'LINC01003', 'TCAF1', 'PRKCI', 'RBMS1', 'PARP4', 'RIOK2', 'GLOD4',
       'FARSA', 
       ###'RP11-386I14.4', # missing from palantir h5ad
       'KAT6B', 'FAM126B', 'CARNMT1', 'IQGAP2',
       'HNRNPUL1', 'VAMP5', 'PIEZO1', 'OGFOD1', 'RNF139', 'HERPUD1',
       'KIF9', 'MCPH1', 'XPC', 'XPO1', 'PPP4R3A',
]
BM_CORRELATED_WITH_MKI67_AND_TOP2A_IN_BOTH_MPP_AND_MEBEMP_L_GENE_NAMES = [
    # in palantir BM CD34+ calculated corrs with ['MKI67', 'TOP2A'] both in MPP and in MEBEMP-L, and then took all genes with corr > 0.75 in both cases.
    'ESCO2', 'KIF20B', 
    ###'DTYMK', # a bit high basal level
    ###'CKS1B', # a bit high basal level
    ###'NUCKS1', # a bit high basal level
    ###'SMC2', # a bit high basal level
    'CENPN',
       'KIF2C', 'CCNB1', 'SPC25', 'KIFC1', 'NCAPG', 'CKAP2L', 'CENPE',
       'CDCA3', 'KIF15', 'NUF2', 'CDK1', 'CDCA8', 'HMMR', 'ASPM',
       'DLGAP5', 'GTSE1', 'CCNA2', 'CDKN3', 'NUSAP1', 'TPX2', 'CENPF',
       'AURKB', 'UBE2C', 'BIRC5', 'MKI67', 'TOP2A',
]
BM_HIGHER_IN_GMP_L_THAN_CLP_E_CANDIDATES_GENE_NAMES = [
    # >2.5
    'MPO', 'AZU1', 'ELANE', 'CTSG', 'PRTN3', 'RNASE2', 'CFD', 'LYZ',
       'SRGN',
       'CALR', 'HSP90B1', 'IGLL1', 'NPW', 'RAB32', 'LGALS1',
       'CSTA', 'CST7', 'MS4A3', 'CEBPD', 'SDF2L1', 'MGST1',
]
BM_HIGHER_IN_B_AND_PRE_B_THAN_CLP_E_CANDIDATES_GENE_NAMES = [
    # >2.5
    'VPREB3', 'VPREB1', 'CD24', 'IGLL1', 'CD79A', 'LEF1', 'CXCR4',
       'ODC1', 'IGHM', 'CCND3', 'CD79B', 'PAX5', 'SRM', 'AKAP12', 'MS4A1',
       ###'NME1', 
       'PPP1R14B', 'UHRF1', 'EBF1', 'CD81', 'CD19', 
       ###'FABP5',
       'MRPL14',
         'MAD2L2', 
         ###'CYCS', 
         ###'PCNA',
           'ZCCHC7', 'CD22', 'PGAM1',
]
BM_HIGHER_IN_SOME_LARGE_PRE_B_THAN_CLP_E_CANDIDATES_GENE_NAMES = [
    # >4
    'LTB', 'DNTT', 'VPREB1', 'CD79A', 'CYGB', 'CD79B', 'UHRF1',
       'IGLL1', 'VPREB3',
]
BM_HIGHER_IN_MKP_POST_MKP_THAN_CLP_E_CANDIDATES_GENE_NAMES = [
    # >4
    'ITGA2B', 'PLEK', 'PF4', 'PDLIM1', 'TPM1', 'LTBP1', 'RGS18',
       'AL157895.1', 'PPBP', 'FERMT3', 'RAP1B',
]
# BM_HIGHER_IN_PRO_B_THAN_CLP_E_CANDIDATES_GENE_NAMES = [ # no need. almost identical to the signature i constructed for PB.
#     # >2.5
#     'DNTT', 'VPREB1', 'LTB', 'IGLL1', 'CD79A', 'UHRF1', 'ADA', 'CYGB',
#        'JCHAIN', 'ZCCHC7',
# ]

GENES_LOWER_IN_CLP_M_PRSS2_COMPARED_TO_OTHER_CLP_M = [
    'ITGA4',
    'RUNX2',
    'TMSB10',
    'AL096865.1',
    'LY86',
    'ACY3',
    'BLNK',
    'CRYBG3',
    'ID2',
    'NPTX2',
    'LSP1',
    # 'IL7R', # not a great separator
    # 'LGALS1', # not a great separator
]

NUCLEUS_RELATED_GENES = [
    'PTPRC',
    'NRIP1',
    
    # 'DEK',
    'ATRX',
    'ARGLU1',
    'MIS18BP1',
    'KTN1',

    'TXNIP',
    'PNRC1',
    'PNISR',
    'ITGA4',
    'PRRC2C',
    'N4BP2L2',
    'SON',
    'GABPB1-AS1',
    'HNRNPU',
    # 'NCL',
    # 'HSP90AA1',
    'ANKRD12',
    'LUC7L3',
    'BPTF',
    'CDK6',
    'TLN1',
    'LRRFIP1',
    'SCAF11',
    # 'HMGB1',
    'TPR',
    'PCM1',
    # 'NEAT1',
    'PAK2',
    'DDX17',
    'PNN',

    # 'ZEB2',
    'AKAP9',
    'MBNL1',
    'ZRANB2',
    'EIF3A',
    'RBM25',
    'GOLGB1',
    'JAK1',
    'NAP1L1',
    'CELF2',
    # 'SRRM1',
    # 'IQGAP2',
    # 'LCP1',
    # 'SEPT7',
    # 'MTDH',
    'BOD1L1',
    'CCDC88A',
    # 'TPM3',
    # 'PSIP1',
    # 'TCEA1',
    'PHF3',
    # 'CCNI',
    # 'SET',
    'PRPF38B',
    'XRCC5',
    'TTC3',
    # 'CAST',
    'MTPN',
    # 'SERBP1',

    # 'SEC62',
    'SRSF11',
    'RBM39',
    'HNRNPA3',
    'DDX5',
    # 'NUCKS1',
    'VPS13C',
    'HP1BP3',
    'HADHA',
    'ACTR2',
    # 'HSP90AB1',
    # 'SRRM2',
    'TNRC6B',
    'CHD9',
    'FNBP1',
    'NKTR',
    'TAOK3',
    'ZBTB20',
    'PNRC2',
    # 'IFI16',
    'ROCK1',
    'KMT2A',
    'NIPBL',
    'MACF1',
    'ZC3H13',
    'PPIG',
    'GOLGA4',
    # 'HSP90B1',
    # 'XRN2',

    'SERBP1',
    'HSP90AB1',
    'STAG2',
    'IFI16',
    'SUPT16H',
    'SPN',
    'REL',
    'SEPT7',
    'USP15',
    'SEC62',
    'BCLAF1',
    'SEPT6',
    'FTX',
    'YWHAZ',
    'ANKRD28',
    # 'MEF2C',
    'SRRM1',
    'SP100',
    'DDX24',
    'ZNF638',
    'SRRM2',
    'MTDH',
    'RIF1',
    'TCEA1',
    'FUS',
    'SLTM',
    'FOXP1',
    'TOP1',
    'KMT2E',
    'ZEB2',
    'DDX46',
    'NUCKS1',
    'APPL1',
    'ARID4B',
    'HSP90B1',
    'RSL1D1',
    'CD46',
    'NEAT1',
    'SENP6',
    'KIF2A',
    'ITPR2',
    'OSBPL8',
    'PSIP1',
    'ANKRD11',
    'SET',
    # 'SLC40A1',
    'HSP90AA1',
    'CHD4',
    'GPBP1',
    'JMJD1C',
    'YBX1',
    'UPF2',
    'MPHOSPH8',
    'GTF2I',
    'PRPF40A',
    'IQGAP2',
    'PRPF4B',
    # 'SMAP2',
    'REST',
    'ATXN7L3B',
    'U2SURP',
    'HMGB1',
    'TUT4',
    'TAX1BP1',
    'ITGB1',
    'AHNAK',
    'XRN2',
    'GCC2',
    'LCP1',
    'CLINT1',
    'ELF1',
    'MYCBP2',
    'NCL',
    'SMCHD1',
    'APLP2',
    'TPM3',
    'UBXN4',
    'CAST',
    'DHX36',
    'STK4',
    'CCNI',
    'DEK',
    # 'BTG1',
]


BATCHY_NEUT_STATE_SPECIFIC_GENE_NAMES = [
    'ANXA3',
    'ZDHHC19',
    # 'CMTM2', # not specific enough when considering monocytes
    'IL1R2',
    # 'S100P', # not specific enough when considering monocytes
    'FCGR3B',
    'CAMP',
    'LCN2',
    'ARG1',
    'CD177',
    'LTF',
    'PGLYRP1',
    'CRISP3',
    'MMP8',
    'ALPL',
    'MGAM',
    'MMP9',
    'CYP4F3',
]
GENES_HIGHER_IN_BATCHY_MONOCYTE_MEBEMP_DOUBLET_THAN_DC_AND_LOW_IN_DC_NAMES = [
    'S100A12', 'VCAN', 'FCN1', 'CD14', 'SERPINA1', 'CSF3R', 'GIMAP4',
    'ANKRD28', 'SPINK2', 'DUSP6', 'GIMAP7', 'FPR1', 'SOD2', 'SLC40A1',
    'IRAK3', 'MARCKS', 'CEBPB', 'CYP1B1', 'CLEC7A', 'FAM117A', 'NPR3',
    'SMIM24', 'SORL1', 'C5AR1', 'CEBPD', 'LRRK2', 'SLC11A1', 'ZFP36L1',
    'FAM30A',
    # 'S100A8', 'S100A9',
]
GENES_HIGHER_IN_BATCHY_MONOCYTE_MEBEMP_DOUBLET_THAN_MONOCYTES_NAMES = [
    'CDK6', 'STMN1', 'AC084033.3', 'ANKRD28', 'SOX4', 'SPINK2',
    'FAM30A', 'SMIM24', 'MSI2', 'NUCB2', 'SNHG8', 'PRSS57', 'TFDP2',
    'NPR3', 'MYB', 'LDHB', 'HINT1', 'CD34', 'FAM117A', 'CPA3', 'NPM1',
    'BEX3', 'C12orf57', 'PCDH9', 'SLC40A1', 'TFPI', 'PRKACB', 'TESPA1',
    'LAPTM4B', 'SEPT6', 'BEX2',

    'CYTL1', 'AC002454.1', 'GUCY1A1',
    'ZFAS1', 'TSC22D1', 'RPS5', 'AC011139.1', 'GABPB1-AS1',
    'LINC02573', 'EEF1B2', 'PBX1', 'GATA2', 'EBPL', 'AKR1C3', 'FCER1A',
    'EGFL7', 'ITM2A', 'HSP90AB1', 'ANGPT1', 'HEMGN', 'MSRB3', 'TSTD1',
    'MLLT3',

    'MTURN', 'FHL1', 'PEBP1', 'RPS21', 'MAP7', 'GYPC', 'NREP',
    'NRIP1', 'RPS18', 'DPPA4', 'CCND2', 'RPSA', 'SOCS2', 'SNHG7',
    'SERPINB1', 'SEPT11', 'SPN', 'SNHG19', 'KDM5B', 'ITM2C', 'HTR1F',
    'IGHM', 'AL157895.1', 'HMGN1', 'HOPX', 'TCF4', 'KIT',
]
GENES_LOWER_MONOCYTE_THAN_NKTDP_L_NAMES = [
    'SPINK2', 'SOX4', 'STMN1', 'JCHAIN', 'ACY3', 'RUNX2', 'ITM2C',
    'FAM30A', 'TSC22D1', 'SLC38A1', 'IGHM', 'SATB1', 'DDIT4', 'SMIM24',
    'HMGN1', 'CTSW', 'CDK6', 'ID2', 'TCF4', 'NUCB2', 'NPM1', 'PEBP1',
    'LDHB', 'MAF', 'TRBC2', 'TTN', 'SEPT6', 'KLRB1', 'CLNK',
    'AL096865.1', 'OCIAD2', 'ODF2L', 'HOPX', 'TUT4', 'GUCY1A1',
    'HOXA9', 'LINC00865', 'STK17A', 'TFDP2', 'CCR7', 'STRBP', 'NREP',
    'ADAM28', 'TRIB2', 'DPPA4', 'PSIP1', 'ANKRD28', 'ARMH1', 'RUNX3',
    'RERE', 'PHACTR1', 'GPR183', 'NEGR1', 'BCL11A', 'N4BP2', 'KIT',
    'AFF3', 'ZFAS1', 'TNFRSF18', 'NYNRIN', 'NKG7', 'ARHGAP5',
    'AC084033.3', 'C12orf57', 'PCDH9', 'BEX2', 'HOXA10', 'CXXC5',
    'HINT1', 'PRKACB', 'BEX3', 'PNN', 'CIRBP', 'KAT6B', 'ADA', 'RPSA',
    'SNHG8', 'GABPB1-AS1', 'LTB', 'MYL6B', 'PAWR', 'HNRNPA1', 'APEX1',
    'FLT3', 'TTC3', 'ZNF22', 'TSTD1', 'PNISR', 'KIAA0087', 'KDM5B',
    'MACF1', 'HSP90AB1', 'RPS5', 'AC243960.1', 'SNHG7', 'MSI2', 'CCT2',
    'HIST1H1D', 'ITGB7', 'GSN', 'PARP1', 'C6orf48', 'BEX4', 'PDE7A',
    'AKR1C3', 'RNASEH2B', 'HADH', 'RPL3', 'FNBP1L', 'CYFIP2', 'RPAP2',
    'DDAH2', 'C12orf75', 'OFD1', 'PCM1', 'RASD1', 'NUCKS1', 'AUTS2',
    'BLNK', 'SSBP2', 'RPS3', 'MME', 'MLC1', 'RPLP0', 'CCNB1IP1',
    'PRSS57', 'IGHD', 'PRDX2', 'BCL2', 'WDR60', 'ZRANB2', 'BCL7A',
    'TNFRSF4', 'ZSCAN18', 'ACAP1', 'ZNF428', 'RHOH', 'IL18', 'TARBP1',
    'ATXN7L3B', 'S100PBP',
]
GENES_LOWER_IN_BATCHY_MONOCYTE_MEBEMP_DOUBLET_THAN_NKTDP_L_NAMES = [
    'JCHAIN', 'ID2', 'ACY3', 'RUNX2', 'SPINK2', 'DDIT4', 'MAF',
    'SLC38A1', 'KLRB1', 'ITM2C', 'RUNX3', 'SATB1', 'TTN', 'CLNK',
    'IGHM', 'CCR7', 'AL096865.1', 'TSC22D1', 'LINC00865', 'STK17A',
    'TRIB2', 'TNFRSF18', 'GPR183', 'CTSW', 'NKG7', 'GSN', 'ITGB7',
    'TRBC2', 'MEF2A', 'KIAA0087', 'MDFIC', 'ADA', 'NYNRIN', 'CYSLTR1',
    'RASD1', 'MEF2C', 'ITGA4', 'NEGR1', 'EVL', 'HMGN1', 'BLNK', 'SOX4',
]
GENES_HIGHER_IN_BATCHY_MONOCYTE_MEBEMP_DOUBLET_THAN_MEBEMP_NAMES = [
    'FCN1', 'FGL2', 'LYZ', 'CST3', 'CYBB', 'DUSP6', 'CTSS', 'CPVL',
    'GIMAP4', 'COTL1', 'VCAN', 'PSAP', 'TYMP', 'SAMHD1', 'LGALS3',
    'IGSF6', 'MPEG1', 'CLEC7A', 'AC020656.1', 'CD68',

    'FCER1G', 'MARCH1', 'MS4A6A', 'ANXA2', 'CFD', 'MAFB', 'CASP1', 'TNFRSF1B',
    'TYROBP', 'CLEC12A', 'KCTD12', 'CD48', 'MS4A7', 'CD300E',
    'TNFSF10', 'GRN', 'C1orf162', 'NCF2', 'NPC2', 'LGALS2', 'ITGB2',
    'SERPINA1', 'FCGRT', 'CD14',
     
    'CTSZ', 'P2RY13', 'CEBPD', 'CSTA',
    'S100A10', 'CFP', 'CD36', 'DMXL2', 'MARCKS', 'FYB1', 'SLC7A7',
    'APOBEC3A', 'ANXA5', 'TNFSF13B', 'LILRB2', 'GIMAP8', 'DOK2',
    'AP1S2', 'IL10RA', 'TMEM176B', 'CTSH', 'VSIR', 'CCR1', 'BRI3',
    'LYST', 'GIMAP7', 'LPAR6', 'EFHD2', 'CSF1R', 'CARD16', 'ADA2',
    'SCIMP', 'TLR2', 'PRELID1', 'IFITM3', 'AOAH', 'FTL', 'TGFBI',
    'GIMAP1', 'CTSB', 'HCK', 'WARS', 'VMP1', 'BLVRB', 'LGALS1', 'TNFAIP2',

    'HNMT', 'GPR65', 'PLBD1', 'TSPO', 'OAS1', 'RXRA',
    'S100A4', 'ARPC1B', 'GNS', 'SLFN5', 'ATP6V1B2', 'LILRB1', 'RNF213',
    'LST1', 'OGFRL1', 'FGD2', 'REEP5', 'PLEKHO1', 'HLA-DQB1', 'CX3CR1',
    'S100A11', 'HLA-DMB', 'IFNGR1', 'CEBPB', 'IQGAP1', 'LRRK2', 'TLR4',
    'CD93', 'SNX10', 'CD33', 'HLA-DRB1', 'ASGR1', 'FCGR2A', 'SAT1',
    'NAGK', 'OAZ1', 'CD86', 'STXBP2', 'GBP1', 'KLF2', 'MX1', 'IFI6',
    'CYBA', 'LRP1', 'LY86', 'CAST', 'PSMB9', 'MNDA', 'HLA-B', 'JAML',
    'RNASE6', 'TMEM176A', 'FLNA', 'GABARAP', 'HLA-DRA', 'GLRX', 'HRH2',
    'CD4', 'ABI3', 'NAAA', 'THEMIS2',
]
GENES_LOWER_IN_MONOCYTE_THAN_IN_BATCHY_CAMP_MONOCYTES = [
    # use both high 'S100A8' and 'S100A9' and also this?
    
    'MMP9', 
    'S100P',
    #  'S100A12', # hmmm
     'LTF',
     'AC084033.3',
     'LCN2',
     'CDK6',

       'CD177',
        'PGLYRP1',
        'SPINK2',
        'STMN1',
        'ANKRD28',
        'SOX4',
    

    # 'SNHG8', 'SMIM24', 'FAM30A', 'NPR3', 'MMP8', 'TFDP2', 'FCGR3B',
    #    'MSI2', 'CYP4F3', 'MYB', 'NUCB2', 'CAMP', 'CMTM2', 
    # #    'S100A8', # hmmm
    #    'PRSS57', 'C12orf57', 'CRISP3', 'SLC40A1', 'HINT1', 'ALOX5AP',
    #    'CPA3', 'CYTL1', 'AC002454.1', 'TFPI', 'FAM117A', 'ZFAS1', 'IGHM',
    #    'AVP', 'SLC2A3', 'RPL23A',

    # 'NAMPT', 'HEMGN', 'LDHB', 'SNHG19', 'PRKACB', 'HIST1H1D',
    #    'LINC02573', 'BEX3', 'AC011139.1', 'RPS29', 'BEX2', 'CD34',
    #    'PCDH9', 'ANGPT1', 'EGFL7', 'CST7', 'ANXA3', 'BASP1', 'CLC',
    #    'DPPA4', 'RPS21', 'TSTD1', 'LAPTM4B', 'RPS5', 'SEPT6', 'RPS18',
    #    'SNHG25', 'MSRB3', 'GUCY1A1', 'EBPL', 'TESPA1', 'ARG1', 'CD24',
    #    'NPM1', 'IL1R2', 'AKR1C3', 'MAP7', 'GABPB1-AS1', 'MCTP2',
    #    'TSC22D1', 'MS4A1', 'ZDHHC19', 'GATA2', 'HSPB1', 'ALPL', 'RPL31',
    #    'PTMS', 'SNHG7', 'MME', 'KIT', 'GYPC', 'ETS1', 'SVIP', 'MGAM',
    #    'CD82', 'SOCS2', 'SPTBN1', 'ZNF385D', 'EEF1B2', 'CCDC171', 'MLLT3',
    #    'FHL1', 'HSP90AB1', 'TCF4', 'ARMH1', 'IGLL1', 'RPS28', 'TRBC2',
    #    'MZB1', 'VNN2', 'PROM1', 'ODF2L', 'AC245060.5', 'PBX1', 'TSPAN2',
    #    'ITM2A', 'HIST1H4C', 'RGL4', 'GPSM2', 'NREP', 'RFLNB', 'IGKC',
    #    'NRIP1', 'SYNE1', 'TFRC', 'SPN', 'RPL22L1', 'SEPT11', 'TXN',
    #    'ST13', 'AC103591.3', 'MDK', 'RPL37A', 'HACD1', 'KDM5B', 'ADAM28',
    #    'GUCY1B1', 'CCND2', 'IGSF10', 'HIST1H1E', 'MAP7D3', 'SYNE2',
    #    'CYSTM1', 'AL157895.1', 'ST8SIA6', 'SMIM27', 'ZNF521', 'HMGN1',
    #    'ERG', 'C1orf21', 'HOPX', 'LRRC75A', 'PKIG', 'C1QTNF4', 'IFITM2',
    #    'HTR1F', 'RALGPS2', 'SCN9A', 'BCL7A', 'MED12L', 'RPSA', 'RPS8',
    #    'PRDX2', 'RPS27', 'PEBP1', 'BCL11A', 'EGR1', 'MTURN', 'PSIP1',
    #    'RPL39', 'STRBP', 'ATF7IP2', 'BANK1', 'ZMYND8', 'EREG', 'TCEAL4',
    #    'PNN', 'LSM7', 'CNRIP1', 'RPL36', 'RHOH', 'RPL37', 'RPL38', 'RHEX',
    #    'HMGN5', 'PADI4', 'RPL12', 'MEIS1', 'BAALC', 'FOXP1', 'ITM2C',
    #    'PLCG2', 'KIFAP3', 'SNRPD2', 'TFF3', 'ZBTB20', 'ANKRD36C',
    #    'AC137767.1', 'LTB', 'CD109', 'SRPK1', 'ARGLU1', 'RPL5', 'PCLAF',
    #    'FABP5', 'LRBA', 'FCER1A', 'RPS15A', 'ADGRG3', 'NCOA7', 'IGF2BP2',
    #    'HOXA9', 'IGLC2', 'TCEAL9', 'TUT4', 'HMGB1', 'CXCR2', 'RBMX',
    #    'PBXIP1', 'MMP25', 'CBX5', 'IGHD', 'HBD', 'HMGA1', 'CCDC34',
    #    'ATP8B4', 'SSBP2', 'AC092683.1', 'RPS12', 'TEC', 'ATP6V0A2',
    #    'VPS13A', 'KCNQ1OT1', 'IKZF2', 'AL360012.1', 'UPF2', 'LIMA1',
    #    'GCSAML', 'ZNF22', 'CEP126', 'NASP', 'AC243960.1', 'CASC15',
    #    'R3HDM4', 'PHACTR4', 'TPM1', 'ZNF439', 'PAICS', 'SELL', 'USP11',
    #    'RPL23', 'MZT2A', 'RPL34', 'PI3', 'RPS3', 'Z93241.1', 'CLEC11A',
    #    'OCIAD2', 'CKS2', 'C1QBP', 'RPS23', 'RPL32', 'PDZD8', 'MYCN',
    #    'KRT8',
]
PRESUMABLY_NEUTROPHIL_SPECIFIC_GENES = ['CAMP', 'LTF', 'LCN2']
RIBO_PROT_GENES_SILENT_IN_PRESUMABLY_NEUTROPHILS = [
    'RPL26', 'RPL30', 'RPL18A', 'RPL6', 'RPL37A', 'RPL18', 'RPS25',
       'RPL23A', 'RPL22', 'RPL4', 'RPSA', 'RPS10',
]

GENES_LOWER_IN_MONOCYTE_THAN_IN_BATCHY_MONOCYTES = [
    'CDK6', 'AC084033.3', 'STMN1', 'SPINK2', 'ANKRD28', 'SOX4', 'MMP9',
       'S100P', 'FAM30A', 'SNHG8', 'SMIM24', 'S100A12', 'CD177', 'NPR3',
       'MSI2', 'TFDP2', 'LTF', 'NUCB2', 'LCN2', 'MYB', 'PRSS57', 'HINT1',
       'C12orf57', 'PGLYRP1', 'SLC40A1', 'CPA3', 'FAM117A', 'TFPI',
       'LDHB', 'CYTL1', 'AC002454.1', 'CD34', 'FCGR3B', 'ZFAS1', 'BEX3',
       'PRKACB', 'PCDH9', 'CMTM2', 'BEX2', 'IGHM', 'HEMGN', 'LINC02573',
       'AC011139.1', 'AVP', 'NPM1', 'EGFL7', 'LAPTM4B', 'ANGPT1',
       'CYP4F3', 
       'SNHG19', 'SEPT6', 'TESPA1', 
       'TSTD1',
       'DPPA4', 'GUCY1A1', 
       'EBPL', 'AKR1C3', 'MSRB3', 'TSC22D1',
       'GABPB1-AS1', 
       'GATA2', 
       'MMP8', 'HIST1H1D',
       'MAP7', 'EEF1B2', 'CAMP', 'HSP90AB1', 'PBX1',
       'GYPC', 'SNHG7',
       'MLLT3', 'ITM2A', 'SOCS2', 'FHL1', 'ALOX5AP', 'KIT', 'CST7',
       
        'MCTP2', 'CRISP3', 'TCF4', 'NAMPT', 'NREP', 'ZNF385D',
       'NRIP1', 'ARMH1', 'PTMS', 'BASP1', 'CD82', 'SPN', 'SEPT11',
       'SNHG25', 'SLC2A3', 'CCND2', 'PROM1', 'KDM5B', 'SPTBN1', 
    #    'S100A8',
       'HSPB1', 'AL157895.1', 'CCDC171', 'PEBP1', 'ODF2L', 'IGLL1',
       'HMGN1', 
        'MTURN', 'FCER1A', 'HTR1F', 'SVIP', 'ST13',
       'HOPX', 'MDK', 'AC245060.5', 'ZNF521', 'HIST1H4C', 'ERG', 'IGSF10',
       'TSPAN2', 
        'MS4A1',
        'C1orf21', 'MME', 'ITM2C', 'ST8SIA6',
       'ETS1', 'MED12L', 'PSIP1', 'BCL11A', 'MZB1', 'PNN', 

       'CLC', 'TFRC', 'BCL7A', 'TXN', 'TRBC2', 'ADAM28', 'ANXA3',
       'ZDHHC19', 'HACD1', 'LRRC75A', 'GUCY1B1', 'C1QTNF4', 'ATP6V0A2',
       'SSBP2', 'MAP7D3', 'SCN9A', 'STRBP', 'GPSM2', 'PRDX2', 'ANKRD36C',
       'TCEAL4', 'ZMYND8', 'EREG', 'IL1R2', 'ALPL', 'RHEX', 

       
        'KCNQ1OT1', 'ZBTB20', 'SMIM27', 'LRBA', 
        'RALGPS2', 

       
        'IGF2BP2', 'PKIG', 'IGKC', 'HIST1H1E', 'MEIS1', 'SNRPD2',
       'KIFAP3', 
        'HMGA1', 'LSM7', 'HOXA9', 
        'ZNF22',
       'BAALC', 'CD109', 'NCOA7', 'RGL4', 'HBD', 'LIMA1', 

       'OCIAD2', 
        'IMPDH2', 'MZT2A', 'FOXP1', 'HLTF', 'SERPINB1',
       'CRHBP', 'PCLAF', 
        'TUT4', 'IKZF2', 
        'MGAM',
       'HMGB1', 
        'CNRIP1', 
        'CD24', 
        'CASC15', 
        'SLC38A1',
       'NASP', 
        'PBXIP1', 'ARGLU1', 'AC103591.3', 
        'ARG1',
       'RBMX', 'FABP5', 'VPS13A', 'PDZD8', 
        'ATF7IP2'
    
    #    'RPL37A',
    #    'RPL23A', 
    #    'RPS5', 
    #    'RPS21', 
    #    'RPS18', 
    #    'RPS29', 
    #    'RPL31',
    #    'RPSA',
    #    'RPS28',
    #    'RPS8',
    #    'RPL22L1',
    #    'RPL36',
    #    'RPL5',
    #    'RPS27',
    #    'RPS3',
    #    'RPL12',
    #    'RPL39',
    #    'RPS15A',
    #    'RPS12',
    #    'RPL37',
    #    'RPS6',
    #    'RPL38',
    #    'RPS23',
    #    'RPL10A',
]
GENES_HIGHER_IN_MONOCYTE_PROG_THAN_IN_MONOCYTES = [
    'CDK6',
 'AC084033.3',
 'STMN1',
 'SOX4',
 'SPINK2',
 'ANKRD28',
 'SNHG8',
 'SMIM24',
 'PRSS57',
 'MSI2',
 'MYB',
 'TFDP2',
 'FAM30A',
 'NUCB2',
 'LDHB',
 'NPR3',
 'CYTL1',
 'FAM117A',
 'CD34',
 'HINT1',
 'CPA3',
 'SLC40A1',
 'TFPI',
 'NPM1',
 'SEPT6',
 'BEX3',
 'C12orf57',
 'GATA2',
 'PRKACB',
 'TESPA1',
 'LRRC75A',
 'TSC22D1',
 'EBPL',
 'BEX2',
 'PBX1',
 'LAPTM4B',
 'PCDH9',
 'ITM2A',
 'HSP90AB1',
 'CCND2',
 'ZFAS1',
 'ANGPT1',
 'GUCY1A1',
 'HEMGN',
 'GABPB1-AS1',
 'AKR1C3',
 'FHL1',
 'MAP7',
 'AC002454.1',
 'EGFL7',
 'AC011139.1',
 'SNHG7',
 'MLLT3',
 'NREP',
 'RPS5',
 'FCER1A',
 'GYPC',
 'NRIP1',
 'HMGN1',
 'MSRB3',
 'SEPT11',
 'MTURN',
 'KDM5B',
 'KIT',
 'LINC02573',
 'DPPA4',
 'ZNF385D',
 'SPN',
 'RPL23A',
 'EEF1B2',
 'RPSA',
 'PEBP1',
 'HOPX',
 'AL157895.1',
 'RPS28',
 'AVP',
 'ARMH1',
 'SNHG19',
 'PNN',
 'SOCS2',
 'RPS29',
 'IGLL1',
 'TSTD1',
 'BCL11A',
 'HIST1H4C',
 'RPS18',
 'ST8SIA6',
 'MED12L',
 'CCDC171',
 'ZNF521',
 'TCF4',
 'KCNQ1OT1',
 'C1orf21',
 'HMGA1',
 'ITM2C',
 'RPS21',
 'ATP6V0A2',
 'OCIAD2',
 'PTMS',
 'PSIP1',
 'CCNB1IP1',
 'HTR1F',
 'ERG',
 'IMPDH2',
 'ZBTB20',
 'ODF2L',
 'RHEX',
 'IGSF10',
 'MYL6B',
 'ZNRF1',
 'ZMYND8',
 'HACD1',
 'HSPB1',
 'ST13',
 'HNRNPA1',
 'SPTBN1',
 'BCL7A',
 'HLTF',
 'LRBA',
]

N227_UPREGULATED_GENES = [
    'S100A4', # universal?

    'HBA1', 'HBG2', # specifically in MKP? (didn't check thoroughly other cell states)
]

# remember: the color doesn't really matter. it is for our general intuition only.
GENES_HIGHER_IN_LESS_DIFFERENTIATED_EP_OF_NIMROD = ['PBX1', 'BIN2', 'SOD2', 'VCL', 'HEMGN']
GENES_HIGHER_IN_MORE_DIFFERENTIATED_EP_OF_NIMROD = ['APOC1', 'BLVRB', 'MPC2']



NKTDP_L_SOMEWHAT_SPECIFIC_GENE_NAMES = [
    'ID2', 'JCHAIN', 'ACY3', 'RUNX2', 'MAF', 'DDIT4', 'LCP1', 'SAMHD1',
    'KLRB1', 'SPINK2', 'SATB1', 'MEF2A', 'RUNX3', 'SLC38A1', 'LSP1',
    'AL096865.1', 'TTN', 'NKG7', 'CLNK', 'LINC00865', 'LGALS1',
    'CORO1A', 'CYTH4', 'IGHM', 'MEF2C', 'GPR183', 'GSN', 'RAB11FIP1',
    'TYROBP', 'ITM2C', 'STK17A', 'ITGB7', 'SPECC1', 'AHR', 
    # 'SOX4',
    'CAPG', 'CCR7', 'TNFRSF18', 'CYSLTR1', 'MDFIC', 'NEGR1', 'TRIB2',
    'NFIL3', 'NYNRIN', 'PAK1', 'SCPEP1', 'HDAC9', 'IRF8', 'TRAF3IP3',
    'ATP10A', 'STK17B', 
    # 'GSTP1', 
    'ITGAL', 'CTSW', 'LTB', 
    # 'CALM2',
    # 'TMSB4X', 
    'TCF4', 'RASD1', 'ZNF217',
]

B_SOMEWHAT_SPECIFIC_GENE_NAMES = [
    # NOTE: would probably be better to use higher_in_b_than_all_other_hspcs, which seems to distinguish better from lymphoid progenitors, especially pro-B and pre-B
    'MS4A1', 'IGKC', 'CD79A', 'LTB', 'BANK1', 'IGHM', 'RALGPS2',
    'CD79B', 'KLF2', 'IKZF3', 
    # 'HLA-DQB1', 
    'BIRC3', 'CD48', 'MARCH1',
    # 'HLA-DQA1', 
    'BTG1', 
    # 'CD52', 
    'CD24', 'LBH', 'LINC00926', 'ARHGAP24',
    'CCDC50', 'SP140', 'FCRL2', 'SPIB', 'TNFRSF13C', 'LINC01857',
    'FAM107B', 
    # 'CD74', 
    'AIM2', 'NCF1', 
    # 'HLA-DMB', 'HLA-DRB1', 
    'ISG20',
    'IGHD', 'CD22', 'POU2F2', 'FAM129C', 
    # 'HLA-DPA1', 
    'FCRL5',
    # 'HLA-DPB1', 
    # 'RIPOR2', 
    # 'HLA-DRA', 
    # 'EVI2B', 
    'TMEM154', 'SCIMP',
    'JCHAIN', 'FCMR', 'CYBB', 'CPNE5', 'HVCN1', 'FCRL3', 
    ###'VPREB3', # higher in pre-B
    'EZR', 'ID3', 'FCRLA', 'ETS1', 'GPR183', 'PAX5', 'SH3BP5',
    'LINC02397', 'CLECL1', 'TNFRSF13B', 'TLR10', 'MARCKS', 'OSBPL10',
    'BLK', 'LY86', 'SWAP70', 
    # 'CD37', 
    'COBLL1', 'CTSH', 'SMIM14',
    'KIAA0040', 
    # 'KIAA1551', 
    # 'STX7', 
    'CXCR4', 
    ###'BLNK', # higher in PB (i.e., late) pro-B
    'POU2AF1',
    'MTSS1', 'IGHA1', 'IGLC2', 'AFF3', 'COTL1', 'ATM', 'P2RX5',
    'CLEC2D', 'IRF8', 
    ###'EBF1', # higher in pre-B and PB (i.e., late) pro-B
    'CD27', 'PARP15', 'IGLC3', 'TBC1D9',
    'CTSZ', 'FCGR2B', 'LINC01781', 
    # 'SP110', 
    'P2RY10', 
    # 'TMSB4X',
    'FCRL1', 'FGD2', 'CYB561A3', 'TTN', 'HERPUD1', 'NCOA3', 
    # 'HLA-DOB',
    'RAB30', 'GGA2', 'TRAF3IP3', 'CHCHD10', 'SEL1L3', 
    # 'SMAP2',
    'RNASET2', 
    # 'CALM1', 
    # 'CORO1A', 
    'CD180', 'ZFP36L1', 'BCL2', 
    # 'SP100',
    'PRKCB', 'CD40', 'CHD7', 
    'MNDA', 
    'CD19', 'CCR6', 'TPD52', 'TRIM38',
    # 'CYBA', 
    'DAPP1', 'BTG2',
]








STRONG_MONOCYTE_SPECIFIC_GENE_NAMES = [
    'S100A9', 'S100A8', 
    'FCN1', 
    'VCAN', 
]

STRONG_MONOCYTE_GENE_NAMES = [
    *STRONG_MONOCYTE_SPECIFIC_GENE_NAMES, 
    'MNDA', # around -14 in B cells
    'CYBB', # around -13.3 in B and NKTDP
    'FGL2', # around -11 in DC
]
MONOCYTE_GENE_NAMES = [
    *STRONG_MONOCYTE_GENE_NAMES,
    'S100A12', 
    'TNFRSF1B',
    'FGR',
    'FCGR2A',
    'NCF2',
    'SLC11A1',
    'CSTA', 'P2RY13',
    'CD14',
    'GIMAP4', 'CFP', 
    # 'CEBPD', # a bit high in GMP
    'MPEG1', 'MS4A6A', 'MS4A7',
    'CLEC7A', 'LRRK2',
    'KCTD12',
    'SLC7A7', 'LGALS3',
    'SERPINA1',
    'IGSF6',
    'CD68',
    'CD300E', 'RAB31', 'CD93',
    'FPR1', 'LILRB3', 'LILRB2', 'LILRB1', 'IL17RA', 'ADA2',
    'LGALS2',
    'RGS2', 'TGFBI',
    'MAFB',
]

HIGHER_IN_PRE_B_THAN_CLP_M = [
    # >3
    'IGLL1', 'DNTT', 'VPREB3', 'VPREB1', 'AKAP12', 'ARPP21', 'CD79B',
       'SOCS2', 'EBF1', 'RAG1', 'TLE4', 'CD24', 'MYB',
       # 2.5-3
       'UHRF1',
       'CD79A', 'CD9', 'ZCCHC7', 'HPS4', 'XBP1', 'IL7R', 'RAG2', 'BLNK',
       'SYNE3', 'LINC01013',
       # 2.2-2.5
       'CD38', 'CMTM8', 'NEIL1', 'TP53INP1', 'ID3',
       'CD81', 'CYGB', 'MT1X',
]
HIGHER_IN_PRE_B_THAN_CLP_M_WITHOU_DNTT = [x for x in HIGHER_IN_PRE_B_THAN_CLP_M if x != 'DNTT']

MOST_ANTI_CORR_WITH_DNTT_IN_HSC_TO_CLP_M = [
    'RPL12', 'AVP', 'RPL5', 'RPS6', 'RPL10A', 'RPS18', 'RPS3A', 'RPS8', 'RPL13', 'HSP90AB1', 'RPS4X', 'HINT1', 'RPL37A', 'RPL3', 'EEF1B2', 'RPL32', 'RPS23', 'RPS2', 'RPS14', 'RPS12', 'RPL7A', 'RPL15', 'RPS24',
]
FILTERED_RP_MOST_ANTI_CORR_WITH_DNTT_IN_HSC_TO_CLP_M = [
    x for x in MOST_ANTI_CORR_WITH_DNTT_IN_HSC_TO_CLP_M if x not in ('AVP', 'HINT1', 'RPS4X', 'EEF1B2', 'HSP90AB1')]
HIGHER_IN_BEMP_THAN_MEBEMP_L_GENES = [
    'HDC', # best?
    'ALOX5AP', # very nice
    'LMO4', # nice

    'FAM83F',
    'BACE2',
    'CNRIP1',
    'HPGDS',
    #  'CASP3',
    'SRGN',
    'CD63',
    'CNKSR3',
    'TPSB2',
    'TPSAB1',
    'IKZF2',
    'GPR183',
    'TENT5A',
    'KIT',
    'GRAP2',

    # a bit noisy inside BEMP?
    'MS4A2', 'MS4A3', 'CLC', 'HPGD',
]
NKT_SOMEWHAT_SPECIFIC_GENE_NAMES = [
    'IL32', 'NKG7', 'GNLY', 'GZMA', 'CCL5', 'CD2', 'GIMAP4', 'GIMAP7',
    'SYNE2', 'CD3G', 'CD3D', 'CST7', 'ARL4C', 'KLF2', 'CD3E', 'KLRD1',
    'SLFN5', 'IL7R', 'TRBC1', 'LINC00861', 
    # 'CALM1', 
    'PRF1', 'FYB1',
    'CD247', 'ETS1', 'GZMB', 'CD48', 'TRAC', 'KLRB1', 'GZMH', 'PYHIN1',
    'RORA', 
    # 'HCST', 
    'SAMD3', 'FGFBP2', 'FAM107B', 'GZMM', 'CCL4',
    'LCK', 'TRBC2', 'IKZF3', 'ID2', 
    # 'BTG1', 
    'BCL11B', 
    # 'S100A10',
    # 'PTPRC',
    'MYBL1', 'CX3CR1', 'LBH', 'SKAP1', 'CLEC2D', 'ITGB2',
    # 'MYL12A', 
    'SAMHD1', 
    # 'RARRES3', 
    # 'TMSB4X', 
    'GBP5', 'FYN', 'GPR65',
    'PPP2R5C', 'FCGR3A', 'IL10RA', 'AAK1', 
    # 'RNF213', 
    'RASGRP1', 'TRDC',
    'CD8A', 'KLRF1', 
    # 'ATM', 
    # 'TRAF3IP3', 
    'PRDM1', 
    # 'HLA-C', 
    # 'KIAA1551',
    'CD7', 'STK17A', 'KLRG1', 'IL2RG', 
    # 'HLA-F',
    'GIMAP1', 'DENND2D',
    
    'CD96', 'C12orf75', 
    # 'CORO1A', 
    'LAT', 'LITAF', 
    # 'CTSW', 
    'EVL',

    'SYNE1', 
    # 'S100A4', 
    # 'S100A11', 
    'SYTL3', 'ITGAL', 
    # 'HLA-B', 
    'SH2D1A',
    'EFHD2', 'SLFN12L', 'IL2RB', 'PTPN4', 'SPOCK2', 'GNG2', 'TRG-AS1',
    'CD226', 
    # 'B2M',
    # 'TXNIP', 
    'SPON2', 'DOK2', 'LYAR', 'GZMK', 'STK17B',
    'PTGDR', 'PRKCH', 'ISG20', 
]

HIGHER_IN_CLP_THAN_PRO_B = [
    'PRSS2', 
    ###'SPINK2', 
    'CD164', 'LMO4', 'HOPX', 'KIAA0087', 
    'IRF9',
    'HOXA9',

    #     'RPS6KB1', 'ARID5B', 'PYGL', 'AC090673.1', 'ZNF503',
    #    'GCSAML', 'RXFP1', 'PTMS', 'TFPI', 'MYCN', 'ERMP1', 'LYN', 'PDIA4',

    #    'NPR3', 'MAF', 'HOXA6', 'IFNG-AS1', 'ARHGAP22', 'AGPAT5', 'CLSTN3',
    #    'CNN2', 'CHMP1B', 'ZNF514',
]

DISEASE_PB_CELL_SCORE_DESC_TO_GENES = {
    'higher_in_n257_c_unassigned_mebemps': [
        # diff exp when comparing c_state=state_unassigned & state.isin(MEBEMP-E, MEBEMP-L), N257 vs rest.
        # log_ratio>2, rest_expr<-13
        'CDKN1C', 'CD84', 'MAN2A2', 'TIMP3', 
        ###'AL034397.3',
        'CAVIN1', 'KCNQ1OT1', 'CD9', 'MDK', 'NPR3', 'TP53INP1',
        'MIR181A1HG', 'CSF2RB', 'CD7', 'ZMYND8', 'ANGPT1', 'ABHD2',
        'ZBTB20', 'FCER1A', 'PTMS', 
        ###'POLR2J3---1', 
        'TAX1BP1', 'ORAI2',
        'FES', 'CITED2', 'GLUL', 
        'RGS18', 
        'GCC2', 'TLE4', 'FADS1', 'CD47',
        'PDZD8', 
        ###'AC008105.3', 
        ###'AHNAK', 
        'PLXNC1', 'EIF4G2', 'NFIX',
        'SLC39A10',
    ],
}

PB_SIG_NAME_TO_GENES = {
    'higher_in_clp_m_and_pro_b_than_mpp': [
        # >2
        'LTB', 'MME', 'JCHAIN', 'DNTT', 'IGHM', 'SCN3A', 'TCF4', 'MZB1',
       'MEF2C', 'LST1', 'DDIT4', 'MEF2A', 'AFF1', 'VPREB1', 'LINC00865',
       'AC245060.5', 'BCL11A', 'CYFIP2', 'ADA',
       # 1.5-2
       'CD22', 
       ###'TMSB4X', 
       'DCK',
       'CORO1A', 'ZFP36L1', 'EVL', 'TRBC2', 'PRKCB', 'RASD1', 'NEGR1',
       'RNASEH2B', 'STK17B', 'LAIR1', 'PAG1', 'CD53', 'MED13L', 'RNF213',
       'SLC38A1', 'KLF6', 'IRF1', 'PPP1R18', 'LAPTM5', 'SEPT9', 'SLC2A5',
       'AC243960.1', 'ETS1', 'CD74', 'IQGAP1', 'IER2', 'ARL4C', 'CALM1',
    ],
    # 'higher_in_dc_than_monocyte': [ # 240619: currently unused
    #     # >3
    #     'TCF4', 'IRF8', 'C12orf75', 'SPIB', 'UGCG', 'PLD4', 'CCDC50',
    #    'LILRA4', 'ITM2C', 
    #    'IGKC', 
    #    'JCHAIN', 
    #    'PPP1R14A', 'PPP1R14B',
    #    'PLAC8', 'BCL11A', 
    #    ###'STMN1', 
    #    'MZB1', 'SCN9A', 'SCT', 'CLIC3',
    #    'ALOX5AP', 'CLEC4C', 'IL3RA', 'ARL4C', 'APP', 'IRF7', 'SLAMF7',
    # ],
    'higher_in_dc_than_nkt': [
        # >2
        'CST3', 'HLA-DRA', 'IRF8', 'TCF4', 'CD74', 'PLD4', 'HLA-DPA1',
       'CCDC50', 'BCL11A', 'CCDC88A', 'LYZ', 'APP', 'SPIB', 'HLA-DPB1',
       'FCER1G', 'CTSZ', 'HLA-DRB1', 'LILRA4', 'RNASE6', 'GAPT',
       'HLA-DMA', 'MEF2C', 'UGCG', 'MPEG1', 'GRN', 'HLA-DMB', 'PTPRE',
       'IGKC', 'PPP1R14A', 'FCGRT', 'JCHAIN', 'TGFBI', 'ITM2C', 'CYBB',
       'IL3RA', 'PPP1R14B', 'IRF7', 'SCT', 'SCN9A', 'TYROBP', 'CLEC4C',
       'HLA-DQA1', 'MZB1', 'CIITA', 'JAML', 'ALDH2', 'SOX4', 'SULF2',
       'SPI1', 'PLAC8', 'HERPUD1', 'LILRB4', 'RNF130', 'FGD2', 'AFF3',
       'BLNK', 'SAMHD1', 'KCTD12', 'DAB2', 'LYN', 'CAPG', 'SERPINF1',
       'VIM', 'IGSF6', 'SHTN1', 'BASP1', 'LST1', 'SYK', 'FGL2',
       'CYB561A3', 'DPYSL2', 'DNASE1L3', 'MYO1E', 'TNNI2', 'PHACTR1',
       'CTSH', 'PYCARD', 'HLA-DQB1', 'NCF1', 'H2AFY', 'MARCH1', 'NPC2',
       'GSN', 'AIF1', 'SEC61B', 'PLP2', 'HLA-DRB5', 'CBFA2T3', 'C12orf45',
       'RUNX2', 'SPINT2', 'FAM129C', 'UNC93B1', 'FTH1', 'TSPAN13',
       'WDFY4', 'CLEC12A', 'CSF2RA', 'FCER1A', 'FLT3', 'LGMN', 'FABP5',
       'CNPY3', 'GSTP1', 'STX7', 'TXN', 'APEX1', 'TYMP', 'MGST2', 'ALCAM',
       'ANXA2', 'P2RY14', 'CD68', 'GAS6', 'BTK', 'AP1S2', 'CAT', 'AXL',
       'TBC1D9', 'VAMP8', 'TNFSF13B', 'GNG7', 'RNASET2', 'BCL7A',
       'SNRNP25', 'NREP', 'CLEC4A', 'CRYBG3', 'TNFAIP2', 'TPM2', 'CALM2',
       'IDH3A', 'CORO1C', 'LGALS1', 'KLF4', 'AC119428.2', 'CTSB', 'ZEB2',
       'NAGK', 'TAGLN2', 'OTULINL', 'UPK3A', 'RAB11FIP1', 'RHEX', 'MGLL',
    ],
    'higher_in_pre_b_than_pro_b': [
        # >2
        'AKAP12', 'ARPP21', 'CD9', 'VPREB3', 'DEPP1', 'CD37', 'LINC01013',
        'HMHB1', 'SMAD1', 'CMTM8', 'PTPRE', 'LHPP', 'CXCR4', 'LEF1',
        'CD79B', 'FOXO1', 'BACH2', 'LIG4', 'CD24', 'MLXIP', 'ELK3', 'YBX3',
        'TOP2B', 'AP001059.2', 'LIMD2', 'GAB1', 'TCF3', 'TNK2', 'PHACTR1',
        'MPP1', 'P4HA2', 'CDK13', 
        ###'SOX4', 
        'SYNE3', 'SNX2', 
        'HIST1H1C',
        'MND1', 'CLEC2D', 'BTG2', 
        ###'TRGV9', 
        'REL', 'RAG1', 
        ###'RPS4Y2',
        'ALDH5A1', 'TERF2', 'MT1X', 'CTBP2', 'PSD3', 'PXDN', 'TMEM243',
        'UTRN', 'MDM2', 'VPREB1', 'LINC02227', 'LRIG1', 'LMO4', 'HPS4',
        'TAPT1', 'DGKD', 'JMJD1C',
    ],
    # 'higher_in_pro_b_than_clp': [ # 240619: currently unused
    #     # copied higher_in_pre_b_than_all_hspcs, filtered it, and added some.
    #     'ARPP21', 'RAG1', 'SOCS2', 'RAG2', 
    #     ###'VPREB1', 
    #     'CMTM8', 'HMHB1',
    #     'UHRF1', 'HPS4', 'LCN6', 'TLE4', 'IGLL1', 'AKAP12', 'EBF1', 'CYGB',
    #     ###'TOP2B', 
    #     'PAG1', 'CMTM7', 
    #     'HHIP', 'RNF152',
    #     ###'PPP1R14B', 
    #     'SLC8A1-AS1', 'AGPS', 'PTPRE', 'MYLK', 'CTGF',
    #     'BAHCC1', 
    #     ###'RCSD1', 
    #     'LINC00426', 'DYRK3', 'SH2D4B', 'SMARCD2',
    #     'ZCCHC7', 'NRXN2', 'MT1X', 'P4HA2', 
    #     ###'AL713998.1', 
    #     'IRX1',

    #     # added according to diff exp between pro-B and CLP and filtering
    #     'IL7R',
    #     ###'JCHAIN',
    #     'CD79A',
    #     'BLNK',
    #     'VPREB3',
    #     'CD79B',
    # ],
    'higher_in_pro_b_than_clp_m': [
        # >2
        'IGLL1', 'IL7R', 'DNTT', 'UHRF1', 'SOCS2', 'VPREB3', 'VPREB1',
       'EBF1', 'ZCCHC7', 'TLE4', 'BLNK', 'CD79A', 'CD79B', 'AC002454.1',
       'CYGB', 'AGPS', 'FBXW7', 'MYC', 'CD81', 'MIR181A1HG', 'MYB',
       'RAG1', 'BAHCC1', 'TP53INP1', 'RAG2', 'XBP1', 'HPS4', 'GAPDH',
       'FAM129C', 'CALCOCO2', 'GLRX', 'H2AFY', 'AKAP12', 'SYNGR1',
       'SLC8A1-AS1', 'VIM', 'IGFBP7', 'CD24', 'TOP2B', 'CD38', 'ARPP21',
       'PAG1', 'RGL4', 'RCSD1', 'CTGF',
    ],
    # 'higher_in_pre_b_and_b_than_clp_m': [ # 240619: currently unused
    #     # >1.5
    #     'CD79B', 'CD79A', 'CD24', 'PAX5', 'VPREB3', 'ISG20', 'FAM107B',
    #    'ID3', 'POU2AF1', 'CD72', 'TNFRSF14', 'KIAA0040', 'EZR', 'CXCR4',
    #    'EBF1', 'SNX2', 'RUBCNL', 'UCP2', 'CD27', 'MARCKS', 'FAM129C',
    #    'P2RX5', 'QRSL1', 
    #    ###'HLA-DRA', 
    #    'MACROD2', 'BLNK', 'BACH2', 'GNG7',

    #    # 1.4-1.5
    #    'PXK', 'BTG1', 'CLEC2D', 'FAM3C', 'HLA-DQB1', 'SMAP2', 'BANK1',

    #    # 1.3-1.4
    #    'SMIM14', 'DGKD', 'CHCHD10', 
    #    ###'ZEB2',
       
    #    # 1.2-1.3
    #    'ODC1', 
    #    ###'CD74', 
    #    'CSGALNACT1',
    #    'REEP5', 'RASGRP2',
    # ],
    # 'higher_in_pre_b_and_b_than_pro_b': [ # 240619: currently unused
    #     # >1.5
    #     'CD37', 'LIMD2', 'CLEC2D', 'REL', 'HLA-DMB', 'AFF3', 'ISG20',
    #    'DGKD', 'TNFRSF14', 'TBC1D10C', 'HLA-DQB1', 'HLA-DMA', 
    #    ###'CD74',
    #    'BTG1', 'CD27', 'RUBCNL', 'HLA-A', 'FCMR', 'PNPLA8', 'IGHD',
    #    'CD72', 'TRIM38', 'BTG2', 'P2RX5', 'PAX5',

    #    # 1.3-1.5
    #    'GSAP', 'SESN3',
    #    ###'HLA-DRA', 
    #    'SNX2', 'FOXP1', 'ARID5B', 'CXCR4', 'STIM2', 'MACROD2',
    #    'ATM', 'KIAA0040', 'BTN3A1', 'UCP2', 'RHOQ', 'NLRC5',

    #    # 1.2-1.3
    #    'CD79B',
    #    'ABLIM1', 'BANK1', 'HIVEP1', 'EZR', 'LMO4', 'STX12', 'SMAP2',
    # ],
    'higher_in_b_than_pre_b': [
        # >3
        ###'MS4A1', 
        ###'IGKC', 
        'IKZF3', 'RALGPS2', 'KLF2', 'BIRC3', 'MARCH1',
        'IGLC2', 'FCRL5', 'CD48', 'CD52', 'SCIMP', 'IGLC3', 'SPIB',
        'ARHGAP24', 'BANK1', 'LINC01857', 'POU2F2', 'EMP3', 'TMEM154',
        'LINC00926', 'LINC02397', 'SP140', 'LYN', 'TNFRSF13B', 'FCRL3',
        'MTSS1', 'CYBB', 'TNFRSF13C', 'PTPRC', 'CTSZ',

        # 2-3
        'HVCN1', 'FCRL2',
        'COTL1', 'LSP1', 'MAP3K1', 'SEL1L3', 'CTSH', 'HLA-DQA1', 'CLECL1',
        'AIM2', 'FCRLA', 'AHNAK', 'GPR183', 'LGALS1', 'TTN', 'PLEK',
        'SH3BP5', 'IGHA1', 'TLR10', 'TBC1D9', 'ADAM28', 'CPNE5', 'CD37',
        'DAPP1', 'LBH', 'OSBPL10', 'NCF1', 'SWAP70', 'FGD2', 'HLA-DRB1',
        'TNFAIP8', 'CCDC50', 'CD40', 'MPEG1', 'RIPOR2', 'BCL2', 'BLK',
        'ATP2B1', 'FCGR2B', 'COBLL1', 'CNN2', 'FCER2', 'HERPUD1', 'PTPN6',
        'IGHM', 'SETBP1', 'MNDA', 'IGHG3', 'CTSS', 'CIB1', 'HLA-DOB',
        'EVI2B', 'SYNGR2', 'CD84', 'CD180', 'S100A6', 'P2RY10', 'PARP15',
        'SLAMF6', 'CYB561A3', 'FTL', 'LRRK2', 'RASGRP3', 'CCR6', 'LTB',
        'KIAA1551', 'RAB30', 'ETS1', 'LINC01781', 'FCRL1', 'TFEC',
    ],
    'higher_in_nktdp_than_clp_m_and_pro_b': [
        # >2
        'SAMHD1', 'LGALS1', 'ID2', 'ACY3', 'NFKB1', 'TTN', 'ITGB7',
       'NFIL3', 'GPR183', 'FYB1', 'CYBB', 'GNLY', 
       ###'CALM2', 
       'COTL1',
       'RUNX3', 'SCPEP1', 'PLD4', 'RAB11FIP1', 'SYT2', 'SEMA4D', 'AHR',
       'SEL1L3', 'KIT', 'ATP10A', 'RIPOR2',
    ],
    # 'clp_abnormal_sig': [ # 240619: currently unused
    #     # NOTE: this seems like CLP-E... leave it for now...

    #     # iter3: filter iter2 by heatmap across CLPs (MCs) and MPPs similar to some CLPs in mds ult.
    #     'CRYGD','PCBD1','SERPINE2','EREG','SKAP1','ZNF385D','MYCT1','MYC','PDLIM1','TNFRSF14',
    #     ###'VAMP8',
    #     'CD84','CD9','TNFSF13B','MT2A',
    #     ###'HSPD1',
    #     ###'HSP90AB1',
    #     ###'RSL1D1',
    #     'TUBA1B',
    #     ###'PLEK',
    #     ###'CKS2',
    #     'S100A10',
    #     ###'AC084033.3',
    #     ###'AC011139.1',
    #     'HSPB1','RHEX','CAVIN1','BEX1',
    #     ###'AC026369.3',
    #     ###'KIT',
    #     ###'MDK',
    #     'ZBTB16','MEG3',
    #     ###'SERPINB1',
    #     'ANGPT1','AVP','IL1B','PPA1',
    #     'CSF3R','MGST1','MYB','C1QTNF4','PRSS57',
        
    #     # # iter2: corr>0.5 (with ['CSF3R', 'MGST1', 'IL1B', 'AVP']) in mds illu:
    #     # 'MGST1', 'CSF3R', 'AVP', 'IL1B', 'PPA1', 'RHEX', 'C1QTNF4',
    #     # 'RAB32', 'PRSS57', 'MYB', 'UCHL1', 'CAVIN1', 'HSPB1', 'SKAP1',
    #     # 'IMPDH2', 'PRSS21', 'PTMS', 'MYC', 'CCDC42', 'BEX1', 'EEF1B2',
    #     # 'MDK', 'AC026369.3', 'MEG3', 'MT2A', 'TMEM98', 'MYCT1', 'GATA2',
    #     # 'TNFRSF14', 'CD9', 'CRHBP', 'CD84', 'RAB7B', 'NOG', 'MT1X',
    #     # 'PDLIM1', 'ZBTB16', 'ZNF385D', 'CFD', 'KIT', 'TM4SF1', 'NPW',
    #     # 'NFE2', 'CPXM1', 'TNFSF13B', 'ANPEP', 'MSRB3', 'ANGPT1',
    #     # 'AL590226.1', 'IFITM3', 'SERPINB1', 'PCBD1', 'FSCN1', 'TLE4',
    #     # 'NPDC1', 'SORL1', 'TALDO1', 'RPL10A', 'EREG', 'ADGRG1', 'HSP90AB1',
    #     # 'AC011139.1', 'MT1F', 'NT5M', 'RPL36A', 'C4orf48', 'ERG', 'CRYGD',
    #     # 'CPVL', 'YBX3', 'HSPD1', 'HINT1', 'SNHG8', 'GNPDA1', 'PLEK',
    #     # 'PYGL', 'CD151', 'HOXB5', 'ACSM3', 'PLAC8', 'CPA3', 'CKS2', 'LAG3',
    #     # 'S100A10', 'RSL1D1', 'TCEAL4', 'CENPV', 'STON2', 'TUBA1B',
    #     # 'ATP1B1', 'TUBA1C', 'CEBPB', 'CHST12', 'TFEC', 'FAH', 'APLP2',
    #     # 'VAMP8', 'S100A6', 'ITGA6', 'KRT8', 'ANXA2', 'SELENOP',
    #     # 'PRKAG2-AS1', 'CREG1', 'SYNGR1', 'EBPL', 'GSTM3', 'GUCY1B1',
    #     # 'AC084033.3', 'EEF2', 'IRAK3', 'SMIM10', 'TCTEX1D1', 'FXN', 'CDK6',
    #     # 'LYL1', 'ETV6', 'AJ009632.2', 'TAL1', 'SOCS3', 'RPL4', 'KRT18',
    #     # 'SERPINE2', 'LINC00891', 'RAB13', 'NENF', 'DRAM1',

    
    #     # 'CSF3R', 'MGST1', 'IL1B', 'AVP',
    # ],
    'prss2_sig': [
        # iter4: filtered iter3 according to heatmap across CLP MCs.
        'PRSS2',
        ###'CD34', # nicely correlated with the rest, but when the module is high, CD34 can increase without it...
        'RXFP1','ADGRG1','ARID5B','TCTEX1D1','PYGL','GCSAML',
        ###'NCOA7', # somewhat strong
        'HEMGN',
        'SLC2A5', # comment out?
        ###'LTB', # comment out?
        'CYTL1',
        'CD52', # comment out? strong.
        ###'SNHG7', # comment out? strong. 240424: not good enough given its strength relative to others, i think. commenting out a bit late.
        'LIMS1','EGFL7','LYL1','RUNX1','CALCRL','TFPI',
        ###'PTMS', # looked somewhat bad in marker heatmap

        # # iter3: filtered iter2 according to heatmap across CLP MCs.
        # 'CD82','S100Z','IL3RA','AKR1C3','FHL3','GNAI1','LY6E','EFNA1','BIN1','SH3TC1','ACOT11','DAPK1','FAM69B','ERGIC1','REL','PRKCB','QPRT','RCSD1','HMGA1','DNAJC1','SLC2A5','LTB','PADI4','STK32B','LRRC26','TRBVB','VCL','LILRB2','COBLL1','PACSIN1','NUP214','CST3','MAP1A','GYPC','MPG','HCST',
        # ###'FAM30A',
        # ###'HOPX',
        # 'PRSS2', # maybe should stay
        # ###'CD99',
        # 'RXFP1','ADGRG1','ARID5B','TCTEX1D1','PYGL','GCSAML','NCOA7','HEMGN','CYTL1','CD52','SNHG7','LIMS1','EGFL7','LYL1','RUNX1','CALCRL','TFPI','PTMS',
        
        # # iter2: corr>0.4 with first iter sum
        # 'CD52', 'EGFL7', 
        # 'PRSS2', 
        # 'CYTL1', 'SNHG7', 'LIMS1', 'PTMS',
        # 'CD34', 'HEMGN', 'NCOA7', 'TFPI', 'RUNX1', 'HCST', 'LYL1', 'MPG',
        # 'ZFAS1', 'RGS10', 'SLC2A5', 'LTB', 'GCSAML', 'INKA1', 'PYGL',
        # 'RXFP1', 'COBLL1', 'MAP1A', 'NUP214', 'LILRB2', 'CALCRL', 'GIMAP7',
        # 'MTURN', 'ARID5B', 'ADGRG1', 'HMGA1', 'TRBVB', 'CST3', 'PACSIN1',
        # 'HLA-DPB1', 'TIMP3', 'TCTEX1D1', 'MZB1', 'NFE2', 'SVIL', 'XIRP2',
        # 'GYPC', 'SPTBN1', 'GPR146', 'AJ009632.2', 'LMO2', 'RPL31',
        # 'CHST12', 'ANP32B', 'RBPMS', 'PADI4', 'TALDO1', 'LY6E', 'ERGIC1',
        # 'RPLP0', 'IL3RA', 'CD82', 'LYN', 'DNAJC1', 'AL356275.1', 'CD99',
        # 'TKT', 'CD22', 'XXYLT1-AS2', 'IDS', 'ACOT11', 'S100Z', 'NRIP1',
        # 'VCL', 'NPR3', 'HLF', 'LINC02573', 'EIF3E', 'ZBTB8A', 'FXYD5',
        # 'FHL3', 'ITM2A', 'RAB13', 'S100A11', 'SH3TC1', 'HOPX', 'SERPINB6',
        # 'SYTL4', 'STK32B', 'NHSL2', 'CDH26', 'LAPTM5', 'GNAI1', 'PRKCB',
        # 'AL590226.1', 'HLA-DPA1', 'IGF2BP2', 'BIN1', 'AKR1C3', 'LRRC26',
        # 'CNST', 'NPDC1', 'GIMAP1', 'PBXIP1', 'SETBP1', 'LAPTM4B', 'FAM30A',
        # 'QPRT', 'GATA2', 'SMC4', 'STARD9', 'CRHBP', 'FAM69B', 'MSI2',
        # 'DAPK1', 'TAGLN2', 'RCSD1', 'ARHGAP22', 'REL', 'FBXW9', 'ABCA13',
        # 'AC016735.1', 'EFNA1',
        
        # # actually started by looking for genes correlated with HOPX (corr>0.3), then plotted in a heatmap across CLP, and found that HOPX is relatively lonely, so switched to focus on the next most famous gene, RUNX1. took from the heatmap the genes most correlated with RUNX1. similarity was mostly in MCs low in RUNX1.
        # 'RXFP1','ADGRG1','ARID5B','TCTEX1D1','PYGL','GCSAML','NCOA7','HEMGN','RGS10','LYL1','CD52','EGFL7','HCST','LIMS1','MPG','CYTL1','CALCRL','RUNX1','TFPI','PTMS',
    ],
    # 'higher_in_nktdp_than_bemp': [ # 240619: currently unused
    #     # >3
    #     'ACY3', 'SPINK2', 'CD74', 'ID2', 'SAMHD1', 'LSP1', 'HLA-DRA',
    #    'JCHAIN', 'MEF2C', 'HLA-DPB1', 'CORO1A', 'RUNX3', 'HLA-DPA1',
    #    'HLA-DRB1', 'CYTH4', 'RUNX2', 'DDIT4', 'IGHM', 'SLC38A1', 'NFIL3',
    #    'GSTP1', 'GSN', 'RIPOR2', 'NFKB1', 'SCPEP1', 'TTN', 'MEF2A',
    #    'PLD4', 'NAP1L1', 'ITM2C', 'AHR', 'RAB11FIP1', 'TNFRSF18', 'MAF',
    #    'BLNK', 'CALM2', 'LTB', 'MDFIC', 'ATP10A', 'AFF3', 'GLIPR1',
    #    'PLEKHO1',
    # ],
    # 'higher_in_nktdp_than_pre_b_and_b': [ # 240619: currently unused
    #     # >2
    #     'ID2', 'ACY3', 'SAMHD1', 'SPINK2', 'RUNX2', 'CTSW', 'NFIL3',
    #    'TYROBP', 'CALM2', 'FYB1', 'RUNX3', 
    #    ###'LGALS1', 
    #    'AL096865.1',
    #    'TNFRSF18', 'NKG7', 'MAF', 'AIF1', 'NYNRIN', 'AHR', 'DDIT4',
    #    'GNLY', 'GSTP1', 'TESC', 'LSP1', 'HCST', 'SPECC1', 'NFKB1', 'IL18',
    #    'CD7', 'TNFRSF4', 'SEMA4D', 'ATP10A', 'GSN', 'IRF2BP2', 'SYT2',
    #    'FOS', 
    #    ###'JCHAIN', 
    #    'LCP1', 'CAPG', 'KIT', 'MLC1', 'PAK1', 'CYTH4',
    #    'SPON2', 'ARMH1', 'FXYD5', 'CD3E', 'ITGB7',
    # ],
    'higher_in_nktdp_than_all_myeloid_hspcs': [
        # >2
        'SAMHD1', 
        ###'LGALS1', # high in BEMPs, especially in later ones
        'ID2', 'ACY3', 'NFKB1', 'TTN', 'ITGB7',
       'NFIL3', 'TESC', 'SCPEP1', 'CYBB', 
       ###'GPR183', # high in late BEMPs
       'GNLY', 'FYB1',
       'CALM2',
       'COTL1', 'RUNX3', 'PLD4', 'IGFBP7', 'RAB11FIP1', 'SYT2',
       'MS4A1', 'SEMA4D', 'AHR', 'SEL1L3', 'RNASE6', 
       ###'KIT', # higher in BEMP than in NKTDP
       'ATP10A',
       'TLR10', 'RIPOR2', 'LY86',
    ],
    'higher_in_mebemp_l_than_mpp': [
        # 240107: mc_utils.get_diff_exp_df(n_mc_ad, n_mc_ad.obs['state'] == 'MEBEMP-L', n_mc_ad.obs['state'] == 'MPP')
        *BEST_HIGHER_IN_MEBEMP_L_THAN_MPP_GENES, 

        # log_ratio 0.5-1.5 and MPP expression <-15, filtered manually according to marker heatmap
        'STXBP5', 'TFR2', 'PLXNC1', 'GDF11', 'ITGA2B', 'HPGDS', 'DNAJC9',
        'LDB1', 'TAL1', 'STAM', 'NLK', 'PLIN2', 'DEPTOR', 'MIR4435-2HG',
        'SEMA3C', 'DLC1', 'FYB1', 'STXBP6', 'C2orf88', 'RGS18', 'CD38',
        'ZFP36L1', 'PPP3R1', 'ACSM3', 'HBS1L', 
        ###'PCLAF', 
        'LXN', 'TPSB2',
        'TIMP3', 'MTHFD2', 'BLVRB', 'APOC1', 'KIFAP3', 'ARMC8', 'PDLIM1',
        'FADS1', 'TIMP1', 'TRAPPC6B', 'MYH10', 
        ###'CENPU', 
        'KCNH2', 'PIM1',
        'ABCC5', 'YARS', 'ANKRD27', 'SMIM1', 'MAN2A2', 'ICAM4', 'P2RX1',
        'FAM45A', 'VKORC1L1', 'SLC27A2', 'EPOR', 'GMPR', 'FADS2', 'MPST',
        'TPSAB1', 'TRIB2', 'FREM1', 'CBFA2T3', 'CCDC28B', 
        ###'MS4A2', 
        'CD7',
        'STON2', 'ERBIN', 'MICAL2', 'TST', 'NECAB1', 'COTL1', 'RFX7',
        'EMID1', 'UROS', 'BOLA3-AS1', 'AGPAT5', 'ZBTB16', 'RYR3', 'ITGAE',
        'NMNAT3', 'ARHGAP22', 'TUBA1C', 'RAP2B', 'PABPC4', 'FLNA',
        'ZDHHC2', 'LYAR', 'NBAS', 'FES', 'MBOAT2', 'ZFPM1', 'ARHGEF12',
        'PLCXD1', 'RIPOR3', 'ANK1', 'DHRS11', 'ASRGL1', 'MPIG6B', 'MEX3B',
        'CMAS', 'RUVBL1', 'WEE1', 'NUP98', 'CD36', 'P2RX5', 'GRK4',
        'SLC39A4', 'NOLC1', 'AOAH', 'ODC1', 'GPATCH4', 'ASAP1', 'SLBP',
        ###'IVNS1ABP', 
        'MRPL14', 'PIK3CB', 'MRPS26', 'SLC14A1', 'VASP', 'FXN',
    ],
    'higher_in_mkp_than_mebemp_l': [
        # >2, filtered manually according to marker heatmap
        'MT2A', 'IFITM3', 'LTBP1', 'RAB27B', 'LAT', 'ITGA2B',
        ###'PFN1', 
        'CD36', 
        # 'ACTB', 
        'ACTN1', 'EIF5A', 'PPP1R14A', 'PDLIM1', 'PPIB', 'PDIA6', 'SNRPG', 'CD63', 'HSPA5', 'POLR2L', 'CYCS', 'HMGB2', 'MANF', 'HSP90B1', 'LY6E', 'TUBB4B', 
        ###'DAD1', 
        'CMTM5', 'ANAPC11', 'GP9', 'MIF', 'TUBB', 
        ###'TPM1', 
        'MT1E', 'PRDX1', 'TPI1', 'PPP1R14B', 'CALM1', 'C1QBP', 'LSM4', 'NAA38', 'SLIRP', 
        ###'PRKAR2B', 
        ###'HNRNPAB', 
        'CAMTA1', 'CENPW', 'PDIA3', 'NDUFAF8', 'WDR1', 'UQCRQ', 'SMC2', 'ARPC4', 'CAMK1', 'FCGR2A', 'SNRPB', 'PA2G4', 'ATAD2', 'COX17', 'PRDX4',
    ],
    'higher_in_mkp_than_mebemp_l_and_eryp': [
        # 240622: ah? why isn't PF4 here??
        ###'H2AFZ', 
        ###'TUBA1B', 
        'PKM', 
        ###'HMGB2', 
        'HSPA5', 'ITGA2B', 
        ###'RAP1B',
       ###'HSP90B1', 
       'PDLIM1', 
       ###'DAD1', 
       'FERMT3', 'TUBB4B', 'MEF2C', 'LY6E',
       'MKI67', 'LTBP1', 'ANAPC11', 'NME1', 'CENPF', 'CACYBP', 
       ###'CD63',
       'NRGN', 'TYMS', 'PPIF', 
       ###'EIF5A', 
       'SNRPG', 
       ###'PLEK', 
       'UBE2C', 
       ###'PRDX1',
       'TARS', 
       ###'TMSB4X', 
       'DTYMK', 'CAVIN2', 'PCLAF', 'MMRN1', 'CKS1B',
       'PRDX4', 'TPM1', 'ASPM', 
       ###'PDIA3', 
       'CENPX', 'RPA3', 'POLR2L',
       ###'HMGN2', 
       'CCND3', 'CDKN3', 'KIF11', 'PRKAR2B', 
       ###'PA2G4', 
       ###'SNRPB',
       'PCNA', 
       ###'ACTB',
        'FKBP3', 'CYCS', 
        ###'HSPE1', 
        'SMC2', 'TGFB1',
       'MPIG6B', 'LDHA', 
       ###'ERH', 
       'SGO1', 'MANF', 'MT2A', 'TTC27', 'CHEK1',
       'TPM4', 'PDIA6', 'MINOS1', 
       ###'TPI1', 
       ###'DUT', 
       'SDF2L1', 
       ###'PFN1',
       'EIF4EBP1', 'KIF2A', 
       ###'SNRPD1', 
       ###'PSMA7', 
       'ATP6V0B', 'TMEM230',
       'TUBA1C', 'CCNB2', 
       ###'RAN', 
       'POLR2J', 'TOP2A', 'KIF15', 'BIRC5',
       'PTTG1', 
       ###'HNRNPAB', 
       ###'DEK', 
       'DCXR',
    ],
    # 'higher_in_mkp_than_ep': [ # 240619: currently unused
    #     ###'PLEK', 
    #     ###'TMSB4X', 
    #     'RAB27B', 'LTBP1', 'CD9', 'TPM4', 
    #     ###'PPIF', 
    #     'LAT',
    #     'PKM', 'GSN', 'FERMT3', 'MMRN1', 'SLC44A1', 
    #     ###'IFITM3', 
    #     ###'PTTG1',
    #     'STOM', 
    #     ###'MT2A', 
    #     ###'UBE2C', 
    #     ###'ACTN1', 
    #     'PDLIM1', 
    #     ###'CD63', 
    #     ###'BIRC5',
    #     'CAP1', 'ZYX', 
    #     ###'DAD1', 
    #     ###'TYMS', 
    #     ###'MKI67', 
    #     'GP9', 'VCL', 'MEF2C',
    #     ###'ACTB',
    #     'CMTM6', 'RGS18', 'CMTM5', 
    #     ###'ANAPC11', 
    #     ###'H2AFZ', 
    #     'CAVIN2',
    #     ###'PFN1',
    #     ###'ITGA2B',
    #     ###'CDKN3', 
    #     ###'HSPA5', 
    #     ###'MANF', 
    #     'CCND3', 'CDKN2D',
    #     'RAP1B', 'VASP', 
    #     ###'MT1E', 
    #     'RGS10', 
    #     ###'PDIA6', 
    #     'PF4', 
    #     ###'TUBA1B', 
    #     ###'PPIB',
    #     ###'HMGN2', 
    #     ###'SH3BGRL3', 
    #     ###'HSP90B1', 
    #     ###'TOP2A', 
    #     'FYB1', 'LIMS1', 'MPIG6B',
    #     'MYH9', 'ANGPT1', 'AL157895.1', 
    #     ###'SNRPG', 
    #     ###'CKS1B',
    # ],
    'higher_in_bemp_than_mebemp_l_and_eryp': [
        'LMO4', 'MS4A2', 'HDC', 'TPSB2', 'MS4A3', 'TPSAB1', 'ALOX5AP',
        'CLC', 
        'HPGD', 'BACE2', 'TENT5A', 'GPR183', 'GRAP2', 'SRGN',
        'FAM83F', 
        # 'S100A11', 
        'CD63', 'LGALS1', 'IKZF2', 'HPGDS', 
        ###'ANXA1',
        'CNKSR3', 'KIT', 'VWA5A', 
    ],
    'mebemp_l_hla_sig': [
        # NOTE: a better name would have been mebemp_l_mhc_ii_sig
        # in nimrod_148_atlas, this sig did not look good across MEBEMP-L metacells, but across MEBEMP-E metacells it looked good enough.
        'CD74','HLA-DPB1',
        'HLA-DMA', # maybe would have been better without it? nevermind though
        'HLA-DPA1','HLA-DRB1',
        'HLA-DRA', # correlated with the others above (in c_expr vs proj_expr log_ratio, across donors)
    ],
    'bemp_cd74_sig': [
        
        'LMO2','BLVRB','CD74','NAP1L1','CD34','HLA-DRA','MYC','KLF1','IGLL1','PLAC8','MAP7','CAT',
        # corr>0.7 with iter2 above genes, and not in s_phase_sig
        ###'HSP90AB1', # too strong
        'H2AFY', # strong. histone, though replication independent. but looks nice in MC leave one out gene gene scatter
        ###'NPM1', # too strong
        'HSPD1', # strong. heat shock. but looks nice in MC leave one out gene gene scatter
        'HBD', 'GYPC', 'NAA38',
        'IMPDH2', 
        ###'PABPC1', # strong, and looks somewhat meh in MC leave one out gene gene scatter
        ###'MEF2C', # looked meh in leave one out analysis for cells, and looks meh in MC leave one out gene gene scatter
        ###'PRKDC', 'BEX3', # looked meh in leave one out analysis for cells.
        ###'EIF3E', # strong, and looks meh in MC leave one out gene gene scatter
        ###'ANP32B', # strong, and looks meh in MC leave one out gene gene scatter
        'WDR43', 'PA2G4', 'NPR3',

        
        # 'CD34', 'HLA-DRA', 'HBD', 'LMO2', 'H2AFY', 'NAP1L1', 'BLVRB',
        # 'NPM1', 'MYC', 'CAT', 'HSP90AB1', 'NPR3', 'HSPD1', 'PRKDC',
        # 'NAA38', 'MED12L', 'LARS', 'BEX3', 'NOP56', 'CD74', 'FAM117A',
        # 'IGLL1', 'HLA-DRB1', 'C1QBP', 'CPVL', 'MEF2C', 'SLC39A3', 'PLAC8',
        # 'PA2G4', 'IMPDH2', 'SNHG19', 'PPA1', 'REXO2', 'CMBL', 'BZW2',
        # 'MYH10', 'KLF1', 'NPW', 'WDR43', 'LSM5', 'PABPC1', 'TIMM13',
        # 'WEE1', 'LINC02573', 'CIITA', 'FAM69B', 'SLC25A6', 'MRTO4',
        # 'FKBP4', 'RAN', 'SORD', 'EIF4EBP1', 'METAP2', 'RNF24', 'C2orf88',
        # 'ENO1', 'NASP', 'CD164', 'HNRNPAB', 'NCL', 'HLA-DPB1', 'LSM7',
        # 'SNRPD1', 'EIF3E', 'RPL22L1', 'MAP7', 'PVT1', 'SNHG8', 'CCDC34',
        # 'NME4', 'ANP32B', 'GCSH', 'PLXNC1', 'H2AFZ', 'HAGHL', 'SLC27A5',
        # 'CYTL1', 'GPATCH4', 'MPC2', 'NME1', 'ZBTB16', 'TMEM65', 'FAM178B',
        # 'MCM4', 'SERBP1', 'IGSF10', 'DAAM1', 'PRDX2', 'LAPTM4B', 'ATIC',
        # 'TRBC2', 'AVEN', 'DNPH1',
    ],
    'bemp_ms4a2_sig': [
        # top genes in higher_in_bemp_than_mebemp_l_and_eryp missing from higher_in_ms4a2_high_bemp_than_other_bemp: LMO4 and HDC. these two seem to increase earlier in the trajectory...
        # second iteration - comparing BEMPs high and low in the first iteration gene module.
        'MS4A2', 
        'HPGD', 'CLC', 'MS4A3', 'TPSB2', 'CPA3', 'TPSAB1',
        'GPR183', 'SLC18A2', 'HPGDS', 'VWA5A', 'RPS6KA5', 'GRAP2', 'BACE2',
        'CLEC12A', 'ACSL4', 'ALOX5AP', 'SLC45A3', 'SLFN5', 'CLU', 'POU2F2',
        'BTG1', 'IKZF2', 'C12orf75', 'KIT', 'FAM83F', 'LGALS1',
        ###'AL157895.1', 
        'LTC4S', 'ITGB8', 'ERN1', 'TENT5A', 'FCER1G', 'KRT1',
        'RAB27B', 'CDK15', 'NMT2', 'STAP1', 'DENND1B', 'ATP6V0A2', 'ARSG',
        'TPSD1', 'GCSAML', 'FCER1A', 'LGALS8', 'IL1RL1', 'EXD3', 'C4orf48',
        'ID2',
        
        # calculated using BEMP MCs of the MDS illu MC model excluding the cluster of MCs of which most MCs were (relatively) dominated by a single donor
    #     'MS4A2', 'MS4A3', 'HPGD', 'TPSB2', 'CLC', 'GPR183', 'SLC18A2',
    #    'VWA5A', 'POU2F2', 'CPA3', 'RPS6KA5', 'FAM83F', 'CLEC12A', 'IKZF2',
    #    'SLFN5', 'GRAP2', 'ACSL4', 'CLU', 'BACE2', 'C12orf75',
    #    ###'AL157895.1', 
    #    'ITGB8', 'BTG1', 'SLC45A3', 'ERN1', 
       ###'DTWD2', # looked meh in leave_one_out.
       
       ###'HDC',
       ###'LMO4',
    ],
    # 'higher_in_ms4a2_low_bemp_than_other_bemp': [
    #     # 'HLA-DRA', 'LGALS1', 'CD74', 'HLA-DPA1', 'LST1', 'CPPED1', 'PLEK',
    #     # 'HLA-DRB1', 'CYTL1', 'HLA-DPB1', 'CPA3', 'MSI2', 'STXBP6', 'CD34',
    #     # 'MGST3', 'MIR181A1HG', 'HBD', 'HLA-DRB5', 'PIM1', 'LY6E', 'MPP1',
    #     # 'FREM1', 'IFI16', 'STXBP5', 'TYROBP', 'CDKN1C', 'ANGPT2', 'CRYGD',
    #     # 'GRAP2', 'IL3RA', 'STK4', 'COTL1', 'TESPA1', 'TESC', 'FKBP5',


    #     # the following is before filtering by leave_one_out results
    #     # mc_mask here is to exclude the cluster of MCs of which most MCs were (relatively) dominated by a single donor
    #     # mc_mask = (mc_ad.obs['state'] == 'BEMP') & (~mc_ad.obs_names.isin([
    #     #     'M2003.46','M2023.51','M1984.39','M2050.37','M1990.01','M2049.95','M1994.51','M2035.26','M2008.39','M2030.18','M2012.07','M2041.11','M2025.66','M2037.87','M2051.40','M2021.08','M1995.50','M2009.06']))
    #     # diff_exp_df = mc_utils.get_diff_exp_df(
    #     #     curr_ad,
    #     #     mc_utils.mc_mask_to_c_mask(mc_ad, c_ad, mc_mask) & c_ad.obs['c_state'].isin(['BEMP']) & (c_ad.obs['bemp_ms4a2_sig'] > -7.8) & (c_ad.obs['bemp_ms4a2_sig'] < -6.8) & ~c_ad.obs['diagnosis'].isin(['normal']),
    #     #     mc_utils.mc_mask_to_c_mask(mc_ad, c_ad, mc_mask) & c_ad.obs['c_state'].isin(['BEMP']) & (c_ad.obs['bemp_ms4a2_sig'] > -7.8) & (c_ad.obs['bemp_ms4a2_sig'] < -6.8) & c_ad.obs['diagnosis'].isin(['normal']),
    #     # )
    #     'LGALS1', 
    #     'HLA-DRA', 
    #     'LST1', 
    #     'HLA-DRB1', 
    #     'CPPED1', 
    #     'HLA-DPB1',
    #     'CYTL1', 'CD74', 'MIF', 'PIM1', 'CRYGD', 'STXBP6', 'RABAC1',
    #     'PLEK', 'TYROBP', 'ETHE1', 'LY6E', 'MGST3', 'GRAP2', 'SMDT1',
    #     'STK4', 'PNMT', 'COTL1', 'MSI2', 'PCMT1', 'TESC', 'RPS19BP1',
    #     'ADCYAP1', 'MIR181A1HG', 'FREM1', 
    #     'HLA-DPA1', 
    #     'TRBC2', 'MPP1',
    #     'HBD', 'CPA3', 'FKBP5', 'CD38', 'TTC5', 'C11orf74', 'STXBP5',
    #     'BAX', 'ANGPT2', 'CD34', 'MARCKSL1', 'IL3RA', 'CTSC', 'IFI16',
    #     'TIMP3', 'LTC4S', 'DBI', 'SELENOW', 'IL18', 'TRMT1', 'TAF9',
    #     'ATP6V1G1', 'FAM114A2', 'ZNRF1', 'POLR2K', 'FABP5', 'PCID2',
    #     'ANXA11', 'CXXC5', 'SPATC1L', 
    #     'HLA-DRB5', 
    #     'CDKN1C', 'TESPA1',
    #     'OXLD1', 'DECR1', 'PMVK', 'ITM2A', 'CYB5B', 'PFDN2', 'LAPTM4B',
    #     'GTF3C6', 'EMID1', 'IER2', 'LSP1',
    # ],
    # 'higher_in_bemp_than_ep': [ # 240619: currently unused
    #     # >2
    #     'LMO4', 'HDC', 'MS4A2', 'TPSB2', 'MS4A3', 'TPSAB1', 
    #     ###'VIM',
    #     'ALOX5AP', 'HPGD', 
    #     ###'ANXA1', 
    #     'HPGDS', 
    #     ###'AHNAK', 
    #     ###'LGALS1', 
    #     'IKZF2',
    #     'CLC', 'SRGN', 'BACE2', 'GATA2', 'TENT5A', 'GPR183', 'GRAP2',
    #     ###'S100A11',
    #     'POU2F2', 'FAM83F', 'CD63', 'CPA3', 'LTC4S', 'VWA5A',
    #     'CLEC12A',
    # ],
    'higher_in_mpp_than_mebemp_l': [
        # >1.5
        'SPINK2', 
        'AVP', 
        'CD52', 'HOPX', 'SELL', 'C1QTNF4', 
        ### 'VIM', 
        'CSF3R',
        'GLIPR1', 
        ### 'ANXA1',

        # added manually
        'FAM30A',

        # log_ratio 0.5-1.5 and MEBEMP_L expression <-15, filtered manually according to marker heatmap
        'MZB1', 'IGHM', 'MDK', 'ITM2C', 'SORL1', 'ADGRG6', 'PROM1',
        'ATP8B4', 'GIMAP7', 'IFITM3', 'NPDC1', 
        ###'EMP1', 
        'CRHBP', 'SPTBN1',
        'FLT3', 
        ###'TSPAN2', 
        ###'S100A10', 
        'BAALC', 'MBOAT7', 'GBP2', 'BIN1',
        'SOCS3', 'AJ009632.2', 'RBPMS', 
        ###'S100A11', 
        ###'IGLL1', 
        'GBP4',
        'PRAM1', 'SELENOP', 'NOG', 'RNF125', 'PLCB1', 'EVI2B', 
        ###'WDR49',
        'FOSB', 'AC020916.1', 'SCN9A', 'ANXA6', 'BEX1', 'DUSP6',
        'C9orf139', 'HLA-DQB1', 
        ###'DPYSL3', 
        'KRT18', 'TRIM22', 'IL6R',
        'PLEKHO1', 'ANTXR2', 'MYO5C', 'SPNS3', 'DUSP1', 'KRT8', 'RGS10',
        'BTG2', 'TLE4', 'TAGAP', 'DST', 'APP', 'TRPS1', 'RARRES3',
        'CEP126', 'AC026369.3', 'MGST1', 'NFKBIZ', 'CBL', 'POU2F2',
        'ADAM28', 'TCP11L2', 'LIMCH1', 'SPTAN1', 'TCF7L2', 'ADGRE5',
        'RNASET2', 'CALN1', 'MDFIC', 'TM4SF1', 'CD9', 'OPTN', 'HLF',
        'RFLNB', 'PGM2L1', 'RAB32', 'FRMD4B', 'DRAM1', 
        ###'LMNA', 
        'PYGL',
        'CALCRL', 'SULT1C4',
    ],
    'higher_in_mpp_than_bemp': [
        # >2
        'SPINK2', 'HLA-DRA', 'CD74', 'NAP1L1', 'NPR3', 'AVP', 'NRIP1',
        'HLA-DPB1', 'EGFL7', 'HLA-DPA1', 
        ###'PCDH9',
        'HOPX', 'AIF1',
        'LAPTM4B', 'CSF3R', 'CRHBP', 'C1QTNF4', 
        ###'FAM30A', 
        'IFI16',
        'HLA-DRB1', 
        ###'AC011139.1', 
        'HEMGN', 
        ###'GSTP1',
        
        # 1.5-2
        'MDK', 'GLIPR1',
        'VAMP5', 'GUCY1A1', 'ADGRG6', 
        ###'LMO2', ###'SELL', 
        'CD44', 
        ###'PLAC8',
        'LYSMD2', 
        ###'CD52', 
        'FHL1', 'MLLT3', 'PROM1', 'CD34', 'SSBP2',
        'SMIM24', 'GBP4', 'RGS10', 'MGST3', 'SORL1', 
        ###'MEF2C', 
        'BAALC',
        'LST1', 
        ###'PRDX1'
    ],
    'higher_in_mpp_than_mkp': [
        'SPINK2', 'CD52', 'CD74', 'SELL', 
        ###'AHNAK', 
        'PNRC1', 
        ###'ANXA1', 
        'AVP',
        ###'ANKRD28', 
        'HLA-DPB1', # uncomment?
        'SMIM24', 'C1QTNF4', 
        ###'IDS', 
        'CSF3R',
        ###'PBXIP1', 
        'CD99', # uncomment?
        'BST2', 'HOPX', 
        ###'KLF6', 
        ###'ATP8B4', 
        ###'MYADM',
        'SORL1', 'PTPRC', 'GLIPR1', 'HLA-DRA', 'CRHBP', 
        ###'AC011139.1',
        ###'CALCOCO2', 
        ###'TSC22D3', 
        'CD44', 'ICAM3',
    ],
    'higher_in_mpp_than_ep': [
        # >1.5
        'SPINK2', 'SELL', 'CD52', 
        ###'VIM',
        'AVP', 'SMIM24', 
        ###'ANXA1', 
        'NRIP1',
        'HOPX', 
        ###'AIF1', 
        'CSF3R', 'C1QTNF4', 'CRHBP', 
        'IDS', 
        ###'PCDH9', 
        'FOS',
        ###'AC011139.1', 
        ###'AHNAK', 
        'MDK', 'SORL1', 'GLIPR1', 'ATP8B4', 'CD74',
        'SPTBN1', 'PBXIP1', 'VAMP5', 'HLA-DPB1', 'ADGRG6', 
        ###'KLF6', 
        'ITM2C',
        'MZB1', 'TSC22D3', 
        ###'HEMGN', 
        'BST2', 'PTPRC', 'HLA-DPA1', 
        ###'TUBA1A',
        'FHL1', 'MYADM', 
        ###'GAPT',
        'PROM1', 'TSPAN2', 'HLA-E', 'HLA-DRA',
    ],
    'higher_in_ep_than_bemp': [
        # log_ratio > 1.5
        'HBD', 'BLVRB', 'MYC', 'LMO2', 'PRKAR2B', 'HNRNPAB', 'HBB', 'UROD',
        'MEF2C', 
        ###'CYTL1', 
        'PPA1', 'NAA38', 'TPM1', 'MPC2', 
        ###'NAP1L1',
        'ANK1', 'CAT', 'REXO2', 'TMEM14B', 'TMEM14C', 
        ###'LAPTM4B',

        # log_ratio 0.5-1.5 and BEMP expression <-15, filtered manually according to marker heatmap
        'PLAC8', 'PLXNC1', 'NECAB1', 'EPSTI1',
        'PKIG', 'DNPH1', 'CA1', 'CMBL', 'NPW', 'TST', 'PPP1R14A', 'PVT1',
        'ACSM3', 'APOE', 'LY6E', 'DAAM1', 'NOLC1', 'SPTA1', 'UBAC1',
        'B3GNT2', 'FAM210B', 'NDUFAF2', 'DDI2', 'SLC39A8', 'CFAP97',
        'MCM4', 'FAM162A', 'HBQ1', 'CENPF', 'CD36', 'KIFAP3', 'RAP2B',
        'AKR7A2', 'FXN', 'FAM69B', 'MRPL12', 'CHCHD10', 'KLHL23', 'NET1',
        'NDFIP1', 'SNCA', 'SMIM10', 'FHL2', 'MPP6', 'MRTO4', 'FSCN1',
        'GPATCH4', 'SMC2', 'SLC25A37', 'ASAP1', 'SORD', 'PDLIM1', 'PCNA',
        'MCM7', 'TOMM40', 'MRPS26', 'BID', 'FAM178B', 'PSMG4', 'CTNNAL1',
        'GAR1', 'GART', 'NOC3L', 'ATIC', 'STRADB', 'PUM3', 'ELOVL6',
        'MCM5', 'MEST', 'CENPX', 'ODC1', 'CMSS1', 'STOM', 'ALAD', 'PXMP2',
        'TYMS', 'MACROD1', 'PUS7', 'CYC1', 'BOLA3', 'MTHFD1', 'AFF3',
        'SERPINE2', 'MCM6', 'ISOC2', 'IFRD2', 'MRPL4', 'LINC00891',
        'MYH10', 'NOP14', 'WDR12', 'EPCAM', 'AHSP', 'TRAP1', 'BRIX1',
        'SMIM37', 'CDCA7', 'PAK1IP1', 'C20orf27', 'CKAP2', 'POLD2',
        'PSMG1', 'NPM3', 'NOP16', 'CPNE3', 'MCM2', 'RPF2', 'MATK',
        'DCTPP1', 'ZFPM1', 'EEF1E1', 'TTLL12', 'AMMECR1', 'METTL5',
        'EXOSC5', 'ACAT1', 'CHEK1', 'ICAM4', 'GADD45A', 'AKR1C1', 'VRK1',
        'GBGT1', 'HADH', 'BCL7A', 'EBNA1BP2', 'GTPBP4', 'RCL1', 'PCCB',
        'FAM89A', 'MIF', 'SYNGR1', 'CDT1', 'CENPH', 'GMNN', 'MRPL14',
        'WDR3', 'GINS2', 'TFDP1', 'CDH1', 'CPVL', 'GLRX3', 'KLF3', 'HELLS',
        'MICAL2', 'MRPS35', 'HIGD1A', 'RNF169', 'FAM118A', 'RHAG', 'CD320',
        'UTP18', 'TRIP6', 'TTF2', 'PITPNB', 'NAE1', 'TRMT10C', 'COPS3',
        'HDDC2', 'KLHL5', 'ABO', 'MFNG', 'FN3K', 'PPA2', 'HMGN5', 'PSMD6',
        'RRS1', 'BRCA1', 'NTHL1', 'CENPK', 'UBE2V2', 'RPA2', 'CMAS',
        'SURF6', 'CLSPN', 'UTP14A', 'TSR1', 'DTL', 'IMP4', 'ALYREF',
        'SNRPA1', 'CTPS1', 'TGFBRAP1', 'CCDC28B', 'CENPV', 'LAP3',
        'CD3EAP', 'UBE2F', 'RFK', 'NDUFAF4', 'LRRCC1', 'RRP1', 'TRIM58',
        'ZNF385A', 'TNFRSF25', 'SPRY1', 'NIP7', 'MRPL17', 'CISD1', 'GALK1',
        'RPA3', 'HOXA10',
    ],
    'higher_in_mebemp_l_than_gmp_l': [
        # >2
        'SLC40A1', 'FCER1A', 'CYTL1', 'AL157895.1', 'ZNF385D', 'CNRIP1',
       'HBD', 'TPM1', 'GATA2', 'TESPA1', 'PDZD8', 'HEMGN', 'PBX1', 'IL1B',
       'PCDH9', 'CTNNBL1', 'KLF1', 'MINPP1', 'SOD2', 'CPPED1', 'ZMYND8',
       
       # 1.7-2
       'MARCKSL1', 'TLN1', 'ATP6V0A2', 'ZBTB20', 'MED12L', 'CYTOR',
       'CSF2RB', 
       ###'AC011139.1', 
       'GATA1', 'PLXDC2', 'ALDH1A1',
    ],
    'higher_in_eryp_than_gmp_l': [
        'SLC40A1', 'CNRIP1', 
        ###'S100A6', 
        'KLF1', 'HBD', 'APOC1', 'TPM1',
        'PDZD8', 'BLVRB', 'CYTOR', 'CSF2RB', 'GATA1', 'CYTL1', 'CTNNBL1',
        'TFRC', 'CASP3', 'ZNF385D', 'MARCKSL1', 'EMP3', 'IL1B', 'MINPP1',
        'PRKAR2B', 'ZMYND8', 'FCER1A', 'GATA2', 'ITGA2B', 'FBXO7',
    ],
    # 'higher_in_eryp_than_mkp_and_mebemp_l': [ # 240619: currently unused
    #     'APOC1', 
    #     ###'S100A6', 
    #     'BLVRB', 'MYC', 
    #     ###'S100A4', 
    #     'CNRIP1', 'CASP3',
    #    'CSF2RB', 'TMEM14C', 'MPC2', 'PDCD4', 'TRIB2', 'MPST', 'NFIA',
    #    'FAM45A', 'FBXO7', 'HDAC7', 'MYB', 'CYTOR',
    # ],
    # 'higher_in_eryp_than_mkp': [ # 240619: currently unused
    #     'PNRC1', 'FAM117A', 
    #     ###'SLC40A1', 
    #     'ZFP36L2', 'CNRIP1', 'CD74',
    #     'CYTOR', 'PDCD4', 'ANKRD28', 'TXNIP', 'CSF2RB', 'BLVRB',
    #     ###'AC084033.3', 
    #     'MYB', 
    #     ###'AHNAK', 
    #     'SLC12A6', 
    #     'EIF3F', 
    #     'FBXO7', 
    #     'KLF1',
    #     'FOXP1', 'ZNF385D', 
    #     ###'SOX4', 
    #     'CFLAR', 'CYTL1', 'ITGA4', 'BAZ2B',
    #     'SBNO1', 'C6orf48', 'HDAC7', 'ICAM3', 'SESN3', 
    #     ###'ZFAS1', 
    #     'ACSM3',
    #     'RNF24', 'TSC22D1', 'BST2',
    # ],
    'higher_in_ep_than_mpp': [
        # >2
        'APOC1', 
        ###'S100A6', 
        'KLF1', 'MYC', 
        ###'HBD', 
        'CNRIP1', 
        ###'S100A4',
        'BLVRB', 'TFRC', 'GATA1', 'CSF2RB', 'CASP3', 'SLC40A1', 'PRKAR2B',
        'CYTOR', 'ITGA2B', 
        'PDZD8', # is this ok?
        'HBB', 'MPC2', 'FBXO7',

        # 1.5-2
        'HNRNPAB',
        'TMEM14C', 'PCLAF', 
        ###'EIF3A', 
        'PDCD4', # is this ok?
        'EIF5A', 'RIF1', 'TST',
        'MPST', 'CTNNBL1', 'CD36', 'FAM45A', 'ANK1', 
        ###'NCL', 
        'PLIN2',
        'DNAJC9', 'KCNH2', 'PHF6', 'SMIM1', 'NAA38', 'HBS1L', 'TFR2',
        'TRIB2', 'FAM178B', 'EMP3', 
        ###'HSPD1', 
        'PVT1', 'CST3', 
        ###'SERBP1',
        'H2AFZ', 'LDB1', 
        ###'YBX1', 
        ###'MYB', 
        'MINPP1', 'SUPT16H', 'C1QBP',
        'RANBP1', 'TYMS',
    ],
    # 'higher_in_gmp_l_than_bemp': [ # 240619: currently unused
    #     # >2
    #     'MPO', 'C1QTNF4', 'AZU1', 'SPINK2', 'IGLL1', 'NAP1L1', 'ELANE',
    #     'CFD', 'CSF3R', 'PLAC8', 'PRTN3', 'CD44', 'CD74', 'CTSG', 'GSTP1',
    #     'HLA-DRA', 'VAMP5', 'RNASE2', 'PKM', 'PPA1', 'ENO1', 'SMIM24',
    #     'NPW', 'PRAM1', 'MGST1', 'HLA-DPA1', 'GLIPR1', 'RFLNB', 'HLA-DPB1',
    #     'HSPB1', 'APLP2', 'CAT', 'TNFSF13B', 'GYPC', 'HLA-DRB1', 'CD99',
    # ],
    # 'higher_in_gmp_l_than_ep': [ # 240619: currently unused
    #     # >2
    #     'MPO', 'C1QTNF4', 'AZU1', 'SPINK2', 
    #     ###'VIM', 
    #     ###'IGLL1', 
    #     'ELANE',
    #    'CSF3R', 'SMIM24', 'SELL', 'LGALS1', 
    #    ###'CFD', 
    #    'PRTN3', 
    #    ###'CD44',
    #    'CTSG', 'RNASE2', 'ATP8B4', 
    #    ###'ANXA1', 
    #    'RAB32', 'VAMP5', 'PRAM1',
    #    ###'MGST1', 
    #    'PKM', 'CALR', 'RFLNB', 
    #    ###'NAP1L1', 
    #    ###'GLIPR1', 
    #    'FOS',
    #    'CEBPA', 
    #    ###'AHNAK', 
    #    'CLEC12A', 'SORL1', 'ITM2C', 
    #    ###'CPA3', 
    #    'HGF',
    #    'XBP1', 
    #    ###'FLT3',
    # ],
    # 'higher_in_pre_b_than_gmp_l': [ # 240619: currently unused
    #     # >3
    #     'DNTT', 'VPREB1', 'LTB', 'VPREB3', 'MME', 'CD79B', 'SOCS2', 'IGHM',
    #     'ARPP21', 'AKAP12', 'CD79A', 'MZB1', 'EBF1', 'TCF4', 'TP53INP1',
    #     'MEF2C', 'BLNK', 'IL7R', 'ZCCHC7', 'PAG1', 'CD24', 'BTG2', 'RAG1',
    #     ###'FAM30A', 
    #     'JCHAIN', 'LAPTM5', 'RCSD1', 'MYLK',
    # ],
    # 'higher_in_nktdp_than_gmp_l': [ # 240619: currently unused
    #     # log_ratio > 3
    #     'ACY3', 'ID2', 'JCHAIN', 'SAMHD1', 'RUNX2', 'DDIT4', 'TTN',
    #     'MEF2C', 'RUNX3', 'CYTH4', 'CORO1A', 'CALM2', 'ITGB7', 'NFIL3',
    #     'RIPOR2', 'LSP1', 'MEF2A', 'GPR183', 'SCPEP1', 'TNFRSF18', 'MAF',
    #     'SEMA4D', 'NCF1', 'LTB', 'NFKB1', 'CYSLTR1', 'RAB11FIP1',
    #     'SLC38A1', 'BLNK', 'ATP10A', 'STK17A',
    # ],
    'higher_in_gmp_l_than_mpp': [
        # >2
        'MPO', 'AZU1', 'IGLL1', 'ELANE', 'CFD', 'PRTN3', 'C1QTNF4',
        'RNASE2', 'CTSG', 'CALR', 'LGALS1', 
        'SRGN', 'NPW', 'MGST1',
        'CLEC11A', 'CPA3',
        # 1.7-2 - better without, i think (for better separation between MPP and GMP-L)
        'CEBPD', 'P4HB', 'NCOA4', 'LYST', 'CLEC12A',
        'TUBA1B', 
        'HSP90B1', # comment? sounds batchy...
        'PRSS57', 'TNFSF13B', 'CEBPA', 'RAB32',
        'EIF3A', 'HGF', 'MYC', 'LYZ', 'ZEB2',
    ],
    # 'higher_in_bemp_than_unassigned_mebemp_l': [ # 240619: currently unused
    #     'LMO4', 'HDC', 'MS4A2', 'TPSB2', 'MS4A3', 'ALOX5AP', 'CNRIP1',
    #     'TPSAB1', 
    #     ###'S100A6', 
    #     'APOC1', 'CLC',
    # ],
    # 'higher_in_gmp_l_than_unassigned_mebemp_l': [ # 240619: currently unused
    #     # >1.5
    #     'MPO', 'IGLL1', 'C1QTNF4', 'SPINK2', 'AZU1', 'CSF3R', 'MGST1',
    #    'CD44', 'CD99', 'CFD', 'TNFSF13B', 'CTSG', 'PLAC8', 'SMIM24',
    #    'RAB32', 'CLEC11A', 'PRAM1', 'ELANE',
    #    # 1.3-1.5
    #    'HSPB1', 'PRTN3', 'SELL',
    #    'LGALS1', 'RNASE2', 'GLIPR1', 'CEBPA', 'PPA1',
    # ],
    # 'higher_in_unassigned_mebemp_l_than_mebemp_l_and_eryp': [ # 240619: currently unused
    #     'LGALS1', 'CALR', 
    #     ###'S100A10', 
    #     'CPA3', 
    #     'HSP90B1', 
    #     ###'ATP8B4', 
    #     'XBP1',
    #     ###'AHNAK', 
    #     'SRGN', 'IGLL1', 'TPSAB1', 'CD99', 'FOS', 'TPM4', 'MGST1',
    #     'RNF168', 
    #     ###'CALCOCO2', 
    #     'RHEX', 'TPSB2', 'C1QTNF4', 'IKZF2', 'ZEB2',
    #     'SORL1', 'CD44', 'PKM', 'LRRC75A', 'RAB32', 'CLEC12A',
    # ],
    # 'higher_in_gmp_l_than_mebemp_l': [ # 240619: currently unused
    #     # >1.5
    #     'MPO', 'IGLL1', 'AZU1', 'C1QTNF4', 'ELANE', 'CFD', 'PRTN3',
    #    'SPINK2', 
    #    ###'VIM', 
    #    'CSF3R', 'CTSG', 'MGST1', 'RNASE2', 'LGALS1',
    #    'CD44', 'CALR', 'PRAM1', 'RAB32', 'CD99', 'NPW', 'TNFSF13B',
    #    'CLEC11A', 'CEBPA', 
    #    ###'NAP1L1', 
    #    'SMIM24', 'PLAC8', 'CLEC12A', 'SELL',
    #    'CEBPD', 'HSPB1', 'PKM', 'HSP90B1', 'HGF', 'GLIPR1', 'ATP8B4',
    #    'LYST', 'VAMP5', 'TUBA1B', 
    #    'NUCB2', 
    #    'FLT3', 
    #    ###'ANXA1', 
    #    'XBP1',
    #    'P4HB', 
    #    ###'PRSS57', 
    #    'RFLNB', 'LYZ', 'PLPPR3', 
    #    ###'ENO1', 
    #    'PPA1', 'MZB1',
    #    'IL6R', 'CST7', 'NCOA4', 'ITM2C', 'KBTBD11', 'S100A11', 
       
    #    # 1-1.5
    #    'FNDC3B',
    #    'PPIB', 'FAM107B', 'ADGRE5', 'C9orf139', 'GYPC', 'SPNS3', 
    #    ###'GSTP1',
    #    'S100A10', 'FOS', 'APLP2', 'SUCNR1', 'DBI', 'CANX', 'CD48',
    #    'ANKRD28', 
    #    ###'SRGN', 
    #    ###'GAPDH', 
    #    'CRYBG1', 'HSPA5', 'CSTA', 'RASGRP2',
    #    'SYNGR1', 'SORL1', 'PYGL', 'CITED4', 'IGFBP7', 'TOP1MT', 'ARMH1',
    #    ###'CD74', 
    #    'CALCOCO2', 
    #    ###'HLA-A', 
    #    'ANXA2', 'TRIM14', 'KLF6', 'DUSP6',
    #    'F13A1', 'IFITM3', 'IGFBP2', 'EMILIN2', 
    #    ###'AHNAK', 
    #    'BCL2', 'EVI2B',
    #    ###'VAMP8', 
    #    'MYADM', 'PTPRE', 'P2RY8', 'RAB27A', 'PLCB1', 'NAALADL1',
    #    'CYBA', 'RNASET2', 'TENT5A', 'HCST', 'RNASE3', 'LCP1', 'LSP1',
    #    'EMB', 'MS4A3', 'NPDC1', 'AGTRAP', 'PROM1', 'ATP2B1', 'UBE2J1',
    #    'SPI1', 'MCL1', 'LAIR1', 'SPARC', 'KCNAB2', 'MAMDC2', 
    #    'HOPX',
    #    'APPL1', 'SAMHD1', 'CEACAM4', 'VAT1', 'PIK3R1', 'RAB7B', 
    #    ###'HNRNPU',
    #    'TMED10', 'ERLIN1', 'GPI', 'BRI3BP', 'PDIA6', 'MAP3K1', 'MRPL57',
    #    'TRAM1', 'CDCA7', 'SLC2A4RG', 'APP', 'SERPINB8',
    # ],
    'higher_in_bemp_than_gmp_l': [
        # >2
        'LMO4', 'HDC', 'MS4A2', 'CNRIP1', 'TPSB2', 'TPSAB1', 'GATA2',
        'S100A6', 'ALOX5AP', 'SLC40A1', 'CSF2RB', 'CASP3', 'FCER1A',
        'APOC1', 'HPGDS', 'S100A4', 'AL157895.1', 'HPGD', 'HBS1L', 'CLC',
        'CD63', 'ALDH1A1', 'IKZF2', 'SOX4', 'CTNNBL1', 'FBXO7', 'ST8SIA6',
        'NLK', 'MS4A3', 'FTH1', 'BACE2', 'EMP3', 'PDZD8', 'ATP6V0A2',
        'TIMP1', 'TP53INP1', 'P2RX1', 'GDF11', 'GRAP2', 'TFRC', 'CNKSR3',
        'TSC22D1', 'TAL1', 'SNX5', 'TESPA1', 'GPR183', 'FAM83F', 'KCNH2',
        'PBX1', 'LAPTM5', 'MARCKSL1', 'KLF1', 'SUPT20H', 'MINPP1', 'CD82',
        'LTC4S', 'RYR3', 'ACSL4', 'GATA1', 'ERN1', 'TRIB2', 'ITGA2B',
        'NMT2', 'EXOSC8', 'EXD3', 'ZFP36L1', 'SEMA3C',
    ],
    # 'higher_in_mbemep_l_than_gmp_e': [ # 240619: currently unused
    #     # >1.5
    #     'SLC40A1', 'HBD', 'CNRIP1', 'KLF1', 'FCER1A', 'CYTL1', 'PDZD8',
    #     'TFRC', 'CSF2RB', 'MINPP1', 'GATA1', 
    #     ###'S100A4', 
    #     'CTNNBL1', 
    #     ###'S100A6',
    #     'ZMYND8', 'PLEK', 'BMP2K',
    #     # 1.3-1.5
    #     'GATA2', 'CYTOR', 'AL157895.1', 'TLN1',
    #     'FTH1', 'PRKAR2B', 'ABCC4', 'ZEB2', 'PSTPIP2', 'ATP6V0A2', 'TPM1',
    #     'CD84',
    # ],
    'higher_in_bemp_than_mebemp_l': HIGHER_IN_BEMP_THAN_MEBEMP_L_GENES,
    'higher_in_ep_than_mebemp_l': [
        # >1.5
        'APOC1', 'BLVRB', 
        ###'S100A6', 
        'MYC', 'HBB', 'MPC2',
        
        # log_ratio 0.5-1.5 and MEBEMP-L expression <-15, filtered manually according to marker heatmap
        'FAM178B', 'PVT1', 'CD36', 'ANK1', 'TST', 'CA1', 'APOE', 'NFIA', 'FAM45A', 'ITGA2B', 'SMIM10', 'TRIB2', 'DDI2', 'PPP1R14A', 'KCNH2', 'SMIM1', 'LYAR', 'TMEM141', 'ELOVL6', 'SYNGR1', 'EPCAM', 'NECAB1', 'PLIN2', 'MTSS1', 'DAAM1', 'CMBL', 'RREB1', 'TESC', 'PCAT18', 'UBAC1', 'TFDP1', 'SLC25A37', 'CDC42BPA', 'PNMT', 'SLC26A2', 'AHSP', 'SPTA1', 'GBGT1', 'FAM89A', 'ATP13A3', 'RYR3', 'GAR1', 'FAM118A', 'ACAT1', 'ISOC2', 'FAM83D', 'NDUFAF2', 'CBLL1', 'UCA1', 'UGGT1', 
        ###'TOMM40', 
        'MACROD1', 'CSF1', 
        ###'EBNA1BP2', 
        ###'MTHFD1', 
        'HBQ1', 'NMNAT3', 'TNFRSF25', 'SORD', 'TTLL12', 'HMBS', 'CTNNAL1', 'AKR1C1', 'SLC39A8', 'RNH1', 'CDH1',

        # added later
        'REXO2', 'TMEM14C', 'TMEM14B',
    ],
    # 'filtered_rp_most_anti_corr_with_dntt_in_hsc_to_clp_m': FILTERED_RP_MOST_ANTI_CORR_WITH_DNTT_IN_HSC_TO_CLP_M, # 240619: currently unused
    'higher_in_pre_b_than_all_hspcs': [
        # >0.8
        'ARPP21', 'RAG1', 'SOCS2', 'RAG2', 'VPREB1', 'CMTM8', 'HMHB1',
        'UHRF1', 'HPS4', 'LCN6', 'TLE4', 'IGLL1', 'AKAP12', 'EBF1', 'CYGB',
        ###'AL161912.2', # 240129: was missing from /dummy/dummy/dummy/raid/mds/new_N280_exps/intermediate_output/raw_without_excluded_genes_and_cells.h5ad), and quite weak anyway.
        'TOP2B', 'PAG1', 'CMTM7', 
        ###'DNTT', 
        'HHIP', 'RNF152',
        'PPP1R14B', 'SLC8A1-AS1', 'AGPS', 'PTPRE', 'MYLK', 'CTGF',
        'BAHCC1', 'RCSD1', 'LINC00426', 'DYRK3', 'SH2D4B', 'SMARCD2',
        'ZCCHC7', 'NRXN2', 'MT1X', 'P4HA2', 'AL713998.1', 'IRX1',
    ],
    # 'higher_in_nkt_nktdp_than_dc_nktdp': [ # 240619: currently unused
    #     # TODO: maybe improve by another iteration based on these genes and cells? also leave one out?
    #     'MAF', 'SATB1', 'KLRB1', 
    #     ###'FAM30A', 
    #     ###'SOX4', 
    #     'TRIB2', 'TOX2',
    #    'TRBC2', 'TSC22D1', 'LMO4', 'TNFRSF18', 'ETS2', 'IL7R', 
    #    ###'IGHM',
    #    'PCDH9', 'DNMT3B', 
    #    ###'HOPX', 
    #    'CLEC2D', 'IL17RE', 'LTB', 'THRA',
    #    'ST3GAL1', 'RORC', 
    #    ###'TRG-AS1', 
    #    ###'JCHAIN', 
    #    'NUCB2', 'RFLNB',

    #    'TCF7',
    # ],
    # 'higher_in_dc_nktdp_than_nkt_nktdp': [ # 240619: currently unused
    #     # TODO: maybe improve by another iteration based on these genes and cells? also leave one out?
    #     'SAMHD1', 'NFKB1', 
    #     ###'CYBA', # not weak in other nktdp
    #     'SPI1', 'COTL1', 'TLR10', 
    #     ###'CD74', # not weak in other nktdp
    #    'PLD4', 
    #    ###'HLA-DPA1', # not weak in other nktdp
    #    ###'LGALS1', # not weak in other nktdp
    #    ###'HLA-DRA', # not weak in other nktdp
    #    ###'VIM', 
    #    ###'HLA-DPB1', # not weak in other nktdp
    #      'DBI',
    #    'IRF8', 'NDRG2', 'SEL1L3', 'FGL2', 'GRN', 
    #    ###'HLA-DRB1', # not weak in other nktdp
    #    'HLA-DQA1',
    #    'CST3', 'CTSZ', 'HLA-DQB1',
    # ],
    'higher_in_dc_than_all_hspcs': [
        # >2
          'CST3', 'FCER1G', 'PTPRE', 'IRF8', 'S100A10', 'TGFBI', 'LILRA4',
       'CLEC4C', 'SLAMF7', 'IRF7', 'SCT', 'UGCG', 'APP', 'ANXA2', 'LYZ',
       'CLIC3', 'C12orf75', 'SULF2', 'CD4', 'OTULINL', 'PTGDS', 'IGSF6',
       'LILRB4', 'FCGRT', 'GRN', 'GZMB', 'DAB2', 'SHTN1', 'AXL', 'RNASE6',
       'ANXA5', 'S100A11', 
       ###'PLEK', 
       ###'PLAC8', 
       'HERPUD1', 'FGL2', 'CTSZ',
       'GAS6', 'TYROBP',
    ],
    'higher_in_p_dc_than_as_and_c_dc': [
        'GZMB', 'JCHAIN', 'MZB1', 'NUCB2', 'DERL3', 'CCDC186', 'SMPD3',
        'FAM129C', 'TPM2', 'STMN1', 'PPP1R14B', 'SELENOS', 'ITM2C',
        'COBLL1', 'ZFAT', 'SEC61G', 'MAP1A', 'CYP46A1', 'MAPKAPK2', 'ETS1',
        'BCL11A',  # comment? a bit strong in pDC
        'CDK2AP2', 'CYSLTR1', 
        'UGCG',  # comment? a bit strong in pDC
        'ERN1', 'SLC15A4', 'KCTD5',
        'PACSIN1', 'SEC11C', 'LMAN1', 
        'SPCS1',  # comment? a bit strong in pDC
        'FCRLA', 
        'IRF7',  # comment? a bit strong in pDC
        'RASD1',
        'HM13',
    ],
    'higher_in_as_dc_than_p_dc': [
        'S100A10', 
        'LYZ', 'FGL2', 
        'ANXA1', 
        'COTL1', 'AXL', 'CTSH', 'ANXA2',
        'KLF4', 'PPP1R14A', 'C20orf27', 'SPI1', 'LST1', 'FAM129A',
        'CX3CR1', 'CLEC10A', 'PPA1', 'CD22', 'FCGRT', 'TNFAIP2', 'PLA2G16',
        'LIMD2', 'PTGDS', 'CFD', 'CFP', 'CD2',
    ],
    'higher_in_c_dc_and_c_monocyte_than_as_dc': [
        # generated the list by comparing cDC and AS-DC
        'FCN1', 'CD52', 'LGALS2', 'CPVL', 'MARCKS', 'MNDA', 'IFITM3',
       'S100A9', 'CAST', 'ZFP36L1', 'DUSP6', 'NCF2', 'SERPINA1', 'MS4A6A',
       'KLF2', 'CLEC7A', 'SCIMP', 'LGALS3', 'CD1C', 'MAFB', 'FYB1',
       'GIMAP7', 'RBPJ', 'GIMAP4', 'CEBPD', 'C1QA', 'GBP1',

        # # according to onek1k
        # 'S100A9', 
        # 'FCER1A', 'FCN1', 'LGALS2', 'CD1C', 'CSTA', 'CD52',
        # 'CPVL', 
        # 'S100A8', 
        # 'MNDA', 'LGALS3', 'DOK2', 'CEBPD', 'AIF1',
        # 'CEBPB', 'ID2', 'VCAN', 'CAST', 'PLBD1', 'NR4A2',
    ],
    'higher_in_monocyte_than_c_dc1': [
        # i.e., nc_c_interm_monocyte
        'S100A4', 'TYROBP',
    ],
    # 'higher_in_c_dc1_than_monocyte': [
    #     # i.e., nc_c_interm_monocyte # didn't use this in the end
    #     'C1orf54', 'CLEC9A', 'DNASE1L3', 'CPNE3', 'IRF8', 'XCR1', 'WDFY4',
    #     'SNX3', 'NET1', 'RAB7B', 'BASP1', 'RGCC', 'STMN1', 'LIMA1', 'PTMS',
    #     'RGS10', 'CLNK', 'SLAMF7', 'CCND1', 'P2RY14', 'HLA-DOB', 'BTLA',
    #     'SLC38A1',
    # ],
    
    # 'higher_in_c_dc1_than_nktdp': [
    #     'LYZ', 'S100A10', 'CST3', 'MNDA', 'CPVL', 'LGALS2', 'C1orf54',
    #    'PTPRE', 'SHTN1', 'CLEC9A', 'DNASE1L3', 'RGS10', 'ANXA2',
    # ],
    # 'higher_in_nktdp_than_c_dc1': [
    #     'ACY3', 'SPINK2', 'JCHAIN', 'CYSLTR1', 'RUNX2', 'DDIT4', 'ARMH1',
    #    'CTSW', 'TTN', 'SEMA4D', 'SCPEP1', 'RUNX3', 'TUT4', 'TESC',
    #    'NUCB2', 'ITM2C', 'IQGAP2', 'GSN', 'SATB1', 'CAPG', 'NFIL3',
    #    'PLP2',
    #    'RNASEH2B', 'PLD4', 
    #    ###'AL096865.1', 
    #    'CYTH4', 'LAT2', 'MAF',
    #    'IGHM', 'RASD1', 'TNFRSF18', 'ITGB7', 'NKG7', 'CYBB', 'POLB',
    #    'MED13L', 'ATP10A', 'MACF1', 'TRAF3IP3', 
    #    ###'POLR2J3---1', 
    #    'AKAP9',
    #    'CXorf21', 'BLNK',
    # ],
    
    'higher_in_interm_monocyte_than_c_dc2': [
        # TODO: overwrite in c_ads
        'FCN1', 'CFD', 'CYBB', 'CD14', 'DUSP6', 'GIMAP4', 'SERPINA1',
       'MARCKS', 'GIMAP7', 'VCAN', 'TMEM176B', 'KLF2', 'MAFB',
        
        # 'FCN1', 'CFD', 'CYBB', 'CD14', 'DUSP6', 'GIMAP4', 'SERPINA1',
        # 'MARCKS', 'GIMAP7', 'TMEM176B', 'KLF2', 'MAFB', 'BRI3', 'CEBPB',
        # 'APOBEC3A', 'NCF1', 'CD68', 'CD300E', 'DMXL2', 'SLC7A7', 'S100A12',
        # 'LILRB2',
        
        # 'S100A8', 'VCAN', 'FCN1', 'CD14', 'DUSP6', 'SERPINA1', 'CFD',
        # 'CD300E', 'ASAH1', 'STXBP2', 'S100A12', 'CEBPB', 'TNFRSF1B',
        # 'SLC7A7', 'RXRA', 'CTSB', 'TNFSF10', 'CD36',
    ],
    'higher_in_as_dc_than_c_dc2': [
        'IRF8', 'LILRA4', 'ALOX5AP', 'TCF4', 'SCT', 'DAB2', 'APP', 'UGCG',
        'PPP1R14A', 'PTGDS', 'CLIC3', 'C12orf75', 'IL3RA', 'SPIB', 'IGKC',
        'HERPUD1', 
        'SOX4', # comment?
        'CYB561A3', 'CLEC4C', 'PPP1R14B', 'CCDC50',
        'PLP2', 'AXL',
    ],
    # 'higher_in_c_dc2_than_as_dc': [
    #     'IFITM3', 'MNDA', 'LGALS2', 'MS4A6A', 'FCER1A', 'CSTA', 'S100A9',
    #     'CD52', 'LGALS3', 'GBP1', 'FCN1', 'CAST', 'RBPJ', 'RETN', 'GPAT3',
    #     'CEBPD', 'SCIMP',
    # ],
    'higher_in_c_dc2_than_interm_monocyte': [
        ###'HLA-DQA1', # not weak in intermMonocyte
        ###'HLA-DQB1', # not weak in intermMonocyte
        'FCER1A', 'SLC38A1', 'BASP1', 'CLEC10A', 'C12orf45', 'ENHO',
        'SNHG8', 'TCEA3', 'SNHG25', 'NREP', 'AFF3', 'CCND2', 'NDRG2',
        'C12orf75', 'PDLIM1', 'CD2', 'HMGN1', 'SPATS2L', 'ADAM28', 'ALCAM',
        'EPB41L2', 'PPP1R14A', 'ATP1B1', 'LIMA1', 'PKIB', 'CD1C', 'MZT2A',
        'RUNX2', 'MRC1', 'CLIC2', 'DHRS9', 'FLT3', 'PTMS', 'CST7',
        'BCL11A', 'SEPT6', 'ITM2C', 'CBX5',
        
        'BCL7A', 
        ###'RETN', 
        'GPR183',
        'TCEAL8', 'APEX1', 'P2RY14', 'STMN1', 'NAV1', 'PNN', 'BIN1',
        # ###'AC245014.3', 
        # 'LDHB', 'ZEB1', 'METAP2', 'APPL1', 'PHACTR1',
        # 'SNHG19',
        
        
    #     'CD1C', 'FCER1A', 'ENHO', 'CLEC10A', 'NDRG2', 'BASP1', 'PKIB',
    #     'GSN', 'CLIC2', 'GARS', 'NCOA7', 'HLA-DOA', 'RGS10', 'DNASE1L3',
    #     'SINHCAF', 'AFF3', 'SRI', 'TMEM14C', 'EEF1E1', 'CD2', 'CCDC144A',
    #     'C12orf75', 'PON2', 'FAM43A', 'SLC25A19', 'NCKAP5', 'KHDRBS1',
    #     'SPECC1', 'IL18', 'SLC38A1',
        
    #     'TP53BP1', 'ARL4C', 'BATF3', 'SCD',
    #    'BOLA1', 'GRIP1', 'DIMT1',
    #    'DHRS3', 'TUBA1C', 'HSPD1', 'THRAP3',
    #    'CCND2', 'HMGN1', 'PHACTR1', 'OSTC', 'FAM104B', 'CCSER1', 'PGP',
    #    'SRSF8',
    ],
    'higher_in_c_dc2_than_all_hspcs': [
        # constructed on 240401 by diff exp between cDC2 MCs excluding c_state=monocyte and HSC_MPP, BEMP, CLP, GMP-L, NKTDP
        'MNDA', 'CST3', 'LYZ', 'S100A10', 'MS4A6A', 'CLEC10A', 'FCER1G',
        'IGSF6', 'LGALS2', 'S100A9', 'HLA-DQB1', 'HLA-DQA1', 'S100A11',
        # 'ANXA2', 'ANXA5', 
        # 'C1orf162', 'CSTA', 'CD48', 'LGALS3',
    ],

    'higher_in_c_dc1_than_as_dc_and_c_dc2': [
        'CLEC9A', 
        ### 'THBD', # (AKA CD141) according to https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01088/full, but not in our data??
        'XCR1', # according to https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01088/full and our data
        'BATF3', # according to https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01088/full, but only slightly in our data?
        'ID2', # according to https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01088/full and our data
        'CPNE3', # according to our data? also kind of fits 1k1k
        'C1orf54', # according to https://www.proteinatlas.org/ENSG00000118292-C1orf54 and our data, and 1k1k?
        'CCND1', # according to https://www.proteinatlas.org/ENSG00000110092-CCND1 and our data, and 1k1k?
        'RAB7B', # according to https://www.proteinatlas.org/ENSG00000276600-RAB7B and our data, and 1k1k?
        'SNX3', # according to our data, and 1k1k? https://www.proteinatlas.org/ENSG00000112335-SNX3 actually says it is quite the same in pDC
        'DNASE1L3', # according to our data, and 1k1k (higher than in AS-DC and cDC2, but lower than pDC)? https://www.proteinatlas.org/ENSG00000163687-DNASE1L3 actually says it is higher in pDC
        'CAMK2D', # according to our data and 1k1k, and https://www.proteinatlas.org/ENSG00000145349-CAMK2D?
        'CADM1', # according to our data and 1k1k, and https://www.proteinatlas.org/ENSG00000182985-CADM1
        'IDO1', # according to our data (though only slightly higher than cDC1) and 1k1k, and https://www.proteinatlas.org/ENSG00000131203-IDO1/immune+cell
        'SHTN1', # according to our data and 1k1k, and https://www.proteinatlas.org/ENSG00000187164-SHTN1/immune+cell (which says it is quite the same in pDC)
        'WDFY4', # according to our data and 1k1k, and https://www.proteinatlas.org/ENSG00000128815-WDFY4
        'CYB5R3', # according to our data (not strong, though) and 1k1k. 
        ###'RGS10', # according to our data and 1k1k, and https://www.proteinatlas.org/ENSG00000148908-RGS10/immune+cell. a bit strong in cDC2 (and also AS-DC?)
        'CLNK', # according to our data and 1k1k
        'NET1', # according to our data and 1k1k
        'RAB32', # according to our data and 1k1k, and https://www.proteinatlas.org/ENSG00000118508-RAB32
        'AP3M2', # according to our data and 1k1k
        'LYRM4', # according to our data and 1k1k
        'BASP1', # according to our data and 1k1k
        'LIMA1', # according to our data and 1k1k
        'HLA-DOB', # according to our data and 1k1k
        'BTLA', # according to our data and 1k1k
    ],
    # 'higher_in_nc_monocytes_than_c_dc': [ # 240619: currently unused
    #     # according to onek1k
    #     'FCGR3A', 
    #     'SMIM25', # SMIM25 is PELATON in onek1k,
    #     'MS4A7', 'CFD', 'SERPINA1', 'IFITM3', 'CD68',
    #    'CDKN1C', 'HES4', 'TCF7L2', 'MAFB', 'BCL2A1', 'C5AR1', 'SLC11A1',
    #    'LILRB2', 'PECAM1', 'MTSS1', 'SLC7A7', 'TESC', 'LYPD2', 'TNFRSF1B',
    #    'LILRA5', 'VMO1', 'CKB', 'HMOX1', 'PILRA', 'ASAH1',
    # ],
    # 'higher_in_c_dc_than_nc_monocytes': [ # 240619: currently unused
    #     'HLA-DQA1', 'HLA-DQB1', 'LGALS2', 'HLA-DMA', 'CLEC10A', 'HLA-DRB5',
    #     'CD1C', 'MS4A6A', 'CPVL', 'HLA-DMB', 'CIITA', 'KCNK6', 'RBPJ',
    #     'CLIC2', 'FCER1A', 'RGS10',
        
    #     # # according to onek1k
    #     # 'FCER1A', 'CLEC10A', 'LGALS2', 'HLA-DQA1', 'CD1C', 'CAPG',
    #     # 'HLA-DMA', 'HLA-DMB', 'RNASE6', 'ALDH2', 'ENHO', 'GSN', 'HLA-DQB1',
    #     # 'JAML', 'CD99', 'LY86', 'SELL', 'ARL4C', 'PHACTR1', 'PPA1',
    #     # 'GPR183', 'SPINT2', 'NDRG2', 'LDHB',
    # ],
    # 'higher_in_monocytes_than_c_dc': [ # 240619: currently unused
    #     'S100A8', 'VCAN', 'S100A12', 'CD14', 'CD36', 'CTSD', 'CSF3R',
    #    'CYP1B1', 'FPR1', 'RGS2', 'CFD', 'LRRK2', 'IL17RA', 
    #    ###'CDA', # looked meh in leave_one_out.
    #    'STXBP2', 'RBP7', 'MEGF9', 'CLEC4E',
        
    #     # # according to onek1k
    #     # 'S100A8', # comment?
    #     # 'S100A12', 'CD14', 'CFD', 'VCAN', 'NCF1', 'RBP7', 'CXCL8', 'CDA',
    #     # 'SERPINA1', 'CYBB', 'CTSD', 'TMEM176B', 'SLC11A1', 'SLC7A7',
    #     # 'SMIM25', # SMIM25 is PELATON in onek1k,
    # ],
    # 'p_as_dc_somewhat_specific': [ # too weak...
    #     ###'CST3', 
    #     'IRF7', 'UGCG', 'LILRA4', 'RNASE6', 'IL3RA', 'HERPUD1',
    # ],
    'c_intermediate_nc_monocyte_somewhat_specific': [
        'MNDA', 'PSAP', 'CEBPB', 'SAT1', 'MS4A7', 
        ###'CTSS', 
        'LILRB2',
       'DUSP6', 
       ###'FTL', 
       'BRI3', 'FGL2', 'CYBB', 'STXBP2', 'MS4A6A',
       'CLEC12A', 'TSPO', 'TNFRSF1B', 'TYMP', 'ASAH1',
    ],
    # 'higher_in_c_dc_than_monocytes': [ # 240619: currently unused
    #     'HLA-DQA1', 
    #     'HLA-DPB1', # not weak in monocytes (~-12.1)
    #     'HLA-DPA1', # not weak in monocytes (~-12.1)
    #     'FCER1A', 'HLA-DQB1',
    #    'CLEC10A', 'CD1C', 'ENHO', 'HLA-DRB5', 'TUBA1B', 'HLA-DMA', 'DBI',
    #    'DNASE1L3', 'PEBP1', 'BATF3', 'DDAH2', 'SLC38A1', 'GABARAPL2',
    #    'ARL4C', 'ATP5MC1', 'NPM1', 
    #    ###'RPL36A', 
    #    'CLIC2', 
    #    ###'SELENOW', # looks bad in leave_one_out.
    #    'UVRAG',
    #    'NDRG2', 'MRC1', 
    #    ###'HNRNPR', # looked meh in leave_one_out.
    #    ###'TIMM8B', # looked meh in leave_one_out.
       
        
    #     # # according to onek1k
    #     # 'FCER1A', 
    #     # ###'HLA-DQA1', 
    #     # 'CLEC10A', 'CD1C', 'ENHO', 'PLD4', 
    #     # ###'GSN',
    #     # 'NDRG2', 'ARL4C', 
    #     # ###'PEBP1', 'PPA1', 'SPINT2', 'HMGA1', 'C1QBP',
        
    #     # oops. didn't realize i have both higher_in_c_dc_than_monocyte and higher_in_c_dc_than_monocytes...
    #     # 'higher_in_c_dc_than_monocyte': [
    #     #     'HLA-DQA1', 'HLA-DQB1', 
    #     #     'DNASE1L3', # not higher in cDC than monocytes in 1k1k
    #     #     'BASP1', 'FCER1A', 'CD1C',
    #     #     'NCOA7', # not much higher in cDC than monocytes in 1k1k
    #     #     'NDRG2', 
    #     #     'CLEC10A', 
    #     #     'RGS10', # not much higher in cDC than monocytes in 1k1k
    #     #     'ENHO', 'GSN', 
    #     #     ###'TUBA1B', # same level in cDC and ncMonocytes in 1k1k
    #     #     'IDO1', # not higher in cDC than monocytes in 1k1k
    #     #     'ATP5MC1', 'CLIC2', 
    #     #     ###'BATF3', # lower in cDC than ncMonocytes in 1k1k
    #     # ],
    # ],
    # 'higher_in_b_than_monocytes_and_nc_monocytes': [ # 240619: currently unused
    #     'MS4A1', 'IGKC', 'IGHM', 'RALGPS2', 'CD79A', 'BANK1', 'IGHD',
    #     'IGLC2', 'IKZF3', 'FCRL5', 'CD22', 'FCMR', 'FCRL2', 'TCF4',
    #     'LINC00926',
    # ],
    # 'higher_in_c_dc_than_monocytes_and_nc_monocytes': [ # 240619: currently unused
    #     'FCER1A', 'CLEC10A', 'HLA-DQA1', 'CD1C', 'ENHO', 'HLA-DPB1',
    #    'HLA-DQB1', 'HLA-DRB1', 'HLA-DMA', 'HLA-DMB', 'HLA-DRA', 'GSN',
    #    'NDRG2', 'ARL4C', 'PPA1', 'HMGN1', 
    #    ###'CD74', 
    #    'SPINT2', 'HMGA1',
    #    'PEBP1',
    # ],
    'higher_in_p_and_as_dc_than_monocyte': [
        'IRF8', 'TCF4', 'UGCG', 'CCDC50', 'C12orf75', 'ITM2C', 'PPP1R14B',
        'SPIB', 'ALOX5AP', 'JCHAIN', 'LILRA4', 'MZB1', 'PLD4', 'BCL11A',
        'APP', 'IGKC', 'IRF7', 'STMN1', 'CLEC4C', 'PPP1R14A', 'SCT',
        'SCN9A', 'CLIC3', 'ARL4C', 'IL3RA', 'GZMB', 'HERPUD1', 'DAB2',
        'RNASE6', 'RUNX2', 'SEPT6', 'APEX1', 'HMGN1', 'CD2AP', 'OFD1',
    ],
    # 'higher_in_c_dc_than_p_dc': [ # 240619: currently unused
    #    'LYZ', 'COTL1', 'MNDA', 'FGL2', 'CPVL', 'FCN1', 'LST1', 'LGALS2',
    #    'MARCKS', 'NCF2', 'CTSH', 'IFITM3', 'IGSF6', 'TIMP1', 'TYMP',
    #    'S100A9', 'CSTB', 'ZFP36L1', 'LGALS3', 'TNFAIP2', 'CAST', 'MS4A7',
    #    'TSPO', 'CFP', 'KLF2', 'SCIMP', 'SERPINA1', 'DUSP6', 'SPI1',
    #    'CEBPD', 'RBPJ', 'CLEC7A',

        
    #     # # according to onek1k
    #     # 'LYZ', 'IFI30', 'COTL1', 'ANXA1', 'CFP', 'CSTA', 'CLEC10A', 'SPI1',
    #     # 'LST1', 'LGALS2', 'MNDA', 'TIMP1', 'S100A9', 
    #     # ###'S100A10', 
    #     # ###'FCER1A',
    #     # ###'AIF1', 
    #     # 'CEBPD', 'CPVL', 'CASP1', 'LY86', 'FCN1', 'CD1C', 
    #     # ###'ANXA2',
    # ],
    # 'higher_in_as_dc_than_c_dc': [ # 240619: currently unused
    #     'PPP1R14A', 'PLAC8', 'LILRA4', 'IRF8', 'C12orf75', 'ALOX5AP',
    #    'AXL', 'ITM2C', 'TCF4', 'SCT', 'IRF7', 'SIGLEC6', 'CYB561A3',
    #    'CLIC3', 'DAB2', 'SPIB', 'SOX4', 'LIME1', 'CCDC50', 'APP', 'PTGDS',
    #    'IGKC', 'IL3RA', 'TSPAN13', 'DNASE1L3', 'UGCG', 'PLD4',
    # ],
    # 'higher_in_p_dc_than_as_dc': [ # 240619: currently unused
    #     'JCHAIN', 'GZMB', 'MZB1', 'SMPD3', 'NUCB2', 'CCDC186', 'DERL3',
    #     'COBLL1', 'CYP46A1', 'FAM129C', 'MAP1A', 'ZFAT', 'SELPLG', 'TPM2',
    #     'RASD1', 'PACSIN1', 'SLC15A4', 'SELENOS', 'MAPKAPK2', 'CUEDC1',
    #     'ERN1', 'CYSLTR1', 'HS3ST3B1', 'PDIA4', 'RRBP1', 'PARVB', 'MPEG1',
    #     'P2RY14', 'TCL1A', 'SLC35F3', 'MAN1A1', 'ETS1', 'NCF4', 'SORL1',
    #     'FCRLA', 'CMKLR1', 'NPC1',
    # ],
    # 'higher_in_mpp_than_clp': [ # 240619: currently unused
    #     'SERPINB1', 'PRSS57', 'MYB', 
    #     ###'AC084033.3', 
    #     'AVP', 'IL1B',
    #    ###'AC011139.1', 
    #    'SLC40A1', 
    #    ###'ANXA1', 
    #    'ZNF385D', 'SELL', 'CSF3R',
    #    ###'AL157895.1', 
    #    'CPA3', 'EREG', 
    #    ###'EEF1B2', 
    #    'C1QTNF4',

    #    'MGST1',
    #    ###'HSP90AB1', 
    #    'PLXDC2', 'FCER1A', 'MDK', 'IGLL1', 'GATA2', 'FAM117A',
    #    'ZBTB16', 'PPA1', 'CD84', 'SLC39A3', 
    #    ###'LDHB', 
    #    'PRKG2',
    # ],
    'higher_in_pre_b_than_clp_m_and_nktdp': [
        # >2
        'IGLL1', 'VPREB3', 'AKAP12', 'ARPP21', 'DNTT', 'CD79B', 'VPREB1',
       'CD9', 'CD24', 'RAG1', 'EBF1', 'SOCS2', 'LINC01013', 'HPS4', 'MYB',
       'TOP2B', 'TLE4', 'CMTM8', 'XBP1', 'HMHB1', 'CXCR4', 'RAG2', 'LEF1',
       'NEIL1', 'SYNE3', 'YBX3', 'PAX5', 'CD79A', 'ID3', 'ZCCHC7',
       'TP53INP1', 'MT1X', 'MDM2', 'AP001059.2', 'MLXIP', 'PTPRE',
       'DEPP1', 'CD81', 'SMAD1', 'LHPP', 'SPTBN1', 'MYLK', 'CD38', 'GLRX',
       'GAB1', 'ZEB2', 'PPP1R14B', 'LCN6', 'CMTM7', 'MPP1', 'MND1',
       'ISG20', 'TSPAN14', 'RCSD1', 'UHRF1', 'FOXO1', 'BACH2',
       'LINC00426', 'PAG1', 'ERG', 'PXDN', 'SNX2', 'GNG7', 'HHIP-AS1',
       'QRSL1', 'AP001059.3', 'TCF3', 'TMEM243', 'HHIP', 'SSBP2', 'P4HA2',
       'CALCOCO2', 'AEBP1', 'FAM241A', 'LINC02227', 'PSD3', 'TAPT1',
       'CHST15', 'SMARCD2', 'ALDH5A1', 'CYGB', 'LIG4', 'TMEM263',
       ###'AL161912.2', # 240129: was missing from /dummy/dummy/dummy/raid/mds/new_N280_exps/intermediate_output/raw_without_excluded_genes_and_cells.h5ad), and quite weak anyway.
       'GADD45A', 'RPS4Y2', 'POU2AF1', 
       ###'MME', 
       'HMCES',
       'SPTA1', 
       ###'MZB1', 
       'VDAC1', 'LAPTM5', 'KIAA0040', 'GBP4',
    ],
    'higher_in_bemp_than_clp_m_and_nktdp': [
        # >2
        'HDC', 'MS4A2', 'CNRIP1', 'TPSB2', 
        ###'LMO4', 
        'TPSAB1', 'GATA2', 'MS4A3', 'MYB', 'SLC40A1', 'HPGDS', 'CSF2RB', 'CASP3', 'CPA3', 'ALOX5AP', 'CD84', 'FCER1A', 'APOC1', 
        ###'AC084033.3', 
        'AL157895.1', 'SRGN', 'HPGD', 'PRSS57', 'IKZF2', 'CLC', 'PDZD8', 'ZEB2', 'BACE2', 'CNKSR3', 'HBS1L', 'CTNNBL1', 'ALDH1A1', 'SERPINB1', 'NLK', 'TAL1', 'CLU', 'FBXO7', 'KLF1', 'IL1B', 'FAM83F', 'SNX5', 'PLIN2', 'KCNH2', 'GATA1', 'GRAP2', 'RYR3', 
        ###'CD63',
    ],
    # 240403: the following "unclear" lists were compiled while trying to understand N403 unassigned cells. gave up when it seemed like these cells are specific to N403, and she has both del(5q) and del(7)
    # 'unclear': [
    #     'UGCG', 'LBH', 'TCF4', 'AXL', 'TNS3', 'RP2', 'DNASE1L3', 'POLR2H',
    #    'UMAD1', 'ZNF333', 'GPI', 'CEP295', 'METTL8', 'TSGA10', 'SEC31B',
    #    'ANKRD33B', 'TMEM42', 'RBBP8', 'GCC2', 'JAKMIP2', 'ERGIC1',
    #    ###'AJ009632.2', 
    #    'DHX37', 'MANF',
    # ],
    # 'unclear2': [
    #     'IFI6', 'PPP6R1', 'HIST2H2AA4', 'ITGB1BP1', 'IFIH1', 'THOC3',
    #    'RASSF4', 'GADD45B', 'OSBPL9', 'IRAK3', 'ZBTB33', 'DBT', 'OAF',
    #    'BTK', 
    #    ###'AC118549.1', 
    #    'ATG4C', 'C22orf39', 'ATP6V0A1', 'BANP',
    #    'TOP1', 'TMLHE', 'SIL1', 'TK2', 'ENHO', 'ADAM8', 'IFITM1', 'GSDMD',
    #    'TRAF7', 'MOCS2', 'CYP20A1', 'HIGD1A', 'TNFSF10', 'AOAH', 'DHRS9',
    #    'NAIP', 'RNASE2', 'ARID3A', 'FAM111A', 'CA2',
    # ],
    # 'unclear3': [
    #     'IL1B', 'TMEM219', 'IRAK3', 'BROX', 'BTK', 'OPA1', 'ZDHHC20',
    #     'PANK3', 'HIST2H2AA4', 'ZBTB33', 'KDM2B', 'ARFGEF1', 'SCYL2',
    #     'FNDC3A', 'IER5', 'TMEM9B', 'TBP', 'RAP2A', 'GMPR2', 'CYP20A1',
    #     'COL4A3BP', 'PNPT1', 'IVNS1ABP', 'UHRF1BP1L', 'PSMB4', 'LAS1L',
    #     'CAPRIN1', 'TADA3', 'CDC123', 'DCK', 'AGO2', 'ARL8A', 'MOCS2',
    #     'NAIP', 'RAB1B', 'CSF2RB', 'FLYWCH2', 'PDAP1', 'PPP6R1', 'CRYBG1',
    #     'RFC1', 'TDP2', 'SNX30', 'TP53INP1', 'PLBD1', 'HACD4', 'SRSF10',
    #     'NRBF2', 'ECPAS', 'ZMYM1', 'TMC6', 'PTAR1', 'GSTK1', 'UBA6',
    #     'MGST2', 'BCAT2', 'TOR1A', 'MFSD14B', 'STX10', 'RALGAPA2',
    #     'AC118549.1', 'TM7SF3', 'NRBP1', 'IFIH1', 'ALOX5', 'THOC3', 'BFAR',
    #     'ENTPD4', 'COPS8', 'ERCC6L2', 'DHX38', 'HNRNPUL2', 'PGM2L1',
    #     'POLR3D', 'TK2', 'TOP1', 'JAK3', 'HACD3', 'OAF', 'CUL1', 'NFYC',
    #     'PYROXD1', 'TULP4', 'PARVB', 'MUS81', 'FAM43A', 'FAM172A',
    #     'PRKACA', 'ATP10D', 'DPY19L4', 'BIN1', 'BRD8', 'MED29',
    # ],
    # 'unclear4': [
    #     'DCHS1', 'EXPH5', 'SLC25A17', 'MCF2L', 'SDR42E2', 'DUSP27', 'GPC6',
    #    'EFNB3', 'LINC00900', 
    #    ###'AC006581.2', 
    #    'PTGR1', 'LIPE-AS1', 'CLSTN2',
    #    'CSKMT', 'FASTKD1', 'SOCS4', 'CLK2', 
    #    ###'AL138899.1', 
    #    'TAF1B',
    #    'TBC1D7', 'RRAGD', 'SRF', 
    #    ###'AC105446.1', 
    #    'SMIM11A', 'SCAF4',
    #    'TBC1D23', 
    #    'AL353708.1', 
    #    'TMEM240', 'AACS', 'B3GLCT', 'KAT14',
    #    'NEK9', 
    #    ###'AC245060.2', 
    #    'GDI1', 'FLG-AS1', 'GFM1', 'PSMD5',
    #    ###'AL121832.2', 'AC246817.2', 
    #    'CRYGS',

    #    'COL1A2', 'AK4', 
    #    ###'AC106886.5',
    #    'ZNF585A', 'KCTD2', 'HEATR6', 'AC239868.2', 'TREM2', 'CA13',
    #    'DZIP1', 'RPS6KA6', 'AC006547.3', 'NBPF9', 'MIS12', 'AC093585.1',
    #    'SERPINH1', 'ZMYM3', 'NUP85', 'C1RL', 'NRAS', 'RWDD3',
    #    'AL603839.3', 'ANGPTL4', 'TMPPE', 'TRIM46', 'DDX31', 'SESN2',
    #    'BTBD10', 'OLFM4', 'AC138956.1', 'ARHGAP19-SLIT1', 'SEC31B',
    #    'PCAT5', 'TMEM86A', 'MED12', 'HSPBAP1', 'HIST1H2BL', 'GPRC5A',
    #    'AC079354.3', 'CYP2B6', 'SH3D19', 'AC007009.1', 'AC091053.1',
    #    'C17orf49', 'ALDH9A1', 'FOCAD', 'TNFSF11', 'AC024257.5', 'DWORF',
    #    'BICC1', 'ASMTL-AS1', 'MMGT1', 'GRAMD1B', 'DFFB', 'GPR1', 'CNNM1',
    # ],
    # 'unclear5': [
    #     'RRAGD', 'ERGIC1', 'DCHS1', 'AJ009632.2', 'NBPF9', 'RP2', 'PPP6R2',
    #     'HIF1A', 'PTCH2', 'TPM2', 'RBBP8', 'AEBP2', 'TNRC18', 'HOXA10',
    #     'FAIM', 'CENPK', 'TRRAP', 'TTC28', 'SMC4', 'ANAPC1', 'KIAA1586',
    #     'MAF', 'MAP3K3', 'SLC25A17', 'CEP295', 'TIGD1', 'GORASP2',
    #     'AL360012.1', 'SH3D19', 'NMNAT3', 'GPC6', 'AL138899.1', 'NEIL1',
    #     'AL355816.2', 'GDI1', 'MBNL1-AS1', 'EXPH5', 'GTF3C1', 'DFFB',
    #     'KIAA0586', 'SDR42E2', 'JAG1', 'SNX4', 'MCF2L', 'PDXDC1', 'SCAF4',
    #     'CLK2', 'TMEM42', 'DACH1', 'C1orf54', 'TNS3', 'TBC1D23', 'TMEM69',
    #     'CTNNB1', 'TCF4', 'GSTO2', 'SCAPER', 'TBC1D7', 'SLC16A9', 'ZNF493',
    #     'KCTD2', 'EFNB3', 'HCFC2', 'UBE2G1', 'DUSP27', 'SRF',
    # ],
    # 'unclear6': [
    #     'BTK', 
    #     ###'RPS4Y1', 
    #     'TAX1BP1', 'DCK', 'TMEM219', 'CDC123', 'PSMB4',
    #    'ZCCHC10', 'BIN1', 'CRYBG1', 'S100A8', 'SEC11C', 'HACD4', 'BROX',
    #    'NRBF2', 'TADA3', 'UBE2D2', 'UBA6', 'NRBP1', 'RGS19', 'MPP1',
    #    'DNMT3A', 'NCBP2-AS2', 'PPDPF', 'GSTK1', 'PDAP1', 'COL4A3BP',
    #    'ZNF655', 'DHX38', 'GUK1', 'CUL1', 'ANKRD36C', 'CAPRIN1', 'PANK3',
    #    'PGLS', 'TMED8', 'ARL8A', 'CSTA', 'TMEM9B', 'FES', 'KDM1B',
    #    'PPP6R1', 'HSPA4', 'DESI2', 'IFI6', 'CCNK', 'ARFGEF1',
    #    'HIST2H2AA4', 'MAP4', 'IFI44L', 'IGKC', 'ZCCHC8', 'POLR3D',
    #    'ARHGAP26', 'TYMP', 'XAF1', 'TP53INP1', 'OPA1', 'ISCU', 'GMPR2',
    #    'IRAK3', 'UBE2J2', 'IER5', 'RAB1B', 'ATP5F1D', 'CISD3', 'IFIH1',
    #    'OAS2', 'PTAR1',
    # ],
    # 'unclear7_somewhat_n403_specific': [
    #     'DPYSL2', 'SDR42E2', 'DCHS1', 'LIPE-AS1', 'PTCH2', 'PKIB',
    #     'AC245014.3', 'EXPH5', 'PDXDC1', 'MCF2L', 'TMEM107', 'STK38L',
    #     'GPC6', 'EFNB3', 'FZD2', 'AACS', 'EBF4', 'CSKMT', 'FASTKD1',
    #     'SLC25A17', 'AC006581.2', 'SNHG25', 'AC010632.1', 'CS', 'DDA1',
    #     'SOCS4', 'ZMAT3', 'COL1A2', 'TPM2', 'B3GNT7', 'DDX31', 'PPP6R2',
    #     'AL138899.1', 'P4HA1', 'CNNM4', 'FBXO46', 'TRA2B', 'NRAS',
    #     'LRRC69', 'SEC31B', 'RBMS2', 'KAT14', 'AC105446.1', 'AC103591.3',
    #     'DHX37', 'SMIM11A', 'GFM1', 'RIBC2', 'AL353708.1', 'PSMD5',
    #     'TMEM240', 'PTPN2', 'MED4', 'STRN3', 'ZNF197', 'MAFIP',
    # ],
    # 'higher_in_n403_unassigned': [
    #     'SNHG25', 
    #     ###'XIST', 
    #     'TCEAL8', 'BEX1', 'HIF1A', 'TOP2B', 'HNRNPH3',
    #     'SPINK2', 'AMD1', 'ERGIC1', 'C8orf33', 'UBAP2L', 'NUCB2', 'RUNX1',
    #     'SLC25A36', 'LLPH', 'NDUFAF8', 'TMEM107', 'GALNT1', 'CLEC11A',
    #     'LYRM2', 'CIR1', 'RABGAP1', 'CD1E', 'IL6R', 'YBX3', 'WNK1',
    #     'PPP6R2', 'PTCH2', 'CORO1C', 'GALNT7', 'BEX3', 'LRPPRC', 'ATP2B4',
    #     'SMIM24', 'GPBP1L1', 'MYC', 'RC3H1',
    # ],
    # 'higher_in_spink2_bemp_than_bemp': [
    #     # 240403: only a single MC now, and the diff exp between it and other BEMPs is not very strong, so moving on.

    #     # 'SPINK2', 'PLAC8', 

    #     'SPINK2', 'HLA-DRA', 'PLAC8', 'HLA-DPA1', 'RNASE2', 'PPA1',
    #     'NAP1L1', 'PCLAF', 'CD74', 'CEBPA', 'RUNX1', 'CSF3R', 'EGFL7',
    #     'H2AFY', 'GABPB1-AS1', 'SNRPD1', 'EID1', 'SLC26A2', 'SLC39A3',
    #     'MIA3', 'C1orf35', 'LMO2', 'MTF1', 'PDCD5', 'NPR3',
    # ],
    # 'lower_in_spink2_bemp_than_bemp': [
    #     'MS4A2', 'TPSAB1', 'TPSB2', 'HPGDS', 'KIT', 'VWA5A', 'FCER1G',
    #     'ACSL4', 'HPGD', 'CNRIP1', 'PAG1', 'ALOX5AP', 'AL157895.1', 
    #     ###'CPA3',
    #     'SLFN5', 'FKBP5', 'NMT2', 'FYB1', 'GRAP2', 
    #     ###'ITM2B', 
    #     'BACH1',
    #     'PLIN2', 'UFL1', 'MIDN', 'CLU', 'YTHDC1', 'HBS1L',
    # ],
    'higher_in_clp_m_and_nktdp_than_mpp': [
        # >2
        
        ###'JCHAIN', 
        'MEF2C', 'ACY3', 'LTB', 'CYTH4', 'RUNX2', 'CLNK',
        'CORO1A', 'ID2', 'LCP1', 'NKG7', 
        ###'IGHM', 
        'TRIB2', 'CCR7', 'TYROBP',
        'DDIT4', 'TRAF3IP3', 'MEF2A', 'ITGAL', 'SATB1', 'HCST', 'CD53',
        'AL096865.1', 'SPECC1', 'LMO4', 'SLC38A1', 'HLA-DMB', 'LSP1',
        'HLX', 'TCF4', 'STK17A', 'AFF3', 'CD74',
        # # 1.8-2
        # 'TNFRSF4', 'ITM2C',
        # 'TMSB4X', 'MAF', 'EVL', 'LINC00865', 'INTS6', 'ANKRD44',
        # 'TNFRSF18', 'RAB31', 'ZFP36L1', 'RASD1', 'PIK3AP1', 'NEGR1',
        # # 1.5-1.8
        # 'STK17B', 'ATM', 'CYFIP2', 'PDE7A', 'GNAS', 'AKNA', 'PALLD',
        # 'RNASET2', 'PLP2', 'PREX1', 'MED13L', 'CTSZ', 'REL', 'TRG-AS1',
        # 'ANKRD11', 'SPNS3', 'RUNX3', 'HLA-DRB1', 'NUDT8', 'IGHD', 'CARD11',
        # 'SEPT9', 'TXNIP', 'JAML', 'IQGAP1', 'MDFIC', 'SP110', 'RASSF2',
        # 'CALM1', 'GSN',
    ],
    'higher_in_monocyte_than_gmp_l': [
        'S100A9', 'S100A8', 'FCN1', 'VCAN', 'TYROBP', 'MNDA', 'CYBB',
        'FCER1G', 'FGL2', 'S100A12', 'SERPINA1', 'CD14', 'MPEG1', 'NCF1',
        'MS4A6A', 'COTL1', 'TNFRSF1B', 'GIMAP4', 'TYMP', 'NCF2', 'FYB1',
        'CD36', 'JAML', 
        ###'AC020656.1', 
        'LGALS3', 'IGSF6',
    ],
    # 'higher_in_cfd_tryptase_gmp_l_than_gmp_l': [
    #     'CFD', 'TPM4', 'S100A10', 'TPSB2', 'TPSAB1', 'S100A6', 'CD9',
    #    'TRGC2', 'EMP3', 'KLF2', 'SUPT4H1', 'CD63', 'MYL12A', 'ANXA1',
    #    ###'GSN', # looked meh in leave_one_out.
    #    'ANXA2', 'ATP2B1', 'LAPTM5', 
    #    ###'H1F0', # looked meh in leave_one_out.
    #    ###'HIST1H2AC', # looked meh in leave_one_out.
    #    ###'TPST2', # looked meh in leave_one_out.
    #    'EMP1', 'TAGLN2',
    
    # #     'S100A10', 'S100A8', 'S100A9', 'FCN1', 'TYROBP', 'S100A11',
    # #    'FCER1G', 'TMEM176B', 'KLF2', 'VCAN', 'CST3', 'LGALS3', 'TPSB2',
    # #    'CTSS', 'LGALS2',
    # ],
    # 'lower_in_cfd_tryptase_gmp_l_than_gmp_l': [ # 240619: currently unused
    #     # according to g_12_09_23_1_a cells
    #     'MPO', 'PRTN3', 'IGLL1', 'RBM3', 'NPW', 'TYMS', 'HIGD2A', 
    #     ###'ANP32E', # looked meh in leave_one_out.
    #    'NASP', 'EIF3A', 'CPA3', 'IGFBP2', 'TENT5A', 'PLAC8', 'MRPS26',
    #    ###'SNHG8', # looked meh in leave_one_out.
    # ],
    # 'higher_in_g_12_09_23_1_a_gmp_l_than_gmp_l': [ # 240619: currently unused
    #     # g_12_09_23_1_a_mc_names: ['M452.81', 'M482.75', 'M1396.63']
    #     'TPSB2', 'CFD', 'DLK1', 'APOE', 'TPSAB1', 'GATA2', 'CD9', 'SPP1',
    #     'S100A10', 'CHI3L1', 'CD63', 'MMP9', 'PDLIM1', 'S100A6', 'PRSS21',
    #     'ZNF385D', 'PRSS3',
        
    #     # 'TPSB2', 'CFD', 'ELANE', 'APOE', 'DLK1', 'AZU1', 'GATA2', 'SPP1',
    #     # 'LYZ', 'TPSAB1', 'CD9', 'PRSS21', 'PRTN3', 'CST7', 'CHI3L1',
    #     # 'MMP9', 'S100A10',
    # ],
    # 'higher_in_weird_gmp_l_than_gmp_l': [ # 240619: currently unused
    #     # gu_mc_mask = (mc_ad.obs['state'] == 'GMP-L') & (mc_utils.get_genes_expr(mc_ad, 'MPO') < -11) & (mc_utils.get_genes_expr(mc_ad, 'TPSB2') > -15)
    #     # gu_mc_mask = (mc_ad.obs['state'] == 'GMP-L') & (mc_utils.get_genes_expr(mc_ad, 'MPO') < -11)
    #     # gu_mc_names = mc_ad.obs_names[gu_mc_mask]
    #     # gu_mc_names
    #     'TPSB2', 'GATA2', 'DLK1', 'APOE', 'TPSAB1', 'SPP1', 'CD9',
    #     'ZNF385D', 'PDLIM1', 'CHI3L1', 'UCHL1', 'CYTL1', 'HLA-DQA1',
    #     'PRSS21', 'MMP9', 'S100A10', 'RHEX', 'CD52', 'SLC39A3',
    # ],
    'higher_in_gmp_l_than_all_other_hspcs': [
        'MPO', 'AZU1', 'ELANE', 
        ###'IGLL1', # a bit high in MPP. after removing IGLL1, i realized the gene sig became higher_in_gmp_l_than_all_other_hspcs rather than higher_in_gmp_l_than_all_other_myeloid_hspcs, as now pre-B and pro-B were also very low in it.
        'PRTN3', 'CFD', 
        ###'C1QTNF4', # a bit high in MPP
       'RNASE2', 'CTSG', 'MGST1', 'CEBPD', 'RAB32', 'TNFSF13B', 'CLEC11A',
       'CEBPA', 'LYST', 'PLPPR3', 'HGF',

        # diff exp between cells high and low in genes above. seems like these are not helpful. 
        # 'PLAC8', 'PRAM1', 'PKM', 'NPW', 'XBP1',
    ],
    'mpp_slc40a1_sig': [
        
        'SLC40A1', 
        ###'AL157895.1', 
        'FCER1A', 'ZBTB20', 'PBX1', 'TESPA1',
       'CPPED1', 'ZNF385D', 'GATA2', 'SESN3', 'PRKACB', 'ATP6V0A2',
       'SNX5', 'ZEB2', 'PDZD8', 
       ###'CYTL1', 
       ###'RNF130', # looked somewhat bad in marker heatmap
       'TNFAIP8', 
       ###'LIMS1', # looked somewhat bad in marker heatmap
       ###'MARCKSL1', # looked somewhat bad in marker heatmap
       'MED12L', 'PLEK', 'PLXDC2', 
       ###'ETS2', # looked somewhat bad in marker heatmap
       'BMP2K', 'CTNNBL1',
       'CASP3', 'CD84', 
       ###'HIST1H2AC', # looked somewhat bad in marker heatmap
       ###'SOD2', # looked somewhat bad in marker heatmap
       'NLK', 
       ###'RGS18', # looked somewhat bad in marker heatmap
       'ZMYND8',
       ###'MYCT1', # looked somewhat bad in marker heatmap
       'MYB', 'SLC12A6', 'MINPP1', 
       ###'CLU', # looked somewhat bad in marker heatmap
       'SLC10A5', 'ZBTB16',
       ###'HIST1H2BF', # looked somewhat bad in marker heatmap 
       ###'MPP1', # looked somewhat bad in marker heatmap
       'EIF2AK1', 'TIMP3',
        
        
        # 'SLC40A1', 'FCER1A', 'PDZD8', 'GATA2', 'PLEK', 'MINPP1',
        # 'ATP6V0A2', 'BMP2K', 'CD84', 
        # # 'TLN1', # strong
        # 'CPPED1', 'AL157895.1',
        # 'CNRIP1', 'ZMYND8', 'TESPA1', 'PBX1', 'TIMP3', 'RNF130', 'PSTPIP2',
        # 'CTNNBL1', 'CYTL1', 'ZNF385D', 
        # ###'AL034397.3', 
        # 'RGS18', 'SOD2',
        # 'HBD', 'SNX5', 'MARCKSL1', 'ZEB2', 'GATA1', 'CASP3', 'CSF2RB',
        # 'ZBTB20', 'TNFAIP8', 'SLC12A6', 'FREM1', 'KLF1', 'NLK', 'PKIG',
        # 'ABCC4', 'LIMS1', 'TAL1', 'STAM', 
        # ###'FTH1', # too strong
        # 'PRKACB', 'TFRC',
        # ###'AC011139.1', 
        # 'ZBTB16', 'STXBP5', 'PLXDC2', 'DEPTOR', 'SLC10A5',
        # 'CPA3', 'ETS2', 'GMPR', 'IGSF10', 'GALNT1', 'RIPOR3', 'C2orf88',
        # 'MYB', 'GNAQ', 'HPGDS', 'TFR2', 'MPP1', 'MED12L', 'GDF11', 'MYCT1',
        # 'HIST1H2AC', 'PLXNC1', 'CLU', 'TRIM24', 'SESN3', 'EIF2AK1',
        # 'HIST1H2BF',
    ],
    # 'gmp_l_fam30a_sig': [
    #     # gmp_mask = mc_ad.obs['state'].isin(['GMP-L']) & (mc_ad.obs['max_donor_id_frac'] < 0.35) & ~mc_ad.obs_names.isin(['M382.96'])
    #     # iter2: iter1 filtered by marker heatmap across mds ult gmp_mask, see above.
    #     'HOPX','CD52','FAM30A','IQGAP2','GUCY1B1','ALDH1A1','HTR1F','LAPTM4B','GCSAML','DRAP1','PTMS','TSTD1','ADGRG1','HEMGN','TPM1',

    #     # # gmp_mask = mc_ad.obs['state'].isin(['GMP-L']) & (mc_ad.obs['max_donor_id_frac'] < 0.35) & ~mc_ad.obs_names.isin(['M382.96'])
    #     # # iter1: corr>0.6 with ['FAM30A', 'CD52'] across mds illu gmp_mask, see above, and not in gmp_l_ahnak_sig.
    #     # 'CD52', 'FAM30A', 'ANGPT1', 'FTH1', 'HOPX', 'AKR1C3', 'TSTD1',
    #     # 'ICAM3', 'ITM2B', 'TSC22D1', 'RBPMS', 'YPEL3', 'MDK', 'GIMAP7',
    #     # 'RCSD1', 'FHL1', 'CD99', 'ALDH1A1', 'CAVIN1', 'HEMGN', 'SPINK2',
    #     # 'GUCY1A1', 'INKA1', 'PTMS', 'DNAJB6', 'CLIC1', 'FRMD4B', 'IGHM',
    #     # 'SSBP2', 'LIMD2', 'CD37', 'CLEC2B', 'RAB27B', 'DDIT4', 'CD34',
    #     # 'RHOC', 'EGFL7', 'SYPL1', 'BEX4', 'LPIN2', 'HTR1F', 'GCSAML',
    #     # 'ETS2', 'H2AFY', 'SOCS2', 'TM4SF1', 'FXYD5', 'PCMTD1', 'DRAP1',
    #     # 'CDKN1C', 'ZSCAN18', 'LINC02573', 'CXXC5', 'ADGRG1', 'NFKBIA',
    #     # 'BIN1', 'AC008105.3', 'AC002454.1', 'CIRBP', 'IQGAP2', 'ESD',
    #     # 'TIMP1', 'SARAF', 'HEBP1', 'SPINT2', 'TPM1', 'LAPTM4B', 'TALDO1',
    #     # 'TCTEX1D1', 'IL1B', 'PNRC1', 'GUCY1B1', 'ARF5', 'TFPI',
    # ],
    'gmp_l_ahnak_sig': [
        # gmp_mask = mc_ad.obs['state'].isin(['GMP-L']) & (mc_ad.obs['max_donor_id_frac'] < 0.35) & ~mc_ad.obs_names.isin(['M382.96'])
        # iter2: iter1 filtered by marker heatmap across mds ult gmp_mask, see above.
        'AHNAK', # strong, but looked good in leave_one_out gene gene scatter
        'S100A10',
        'ANXA1', # strong, but looked good enough in leave_one_out gene gene scatter
        ###'VIM', # very strong, and looked meh in leave_one_out gene gene scatter
        'CP','MAMDC2','TSC22D3','LMNA','UCP2',
        ###'SH3BGRL3', # strong, and looked meh in leave_one_out gene gene scatter
        'EMP3','TSPAN2','TUBA1A',
        ###'CLIC1', # strong, and looked meh in leave_one_out gene gene scatter
        ###'CD37', # looked meh in leave_one_out gene gene scatter
        'SCN9A','PRMT2','CA5B','SEPT11',
        ###'DNAJB6', # looked somewhat meh in leave_one_out gene gene scatter, and is a chaperone
        'CAMK1D',
        ###'CD99', # strong, and looked meh in leave_one_out gene gene scatter
        'EMP1','MYADM',
        'TAGLN2', # strong, but looked good in leave_one_out gene gene scatter
        'ATP2B1','IDS','CALCOCO2','PBXIP1','DDAH2','MTURN',
        # added manually later by corr>0.8 with iter2 across mds ult gmp_mask, see above.
        ###'ANKRD28', looked meh in leave_one_out gene gene scatter of mds illu
        'PTPRC',
        
        # # gmp_mask = mc_ad.obs['state'].isin(['GMP-L']) & (mc_ad.obs['max_donor_id_frac'] < 0.35) & ~mc_ad.obs_names.isin(['M382.96'])
        # # iter1: corr>0.6 with ['S100A10', 'ANXA1'] across mds illu gmp_mask, see above.
        # 'ANXA1', 'S100A10', 'CRIP1', 'VIM', 'TAGLN2', 'SH3BGRL3', 'TSPAN2',
        # 'TUBA1A', 'SEPT11', 'MAMDC2', 'DDAH2', 'AHNAK', 'TSC22D3',
        # 'PBXIP1', 'LMNA', 'CP', 'MTURN', 'CD99', 'CAMK1D', 'DNAJB6',
        # 'FXYD5', 'IDS', 'CALCOCO2', 'EMP1', 'SCN9A', 'CLIC1', 'RABAC1',
        # 'ATP2B1', 'CD37', 'HPGD', 'MYADM', 'TPPP3', 'CEP126', 'UCP2',
        # 'PRMT2', 'EMP3', 'CDKN1C', 'SOCS3', 'SARAF', 'ZSCAN18', 'FRMD4B',
        # 'CA5B', 'SNX2', 'EFHC1', 'GLIPR2', 'C9orf16', 'ARF5',
    ],
    'gmp_l_elane_sig': [
        
        'PRTN3', 'ELANE', 'AZU1', 'CST7', 'CFD', 'MS4A3', 'LYZ', 'P4HB',
        'CTSG', 'PDIA3', 'MNDA', 'RNASE2', 
        ###'CALR', # strong, and looked meh in leave_one_out gene gene scatter
        'PLPPR3', 'CANX',
        'PDIA6', 'CSTA', 'GPI', 'SERPINB10', 'RNASE3', 'RPN1', 
        ###'SDF2L1', # didn't look good in marker heatmap
        'FAM45A', 'FAM107B', 'SSR4', 'PPIB', 'HMGB2', 'PDIA4', 
        ###'SEC61B', # didn't look good in marker heatmap
        'CEBPD', 'LYST', 'KBTBD11', 
        ###'POLR2L', # didn't look good in marker heatmap
        'NPW', 'MARCKS', 'EMILIN2',
        ###'FKBP2', # didn't look good in marker heatmap
          'STT3A', 'RFLNB', 'FLNA',
        
        
        # 'ELANE', 'CTSG', 
        # ###'SRGN', # too strong
        # 'MS4A3', 'PRTN3', 'P4HB', 'MNDA', 'PDIA3',
        # 'CANX', 'SSR4', 'CST7', 'RNASE2', 'AZU1', 'CEBPD', 'EMILIN2',
        # 'CLEC11A', 
        # 'RAB32', # strong, and looked meh in leave_one_out gene gene scatter
        # 'CFD', 
        # 'HOMER3', 'SERPINB10', 'CSTA', 
        # ###'MPO', # too strong
        # 'PLPPR3', 'FAM45A', 'SMIM20', 'F13A1', 'LYST', 
        # ###'HSP90B1', # strong, and looked meh in leave_one_out gene gene scatter
        # ###'HSPA5', # strong, and looked meh in leave_one_out gene gene scatter
        # 'LYZ', 'PDIA4', 'RGCC', 'KBTBD11', 'AURKAIP1', 'RNASE3', 'MGST1',
        # 'MARCKS', 'NDRG1', 'CTSD', 
        # 'CALR', # strong, and looked meh in leave_one_out gene gene scatter
        # 'RPS17', 'TESC', 'RFLNB',
        # 'MFSD10', 'CEBPA', 'HMGB2', 'KIAA0930', 'PHPT1', 'AC020656.1',
        # 'RPN1', 'FLNA', 'GPI', 'PIK3AP1', 'LINC00926', 'PPIB', 'FNDC3B',
        # 'PDIA6', 'HGF', 'FKBP2', 'EDEM3', 'IFT57', 'FAM107B', 'POLR2L',
        # 'SND1', 'MIF4GD', 'SDF2L1', 'SEC61B', 'GSDMD', 'MAPK12', 'STT3A',
        # 'FOSB', 'P2RX1', 'NPW', 'POM121', 'LARP1B', 'HRH2',
    ],
    'nktdp_lmna_sig': [
        # iter2: iter1 filtered by heatmap of NKTDP MCs.
        'LMNA','AHNAK','TAGLN2','RGCC','EMP1','MYADM','IDS','ATP2B1',
        ###'S100A10', # leave one out gene gene scatter looks meh.
        ###'ANXA1', # leave one out gene gene scatter looks meh.
        
        
        # 'LMNA', 'EMP1', 'TAGLN2', 'EMP3', 'AHNAK', 'GAS7', 'MYADM',
        # 'CRIP2', 'IQGAP1', 'IDS', 'SFXN3', 'RGCC', 'ANXA1', 'SERPINE1',
        # 'TMEM43', 'CRIP1', 'CDC42', 'LRRFIP1', 'CD9', 'DPF1', 'PLEC',
        # 'AAED1', 'ANTXR2', 'EZR', 'ACOT7', 'TSPAN2', 'JUND', 'CACNB4',
        # 'S100A10', 'ATP2B1', 'CA5B', 'ADGRE5', 'SUN1', 
        
        # 'CALM1', 'RFX2',
        # 'HMGCS1', 'DBI', 'TPPP3', 'CD44', 'KIAA1841', 
        # ###'AC020916.1', 
        # 'CKLF',
        # 'TUBB6', 'PRMT2', 'ARHGEF40', 'CD99', 'RTN4', 'CEP126', 'CD58',
        # 'ADAP1', 'MYCBP2', 'CABIN1', 'MFSD6', 'CAST', 'CAMK1D', 'AGPAT4',
        # 'CSPP1', 'TUBB4B', 'GLUD1', 'JADE2', 'S100A6', 'WDR49', 'SVIL',
        # 'NPTX2', 'DENND4A', 'ARL4A', 'WDR1', 'MTURN', 'KLF6', 'ATP1B1',
        # 'PBXIP1', 'FOSB',
    ],
    'nktdp_maf_sig': [
        # iter3: best of iter2, then manually adding to them genes corr with best of iter2
        'MAF','TNFRSF18','TNFRSF4','TRIB2','RORC','IL17RE',
        
        'PTH1R','STK17A','CLEC2D','DNMT3B','TOX2',
        # manually added
        'TLE3', 'CD7',
        'CASC15',
        'MLC1',
        'THRA',
        'EMID1',
        #'DDIT4',
        'TCF7', # a specific MC is off...
        'MAML3',


        
        # iter2: iter1 filtered by heatmap of NKTDP MCs.
        ###'SATB1', # looks like SATB1 leaves the rest when it is high...
        # 'MAF','TNFRSF18','TNFRSF4','TRIB2','RORC','IL17RE','KLRB1','PTH1R','STK17A','CLEC2D','DNMT3B','TOX2',
        

        # 'ZNF3','LINC00672','ZNF217',
        # ###'C1QTNF3', # leave one out gene gene scatter looks bad.
        # 'AQP3','RASGRP1','ZNF493','METTL7A','SPOCK2','FURIN','SAMD9L','MAML3','RB1CC1','COMMD7',
        # ###'AC243960.1',
        # 'RARRES3','C15orf39','ZNF827','DYNLT3','SSBP3','NFE2L2','AAK1','PHIP','AKT3','FYTTD1','FAM160B1',
        # ###'ARL6IP5', # somewhat strong, and looked meh in leave_one_out
        # 'CSAD','BIN2','ZMIZ1',
        # ###'PNRC1', # somewhat strong, and looked meh in leave_one_out
        # 'CHD4',
        # ###'ABCA2','KCTD9','UNC13C','SAMD4B','CGGBP1','TRIM56','VPS4B', # leave one out gene gene scatter looks meh.
        # ###'VGLL4', # looked meh in leave_one_out
        # ###'IER2', # looked meh in leave_one_out
        
        # ###'H3F3B', # strong, and looked meh in leave_one_out
        # ###'H3F3A', # strong, and looked meh in leave_one_out
        # ###'HMGB1', # strong, and looked meh in leave_one_out
        # 'STK17B','LCK','SPSB3','MAML2','CD46','TNRC6C','TNFSF4','MIR1-1HG-AS1','NOTCH1',
        # ###'XBP1', # leave one out gene gene scatter looks meh.
        # 'IQGAP2','TRG-AS1','FRMD4B','CASC15','GOLGA8B','ABLIM1','DPPA4','MLC1','IL7R','SPRY1','TRBC1','CITED2','SLC5A3','SCML4','MRPS6','TLE3','EMID1','THRA','SLC9A3R1',


        
        # 'KLRB1', 'MAF', 'TRIB2', 'ST3GAL1', 'PRL', 'EMID1', 'DNMT3B',
        # 'CASC15', 'THRA', 'COMMD7', 'TOX2', 'ARL6IP5', 'TNFRSF18', 'SPRY1',
        # 'SCML4', 'PCDH9', 'CLEC2D', 'RFLNB', 'FRMD4B', 'HMGB1', 'DPPA4',
        # 'PNRC1', 'DYRK2', 'NFE2L2', 'CITED2', 'RORC', 'TRBC2', 'ETS2',
        # 'ZNRF1', 'TRBC1', 'AC243960.1', 'TNFSF4', 'IL17RE', 'SLC5A3',
        # 'FAT3', 'ABLIM1', 'RARRES3', 'RB1CC1', 'SPOCK2', 'UNC13C',
        # 'SAMD4B', 'H3F3A', 'H3F3B', 'MRPS6', 'NFIA', 'BIN2', 'CHD4',
        # 'FYTTD1', 'SSBP3', 'STK17B', 'SPSB3', 'C15orf39', 'MIR1-1HG-AS1',
        # 'MLC1', 'MAML3', 'DYNLT3', 'LCK', 'ZNF827', 'METTL7A', 'TLE3',
        # 'PHIP', 'TRIM56', 'GOLGA8B', 'XBP1', 'FURIN', 'SATB1', 'TRG-AS1',
        # 'MAML2', 'PTH1R', 'IQGAP2', 'NUCB2', 'APOL6', 'IL7R', 'AKT3',
        # 'AC245060.5', 'TNRC6C', 'CGGBP1', 'TNFRSF4', 'CSAD', 'TCF7L2',
        # 'SAMD9L', 'RASGRP1', 'LTB', 'NOTCH1', 'VPS4B', 'ZNF217',
        # 'FAM160B1', 'ZNF3', 'HMGB2', 'ABCA2', 'OSBPL3', 'ZNF493',
        # 'SLC9A3R1', 'ZMIZ1', 'STK17A', 'PLCB1', 'KCTD9', 'TSC22D1', 'AQP3',
        # 'IGHM', 'IER2', 'AAK1', 'C1QTNF3', 'VGLL4', 'CD46', 'LINC00672',
    ],
    # 'nktdp_nfil3_sig': [
    #     # NOTE: ugh. plotting a marker heatmap - this isn't really a sig... too heterogenous...
    #     # NOTE: somewhat correlated with nktdp_maf_sig
    #     # started from RUNX3, ACY3, TFRC, CALM2, and found that in the sig i get, TESC is the best in leave one out across cells.
    
    #     'TESC', 
    #     'ATP10A', 
    #     ###'CALM2', # very strong, and leave one out gene gene scatter doesn't look good enough to justify adding it.
    #     'ELL2', 
    #     'NFIL3', # best?
    #     ###'ACY3', # strong, and leave one out gene gene scatter looks meh
    #     'CCDC141',
    #     'SEMA4D', 'FYB1', 'ZDHHC14', 'TFRC', 'PTK2B', 'CYSLTR1', 'CD3E',
    #     'DDIT4', 'SYT2', 'CYTH4', 'CTSW', 
    #     ###'LSP1', # strong, and leave one out gene gene scatter looks meh
    #     'FGR', 'MEF2A',
    #     ###'PIK3AP1', # leave one out gene gene scatter looks meh
    #     'TTN', 'POLB', 
    #     ###'TUT4', # a bit strong, and leave one out gene gene scatter looks meh
    #     ###'DYNC1LI2', # leave one out gene gene scatter looks meh
    #     'FXYD7', 'S100PBP',
    #     ###'PPP1R12B', # leave one out gene gene scatter looks meh
    #     ###'NYNRIN', # pretty correlated with nktdp_maf_sig
    #     'IGFBP7', 'KCNMB4', 'MRC2', 'RUNX3', 
    #     ###'AHR', # leave one out gene gene scatter looks meh
    #     'ADGRD1', 
    #     ###'IL2RG', # leave one out gene gene scatter looks meh 
    #     ###'ADGRG5', # leave one out gene gene scatter looks meh
    #     ###'GPR183', # leave one out gene gene scatter looks meh
    #     ###'CAPG', # leave one out gene gene scatter looks meh
    #     ###'ITGB7', # leave one out gene gene scatter looks meh
    #     # added manually from lower corr with TESC:
    #     'SPON2',
        
    #     # 0.5 < corr < 0.57
    #     # 'RAB11FIP1', 'ADARB1', 'TES', 'SCPEP1', 'TXN', 'SPTY2D1', 
    #     # ###'PTH1R', # in nktdp_maf_sig
    #     # 'TRDC', 'LINC01480', 'MED13L', 'SPON2', 'HERPUD1', 
    #     # # 'FXYD5', # strong
    #     # 'SLC35E3', 'EIF4G2', 'CCDC71L', 'CNOT2', 'CARD11', 'MGAT5',
    #     # 'ARMH1', 'LINC01814', 'TDRKH-AS1', 'CD7', 'ARHGEF10', 'EHBP1L1',
    #     # 'FHIT', 'IER5', 
    #     # ###'IL17RE', # in nktdp_maf_sig
    #     # 'MYO9A', 'PLXNC1', 'LAT2', 'LINGO4',
    #     # 'YAF2', 
    #     # ###'RORC', # in nktdp_maf_sig
    #     # 'NFATC4', 'N4BP2', 'AGFG1', 'MYOM2', 
    #     # ###'TLE3', # in nktdp_maf_sig
    #     # 'G3BP2', 'RBBP6', 'IGFBP6', 'PTEN', 'DRAP1',

    #     # iter1: corr with RUNX3, ACY3, TFRC, CALM2, and filtering by leave one out
    #     # 'CALM2', # strong
    #     # ###'ACY3', # strong, and leave one out gene gene scatter looks somewhat meh
    #     # 'TESC', 'ATP10A', 'TFRC', 
    #     # ###'LSP1', # strong, and leave one out gene gene scatter looks meh
    #     # 'PTK2B',
    #     # 'ADGRD1', 'FYB1', 'SEMA4D', 'RUNX3', 'POLB', 'KIAA0930', 'KYNU',
    #     # 'NFIL3', 'RAB11FIP1', 'SYT2', 'TXN', 'CYTH4', 'AHR', 'PIK3AP1',
    # ],
    'nktdp_samhd1_sig': [
        # iter2: iter1 filtered by heatmap of NKTDP MCs.
        'COTL1',
        'SAMHD1', # strong, but looks nice in leave one out gene gene scatter
        ###'VIM',
        'RGS10','NDRG2','SPI1','GRN','IRF8','CST3','MPEG1','C12orf75','JAK2','CNPY3',
        ###'AL592430.1',
        'PLAC8','HLA-DMA',
        
        'CIITA',
        'HLA-DRB1', # strong, but looks good in leave one out gene gene scatter that doesn't contain any of the strong genes
        'HLA-DQB1','HLA-DMB','LTA4H','PYCARD','TIFAB','FCGRT','SEL1L3','PLD4',
        'HLA-DPA1', # strong, but looks good in leave one out gene gene scatter that doesn't contain any of the strong genes
        'HLA-DRA', # strong, but looks good in leave one out gene gene scatter that doesn't contain any of the strong genes
        'HLA-DPB1', # strong, but looks good in leave one out gene gene scatter that doesn't contain any of the strong genes
        'FGD2','CTSZ','TLR10','HLA-DQA1','FGL2',


    
    #     'IRF8', 'GRN', 'HLA-DPA1', 'TIFAB', 'SLAMF7', 'CTSZ', 'PPM1J',
    #    'MCOLN2', 'IL18R1', 'FCGRT', 'FGL2', 'CIITA', 'FGD2', 'TNFAIP8',
    #    'APEX1', 'CST3', 'HLA-DQB1', 'BRI3BP', 'GDI2', 'C12orf75', 'SHTN1',
    #    'LILRA4', 'BTLA', 'FCER1G', 'CFP', 'HLA-DPB1', 'SAMHD1', 'MPEG1',
    #    'SPI1', 'HLA-DMB', 'H2AFY', 'HLA-DRA', 'IGSF6', 'EEF1B2', 'COTL1',
    #    'UBE2E2', 'BID', 'SEMA4A', 'FUOM', 'FAM81A', 'MARCH1', 'CLCN5',
    #    'ANXA5', 'LRRK2', 'NAGA', 'NDRG2', 'ZNF366', 'CCDC88A', 'EMB',
    #    'DHRS9', 'SCIMP', 'PYCARD', 'RGS10', 'HNRNPC', 'LILRB1', 'CNN2',
    #    'TLR10', 'DPP4', 'SLC12A9', 'FAM129A', 'WDFY2', 'EFHD2', 'TPM3',
    #    'HLA-DMA', 'CNPY3', 'SELPLG', 'PLAC8', 'KCNMB1', 'PLD4',
    #    'CACNA2D3', 'SCN9A', 'RAB32', 'CPPED1', 'JAK2', 'S100B', 'CDV3',
    #    'AL592430.1', 'GCSAM', 'SLC8A1', 'KLF4', 'LTA4H', 'RPL4', 'HCLS1',
    #    'RNASET2', 'HLA-DQA1', 'PAK1', 'SEL1L3', 'LYZ', 'LILRA1', 'TPP1',
    #    'BCL6', 'HLA-DRB1', 'MS4A6A', 'PLEK', 'FLT3', 'PTK2', 'VIM',
    #    'SKAP2',
    ],
    'clp_dntt_sig': [
        # NOTE: there is some overlap with higher_in_pro_b_than_clp, but this one has DNTT, which i think is right. maybe optimally would be better to add DNTT to higher_in_pro_b_than_clp etc, but the priority of this is quite low now, i think.
        
        'DNTT', 
        ###'ACTG1', # very strong, but looks nice in MC gene gene scatter. but does not look good in cell leave one out analysis, when excluding DNTT.
        'VPREB1', 
        'IL7R', # not specific?
        'CD79A', 
        ###'CYGB', # much higher in pro-B than in CLP, so I guess it came up here due to later cells misclassified as CLP? across real CLP cells, i guess it is not very correlated with the rest of the sig. i.e., it increases later in the differentiation process.
        ###'UHRF1', # much higher in pro-B than in CLP, so I guess it came up here due to later cells misclassified as CLP? across real CLP cells, i guess it is not very correlated with the rest of the sig. i.e., it increases later in the differentiation process.
        'ADA',
        'GSDME', 
        ###'TBCD', # didn't look good in mds illu and mds ult.
        ###'AC084033.3', 
        'IGLL1', # might look at first a bit bad in mds illu and mds ult sig heatmaps, but i think it is because of CLP-E MCs.
        ###'CALM1', # strong, and looks meh in MC gene gene scatter
        ###'CDK6', # strong, and looks meh in MC gene gene scatter
        ###'BLNK', # in MC leave one out gene gene scatter, looks better in clp_id2_runx2_sig than in clp_dntt_sig.
        'SPON1', 'BAALC', 
        ###'TMSB10',  # very strong, and looks meh in MC gene gene scatter
        ###'LTB',  # strong, and looks meh in MC gene gene scatter
        'CD48', 
        ###'STMN1', # strong, and looks meh in MC gene gene scatter
        ###'MYLK', # looks meh in MC gene gene scatter
        ###'RRBP1', # looks meh in MC gene gene scatter
        ###'HCST', # strong, and looks bad in MC gene gene scatter
        ###'ARPP21', # much higher in pro-B than in CLP, so I guess it came up here due to later cells misclassified as CLP? across real CLP cells, i guess it is not very correlated with the rest of the sig. i.e., it increases later in the differentiation process.
        'FRMD4A',
        ###'GAPDH',  # strong, and looks meh in MC gene gene scatter
        ###'SMIM24', # strong, and looks meh in MC gene gene scatter
        ### 'EIF1', # very strong, and looks meh in MC gene gene scatter
        ###'GAS7', # didn't look good in mds illu and mds ult.
        ###'OSBPL3', # in MC leave one out gene gene scatter, looks somewhat better in clp_id2_runx2_sig than in clp_dntt_sig.
        ###'CYFIP2', # strong, and looks meh in MC gene gene scatter
        'A1BG',
    ],
    'clp_id2_runx2_sig': [
        # iter2: iter1 filtered by heatmap of CLP MCs and by strength.
        ###'GPR183', # looked somewhat bad in marker heatmap
        ###'AHR', # leave one out gene gene scatter looks meh.
        'TESC','TTN',
        ###'IGFBP7', # looked somewhat bad in marker heatmap
        'NFIL3',
        ###'SPCS1', # looked somewhat bad in marker heatmap
        'CAPG','ITGB7','ACY3','BLNK','TNFRSF18','MARCKSL1','RUNX3','FAT3','TNFRSF4',
        'ITGA4',
        ###'STMN1',
        ###'TMSB10',
        'SATB1','DDIT4','RASD1','NYNRIN','TRIB2','STK17B','CYSLTR1',
        ###'JCHAIN',
        'KLRB1','MAF','RUNX2','CRYBG3','LINC00865','RGS1','ID2',
        'METTL7A',
        ###'DCK', # leave one out gene gene scatter looks meh.
        'STAP1',
        ###'MOB1B', # leave one out gene gene scatter looks meh.
        ###'IER2', # leave one out gene gene scatter looks meh.
        ### 'SLC38A1', # seems to increase before most others, though maybe good enough
        'SPECC1',
        ###'DYNC1LI2', # seems to increase before most others
        ###'ARMH1', # leave one out gene gene scatter looks meh.
        ###'TRAF3IP3', # seems to increase before most others
        ###'TSC22D1',
        ###'MYL12A', # seems to increase before most others
        ###'DRAP1', # seems to increase before most others
        'LINC01374',
        ###'H3F3B','ACTB','HMGB1','PTMA','TMSB4X',
        'ZNF608','STK17A','DYRK2','MEF2A',
        ###'TNIK', # seems to increase before most others
        'RIOX2','CLEC2D','LY86','NEGR1',
        ###'LAT2', # leave one out gene gene scatter looks meh.
        'OSBPL3','APP','REEP5','CD7','ATP10A','ZNF217','CYTIP','EVI2B',
        ###'RUFY3', # leave one out gene gene scatter looks meh.
        ###'SLC38A2', # leave one out gene gene scatter looks meh.
        'NREP',
        ###'TES', # leave one out gene gene scatter looks meh.
        ###'DDX5',
        'CELSR1',
        ###'MAML3', # leave one out gene gene scatter looks meh.
        'HDAC9','FYB1','BCL2','BANK1',
        ###'FOS', # notorious, and leave one out gene gene scatter looks meh-ish.
        'MLC1','NCF1','MDFIC','CLEC4A','NFIA',
        ###'ETS2', # leave one out gene gene scatter looks good. but looks somewhat bad in marker heatmap
        'PRL','IL17RE','CCDC141','TOX2','RORC','MRAS','UNC13C',

        ###'GSN', # seems to increase before most others
        ###'GABARAPL2', # seems to increase before most others
        'MRC2','PTH1R','PIEZO2',
        ###'CALM2',
        'TDRKH-AS1','ARHGAP10','NFATC4','SUCLG2','TDRKH','REN',
        ###'HSP90AA1',
        'CTSC','SCPEP1',
        ###'SNRNP25', # leave one out gene gene scatter looks meh.
        'CTSW',
        'WWOX',
        

        
        # 'RUNX2', 'ID2', 
        # ###'AL096865.1', 
        # 'PTMA', 'CRYBG3', 'MAF', 'KLRB1',
        # 'RGS1', 'CYSLTR1', 'FAT3', 'LINC00865', 
        # ###'JCHAIN', 
        # 'TNFRSF18',
        # 'DDIT4', 'MEF2A', 'TRIB2', 'CAPG', 'TTN', 'STK17A', 'STK17B',
        # 'EVI2B', 'PRL', 'ITGA4', 'NFIL3', 'CCDC141', 'NYNRIN', 'TMSB4X',
        # 'ACY3', 'RUNX3', 'MDFIC', 'NFATC4', 'SPCS1', 'TESC', 'ITGB7',
        # 'PIEZO2', 'TMSB10', 'RASD1', 'SOX4', 'UNC13C', 'ACTB', 'BANK1',
        # 'ETS2', 'PTH1R', 'RORC', 'NCF1', 'TNFRSF4', 'HMGB1', 'CYTIP',
        # 'CELSR1', 'MRAS', 'FAM69C', 'H3F3B', 'CCR1', 'FOS', 'ZNF217',
        # 'DDX5', 'TOX2', 'IGFBP7', 'TSC22D1', 'NEGR1', 'APP', 'CTSW',
        # 'CTSC', 'RIOX2', 'HDAC9', 'ARHGAP10', 'SPON2', 

        # 'MLC1', 'BLNK',
        # 'FGFBP2', 'TDRKH-AS1', 'NFIA', 'SATB1', 'IER2', 'RNASE6', 'CALM2',
        # 'DUSP4', 'MYL12A', 'REEP5', 'SPECC1', 'FYB1', 'GSN', 'SLC38A1',
        # 'MARCKSL1', 'RUFY3', 'BCL2', 'ARMH1', 'MGP', 'NREP', 'WNT5B',
        # ###'AC005307.1', 
        # 'SYT2', 'CLEC4A', 'DCK', 'AL136987.1', 'CLEC2D',
        # 'AHR', 'TES', 'SNRNP25', 'WWOX', 'ATP10A', 'GPR183', 'DYRK2',
        # 'PTPRK', 
        # ###'AC090673.1', 
        # 'LPAR1', 'METTL7A', 'SCPEP1', 'LINC01374',
        # 'OSBPL3', 
        # ###'AC114760.2', 
        # 'MOB1B', 'GABARAPL2', 'LY86', 'TRAF3IP3',
        # 'ZNF608', 'TNIK', 'SUCLG2', 'IL17RE', 'LEF1', 
        # ###'AL096799.1',
        # 'KCNQ3', 'STMN1', 'REN', 'RORB', 'STAP1', 'DYNC1LI2', 'HSP90AA1',
        # 'MRC2', 'SLC38A2', 'TDRKH', 'LAT2', 'MAML3', 'CD7', 'ILDR2',
        # 'DRAP1',
    ],
    'higher_in_nktdp_than_clp': [
        'SAMHD1', 'ID2', 'ACY3', 'LGALS1', 'NFKB1', 'ITGB7', 'TTN',
        'GPR183', 'NFIL3', 'SCPEP1', 'TESC', 'CYBB', 'RUNX3', 'FYB1',
        'GNLY', 'COTL1', 'IGFBP7', 'PLD4', 'AHR', 'KIT', 'SYT2',
        'RAB11FIP1', 'RNASE6', 'ATP10A', 'SEMA4D', 'CCDC50', 'SEL1L3',
        'MS4A1', 'LY86', 'BCL2', 'TLR10', 'SPON2', 'BLNK', 'ZDHHC14',
        'NCF1', 'ITGAX',
    ],
    'eryp_plek_sig': [
        

        
        ###'ZEB2','SLC40A1', # strong, and looked meh in leave one out gene gene scatter.
        'SERPINB6','GATA2','TAOK3','DEPTOR','PBX1',
        ###'PNRC1', # strong, and looked somewhat meh in leave one out gene gene scatter.
        'ST8SIA6','CELF2','PCDH9',
        ###'AC011139.1',
        'NRIP1','ARHGEF6','TESPA1','ATP6V0A2','NPR3','FERMT3','ABCC4','MYO1G','CPPED1','IGSF10','LAPTM5','MAX',
        ###'MACF1', # strong, and looked meh in leave one out gene gene scatter.
        'AKAP13','TMEM65','CAMK1D',
        # added if corr>0.8 with the original iter2
        'PLEK', 'GAPT', 'CPA3', 'MSI2', 'KCNQ1OT1', 'FNBP1', 'ANGPT1',
       'MED12L', 
       ###'HLA-DRA', # strong, and looked meh in leave one out analysis across cells.
       'SMIM24', 'SOD2', 'DDAH2', 'HLA-DPB1',
       'STXBP5', 
       ###'SERPINB1', # too strong
       'HEMGN', 'EREG', 'CRHBP',
        
        
        # 'SLC40A1', 'ZNF385D', 'AFF1', 'KAT6B', 'ARHGEF12', 'ARHGEF6',
        # 'CD84', 'ATXN1', 'ZEB2', 'ARFGAP2', 'H1F0', 'NRIP1', 'PNISR',
        # 'TMEM65', 'GSE1', 'CAMK1D', 'CFLAR', 'TRAF3IP3', 'IKZF1',
        # 'ST8SIA6', 'NLK', 'CPPED1', 'PNRC1', 'NFAT5', 'CITED2', 'FRMD6',
        # 'ARID4B', 'BPTF', 'HIST1H2AC', 'C1orf21', 'TESPA1', 'ZMYM2',
        # 'WDFY2', 'PBX1', 'GOLGB1', 'SESN3', 'ALDH2', 'ARL6IP5', 'NPR3',
        # 'MACF1', 'IGSF10', 'ABCC4', 'PCDH9', 'AKAP13', 'BRPF1', 'CDC42SE2',
        # 'FERMT3', 'LEPROT', 'ATP6V0A1', 'SERPINB6', 'NAGK', 'MYO1G',
        # 'AC011139.1', 'EID2B', 'FAM117A', 'GATA2', 'ZDHHC5', 'IKZF2',
        # 'CCZ1', 'SETX', 'FRY', 'FNBP4', 'STAT5B', 'AC090152.1', 'ERG',
        # 'TRIM24', 'NFE2', 'SFT2D3', 'HIST2H2AC', 'WRNIP1', 'SEMA3C',
        # 'CASD1', 'YPEL3', 'ATP6V0A2', 'DNM1L', 'ARHGAP15', 'CHURC1',
        # 'LAPTM5', 'MAX', 'ZMYM5', 'CREBRF', 'TAOK3', 'MADD', 'DEPTOR',
        # 'ABCB1', 'KIAA1109', 'OAZ2', 'VAV3', 'CELF2', 'PPP1R15A',
    ],
    'eryp_apoc1_sig': [
        
        'ITGA2B','AHSP','FAM89A','PLIN2',
        ###'PHF6', # somewhat strong, and looked somewhat meh in leave one out gene gene scatter.
        'TST','MPST','CMBL',
        'TMEM14B', # somewhat strong
        ###'ITGB1', # somewhat strong, and looked somewhat meh in leave one out gene gene scatter.
        'UROD','SMIM10','TMEM141','SYNGR1','SMIM1','PPP1R14A',
        ###'PRDX2', # somewhat strong, and looked somewhat meh in leave one out gene gene scatter.
        'CA1',
        ###'CST3', # somewhat strong, and looked somewhat meh in leave one out gene gene scatter.
        'ANK1','CAST','CDH1','PVT1','DDI2','NFIA','FAM45A','ELOVL6','EPCAM','REXO2','MPC2','FAM178B','CD36','APOE','APOC1',
        ###'MYC', # somewhat strong, and looked meh in leave one out gene gene scatter.
        'BLVRB',
        ###'TMEM14C', # somewhat strong, and looked meh in leave one out gene gene scatter.
        'HBB',
        
        
        # 'BLVRB', 'APOC1', 'ANK1', 'HBB', 'SMIM1', 'APOE', 'TMEM14B',
        # 'FAM178B', 'TMEM14C', 'EPCAM', 'MPC2', 'CA1', 'CD36', 'AHSP',
        # 'REXO2', 'CAST', 'IGLL1', 'FAM45A', 'NFIA', 'SYNGR1', 'TST',
        # 'TMEM141', 'PPP1R14A', 'CDH1', 'CMBL', 'PVT1', 'RHAG', 'FAM89A',
        # 'GBGT1', 'CST3', 'EIF3A', 'KLF1', 'METAP2', 'KCNH2', 'MPST',
        # 'PHF6', 'ATP5IF1', 'UROD', 'EPSTI1', 'PRDX2', 'SMIM10', 'TNFRSF25',
        # 'RPL22L1', 'UBAC1', 'PNMT', 'STRADB', 'MYC', 'ELOVL6', 'FCGRT',
        # 'PPA1', 'TPGS2', 'ITGA2B', 'YBX1', 'RREB1', 'DDI2', 'PLIN2',
        # 'SMIM37', 'TMOD1', 'LYAR', 'GPX4', 'CASP3', 'CDK4', 'TAF9', 'ATG3',
        # 'MARCKSL1', 'HBS1L', 'SEC31A', 'PTRHD1', 'ESD', 'CASP6', 'CD320',
        # 'ITGB1', 'MTSS1', 'TRIB2', 'SLC26A2', 'NME4', 'PRKAR2B', 'CD99',
        # 'PRDX3', 'NIPSNAP2', 'PSMD12', 'PSMG1', 'PMP22', 'SLC39A8', 'ECH1',
        # 'BLVRA', 'C9orf40', 'CSF1', 'AIG1', 'GLRX', 'PABPC4', 'ISOC2',
        # 'TMEM245', 'NUCB2', 'TESC', 'LMNA', 'HSPD1', 'SASS6', 'NECAB1',
        # 'SUCLG2', 'CENPU', 'MRPS12', 'CBLL1', 'STMN1', 'MRPL15',
        # 'MIR4435-2HG', 'RIF1', 'LPCAT3', 'HSPB1', 'TXN2', 'HOOK1', 'MIF',
        # 'LINC00891', 'ZNF706', 'ALDH1A1', 'MGST2', 'RBBP7', 'UCA1', 'XK',
        # 'CPNE3', 'C19orf48', 'HNRNPC', 'FKBP4', 'ATP13A3', 'C1QBP',
        # 'CCDC34', 'LBR', 'PTP4A2', 'RYR3', 'NDUFAF2', 'CEBPZ', 'EXOSC5',
        # 'FAM118A', 'RBX1', 'PCCB', 'ATP5MC1', 'PNPO', 'HIST1H4C', 'GSPT1',
        # 'NFE2L2', 'ATP5MPL', 'PCAT18', 'DDX21', 'AKR1C1', 'FBL', 'IBTK',
        # 'HEBP1', 'DMAC1', 'ZBTB7A', 'THRAP3', 'BAZ1A', 'ARV1', 'SNRPB',
        # 'ALAD', 'PA2G4', 'VAPA', 'NCL', 'ARL2', 'PRKDC', 'SNRPF', 'GTF2A2',
        # 'CD47', 'BRIX1', 'ECI1', 'CDYL', 'HIST1H2BE', 'ATP5MF', 'NET1',
        # 'BTK', 'HBD', 'AMMECR1', 'LEF1', 'ACAT1', 'MACROD1', 'RMDN1',
        # 'NASP', 'RAD23A', 'AKR1C3', 'SPTA1', 'KIF13A', 'HMGA1', 'SNU13',
        # 'U2SURP', 'HADH', 'GADD45GIP1', 'FECH', 'DLD', 'ANKRD12', 'HMGN5',
        # 'WBP11', 'NDUFAB1', 'DUT', 'HPS1', 'RNH1', 'AK6', 'GLRX3', 'POP7',
        # 'POLR2H', 'CASP8', 'HNRNPAB', 'MXD1', 'HES6', 'RAN', 'ETFA',
        # 'KMT5A', 'TFR2', 'PHB2', 'HAGH', 'MGAT4B', 'SORD', 'SLC25A3',
        # 'ENDOG', 'TIMM10', 'CISD2', 'ASAP1', 'PPP1R14B', 'CCT8', 'BCL7A',
        # 'NPW', 'TCEA1', 'RSRC1',
    ],
    'clp_vim_sig': [
        
        'LMNA','AHNAK','ANXA1','ATP2B1','VIM','MYADM',
        ###'CD44', # leave one out gene gene scatter looks meh.
        'MTURN','IDS','PBXIP1','ETS1','TLE4','AAED1','CEP126','S100A6','UCP2','TAGLN2','DPYSL3',
        ###'TSC22D3', # leave one out gene gene scatter looks meh.
        'TSPAN2',
        ###'S100A10', # leave one out gene gene scatter looks meh.
        'DBI','KLF2','KLF6','ANKRD28','TUBB6','SERPINE1','CRIP2','WDR49','DPF1','RGCC','TUBA1A','SOCS3','EMP1','ADGRE5',
        
        
        # 'VIM', 'ANXA1', 'AHNAK', 'MYADM', 'TAGLN2', 'ATP2B1', 'TSPAN2',
        # 'LMNA', 'DPYSL3', 'CRIP2', 'SERPINE1', 'AC020916.1', 'TUBA1A',
        # 'RGCC', 'UCP2', 'ADGRE5', 'ANKRD28', 'EMP1', 'S100A6', 'SOCS3',
        # 'CEP126', 'MTURN', 'PBXIP1', 'AAED1', 'TUBB6', 'DBI', 'DPF1',
        # 'TPPP3', 'CD44', 'TSC22D3', 'KLF2', 'CRIP1', 'TPPP', 'TLE4',
        # 'GAPDH', 'CNIH2', 'S100A10', 'PRMT2', 'TPM3', 'KLF6', 'ANXA2',
        # 'CACNB2', 'LRRFIP1', 'IDS', 'ETS1', 'DGKG', 'WDR49',
    ],
    'mpp_mebemp_vim_sig': [
        # started by corr with AHNAK, but then concluded it is noisy...
        # 'AHNAK', # strong. correlated with sum of all non-strong ones, but the MC scatter is a bit ugly
        'VIM', # very strong. highly correlated with the sum of all non-strong ones.
        # 'TSPAN2', 
        # 'ANXA1', # very strong. correlated with sum of all non-strong ones, but the MC scatter is a bit ugly. 
        # 'TAGLN2', # strong. correlated with sum of all non-strong ones, but the MC scatter is a bit ugly

        'ATP2B1', 'FHL1',
        'MYADM', 
        'S100A10', 
        'LMNA', 'PBXIP1', 'CALCOCO2', 'TUBA1A',
        'ANKRD28', 'SEPT11', 'ATP8B4', 'IDS', 'RPGR', 'KLF6', 'WDR49',
        'DPYSL3', 'EMP1', 'TFPI', 'MAMDC2', 'ANXA2', 'CEP126', 'ANTXR2',
        'TCP11L2', 'TMOD3', 'PGM2L1',

        # correlated with the sum of the genes above. filtered by marker heatmap, and manually by gene gene scatter:
        ###'CD44', 
        ###'S100A11', 
        ###'CD99', 
        'TSC22D3',
        'CRIP2', 'SSBP2', 'CAST', 'PTMS', 'UCP2', 'ZSCAN18', 'CA5B', 'DBI',
        ###'AC020916.1', 
        ###'FOS', 
        'FRMD4B', 'SCN9A', 'ADGRE5', 
        'MBOAT7', 
        'GBP2',
        #    'CACNB2', 
        #    'PRMT2', 
        ###'NAP1L1', 
        #    'SERPINE1', 
        #    'DPF1', 
        #    'NR3C1',
        #    'ARHGEF40',

        # correlated with the sum of the genes above.
        ###'PTPRC', 
        # 'RFX2',
    ],
    # 'higher_in_clp_than_nktdp': [ # 240619: currently unused
    #     'PRSS2', 'HOPX', 'CD34', 'HOXA9', 'MZB1', 'EGFL7', 'KIAA0087',
    #    'MME', 'HEMGN', 'SLC2A5', 'CYTL1', 'COBLL1', 'DNTT', 'TFPI',
    #    'SCN3A', 'BASP1', 'NPR3',
    #    ###'AC004540.2', 
    #    'MTURN', 'INKA1', 'NCOA7',
    #    'RUNX1', 'ADGRE2', 'PTMS', 'MEIS1', 'PROM1', 'GYPC', 'HOXA10',
    #    'SPON1', 'AFF1', 'LYL1', 'ARID5B', 'RBPMS', 'HOXA6', 'HOXA3',
    #    'TNFAIP2', 'OSBPL3', 'MAP1A', 'GNAI1', 'GIMAP7', 'CD22', 'GCSAML',
    #    'PACSIN1', 'LMO2', 'IFNG-AS1', 'VPREB1', 'ADGRG6', 'CD82', 'VSIR',
    #    'QPRT', 'LIMS1', 'BAALC', 'SEMA3C', 'CCR9', 'TMEM38B', 'LPIN1',
    #    'IGF2BP2',
    # ],
    'higher_in_clp_m_than_all_myeloid_hspcs': [
        # >2.5
        'KIAA0087', 
        'LTB', # remove?
        # 'PRSS2', # remove?
        'MME', 
        # 'IGHM', # remove?
        'JCHAIN', 'CCR7',
        'CLNK', 'SCN3A', 'BASP1', 'HOXA9', 'TSC22D1', 'SLC2A5', 'TCF4',
        'DNTT', 'TRBC2', 'COBLL1', 'HCST', 'FAM30A',
        # 2-2.5
        'MEF2C', 'AFF1',
        'TRAF3IP3', 'MZB1', 'HLX', 'TYROBP', 'S100Z', 'SATB1', 'KCNK17',
        'LINC00865', 'CYTH4', 'ITGAL', 'MEF2A', 'CD53', 'OSBPL3', 'LCP1',
        ###'TRIB2', 
        'RUNX2', 'CORO1A', 'PALLD', 'MOB1B', 'CCR9', 'LST1',
        'NKG7',
        # 1.5-2 # already strong enough even without these
        # 'AC245060.5', 'SPECC1', 'RERE', 'SLC38A1', 'EVL', 'TRBVB',
        #     'TXNIP', 'SPON1', 'RAB31', 'ID2', 'LZTFL1', 'PRKCB', 'ACY3',
        #     'JAML', 'IFNG-AS1', 'HOPX', 'DYRK2', 'MED13L', 'LMO4', 'CXXC5',
        #     'TTC3', 'INTS6', 'SEMA3C', 'ANKRD44', 'DCK', 'ITM2C', 'CYFIP2',
        #     'NUP214', 'VSIR', 'KCTD12', 'ARID5B', 'BCL11A', 
        #     'CD99', # remove??
        #     'HOXA10',
        #     'ADA', 'ZFP36L1', 
        #    #  'CD34', 
        #     'VPREB1', 'ADGRE2', 'RNASEH2B',
        #     'ST3GAL1', 'GNAS', 'KLF3', 'IQGAP1', 'LPIN1', 'REL', 'SPINK2',
        #     'N4BP2', 'CD47', 'STK17A', 'EML4', 'DDIT4', 'NPTX2', 'PIK3AP1',
        #     'CD22', 'SH3BP5', 'ADAM28', 'PAG1', 'SYK', 'TMEM38B', 'VCL',
        #     'AL096865.1', 'KLF6', 'ARHGAP25', 'PLP2',
    ],
    # 'higher_in_pre_b_than_clp_m': HIGHER_IN_PRE_B_THAN_CLP_M, # 240619: currently unused
    # 'higher_in_nktdp_l_than_clp_e_candidates': [ # 240619: currently unused
    #     # >2
    #     'ID2', 'ACY3', 'SAMHD1', 'GNLY', 'LGALS1', 'CALM2', 'RUNX3',
    #     'GPR183', 'SCPEP1', 'NFKB1', 'TTN', 'RUNX2', 'NFIL3', 'DDIT4',
    #     'ITGB7', 
    #     ### 'ACTB', # remove? very high in general...
    #     'LY86', 'CYBB', 'COTL1', 'RAB11FIP1', 'CAPG',
    #     'MAF', 'JCHAIN', 'TESC', 'FYB1', 'TNFRSF18', 'HERPUD1', 
    #     ###'TMSB10', # remove? very high in general...
    #     'PIK3AP1', 'RNASE6', 
    #     ### 'TMSB4X', # remove? very high in general...
    #     'IRF8', 'LSP1', 'ATP10A', 'CCDC50',
    #     'TLR10', 'SYT2', 'SPON2', 'DBI',
    # ],
    'strong_c_monocyte_specific': STRONG_MONOCYTE_SPECIFIC_GENE_NAMES,

    # 'monocyte_and_dc_somewhat_specific': [ # 240619: currently unused
    #     'CST3', 'SAMHD1', 'LYZ', 'IRF8', 'FCER1G', 'PTPRE', 'S100A10',
    #     'TYROBP', 
    #     ### 'LGALS1', # strong in N186-AML
    #     'ANXA2', 'CCDC50', 'UGCG', 'RNASE6',
    #     'C12orf75', 'S100A6', 'FGL2',
    # ],
    # 'nc_monocyte_specific': [ # 240619: currently unused
    #     'CDKN1C', 'LYPD2', 'FCGR3A', 'TCF7L2', 'SMIM25',
    #     'NFKBIZ', 'BCL2A1', 'SIGLEC10',
    # ],
    'higher_in_nc_monocyte_than_c_and_interm_monocyte': [
        'FCGR3A', 'CDKN1C', 'TCF7L2', 'MTSS1', 'RHOC', 'LYPD2', 'HES4',
       'CD79B', 'BCL2A1', 'NFKBIZ', 'LTB', 'VMO1', 'PAG1', 'C1QA', 'EVL',
       'SIGLEC10', 'CKB', 'CYTIP', 'CHST2', 'INSIG1', 'KCNQ1OT1',
       'PLAGL2', 'PPM1N', 'ABI3', 'CTSL', 'RRAS', 'PTP4A3', 'HMOX1',
       'BATF3',
    ],
    # 'higher_in_nc_monocyte_than_b': [
    #     'LST1', 'AIF1', 'CST3', 'FCER1G', 'FCGR3A', 'TYROBP', 'SERPINA1',
    #     'FCN1', 'SMIM25', 'CFD', 'IFITM3', 'S100A9', 'CDKN1C', 'LYZ',
    #     'DUSP6', 'CEBPB', 'TCF7L2', 'LILRB2', 'FYB1', 'FGL2', 'MAFB',
    #     'CLEC12A', 'CLEC7A', 'GIMAP4', 'IGSF6', 'SAMHD1', 'SPN', 'CX3CR1',
    #     'CSF1R', 'TIMP1', 'PECAM1', 'GIMAP7', 'STXBP2', 'CFP',
    # ],
    # 'higher_in_b_than_nc_monocyte': [
    #     'MS4A1', 'IGHM', 'IGKC', 'CD79A', 'BANK1', 'RALGPS2', 'IGHD',
    #     'IGLC2', 'IKZF3', 'FCRL5', 'FCMR', 'FCRL2', 'CD22', 'LINC00926',
    #     'SPIB', 'CD24', 'TTN', 'LINC01857', 'LBH', 'FAM129C', 'LINC02397',
    #     'IGLC3', 'FCRL3', 'TNFRSF13C',
    # ],
    'higher_in_monocyte_than_b': [
        'S100A9', 'LYZ', 'S100A8', 'CST3', 'FCN1', 'VCAN', 'AIF1',
        'TYROBP', 'FCER1G', 'LST1', 'S100A12', 'FYB1', 'FGL2', 'CD14',
        'SERPINA1', 'MS4A6A', 'DUSP6', 'CFD', 'CSTA', 'GIMAP4', 'JAML',
        'ANXA1', 'AC020656.1', 'CLEC12A', 'IGSF6', 'SAMHD1', 'CEBPD',
        'CD36', 'CLEC7A', 'CPVL', 'LGALS2', 'KCTD12', 'TNFAIP2', 'LGALS3',
        'GIMAP7', 'NCF2', 'RGS2', 'CSF3R', 'CFP', 'FOS',
    ],
    'higher_in_b_than_monocyte': [
        'MS4A1', 'IGHM', 'IGKC', 'CD79A', 'BANK1', 'RALGPS2', 'CD79B',
        'IGLC2', 'IGHD', 'IKZF3', 'FCMR', 'FCRL5', 'FCRL2', 'CD22', 'LBH',
        'LINC00926', 'FAM129C', 'SPIB', 'FCRL3', 'TNFRSF13C', 'LINC01857',
        'IGLC3', 'CD24', 'LINC02397', 'ETS1', 'BLK',
    ],
    'higher_in_interm_monocyte_than_c_monocyte': [
        'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 
        ###'HLA-DRB1', # a bit strong in classical monocytes
        'HLA-DQB1',
        'FCGR3A',
        'HLA-DMA',
        'HLA-DRB5',
        
        # 'ABI3', 'MTSS1', 'CSF1R', 'LIPA', 'LY6E', 'RHOB',
    ],
    # 'higher_in_nc_monocyte_than_nkt': [ # 240619: currently unused
    #     'LST1', 'AIF1', 'MS4A7', 'CST3', 'SMIM25', 'SERPINA1',
    # ],
    # 'higher_in_memory_b_than_monocytes_and_nc_monocytes': [ # 240619: currently unused
    #     'MS4A1', 'IGKC', 'IGHM', 'BANK1', 'CD79A', 'RALGPS2', 'IGLC2',
    #    'IGHD', 'IKZF3', 'FCMR', 'FCRL2', 'LINC00926', 'CD24', 'FAM129C',
    #    'IGLC3', 'SPIB', 'CD22', 'TNFRSF13C', 'LINC01857', 'LBH',
    #    'LINC02397', 'FCRL5', 'AFF3', 'BLK', 'FCRL3', 'COBLL1', 'LTB',
    #    'TTN', 'SLC38A1', 'TCF4', 'TNFRSF13B', 'BCL2', 'ETS1', 'VPREB3',
    #    'P2RX5', 'PAX5', 'JCHAIN',
    # ],
    # 'higher_in_monocytes_and_nc_monocytes_than_memory_b': [ # 240619: currently unused
    #     'CST3', 'TYROBP', 'AIF1', 'FCER1G', 'LST1', 'FCN1', 'S100A9',
    #    'S100A4', 'S100A11', 'LYZ', 'FYB1', 'SERPINA1', 'DUSP6', 'PSAP',
    #    'CFD', 'FGL2', 'CLEC12A', 'S100A6', 'SAMHD1', 'CTSS', 'GIMAP4',
    #    'CLEC7A', 'IGSF6', 'SRGN', 'GIMAP7', 'TYMP', 'LGALS3', 'BRI3',
    #    'TSPO', 'CFP', 'VSIR', 'ITGB2', 'LILRB2', 'NCF2', 'LGALS1', 'CSTA',
    #    'S100A10', 'IRAK3', 'RGS2', 'P2RY13', 'DMXL2', 'CD300E', 'S100A8',
    #    'CPVL', 'JAML', 'ANXA1', 'LILRB3', 'SAT1', 'TIMP1', 'TNFAIP2',
    #    'CEBPB', 'IFITM3', 'EFHD2', 'SLC11A1', 'CD68', 'TNFSF10', 'MYO1F',
    #    'ANXA5', 'TKT', 'STXBP2', 'TNFSF13B', 'MAFB', 'MNDA', 'RXRA',
    #    'TNFRSF1B', 'CYBB', 'C5AR1', 'PILRA', 'CPPED1', 'ARRB2', 'LCP1',
    #    'FTL', 'VMP1', 'C1orf162', 'DUSP1', 'CEBPD', 'GSTP1', 'GIMAP8',
    #    'SMIM25', 'PRAM1', 'CD4', 'SPI1',
    # ],
    'plasmablast_somewhat_specific': [
        'IGHA1', 'TNFRSF17', 'SEC11C', 'TENT5C', 'IRF4', 
        ###'JCHAIN',
       'TXNDC5', 'IGHA2', 'ABCB9', 
       ###'SUB1',
       'TXNDC11', 
       # 'ISG20', 'SPCS2', # better without these, as they are a bit strong in others?
    ],
    'higher_in_b_than_plasmablast': [
        'HLA-DRA', 'MS4A1', 'LTB', 'HLA-DPA1', 'HLA-DPB1', 'CD37', 'BANK1',
       'HLA-DQB1', 'HLA-DRB1', 'MARCH1', 'RALGPS2', 
       ###'CD74', 
       'HLA-DMB',
       'CD79B', 
       ###'TXNIP', 
       'FCRL2', 'CD22', 'SCIMP', 'ARHGAP24',

       'SP110',
       'MARCKS', 'HLA-DQA1', 'SPIB', 'PLAC8', 'REL', 'BCL11A', 'SWAP70',
       'CD24', 'LBH', 'LINC00926', 'FOXP1', 'POU2F2', 'ZEB2', 'IFITM2',
       'STX7', 'CIITA',
    ],
    'igha_iglc': [
        'IGHA1', 
        'IGHA2',
        'IGLC2',
        'IGLC3',
    ],
    'ighm_ighg': [
        'IGHM',
        'IGHG1', 
        'IGHG2', 
        'IGHG3', 
        # IGKC also fits here, in the group of differentially expressed genes between the two plasmablast states (presumably)
    ],
    'higher_in_plasmablast_than_clp': [
        'IGHA1', 
        'IGKC', # comment?
        'IGLC2', 'TNFRSF17', 'SEC11C', 'ISG20', 'IGHA2',
       'IGLC3', 
       ###'SPCS1', 
       'DERL3', 
       ###'LGALS1', 
       'KLF2', 
       ###'XBP1', 
       'CD27',
       'IRF4', 'PRDM1', 
       ###'UBE2J1', 
       'HERPUD1', 'TENT5C', 'CD79A', 
       ###'SPCS2',
       ###'S100A10', 
       'CD38', 
       ###'MYDGF', 
       ###'HSPA5', 
       ###'SPCS3', 
       'PDIA4', 'FKBP11',
       ###'S100A6', 
       ###'SSR3',
         'SDF2L1', 'AQP3', 'SLAMF7', 
       ###'SELENOS',
    ],
    'higher_in_b_than_all_other_hspcs': [
        'MS4A1', 
        ###'IGKC', # ugh. extremely bursty?
        'IKZF3', 'RALGPS2', 'KLF2', 
        ###'IGLC2', # ugh. extremely bursty?
        'MARCH1',
        'FCRL5', 
        ###'IGLC3', # ugh. extremely bursty?
        'ARHGAP24', 'BANK1', 'BIRC3', 'LINC00926',
        'SCIMP', 'LINC01857', 'LINC02397', 'FCRL3', 'TNFRSF13C',
        'TNFRSF13B', 'FCRL2', 'SP140', 'AIM2', 'FCRLA', 'MTSS1', 'IGHA1',
        'CPNE5', 'LBH', 'OSBPL10', 'TMEM154', 'CLECL1', 'TBC1D9', 'FCGR2B',
        'CD40', 'CD48', 'BLK', 'MNDA', 'HLA-DQA1', 'SPIB', 
        ###'CD84', # 240622: removed as it seemed like it shouldn't be here, when plotting for nimrod_148_atlas
        'SLAMF6', 'FCER2', 'HLA-DOB', 'IGHG3', 'CD37', 'RAB30',
        'LINC01781', 'FCRL1',
    ],
    # 'higher_in_plasmablast_than_other_B': [
    #     # didn't use in the end
    #     # >3.5
    #     'JCHAIN', 'IGHA1', 'TNFRSF17', 'MZB1', 'PPIB', 'SEC11C', 'HSP90B1',
    #    'XBP1', 'DERL3', 'FKBP11', 'SUB1', 'MANF', 'SDF2L1', 'HSPA5',
    #    'SRGN', 'RRBP1', 'SELENOS', 'PRDM1', 'MYDGF', 'MAN1A1', 'LGALS1',
    #    'PDIA4', 'AQP3', 'CD38', 'CALR', 'SSR3',
    # ],
    # 'higher_in_plasmablast_than_other_mature_b': [
    #     # didn't use in the end
    #     ###'JCHAIN', 
    #     'IGHA1', 'PPIB', 'TNFRSF17', 'MZB1', 'XBP1', 'PRDM1',
    #     'DERL3', 'IGHA2',
    # ],
    'higher_in_b_than_nkt': [
        # >3
        'IGKC', 'HLA-DRA', 'MS4A1', 'IGHM', 'BANK1', 'CD74', 'CD79A',
       'RALGPS2', 'CD79B', 'IGLC2', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1',
       'HLA-DQA1', 'IGHD', 'MARCH1', 'CD22', 'MEF2C', 'BCL11A', 'TCF4',
       'FCRL2', 'IGLC3', 'FCRL5', 'SPIB', 'LINC00926', 'HLA-DQB1',
       'HLA-DMB', 'ARHGAP24', 'FAM129C', 'TNFRSF13C', 'LINC02397',
       'SWAP70', 'CD24', 'GAPT', 'LINC01857', 'AFF3', 'HLA-DMA', 'COBLL1',
       'MARCKS', 'SCIMP', 'BLK', 'CD37', 'CIITA', 'CYBB', 'TNFRSF13B',
       'ADAM28', 'PAX5', 'CCDC50',
    ],
    'higher_in_aged_b_than_memory_b': [
        *GENES_HIGHER_IN_ZEB2_HIGH_B,
    ],
    'higher_in_ccdc88a_high_b_than_other_mature_b': [
        # unfiltered
        # 'TCF4', 'S100A6', 'S100A4', 'ABCA6', 'CCDC88A', 'JCHAIN', 'CD27',
        # 'CLNK', 'TYROBP', 'ITGB1', 'TTN', 'AIM2', 'CLECL1', 'ADTRP',
        # 'LINC01857', 'VOPP1', 'IGHA1', 'AHNAK', 'TFEC', 'TRAC', 'HCST',
        # 'METTL8', 'LEF1', 'IGLL5', 'SYNE2', 'ID3', 'FCMR', 'MCOLN2',
        # 'CCDC50', 'SCIMP', 'AC012368.1', 'ANXA2', 'SGPP1', 'AC008105.3',
        # 'IL10RA',
        
        # filtered by marker heatmap
        'MCOLN2','TCF4','AC008105.3','ID3','METTL8','SYNE2','AC012368.1','CCDC88A','TRAC','IGLL5','LEF1','TYROBP','CLNK','ADTRP','ABCA6',
    ],
    'higher_in_naive_b_than_memory_b': [
        'TCL1A', 
        'IL4R', 
        ###'IGHD', 
        ###'PCDH9', 
        'PLPP5', 
        'AL139020.1', 
        'FCER2', 
        ###'YBX3', 
        ###'IGHM', 
        'APLP2', 
        ###'CLEC2D', 
        'TSPAN13', 
        'CXCR4', 
        'TAPT1', 
        'H1FX', 
        'GCNT1', 
        'HVCN1', 
        ###'NCF1', 
        'BACH2', 
        ###'MEF2C', 
        'PHF21A', 
        ###'BCL7A', 
        ###'STX7',
    ],
    'higher_in_nkt_than_nktdp': [
        # >3
        'IL32', 'GIMAP7', 'CCL5', 'SYNE2', 'GIMAP4', 'KLF2', 'GZMA',
        'CD3D', 'CD3G', 'CST7', 'CD2', 
        ###'S100A10', 
        'IL7R', 'LINC00861',
        'GZMM', 'TRBC1', 'ETS1', 'CD52', 'KLRD1', 'PYHIN1', 'CD48',
        'CD247', 
        ###'S100A6', 
        'IKZF3', 'RORA', 'SLFN5', 'GZMH', 'TRAC',
        'IL10RA', 'SAMD3', 'ARL4C', 'SKAP1', 'CALM1', 'BCL11B', 'ISG20',
        'GBP5', 'LCK', 'MYBL1', 'GIMAP1', 'GZMB', 'LAT', 'KLRG1', 'PRF1',
        'MT2A', 'CCL4', 'GZMK', 'SLFN12L', 'NKG7',
    ],
    'nkt_somewhat_specific': [
        'IL32', 'NKG7', 'GNLY', 'GZMA', 'CCL5', 'CD2', 'GIMAP4', 'GIMAP7',
        'SYNE2', 'CD3G', 'CD3D', 'CST7', 'ARL4C', 'KLF2', 'CD3E', 'KLRD1',
        'SLFN5', 'IL7R', 'TRBC1', 'LINC00861', 
        # 'CALM1', 
        'PRF1', 
        ###'FYB1',
        'CD247', 'ETS1', 'GZMB', 'CD48', 'TRAC', 'KLRB1', 'GZMH', 'PYHIN1',
        'RORA', 
        # 'HCST', 
        'SAMD3', 'FGFBP2', 'FAM107B', 'GZMM', 'CCL4',
        'LCK', 'TRBC2', 'IKZF3', 
        ###'ID2', 
        # 'BTG1', 
        'BCL11B', 
        # 'S100A10',
        # 'PTPRC',
        'MYBL1', 'CX3CR1', 'LBH', 'SKAP1', 'CLEC2D', 'ITGB2',
        # 'MYL12A', 
        ###'SAMHD1', 
        # 'RARRES3', 
        # 'TMSB4X', 
        'GBP5', 'FYN', 'GPR65',
        'PPP2R5C', 'FCGR3A', 'IL10RA', 'AAK1', 
        # 'RNF213', 
        'RASGRP1', 'TRDC',
        'CD8A', 'KLRF1', 
        # 'ATM', 
        # 'TRAF3IP3', 
        'PRDM1', 
        # 'HLA-C', 
        # 'KIAA1551',
        'CD7', 
        ###'STK17A', 
        'KLRG1', 'IL2RG', 
        # 'HLA-F',
        'GIMAP1', 'DENND2D',

        'CD96', 'C12orf75', 
        # 'CORO1A', 
        'LAT', 'LITAF', 
        # 'CTSW', 
        'EVL',

        'SYNE1', 
        # 'S100A4', 
        # 'S100A11', 
        'SYTL3', 'ITGAL', 
        # 'HLA-B', 
        'SH2D1A',
        'EFHD2', 'SLFN12L', 'IL2RB', 'PTPN4', 'SPOCK2', 'GNG2', 'TRG-AS1',
        'CD226', 
        # 'B2M',
        # 'TXNIP', 
        'SPON2', 'DOK2', 'LYAR', 'GZMK', 'STK17B',
        'PTGDR', 'PRKCH', 'ISG20',
    ],
    'endothel_somewhat_specific': [
        # mc_utils.get_mean_expr(n_mc_ad, n_mc_ad.obs['state'] == 'Endothel', ret_df=True, layer_name='expr', background_mask_for_rest_max_expr=n_mc_ad.obs['state'] != 'Unknown')
        # expr_rest_max_expr_log_ratio > 1.5
        # 'MGP', 'SPARCL1', 'CAV1', 'ADIRF', 'IFI27', 'TM4SF1', 'EMCN',
        # 'CALD1', 'GNG11', 'SPARC', 'CRIP2', 'EPAS1', 'CAVIN1', 'COL18A1',
        # 'FLT1', 'CALCRL', 'HSPG2', 'DSTN', 'SERPINE1', 'CLIC4', 'ARL4A',
        # 'PLK2', 'ATP8B1', 'PLS3', 'LAMC1', 'DOCK9', 'FBN1', 'IGFBP4',
        # 'MMRN1',

        'VWF', 'POSTN', 'NR2F2', 
        ###'SYNE2', 
        'RNASE1', 'CAV1', 'WWTR1',
        'TM4SF1', 'CLEC14A', 'COL3A1', 'CAVIN1', 'EMCN', 'IGFBP4', 'MMRN1',
        'SPARCL1', 
        ###'S100A10', 
        'FILIP1L', 'PECAM1', 'CLU', 'ARHGAP29',
        'HSPG2', 'LAMA4', 
        ###'S100A11', 
        'TJP1', 'SPARC', 'CD9', 'PRCP',
        'SKIL', 'ARL4A', 'KLF2', 'BMPR2', 'EFEMP1', 
        ###'ANXA5', 
        'CRIP2',
        'FBN1', 'RAPGEF5', 'LDB2', 'IL6ST', 'MTUS1', 'APP', 'SASH1', 'FN1',
        'DSTN', 'ENTPD1', 'CALD1', 'PDLIM1', 'ECE1', 'BGN', 'PTPRB',
        ###'ANXA2', 
        'AQP1', 
        ###'CD93', 
        'ADGRL4', 'SNTB2', 
        ###'TIMP1', 
        'NFIB',
        ###'MARCKS', 
        ###'IL32', 
        ###'ZFP36L1', 
        'LMNA', 'NRP1', 'FAM198B', 'PODXL',
        'EMP1', 'TGFBR3', 'EFNB2', 
        ###'GIMAP4', 
        ###'ITGB1', 
        'S100A13', 'HYAL2',
        'CD59', 'MMRN2', 
        ###'KCTD12', 
        'RAMP2',
    ],
    'lamp3_outlier_specific': [
        'LAMP3', 'IDO1', 'CCL19', 'CD83', 'LAD1', 'NCCRP1', 'CD1E', 'KDM2B', 'SLCO5A1', 'FSCN1', 'FBXO27', 'NMRK1', 'BCL2L14', 'ETV3', 'RARRES2', 'ARHGAP10', 'CDKN1A', 'TSPAN33', 'CLIC2', 'CD1B',
    ],
    's_phase_sig': [
        # used the filtered list below to split MEBEMP-L cells in nimrod_atlas to highest PCNA sig and lowest PCNA sig. most diff expr genes:
        ###'HSPD1', 
        ###'TPI1', 
        ###'HSPA8', 
        'RAN', 
        ###'HSPE1', 
        ###'H2AFZ', # a bit strong and not clearly division related
        ###'PRKDC', # a bit strong and not clearly division related
        ###'PA2G4', # a bit strong and not clearly division related
        'NHP2', 'NME1', 
        ###'SNRPB', 
        ###'SNRPD1', 
        'NME4', 
        ###'SRM', 
        ###'EIF5A',
        ###'NAA38', 
        'RANBP1', 
        ###'ATP5MC1', 
        'DUT', 
        ###'FABP5', 
        'PAICS', 
        'MCM3',
        ###'C1QBP', 
        ###'CCDC85B', 
        'HELLS', 
        ###'MRPS34', 
        'MCM4', 
        ###'DCTPP1', 
        ###'TUBA1B', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        'DNMT1', 'PCNA', 'DNPH1', 
        ###'PRELID1', 
        ###'NOP58', 
        ###'EIF4EBP1', 
        'TYMS',
        'CDK4', 
        ###'TMPO', 
        'SMC1A', 
        ###'GCSH',
        
        # added manually
        'BRCA1',
        'BRCA2', # looked meh in leave_one_out?
        'MCM5',
        'MCM6',
        'MCM7',
        'MCM2',
        'CHEK1',
        'CLSPN',
        'PCLAF',
        # were higher in HSC_MPPs and also MEBEMP-Ls that were higher in the above genes, and protein atlas seemed to indicate they are plausibly expressed during proliferation.
        'GINS2',
        'ATAD5',
        'CDT1',
        'UNG',
        'FEN1',
        'DCTPP1',
        'CDCA7', # comment out?
        'GMNN',
        'USP1',
        'RPA3',
        'ORC6',
        ###'ERH', # a bit strong and not clearly division related
        'PRIM1', # DNA primase subunit 1
        'DHFR',
        'CENPH',
        'XRCC2',
        'LIG1', # DNA ligase 1
        ### 'MSH6', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        'SMC2',
        'CENPU',
        'CENPK',
        'RRM1',
        'CHAF1A',
        'POLE3', # DNA polymerase epsilon 3, accessory subunit
        'BRIP1',
        'LMNB1', # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4615282/: "lamin B1 (LMNB1) is required for proper chromosome condensation in interphase nuclei" and "LMNB1 repression slowed cellular growth due to S-phase delays and increased genomic instability"
        'DTL',
        'RAD51AP1',
        'NUSAP1',
        'CENPM',
        'SUPT16H',
        'SLBP',
        'UHRF1',
        'MCM10',
        'RNASEH2C',
        'NUDC',
        'FANCI',
        ###'XRCC6', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        'POLD2', # DNA polymerase delta 2, accessory subunit
        'CENPF',
        'CDC6',
        'RFC4',
        'WEE1',
        'RFC3',
        'TIMELESS',
        ###'RMI2', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        ###'RFC1', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        'MAD2L1',
        'ANAPC11',
        ###'KIF20B', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        ###'SSRP1', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        'MTHFD1',
        'ZWINT',
        ###'RFC2', # doesn't look good in a marker heatmap of nimrod_atlas MEBEMP-E and MPP
        'DTYMK',
        'DDX11',
        
        # # filtered the list below by manually examining a clustermap of these genes across MEBEMP-E and MEBEMP-L MCs in nimrod_atlas
        # 'MCM7','ATAD5','MCM6','UNG','NLN','GMNN','HMGN2','BRCA2','CENPK','CDT1','BRCA1','DTL','HSP90AA1','CHEK1','FEN1','CLSPN','GINS2','MCM2','WDR76','CD320','ORC6','ANKLE1','SMC1A','PAICS','TMEM14B','LMNB1','TMPO','DUT','MCM5','H2AFZ','PCNA','MCM3','DNMT1',
        # ###'PFN1',
        # 'PCLAF',
        # ###'FABP5',
        # 'HELLS','MCM4','TYMS','HSPE1','EIF5A','CCDC85B','NAA38','GPX4','TPI1','SRM','NOP58','DNPH1','MRPS34','NHP2','GCSH','PRELID1','SNRPD1','SNRPB','EIF4EBP1','NOLC1','DCTPP1','SLBP','CDK4','PRKDC','C1QBP','PA2G4','HSPD1','NCL','YBX1','ACTB','NME1','RANBP1','RAN','NME4','ODC1','HSPA8','ATP5MC1'
        
        # # corr > 0.35 with PCNA across nimrod_atlas MEBEMP-E and MEBEMP-L.
        # 'PCNA', 'TYMS', 'MCM4', 'GINS2', 'HELLS', 'CLSPN', 'MYC', 'HSPE1',
        # 'RANBP1', 'NME1', 'SNRPB', 'FABP5', 'FEN1', 'PCLAF', 'H2AFZ',
        # 'PRKDC', 'MCM3', 'HSPD1', 'ACTB', 'EIF5A', 'LMNB1', 'YBX1',
        # 'ATP5MC1', 'PFN1', 'MCM2', 'APOC1', 'TFRC', 'PA2G4', 'HNRNPAB',
        # 'DNPH1', 'KLF1', 'TST', 'SRM', 'NCL', 'DUT', 'BRCA1',
        # 'FAM45A', 'RAN', 'DTL', 'HMGN2', 'SRSF7', 'CHEK1', 'EIF4EBP1',
        # 'CDT1', 'DCTPP1', 'NPW', 'CD320', 'MCM5', 'NME4', 'MCM6',
        # 'HSP90AA1', 'MRPS26', 'DNMT1', 'NDUFAF8', 'ITGA2B', 'PRELID1',
        # 'PRMT1', 'ANKLE1', 'ORC6', 'KCNH2', 'GCSH', 'TMPO', 'WDR76',
        # 'PXMP2', 'CCT5', 'NOP58', 'CCDC85B', 'UQCRQ', 'HSPA8', 'FAM178B',
        # 'ODC1', 'UNG', 'MRPS34', 'NOLC1', 'SIVA1', 'TMEM14B', 'NAA38',
        # 'C1QBP', 'GPX4', 'ERH', 'NOP56', 'NDUFS6', 'SNRPD1', 'PAICS',
        # 'CENPK', 'PHPT1', 'TPI1', 'CST3', 'HBQ1', 'SLBP', 'ATAD5', 'NHP2',
        # 'SMC1A', 'DKC1', 'ANK1', 'NDUFB2', 'GMNN', 'BRI3BP', 'LY6E',
        # 'MRPL12', 'SUPT16H', 'NMNAT3', 'BRCA2', 'MCM7', 'COA4', 'MTHFD2',
        # 'MLEC', 'GPATCH4', 'UQCC2', 'NLN', 'CDK4', 'PRKAR2B', 'GSPT1',
        # 'SSB', 'MRTO4',
    ],
}

HIGHER_IN_CLP_E_THAN_MEBEMP_L = [
    # log ratio >1
    'SPINK2', 'HOPX', 'FAM30A', 'IGHM', 'CD52', 'GLIPR1', 'LTB',
    'MZB1', 'PRSS2', 'CD99', 'AVP', 'ITM2C', 'CD74', 'SELL', 'GIMAP7',
    'CD44', 'C1QTNF4', 'FLT3', 'BAALC', 'LST1', 'HLA-DPA1', 'SLC2A5',
    'TRBC2', 'HLA-DPB1', 'CSF3R', 'RCSD1', 'RBPMS', 'HLA-DRB1',
    'PRDX1', 'HCST', 'SORL1', 'HLA-DRA', 'IFITM3', 'NCOA7', 'GBP4',
    'SLC38A1', 'KIAA0087', 'GBP2',
    'FOS', # is this ok?
    'NRIP1', 'TSC22D3', 'BST2',
    'PROM1', 'AJ009632.2', 'KLF6', 'GUCY1A1', 'EGFL7', 'IRF1', 'LSP1',
    'ADGRG6', 'GSTP1', 'SPTBN1', 'PTPRC', 'SMIM24', 'DUSP6', 'CALHM6',
    'RNF125', 'CORO1A', 'IDS', 'PRAM1', 'CYTH4', 'PBXIP1', 'RERE',
    'EMP1', 'IQGAP2', 'VAMP5', 'COBLL1', 'HLA-DQB1', 'EVI2B', 'PLCB1',
    'MDK', 'TRIM22', 'TFPI', 'BIN1', 'SVIL', 'CLNK', 'JAML', 'HOXA9',
    'S100Z', 'GIMAP1', 'PLEKHO1', 'NPDC1', 'BCL11A', 'S100A11', 'TCF4',
    'GNAI1', 'BTG2', 'ICAM3', 
    
    'ADAM28', 'ALDH2', 'MME',
    'DNTT', 'ARHGAP25', 'RFLNB', 'ENO1', 'HLA-E', 'JCHAIN', 'ATP2B1',
    'SSBP2', 'HLA-DMB', 'ADA', 'GSN', 'RNASET2', 'TCF7L2', 'CD34',
    'ELK3', 'AIF1', 'SATB1', 'IFI16', 'TCTEX1D1', 'SPNS3', 'RNF213',
    'LAIR1', 'EVL', 'APP', 'KRT18', 'AC243960.1', 'ZFP36L2',
]
HIGHER_IN_GMP_L_THAN_CLP_E = [
    # log ratio >1
    'MPO', 'AZU1', 'CFD', 'ELANE', 'IGLL1', 'CPA3', 'PRTN3', 'RNASE2',
    'MGST1', 'SRGN', 'CTSG', 'CALR', 'LGALS1', 'ZEB2', 'CLEC12A',
    'RAB32', 'LYZ', 'MYC', 'CAT', 'CEBPD', 'TNFSF13B', 'PRSS57',
    'CLEC11A', 'LYST', 'AC084033.3', 'TUBA1B', 'HGF', 'FABP5',
    'SERPINB1', 'NCOA4', 'MS4A3', 'FNDC3B', 'APLP2', 'PCLAF', 'FAM45A',
    'P4HB', 'EIF3A', 'CST7', 'NPW', 'SUCNR1', 'MYB', 'HSP90B1',
    'CDCA7', 'F13A1', 'HSPA5', 'CSF3R', 'CSTA', 'TENT5A', 'KBTBD11',
    'RNASE3', 'PPIB', 'PPA1', 'IGFBP2', 'CPVL', 'CANX', 
    'VIM', # TODO: remove and update threshold
    'PLAC8',
    'NUCB2', 'CD38', 'MCM4', 'C1QTNF4', 'H2AFZ', 'PLEK', 'HSPD1',
    'XBP1', 'FLNA', 'NASP', 'MLEC', 'AC002454.1', 'CEBPA', 'GPI',
    'PDIA6', 'TYMS', 
    'GAPDH', # TODO: remove (i think) and update threshold
    'BRI3BP', 'HSPB1', 'IGFBP7', 'PKM',
    'CENPF', 'RAN', 'AL355922.1', 'NAA38', 'ACTN1',
]

HIGHER_IN_MEBEMP_L_THAN_CLP_E = [
    # log ratio >1
    'SLC40A1', 'CNRIP1', 'CPA3', 'HBD', 'FCER1A', 'AC084033.3', 'TFRC',
    'PDZD8', 'KLF1', 'CTNNBL1', 'CSF2RB', 'AL157895.1', 'S100A6',
    'MINPP1', 'PLEK', 'ZEB2', 'GATA1', 'GATA2', 'PRKAR2B', 'CASP3',
    'CD84', 'S100A4', 'ZMYND8', 'NAA38', 'IL1B', 'SERPINB1', 'MPP1',
    'PSTPIP2', 'CPPED1', 'ABCC4', 'PKIG', 'PDCD4', 'MYB', 'BMP2K',
    'RNF130', 'HPGDS', 'MYC', 'SNX5', 'TLN1', 'PCLAF', 'ATP6V0A2',
    'CYTOR', 'STXBP5', 'SRGN', 'FBXO7', 'TPM1', 'DNAJC9', 'STAM',
    'ZNF385D', 'MARCKSL1', 'PLIN2', 'CYTL1', 'CENPU', 'NFE2', 'NLK',
    'PBX1', 'DEPTOR', 'HLTF', 'TFR2', 'AL034397.3', 'GNAQ', 'TESPA1',
    'SLC27A2', 'ACSM3', 'ITGA2B', 'CITED2', 'PDLIM1', 'SOD2',
    'C2orf88', 'TAL1', 'GALNT1', 'MIR4435-2HG', 'STXBP6', 'DAD1',
    'EIF2AK1', 'PLXNC1', 'GDF11', 'MRPL52', 'RGS18', 'KIFAP3', 'DLC1',

    # log ratio 0.7-1

    # log ratio 0.6-0.7

    # log ratio 0.5-0.6
]

GENES_HIGHER_IN_HIGH_SUB1_B_COMPARED_TO_B_NAMES = [
    'SUB1', 'UBE2J1', 'SRGN', 'ELL2', 'DAD1', 'LGALS3', 'LGALS1',
    'PPIB', 'TNFRSF17', 'JCHAIN', 'TMEM258', 'CYTOR', 'SSR3', 'SEC61G',
    'AQP3', 'SEC61B', 'CD38', 'OSTC', 'DNAJB11', 'PDIA4', 'HSPA5',
    'P4HB', 'IGHA1', 'MAN1A1', 'MANF', 'HSP90B1', 'COX5A', 'CKAP4',
    'TENT5C', 'SEC14L1', 'MTHFD2', 'DNAJC3', 'SEMA4A', 'VIM', 'TMED9',
    'CALR', 'TXNDC11', 'PLA2G16', 'PRDM1', 'HYOU1',
]

B_GENE_NAMES = [
    'ID3',
    'LINC01781',
    'FCRL5', 'FCRL3', 'FCRL2',
    'AIM2',
    'FCRLA',
    'LINC01857',
    'SP140', 'SH3BP5', 'OSBPL10',
    'TLR10', 'ARHGAP24',
    'CPNE5', 'CD24',
    'PAX5',
    'MS4A1',
    'CD27', 'CLECL1',
    'HVCN1', 'GPR183',
    'LINC00926',
    'IRF8', # though maybe too weak?
    'TNFRSF13B', 'IKZF3',
    'FAM129C', 'CD22', 'CD79A', 'SPIB',
    'VPREB3', 'TNFRSF13C',
]

PALANTIR_GMP_MAX_CD79AB_EXPR = -13.7
PALANTIR_B_MAX_MPO_EXPR = -14

# BEST_SEPARATORS_BETWEEN_MS4A1_HIGH_AND_LOW_CLP_E = [
#     'MS4A1', 'IKZF3', 'AIM2', 'VPREB3', 'IGKC', 'MARCH1', 'CD79A',
#     'BANK1', 'SP140', 'FCRL5', 'ID3', 'CD24', 'CYBB', 'CD48', 'PAX5',
#     'BIRC3', 'RALGPS2', 'CD22', 'SMIM14', 'POU2F2',
# ]
BEST_SEPARATORS_BETWEEN_MS4A1_HIGH_AND_LOW_CLP_E = [
    'MS4A1', 
    'IGKC', # a very problematic one...
    'CD22', 'PAX5', 'CD79A', 'MARCH1', 'CD48',
    'IKZF3', 'LINC00926', 'CPNE5', 'RALGPS2', 'TNFRSF13C', 'LBH',
    'BANK1', 'CD74', 'P2RX5', 'HLA-DMB', 'ARHGAP24', 'AIM2', 'PARP15',
    'CD79B', 'FAM129C', 'FCRL2', 'CD40', 'SYNE2', 'SCIMP', 'OSBPL10',
    'CLEC2D', 'LINC01857', 'FCRL1', 'BLK', 'HVCN1', 'ISG20',
    'LINC02397', 'NCF1', 'TTN', 'TNFRSF13B', 'FCRLA', 'HLA-DQB1',
    'HLA-DQA1', 'CD24', 'CCDC50', 'CTSZ', 'FCRL5', 'RAB30', 'TRIM38',
    'SPIB', 'BIRC3', 'FCRL3', 'POU2AF1', 'ALOX5', 'IRF8', 'PNOC',
    'CCR6', 'ID3', 'EBF1', 'ALOX5AP', 'FGD2', 'POU2F2', 'TMSB4X',
]

GENES_HIGHER_IN_CLP_E_COMPARED_TO_B = [
    'LEPROT', 'SNCA', 'SNHG8', 'MSRB3', 'RNF24', 'ANKRD28', 'KIT',
    'ZFAS1', 
    # 'HSP90AB1', 
    'HOXA7', 'HMGA2', 'FOS', 'MED12L', 'DPPA4',
    'TM7SF3', 'ZNF521', 'HACD1', 'MAP1A', 'AVP', 'HOXA9', 'CTSW',
    'CNTLN', 'CREG1', 'RAB13', 'EVA1B', 'BST2', 'PLCB1', 'MAPKAPK3',
    'CD109', 'AMD1', 'TCEAL9', 'CCDC171', 'MLLT3', 'ALDH2', 'CSF3R',
    'HTR1F', 'CDK6', 'HDGFL3', 'GNA15', 'GSTP1', 'PRAM1', 'IGSF10',
    'TRIM24', 'HEMGN', 'HEBP1', 'BCAT1', 'GUCY1B1', 'FRMD4B', 
    # 'EIF3E',
    'ATP8B4', 'AJ009632.2', 'BZW2', 'ALDH1A1', 'DST', 'ANP32B', 'LDHB',
    'NUCB2', 'PRDX6',
    # 'RPLP0',
]
GENES_HIGHER_IN_B_COMPARED_TO_CLP_E = [
    'FCRL2', 'CD27', 'MAP3K1', 'CD79B', 'SMAP2', 'CYBB', 'SMC6',
    'PSMB9', 'FAM107B', 'CDK14', 'SIT1', 'CD19', 'ACP5', 'BLNK',
    'GRB2', 'PCED1B', 'BHLHE41', 'RUBCNL', 'POU2F2', 'BANK1', 'GGA2',
    'SLAMF6', 'EBF1', 'HLA-DOB', 'LIMD2', 'SPIB', 'ETS1', 'LINC02397',
    'CD40', 'MS4A1', 'ARHGAP24', 'GPR183', 'SP100', 'TNFRSF13C',
    'BTG1', 'PARP15', 'LINC00926', 'TLR10', 'ATM', 'MTSS1', 'HLA-DQA1',
    'TNFRSF13B', 'PAX5', 'COTL1', 'ID3', 'FCRLA', 'FCRL1', 'RHBDF2',
    'CSGALNACT1', 'CTSZ', 'CD180', 'TNFAIP8', 'TMSB4X', 'MYL12A',
    'MBD4', 'BACE2', 'GPR65', 'UGCG', 'CLECL1', 'RAB30',
]

BM_DC_SOMEWHAT_SPECIFIC_GENE_NAMES = [
    # started from IRF8, and added some genes which were very correlated with it.
    'IRF8', 'PLD4', 'SPIB', 'IRF7', 'RUNX2',

    'CCDC50', 'UGCG', 'JCHAIN', 
    ### 'TCF4', # very strong in some HSPCs
    'IGKC', 'C12orf75',
    'RNASE6', 'HERPUD1', 'APP', 'LILRA4', 'GAPT', 
    ###'SAMHD1', 
    ### 'PLEK', 
    'PLAC8', 'ALOX5AP', 'MPEG1',
       'GZMB', 
    ### 'CD74', # very strong in some HSPCs
       ###'LGALS1', 
       'FAM129C',
]
BM_GMP_L_GENE_NAMES = [
    'MPO', 'AZU1', 'PRTN3', 'CTSG', 'CST7', 'ELANE', 'RNASE2',
    'MS4A3', 
    # 'SRGN',

    #    'RNASE3', 'KBTBD11', 'SERPINB10', 'MS4A3', 'TENT5A', 'FUT4', 'HGF',
    #    'SUCNR1', 'CLEC11A', 'AL355922.1', 'CALR', 'SRGN', 'FNDC3B',
]
BM_HIGHER_IN_GMP_SERPINB8_THAN_GMP_L_MONO_DOUBLETS = [
       'SERPINB8', 'RUNX2', 'MKI67', 'HMMR', 'GPR183', 'TOP2A', 'NCAPG',
       'NUF2', 'SLC44A1', 'B3GALNT2', 'GTSE1', 'ATAD5', 'SKA2', 'GAPT',
       'GMNN', 'FKBP5', 'CENPW', 'CXorf21', 'CLSPN', 'SUZ12', 'PPM1G',
       'CENPF', 'ARHGAP11A', 'CXCR4', 'HIST1H1D', 'BRI3BP', 'ANP32E',
       'MAGOH', 'KIF15', 'HMGB2',
]
BM_LOWER_IN_GMP_SERPINB8_THAN_GMP_L_MONO_DOUBLETS = [
    'FCN1', 'CTSS', 'PRTN3', 'S100A8', 'SERPINA1', 'CD14', 'GIMAP4',
       'NEAT1', 'CYBB', 'CLEC7A', 'FYB1', 'FCGR2A', 'VCAN', 'IFITM3',
       'GIMAP7', 'FPR1', 'ELANE', 'NCF1', 'TYROBP', 'FGL2', 'S100A12',
       'CARD16', 'NCF2', 'S100A11', 'CALHM6', 'LILRB3', 'S100A9', 'AIF1',
       'CD300E', 'CX3CR1', 'CST7', 'CYP1B1', 'AZU1', 'DUSP6', 'PSAP',
       'RNF213', 'CD52', 
       ###'CPA3', 
       'MAFB', 'GBP2', 'TNFSF10', 'IRF1',
       'SAT1', 'MARCKS', 'LILRB2', 'TNFRSF1B', 'P2RY13', 'GBP1', 'FCGR3A',
       'VAMP5', 
       ###'SERPINB1', 'S100A6', 'CAMLG', 
       'CPVL', 'GBP5', 'LRRK2',
]



NON_RIBOSOMAL_VERY_HIGHLY_EXPRESSED_GENES = [
    # actually some of them, like FAU and UBA52, encode for fusion proteins which contain ribosomal proteins.
    # 'H3F3A', 'TMSB10', 'EEF1B2', 'PTMA', 'COX7C', 'HINT1', 'RACK1',
    # 'HLA-B', 'HLA-DRA', 'HLA-DPB1', 'EEF1A1', 'FTH1', 'FAU', 'NACA',
    # 'TPT1', 'B2M', 'EIF1', 'ATP5F1E', 'UBA52', 'FTL',


    'SH3BGRL3', 'CD52', 'YBX1', 'UQCRH', 'S100A6', 'S100A4', 'PTPRC',
       'H3F3A', 'OST4', 'CALM2', 'TMSB10', 'MZT2B', 'EEF1B2', 'ARPC2',
       'PTMA', 'TMA7', 'SERP1', 'CCNI', 'SNHG8', 'SUB1', 'BTF3', 'NSA2',
       'COX7C', 'HINT1', 'SKP1', 'NPM1', 'RACK1', 'SERPINB1', 'HLA-A',
       'HLA-E', 'HLA-B', 'CLIC1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1',
       'HLA-DPB1', 'HSP90AB1', 'EEF1A1', 'COX7A2', 'PNRC1', 'SNX3',
       'ACTB', 'NDUFA4', 'TOMM7', 'HNRNPA2B1', 'PPIA', 'CHCHD2',
       'SLC25A6', 'UXT', 'SH3BGRL', 'UQCRB', 'COX6C', 'EIF3E', 'EEF1D',
       'ANP32B', 'EIF3F', 'FTH1', 'FAU', 'CFL1', 'ATP5MG', 'GAPDH',
       'ARHGDIB', 'LDHB', 'PFDN5', 'ATP5MC2', 'HNRNPA1', 'MYL6', 'NACA',
       'ARPC3', 'ERP29', 'DYNLL1', 'UBC', 'HMGB1', 'TPT1', 'ITM2B',
       'COMMD6', 'HSP90AA1', 'SRP14', 'SERF2', 'B2M', 'COX4I1', 'CYBA',
       'PFN1', 'UBB', 'EIF1', 'DDX5', 'SUMO2', 'H3F3B', 'ACTG1', 'MYL12A',
       'MYL12B', 'ZFAS1', 'ATP5F1E', 'CIRBP', 'OAZ1', 'EEF2', 'UBL5',
       'UBA52', 'COX6B1', 'EIF3K', 'GMFG', 'SNRPD2', 'NOP53', 'FTL',
       'CD37', 'UQCR10', 'HMGN1',
]



# POTENTIAL_MKP_GENE_MODULE1 = [
#     'RGS18', 'DNM3', 'ARHGAP6', 'RAB27B', 'CMTM6', 'FYB1', 'KIF2A', 'LIMS1', 'FLI1', 'MMRN1', 'FERMT3', 'RAP1B', 'CAVIN2', 'SLC44A1', 'GUCY1B1', 'F2R', 'PECAM1', 'NRIP1', 'EIF4G3', 'TLN1', 'INPP4B', 'VCL', 'BIN2', 'MTPN', 'MED12L', 'PLXDC2',
# ]
# POTENTIAL_MKP_GENE_MODULE2 = [
#     'SLC9A9', 'ST6GAL2', 'PIM1', 'WASF1', 'PTK2', 'KCNK6', 'FYN', 'ILK', 'CLU', 'ADRA2A', 'SYNE3', 'CRLF3', 'SEPT11', 'TACC1', 'LRRFIP2', 'MGLL', 'LAPTM4B', 'LCP2', 'MAP3K5', 'RASA1', 'SOS1', 'PTBP3', 'CHD4', 'PPP1R12A', 'MSN', 'TPP2', 'PDLIM1', 'DAD1', 'DOK2', 'SPX', 'LGALSL', 'CLEC1B', 'PRICKLE1', 'LTBP1', 'ITGB3', 'RASGRP3', 'SELP', 'NFIB', 'CMAS', 'PRKAR2B', 'TPM1', 'MPIG6B', 'C2orf88', 'PSTPIP2', 'MEF2C', 'RUFY1', 'RRN3', 'SNCA', 'BANK1', 'STOM',
# ]


GENES_HIGHER_IN_MKP_COMPARED_TO_EP_AND_MEBEMP_L_AND_VERY_POS_CORRELATED_WITH_ITGA2B_ACROSS_MKPS_AND_LOW_IN_PLATELETS = [
    # 231004: MKP vs EP and MEBEMP-L log ratio > 0.8, correlation with ITGA2B across MKPs > 0.8, and platelet expression < -15.5.
    # also, for many genes here, high expression in platelets makes no sense, e.g. cell cycle.

    'RRN3', # wow good
    'MKI67', # wow good
    'TOP2A', # wow good
    'TPX2', # good
    # 'HSPA5', # good, but sound stress related?
    # 'DNAJB11', # good, but sound stress related?
    'SMC2', # good
    'HMGB3', # good
    'FAM98A', # good
    'CKAP2L', # good
    'TTC27', # good, but a bit too unknown function?
]
GENES_LOWER_IN_MKP_COMPARED_TO_OTHER_MEBEMPS = [
    'APOC1', # one of our best EP markers?
    'CA1', # one of our best EP markers?
    # 'HBD',
    'CNRIP1', 'CSF2RB', 'MYB', 'PNRC1',
]

# MKP_SPECIFIC_GENES_IF_IGNORING_PLATELETS = [
#     'ITGA2B', 
#     # 'CD9', # seems less precise than ITGA2B, but good enough
#     # 'LTBP1', # seems less precise than ITGA2B, but good enough
#     # 'RAB27B',
#     # 'LAT', 
#     # 'PPIF', 
    
#     # 'PLEK', # seems less specific than ITGA2B - not good enough because not lower in many MEBEMP-L MCs.
#     # 'PF4', # seems less specific than ITGA2B - not good enough because not lower in many MEBEMP-L MCs.
#     # 'TPM4', # seems less specific than ITGA2B - not good enough because not lower in many MEBEMP-L MCs.
#     # 'MMRN1', # seems less specific than ITGA2B - not good enough because not lower in many MEBEMP-L MCs.
#     # 'CD36', # seems less specific than ITGA2B - not good enough because not lower in many MEBEMP-L MCs.
# ]
GENES_AT_LEAST_256_FOLD_LOWER_IN_MKP_THAN_IN_PLATELETS = [
    'PPBP', 'TUBB1', 'GNG11', 'PTCRA', 
]

GENES_AT_LEAST_64_FOLD_LOWER_IN_MKP_THAN_IN_PLATELETS_AND_LOW_IN_MKP = [
    'ENKUR', 
    'PTCRA', # https://www.proteinatlas.org/ENSG00000171611-PTCRA: "Pre T cell antigen receptor alpha" and "Single cell type: expression clusteri: Platelets - Hemostasis (mainly)". i think it really is high in platelets, and even if it is high in T cells in our data, i guess it should still be low in MKP, assuming it can't be too high in ambient noise.
    'CLDN5', 
    'SMIM5', # sometimes higher than -16 in MKPs
    'TMEM40', # sometimes higher than -16 in MKPs
    'MAP3K7CL',
    'GRAP2', # seems correlated with ITGA2B.
    'CD226', # seems correlated with ITGA2B.
    'ACRBP', 
    'MARCH2', # seems correlated with ITGA2B.
    
    
    # 'HIST1H2BJ', 
    # 'PLEKHO1', 'MMD',
    #    'TREML1', 'SH3BP5', 'F13A1', 'SPARC',
]
NON_LATERAL_GENES_AT_LEAST_8_FOLD_HIGHER_IN_MKP_THAN_IN_PLATELETS = [
    'S100A4', 
    'SERBP1', 
    'S100A6', 
    'LSM7', 'PPIB', 'HNRNPAB', 'SNRPD2',
    'SNRPF', 'FABP5', 'NSA2', 'UQCRQ', 'SRSF7', 'VDAC1', 'MZT2B',
    'NDUFB2', 'PDIA6', 'EIF5B', 'CCT2', 'CCT3', 'SSBP1', 'ATP5MC1',
    'SLIRP', 'NHP2', 'AC084033.3', 'NME1', 'ANP32E', 'TIMM13',
    'IMPDH2', 'PPA1', 'GYPC', 'GADD45GIP1', 'PAICS', 'MZT2A', 'IL1B',
    'GGCT', 'TFAM',
]
MKP_SEPARATORS_AT_LEAST_4_TIMES_WEAKER_IN_PLATELETS = [
    'TWISTNB', 'PARP1', 'FAM98A', 'HACD1', 'SMC2', 'TPR', 'LDHB', 'HNRNPD', 'SPTA1', 'RRN3', 'LIMA1', 'ILF3', 'DCXR', 'NDUFAF8', 'SUPT16H', 'UROS', 'TTC27', 'CENPW', 'PRDX4', 'NOL7', 'PPM1G', 'GSPT1', 'PSMD1', 'DECR1', 'BMS1', 'DCUN1D5', 'RNF130', 'FKBP3', 'CYC1', 'BRCA1', 'SELENOH', 'ATP5MC3', 'MED12L', 'ATAD2', 'XRCC2', 'NAXE', 'CENPJ', 'CASP8AP2', 'PSMD7', 'EIF5A', 'NUCKS1', 'YARS', 'ST6GAL2', 'CKAP2L', 'SOS1', 'CLCC1', 'HNRNPC', 'UCHL5', 'CAMK1', 'HMGXB4', 'CENPM', 'RABL6', 'THRAP3', 'YBX1', 'STOML2', 'SMARCA5', 'BANF1', 'ZBTB11', 'SLC25A5', 'ITPKA', 'TIMM8B', 'SPOUT1', 'XRCC5', 'ETF1', 'LRRCC1', 'MRPL51', 'PIM1', 'KNL1', 'IDH2', 'ATP1A1', 'HNRNPA0', 'NFU1', 'COX7A2',  'CHCHD1',
]

MKP_SEPARATORS_AT_LEAST_2_TIMES_WEAKER_IN_PLATELETS = [
    # welch < 1e-6
    'TPR', 'LDHB', 'GSPT1', 'XRCC2', 'CENPJ', 'ST6GAL2', 'COX7A2',
    'HDDC3', 'PGRMC2', 'FBXO9', 'SLC44A1', 'UROD', 'ANP32B',

    # welch < 1e-5
    # 'HACD1', 'TPR', 'LDHB', 'RRN3', 'LIMA1', 'UROS', 'TTC27', 'GSPT1',
    # 'DCUN1D5', 'SELENOH', 'XRCC2', 'CENPJ', 'PSMD7', 'ST6GAL2', 'SOS1',
    # 'HNRNPC', 'CAMK1', 'HMGXB4', 'SMARCA5', 'ZBTB11', 'ETF1', 'PIM1',
    # 'IDH2', 'NFU1', 'COX7A2', 'HDDC3', 'NR2F6', 'DTD1', 'PGRMC2',
    # 'CDH6', 'RASGRP3', 'GMPS', 'ADCY6', 'COX8A', 'FBXO9', 'PTPN11',
    # 'SLC44A1', 'SMIM1', 'PDIA3', 'UROD', 'ANP32B', 'THUMPD2', 'TTK',
    # 'FADS2', 'CDC27',

    # welch < 1e-4
    # 'TWISTNB', 'PARP1', 'FAM98A', 'HACD1', 'SMC2', 'TPR', 'LDHB',
    # 'HNRNPD', 'SPTA1', 'RRN3', 'LIMA1', 'ILF3', 'DCXR', 'NDUFAF8',
    # 'SUPT16H', 'UROS', 'TTC27', 'CENPW', 'PRDX4', 'NOL7', 'PPM1G',
    # 'GSPT1', 'PSMD1', 'DECR1', 'BMS1', 'DCUN1D5', 'RNF130', 'FKBP3',
    # 'CYC1', 'BRCA1', 'SELENOH', 'ATP5MC3', 'MED12L', 'ATAD2', 'XRCC2',
    # 'NAXE', 'CENPJ', 'CASP8AP2', 'PSMD7', 'EIF5A', 'NUCKS1', 'YARS',
    # 'ST6GAL2', 'CKAP2L', 'SOS1', 'CLCC1', 'HNRNPC', 'UCHL5', 'CAMK1',
    # 'HMGXB4', 'CENPM', 'RABL6', 'THRAP3', 'YBX1', 'STOML2', 'SMARCA5',
    # 'BANF1', 'ZBTB11', 'SLC25A5', 'ITPKA', 'TIMM8B', 'SPOUT1', 'XRCC5',
    # 'ETF1', 'LRRCC1', 'MRPL51', 'PIM1', 'KNL1', 'IDH2', 'ATP1A1',
    # 'HNRNPA0', 'NFU1', 'COX7A2', 'CHCHD1', 'RNASEH2C', 'HDDC3',
    # 'PSMD8', 'NDUFS1', 'TCP1', 'NR2F6', 'ISOC1', 'DTD1', 'SUZ12',
    # 'SKA2', 'SLC25A3', 'PGRMC2', 'RTCA', 'RASA1', 'CPSF2', 'PSMC5',
    # 'SART3', 'CDH6', 'RASGRP3', 'MRPS18C', 'GSS', 'MRPS10', 'GMPS',
    # 'TPI1', 'PPP1CC', 'NDUFV2', 'MRPL17', 'ADCY6', 'SNUPN', 'COX7B',
    # 'TM9SF3', 'LAGE3', 'DNMT1', 'COX8A', 'FBXO9', 'PTPN11', 'CENPP',
    # 'COPS3', 'SLC44A1', 'SMIM1', 'COPS4', 'SRSF2', 'PDIA3', 'UROD',
    # 'CISD2', 'SRSF4', 'ANP32B', 'THUMPD2', 'TTK', 'FADS2', 'BAG5',
    # 'MESD', 'APOOL', 'EBAG9', 'CDC27', 'SRSF9'
    

    # welch < 1e-3
    # 'SNRPF', 'NME1', 'TWISTNB', 'SNRPD1', 'PARP1', 'CCT3', 'VDAC1',
    # 'ANP32E', 'PTTG1', 'CBX5', 'RPA3', 'FAM98A', 'HACD1', 'CENPE',
    # 'SMC2', 'TPR', 'LDHB', 'UQCRQ', 'HNRNPD', 'SPTA1', 'RRN3', 'LIMA1',
    # 'ILF3', 'DCXR', 'NDUFAF8', 'GPATCH4', 'SUPT16H', 'FAM136A',
    # 'SRSF7', 'C1QBP', 'UROS', 'HMGB3', 'TTC27', 'UQCC2', 'RAN',
    # 'CENPW', 'PPIH', 'PRDX4', 'NOL7', 'PPM1G', 'GSPT1', 'NUDCD2',
    # 'BCCIP', 'HNRNPAB', 'FAM162A', 'OXCT1', 'BRIX1', 'PSMD1', 'GLO1',
    # 'HNRNPM', 'DECR1', 'BMS1', 'DCUN1D5', 'NAA15', 'DDX46', 'SLIRP',
    # 'RNF130', 'SRPK1', 'DKC1', 'PSMA5', 'FKBP3', 'CYC1', 'SSB',
    # 'BRCA1', 'SF3B3', 'CCT5', 'HNRNPA3', 'KPNB1', 'PSMC3', 'SELENOH',
    # 'ZWINT', 'HMGN2', 'GOLGB1', 'ATP5MC3', 'MED12L', 'TAF15',
    # 'HSP90AA1', 'LTV1', 'SFPQ', 'ATAD2', 'XRCC2', 'PSMA4', 'EIF2S1',
    # 'NAXE', 'SNRPG', 'CENPJ', 'CKS1B', 'CASP8AP2', 'USP14', 'PSMD7',
    # 'SUCLA2', 'EIF5A', 'JPT1', 'LDHA', 'NUCKS1', 'ILF2', 'YARS',
    # 'TIMM17A', 'ST6GAL2', 'CKAP2L', 'SOS1', 'HSPA9', 'RAD23A',
    # 'MRPL27', 'CLCC1', 'HNRNPC', 'GTF2H3', 'NCAPG', 'PIEZO2', 'NDUFB8',
    # 'UCHL5', 'NABP2', 'METTL5', 'CAMK1', 'SNRPB', 'EIF4A1', 'HMGXB4',
    # 'CENPM', 'RABL6', 'RBL1', 'THRAP3', 'CACYBP', 'YBX1', 'WDR18',
    # 'MRPL13', 'STOML2', 'SMARCA5', 'GARS', 'BANF1', 'ZBTB11', 'GGH',
    # 'SLC25A5', 'CDKN3', 'ITPKA', 'SASS6', 'NUP37', 'TIMM8B', 'SPOUT1',
    # 'EIF4G1', 'ATP5F1B', 'CUL4A', 'XRCC5', 'RCC2', 'ETF1', 'LRRCC1',
    # 'GATA2-AS1', 'MRPL51', 'TAF9', 'CCT8', 'STIP1', 'MARS', 'EEF1E1',
    # 'NDUFS2', 'PIM1', 'CRLF3', 'RNASEH1', 'KNL1', 'IDH2', 'NDUFS6',
    # 'ATP1A1', 'HNRNPA0', 'DAZAP1', 'PSMB6', 'XRCC6', 'NFU1', 'AURKA',
    # 'COX7A2', 'CHCHD1', 'YEATS4', 'RNASEH2C', 'RNPS1', 'PARPBP',
    # 'SCCPDH', 'SET', 'COPS6', 'UBR1', 'HDDC3', 'TDG', 'COPB1', 'RCAN1',
    # 'PSMD8', 'NDUFS1', 'POLR2L', 'NUDT15', 'CNOT9', 'MRPL14', 'SRSF10',
    # 'TCP1', 'NR2F6', 'PTGR1', 'ANP32A', 'CSE1L', 'ISOC1', 'DTD1',
    # 'DEGS1', 'UQCRFS1', 'TARS', 'RASAL2', 'SUZ12', 'SKA2', 'SLC25A3',
    # 'CXCL2', 'AASDHPPT', 'PGRMC2', 'RTCA', 'COX6C', 'RASA1', 'CPSF2',
    # 'PSMC5', 'SUMO1', 'SART3', 'PPID', 'XPO1', 'TOE1', 'CDH6', 'H2AFV',
    # 'AIFM1', 'RASGRP3', 'MRPS18C', 'GSS', 'NFATC2IP', 'KAZN',
    # 'PRICKLE1', 'ITPA', 'MRPS10', 'THOC7', 'RPIA', 'GMPS', 'TPI1',
    # 'ANKRD35', 'PPP1CC', 'NDUFV2', 'ENY2', 'EIF4E2', 'MRPL17', 'ADCY6',
    # 'PSMG3', 'CENPV', 'PSMB2', 'SNUPN', 'COX7B', 'PRRC2C', 'TM9SF3',
    # 'LAGE3', 'DNMT1', 'COX8A', 'FBXO9', 'L2HGDH', 'PTPN11', 'CENPP',
    # 'COPS3', 'SLC44A1', 'KPNA3', 'PSMA3', 'SMIM1', 'COPS4', 'ENSA',
    # 'SRSF11', 'SRSF2', 'PDIA3', 'UROD', 'CISD2', 'SRSF4', 'OPA1',
    # 'ANP32B', 'PDIA5', 'THUMPD2', 'MLEC', 'DCTD', 'POLE2', 'PPP5C',
    # 'MMEL1', 'MRPL41', 'TTK', 'VCP', 'FADS2', 'BAG5', 'MESD', 'APOOL',
    # 'EBAG9', 'CDC27', 'C1orf43', 'SRSF9',
]

GENES_AT_LEAST_8_TIMES_STRONGER_IN_PLATELETS = [
    # 8 times stronger in platelets, compared to mean MKP metacell
    'PPBP', 'PTCRA', 'TMEM40', 'GNG11', 'CCL5', 'ACRBP', 'MAP3K7CL',
    'TUBB1', 'TREML1', 'CLDN5', 'GRAP2', 'F13A1', 'HIST1H2AC', 'MMD',
    'TUBA4A', 'ENKUR', 'NRGN', 'SMIM5', 'PDZK1IP1', 'MYL9', 'PDGFA',
    'MARCH2', 'AC147651.1', 'AP001189.1', 'CTTN', 'TMEM140', 'PLEKHO1',
    'PF4', 'CD226', 'FAM110A', 'SH3BP5', 'TNNC2', 'HIST1H2BJ', 'SMOX',
    'SNN', 'CLEC1B', 'VIM-AS1', 'TSPAN33', 'LGALSL', 'LGALS12',
    'PTPRJ', 'TSC22D1', 'DENND2C', 'DAB2', 'HIST2H2AA3', 'SAT1',
    'AQP10', 'DAPP1', 'GP9', 'HIST1H3H', 'FRMD3', 'CALD1', 'SCN1B',
    'C19orf33', 'RGS10', 'SPARC', 'RGCC', 'RGS18', 'TRAPPC3L',
    'TNFSF4', 
    'S100A8', 
    'NEXN', 
    'MFAP3L', 
    # 'HBA1', 
    'CDKN1A', 'CTSA',
    'HGD', 'LCN2', 'SLFN14', 'CXCL5', 'FAM107B', 'ABCC3', 'PARD3',
    'RHOBTB1', 'BEND2', 'SH3BGRL2', 'MISP3', 'FRMD4B', 'YPEL5',
    'IGF2BP3', 'PVALB', 'RAB31', 'TSPAN18', 'AL731557.1', 'RIPOR2',
    'S100A9', 
    'TUBA8', 'ZNF185', 'NCK2', 'CDC14B', 'AC090409.1',
    'CLCN3', 'ESAM', 'NCOA4', 'GAS2L1', 'PGRMC1', 'MFSD1', 'CLU',
    'TSC22D3', 'UBE2H', 'CAVIN2', 
    #    'HBG2', 

    'RGS3', 'CD99', 'PTPN18',
    'TNFSF13B', 'ODC1', 'GNAZ', 'HEXIM2', 'RNF11', 'ANKRD9', 'CCDC92',
    'CD68', 'BMP6', 'RIOK3', 'IFRD1', 'FAXDC2', 'SERPINE1', 'ALOX12',
    'CMTM5', 'CDKN2D', 'MYLK', 'PDE5A', 'YWHAH', 'RGS6', 'SPHK1',
    'IRX3', 'ELOVL7', 'P2RY12', 'RAB32', 'AKIRIN2', 'TSPAN9', 'CARD19',
    'SEC14L1', 'Z82206.1', 'ABLIM3', 'ITGB5', 'LINC00853', 'LDLRAP1',
    'TDRP', 'TMEM158', 'NT5C3A', 'WBP2', 'CD47', 'EHD3', 'AVPR1A',
    'DUSP22', 'JAM3', 'RUFY1', 'ACSBG1', 'LYPLAL1', 'SMIM3', 'THBS1',
    'ITM2B', 'GRK5', 'HIST1H4H', 'TMCC2', 'EPB41L3', 'ANO6', 'PIP4P2',
    'HBE1', 'CABP5', 'PCSK6', 'AP001636.3', 'KIAA0513', 'MYL12A',
    'AP003068.2', 'ZNF438', 'C12orf75', 'CD9', 'HPGD', 'TST', 'AMIGO2',
    'RAB30', 'RBM38', 'CD40LG', 'ATL1', 'MINDY1', 'ST3GAL6', 'ZGLP1',
    'MITF', 'GLA', 'REPS2', 'C15orf54', 'CA2', 'WDR11-AS1', 'HLA-E',
    'EMC3', 'DOK2', 'MAX', 'SMPD1', 'PADI4', 'FCER1G', 'FAM81B',
    'TLK1', 'BBC3', 'DMTN', 'INAFM2', 'HACD4', 'CAV2', 'TAX1BP3',
    'UBXN11', 'PDE4D', 'MTURN', 'DCLRE1A', 'SPNS1', 'ASAP1', 'GSTO1',
    'RASGRP2', 'TMEM70', 'CYB5R1', 'AC010542.5', 'EGLN3', 'NPTN',
    'GRHL1', 'HBQ1', 'AL031005.1', 'SSX2IP', 'RSPH9', 'MOB3C', 'RIT1',
    'SVIP', 'NSMCE3', 'C1orf198', 'H2AFJ', 'PEAR1', 'SWI5', 'MPP1',
    'CMIP', 'SUCNR1', 'KLF3', 'OAZ1', 'MLH3', 'ARHGAP18', 'MYZAP',
    'ADD3', 'XPNPEP1', 'TSPOAP1-AS1', 'SLC44A2', 'RAB11A',
    'AC104794.2', 'R3HDM4', 'VEPH1', 'RABGAP1L', 'AC074327.1',
    'PLA2G12A', 
    #    'HBB', 
    'OST4', 'WIPI1', 'C2orf88',
]

GENES_HIGHER_IN_MKP_COMPARED_TO_PLATELTS_AND_ALSO_EP_AND_MEBEMP_L = [
    # *MKP_SEPARATORS_AT_LEAST_4_TIMES_WEAKER_IN_PLATELETS,
    
    # ones higher in platelets were commented out (without further explanation)
    *[
        'TTC27', # what. but not strong.
        'TWISTNB',
        'CAMK1', # similar expression level in MKPs and monocytes. extremely good separator between MKPs and EPs and MEBEMP-Ls, but not strong.
        
        # added on 230307:
        'PSMD7',
        'LIMA1',
        'GSPT1',
        'HNRNPD', # unsure whether we want it here...
        'PRDX4',
        'SNRPG',
        'PTGR1',
        'ITPKA',
        *[
            # 'FERMT3', 
            # 'ST6GAL2', # not good separators from MEBEMP-Ls and EPs
            # 'MGLL', 
            # 'ADRA2A', # not good separators from MEBEMP-Ls and EPs
            # 'DNM3', 'ARHGAP6',
        ], 
        *[
            # 'C2orf88', 
            # 'PIM1', # not good separators from MEBEMP-Ls and EPs
            # 'CAVIN2', 
            # 'CAVIN1', # not good separators from MEBEMP-Ls and EPs
            # 'RAP1B', 'PSTPIP2',
        ], 
        *[
            # extremely good
            # 'STOM', # 4.5 times higher in platelets
            # 'MEF2C', # 1.25 times lower in platelets
            # 'SLC44A1', # excellent, but only around 2 times weaker in platelets
            # 'RUFY1', # 11.8 times higher in platelets
            'RRN3',
        ], 
    ],
    *[
        *[
            # 'SF3B3', 'SUZ12', 'DHX9', 'ILF3', 'PARP1', 'USP14', 'GTF3C4', # not good separators from MEBEMP-Ls and EPs
            # 'CLTC', # borderline. and 1.8 times lower in platelets
            # 'TCERG1', # 230305: not good separator from MEBEMP-Ls and EPs. i don't know whether i missed it yesterday or in the new model it isn't a good separator...
            # 'HSPA4', # https://www.proteinatlas.org/ENSG00000170606-HSPA4: "Heat shock protein family A (Hsp70) member 4" and "Stress response"
            # 'NSD2', # not good separators from MEBEMP-Ls and EPs
            # 'MTDH', # not separating from most BEMPs
        ],
        *[
            # 'TRA2B', 'TSR1', # not good separators from MEBEMP-Ls and EPs
            # 'PSMD1', # borderline
            # 'CCDC47', 'STIP1', 'HNRNPH3', 'PSMD2', 
            # 'PDIA3', # good, but only around 3 times weaker in platelets
            # 'SF3A3', 'DDX24', 'PKP4', 'RFC1', 'EZH2', 'LMNB2', # not good separators from MEBEMP-Ls and EPs
            # 'EIF5', # 1.05 times higher in platelets
            # 'TTF2', 'RABL6', # not good separators from MEBEMP-Ls and EPs
        ],
    ],
    *[
        *[
            # 'CENPM', 'CDC6', 'CHEK1', 'ZWINT', 'MAD2L1', 'CKS1B', 'CLSPN', 'RAD51AP1', 'CYCS', 'DSCC1', # not good separators from MEBEMP-Ls and EPs
            # 'ACTB', # 2.7 times higher in platelets
            # 'VDAC1', # not good separators from MEBEMP-Ls and EPs
            'TPI1', # a bit strong in platelets (only around 3 times weaker compared to MKPs). also, glycolysis, so feels not very right as a marker...
            'LDHA', 
        ],
    ],
    
    *[
        # all higher in MKPs compared to platelets.
        # borderline? half of the MKPs are higher than everything else.
        *[
            # 'PAXX', 'EIF4A3', 'LAGE3', 'DCXR', # not good separators from MEBEMP-Ls and EPs
            'MANF', 
            # 'NOL7', 'MRPS15', 'PSMG1', 'STOML2', 'PSMB5', 'HSPB11', # not good separators from MEBEMP-Ls and EPs
            # 'MT2A', # # a stress gene? TODO: maybe remove from our MKP signature?? unless it and all its friends are forbidden????
            # 'MRPL36', 'COX14', 'TMEM121', 'THOC7', 'COA4', 'TMEM147', 'SDF2L1', 'TIMM8B', 'TXNDC17', # not good separators from MEBEMP-Ls and EPs
        ],
    ],

    # unsure. similarly strong in a minority of EPs and MEBEMP-Ls with low platelet genes
    # *['CBX5', 'CENPF', 'ANP32E', 'NUCKS1', 'KIF20B', 'BRCA2', 'FANCI', 'XRCC6'],
    # *['NEK2', 'TPX2', 'BUB1', 'CCNB1', 'ESCO2', 'BUB1B', 'KNL1', 'KIF11', 'TROAP', 'KIF4A', 'GTSE1', 'NUF2', 'HJURP', 'FOXM1', 'CEP55', 'UBE2C', 'NCAPH', 'CIT', 'CDCA8', 'KIF2C', 'CENPE', 'CDCA2', 'DEPDC1', 'CKAP2L', 'CIP2A', 'FANCD2', 'SHCBP1', 'RTKN2', 'KIF23', 'ANLN', 'MELK', 'SKA1', 'HIST1H3G', 'KIFC1'],
]




PROBABLY_PLATELET_SPECIFIC_GENES = lateral_and_noisy_genes.PROBABLY_PLATELET_SPECIFIC_GENES


FIFTY_STRONGEST_PLATELET_GENES = lateral_and_noisy_genes.FIFTY_STRONGEST_PLATELET_GENES
    
STRONGEST_PLATELET_GENES = lateral_and_noisy_genes.STRONGEST_PLATELET_GENES

# STRONGEST_AND_ENRICHED_PLATELET_GENES = sorted(set(STRONGEST_PLATELET_GENES) - {'B2M','HLA-B','FTL','HMGB1','PTMA'}) # commented out because i didn't update the other set to keep only enriched ones.



UNCLEAR_N208_MODULES = [
#     'HPGDS',
#  'RHEX',
 'GBP4',
 'DPYD',
 'TNFRSF1A',
 'PPM1F',
 'DDHD1',
 'TAPBP',
 'PDE7A',
 
 'CD164',
 'UPF2',
 'UBTF',
 'CTCF',
 'BRD4',
 'RBBP6',
 'DHX36',
 'NEMF',
 'PPIG',
 'SF3B2',
 'NUB1',
 'STAT6',
 'RICTOR',
 'PIAS1',
 'CREBBP',
 'PRDM2',
 'DMTF1',
 'HELZ',
]


MITOCHONDRIALLY_ENCODED_GENE_NAMES = [
    'MT-ND1',
    'MT-ND2',
    'MT-CO1',
    'MT-CO2',
    'MT-ATP8',
    'MT-ATP6',
    'MT-CO3',
    'MT-ND3',
    'MT-ND4L',
    'MT-ND4',
    'MT-ND5',
    'MT-ND6',
    'MT-CYB',
]



RIBOSOMAL_PROTEIN_GENE_NAMES = lateral_and_noisy_genes.RIBOSOMAL_PROTEIN_GENE_NAMES

RIBOSOMAL_PROTEIN_GENE_NAME_PATTERNS = [
    '^RP[SL][0-9]*.*',
]
RIBOSOMAL_PROTEIN_GENE_NAME_EXCLUDE_PATTERNS = [
    '^RPS6K.*',
    # 'RPL22L1', # https://www.proteinatlas.org/ENSG00000163584-RPL22L1: "Ribosomal protein L22 like 1", and from https://en.wikipedia.org/wiki/RPL22L1 it seems like there is evidence for it having other functions, and it isn't as highly expressed as the classic ribosomal proteins (roughly 1/20 the expression of RPL22 (in partial data examined on 221207)). 230225: but it is salt and pepper in current model, so decided to not exclude it.
    'RPL34-AS1', # https://www.genecards.org/cgi-bin/carddisp.pl?gene=RPL34-DT: "RPL34 Antisense RNA 1". sounds like it isn't an RPL gene, and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221207).
    # 'RPL39L', # https://www.proteinatlas.org/ENSG00000163923-RPL39L: "Ribosomal protein L39 like" and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208). 230225: but it is salt and pepper in current model, so decided to not exclude it.
    # 'RPL26L1', # https://www.proteinatlas.org/ENSG00000037241-RPL26L1: "Ribosomal protein L26 like 1" and "This gene encodes a protein that shares high sequence similarity with ribosomal protein L26." and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208). 230225: but it is salt and pepper in current model, so decided to not exclude it.
    'RPS10-NUDT3', # https://www.proteinatlas.org/ENSG00000270800-RPS10-NUDT3: "RPS10-NUDT3 readthrough" and "This locus represents naturally occurring read-through transcription between the neighboring RPS10 (ribosomal protein S10) and NUDT3 (nudix (nucleoside diphosphate linked moiety X)-type motif 3) genes on chromosome 6. The read-through transcript produces a fusion protein that shares sequence identity with each individual gene product." and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208).
    'RPL36A-HNRNPH2', # https://www.proteinatlas.org/ENSG00000257529-RPL36A-HNRNPH2: "RPL36A-HNRNPH2 readthrough" and "This locus represents naturally occurring read-through transcription between the neighboring ribosomal protein L36a and heterogeneous nuclear ribonucleoprotein H2 (H') genes on chromosome X. The read-through transcript produces a protein with similarity to the protein encoded by the upstream locus, ribosomal protein L36a." and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208).
    'RPL3L', # https://www.proteinatlas.org/ENSG00000140986-RPL3L: "Ribosomal protein L3 like" and "This gene encodes a protein that shares sequence similarity with ribosomal protein L3. The protein belongs to the L3P family of ribosomal proteins. Unlike the ubiquitous expression of ribosomal protein genes, this gene has a tissue-specific pattern of expression, with the highest levels of expression in skeletal muscle and heart. It is not currently known whether the encoded protein is a functional ribosomal protein or whether it has evolved a function that is independent of the ribosome." and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208).
    'RPL10L', # https://www.proteinatlas.org/ENSG00000165496-RPL10L: "Ribosomal protein L10 like" and "This gene encodes a protein sharing sequence similarity with ribosomal protein L10. It is not currently known whether the encoded protein is a functional ribosomal protein or whether it has evolved a function that is independent of the ribosome." and "Testis-specific component of the ribosome, which is required for the transition from prophase to metaphase in male meiosis I (By similarity). Compensates for the inactivated X-linked RPL10 paralog during spermatogenesis" and it isn't expressed at all in partial data examined on 221208.
    'RPL7L1', # https://www.proteinatlas.org/ENSG00000146223-RPL7L1: "Ribosomal protein L7 like 1". and "Enables RNA binding activity. Predicted to be involved in maturation of LSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA). Predicted to act upstream of or within blastocyst formation. Predicted to be located in nucleolus. Predicted to be part of cytosolic large ribosomal subunit." but it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208).
    # 'RPS27L', # https://www.proteinatlas.org/ENSG00000185088-RPS27L: "Ribosomal protein S27 like" and "This gene encodes a protein sharing 96% amino acid similarity with ribosomal protein S27, which suggests the encoded protein may be a component of the 40S ribosomal subunit." but it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208). 230225: but it is salt and pepper in current model, so decided to not exclude it.
    'RPL17-C18orf32', # https://www.proteinatlas.org/ENSG00000215472-RPL17-C18orf32: "RPL17-C18orf32 readthrough" and "This locus represents naturally occurring read-through transcription between the neighboring RPL17 (ribosomal protein L17) and C18orf32 (chromosome 18 open reading frame 32) genes. Alternative splicing results in multiple transcript variants. The encoded isoforms share sequence identity with the RPL17 protein, but they include frameshifted C-terminal regions derived from the downstream gene exons." and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208).
    'RPS19BP1', # https://www.proteinatlas.org/ENSG00000187051-RPS19BP1: "Ribosomal protein S19 binding protein 1" and "Direct regulator of SIRT1. Enhances SIRT1-mediated deacetylation of p53/TP53, thereby participating in inhibition of p53/TP53-mediated transcriptional activity." and it isn't as highly expressed as the classic ribosomal proteins (in partial data examined on 221208).
    'RPL34-DT', # https://www.genecards.org/cgi-bin/carddisp.pl?gene=RPL34-DT: "RPL34-DT (RPL34 Divergent Transcript) is an RNA Gene, and is affiliated with the lncRNA class"
    'RPL37A-DT', # https://www.genecards.org/cgi-bin/carddisp.pl?gene=RPL37A-DT: "RPL37A-DT (RPL37A Divergent Transcript) is an RNA Gene, and is affiliated with the lncRNA class"
    'RPS27AP5', # https://www.genecards.org/cgi-bin/carddisp.pl?gene=RPS27AP5: "RPS27AP5 (Ribosomal Protein S27a Pseudogene 5) is a Pseudogene" though "Predicted to be located in ribosome". but obviously isn't one of the canonical ones. and was around -12.3 in oneK1K (doesn't appear in our ref genome, i think), so decided to not consider it as one of the ribo protein genes.

]
mds_in_out_dir_paths.PREPROCESSING_DIR_PATH


MDS_ANALYSIS_PARAMS = {
    'dataset_name_to_main_paths': {
        'illu_mds': {
            'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240618_pb_illu_mds_cytopenia_normal'),
            'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_illu_mds_cytopenia_normal/c.h5ad'),
            'preprocessing_out_dir_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_illu_mds_cytopenia_normal'),
            'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240522_all_exp_df.csv'),
        },
        'ult_mds': {
            'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240623_pb_ult_mds_cytopenia_normal'),
            'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_ult_mds_cytopenia_normal/c.h5ad'),
            'preprocessing_out_dir_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_ult_mds_cytopenia_normal'),
            'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240522_all_exp_df.csv'),
        },
    },

    # 'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240222_pb_exps_potentially_containing_mpn_illu'),
    # 'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_illu_mpn/c.h5ad'),
    # 'preprocessing_out_dir_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_illu_mpn'),
    # 'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240522_all_exp_df.csv'),
    
    # 'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240222_pb_exps_potentially_containing_mpn_ult'),
    # 'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_ult_mpn/c.h5ad'),
    # 'preprocessing_out_dir_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, 'pb_ult_mpn'),
    # 'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240430_all_exp_df.csv'),

    
    # 'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'all_data'),
    # 'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '230919_combined_experiments.h5ad',
    # 'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240227_all_exp_df.csv',
    
    # 'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'all_ult_ill_tech_reps'),
    # 'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240217_all_ult_ill_tech_reps.h5ad',
    # 'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240217_all_ult_ill_tech_reps_exp_df.csv',
    
    
    # 'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'new_N280_exps'),
    # 'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240129_new_N280_exps.h5ad',
    # 'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240129_new_N280_exp_df.csv',
    
    # 'out_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240129_exps_with_donors_with_any_sample_that_waited'),
    # 'raw_ad_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240129_exps_with_donors_with_any_sample_that_waited.h5ad',
    # 'preprocessing_out_exp_df_csv_file_path': os.path.join(mds_in_out_dir_paths.PREPROCESSING_DIR_PATH, '240129_exps_with_donors_with_any_sample_that_waited_exp_df.csv',
    
    
    'features_with_low_corr_across_bio_reps': [
        # 240617: r_value<0.5 across bio reps (calculated corr only for a subset of all features)
        'log_ratio_c_GMP-L_HSC_MPP',
        'mpp_mebemp_vim_sig_HSC_MPP_0.5',
        'mpp_mebemp_vim_sig_MEBEMP-L_0.5',
        's_phase_sig_BEMP_0.5',
        'log_c_NKTDP',
        'eryp_plek_sig_ERYP_0.5',
        'log_c_high_MPP_sig_ERYP',
        'log_c_pro-B?',
        'clp_vim_sig_CLP_0.5',
        'gmp_l_elane_sig_GMP-L_0.5',
        'nktdp_samhd1_sig_NKTDP_0.5',
        'gmp_l_ahnak_sig_GMP-L_0.5',
        'nktdp_lmna_sig_NKTDP_0.5',
        'nktdp_maf_sig_NKTDP_0.5',
        'log_c_DC-prog?',
        'log_c_CFD_tryptase_GMP-L'
    ],
    
    'pb_sig_default_min_donor_cell_count': 10,
    'pb_sig_c_states_repr_to_min_donor_cell_count': {
        'CLP': 20,
        'BEMP': 20,
        'MEBEMP-L': 100,
        'HSC_MPP': 100,
        'ERYP': 20,
    },
    'pb_weak_sig_min_donor_cell_count': 50,
    'pb_sig_quantiles': [0.05, 0.1, 0.5, 0.9, 0.95],

    'pb_sig_name_to_relevant_states': {
        'c_intermediate_nc_monocyte_somewhat_specific': ['Monocytes', 'all HSPCs'],
        'lamp3_outlier_specific': ['mregDC', 'all HSPCs'],
        'nkt_somewhat_specific': ['NKT', 'all HSPCs'],
        'higher_in_b_than_all_other_hspcs': ['B', 'all HSPCs'],
        'higher_in_gmp_l_than_all_other_hspcs': ['all HSPCs'],
        'higher_in_dc_than_all_hspcs': ['DC', 'all HSPCs'],
        'strong_c_monocyte_specific': ['Monocytes', 'all HSPCs'],
        'higher_in_c_dc2_than_all_hspcs': ['DC', 'all HSPCs'],
        'endothel_somewhat_specific': ['Endothel', 'all HSPCs'],
        'ighm_ighg': ['B', 'all HSPCs'],
        'higher_in_pre_b_than_pro_b': ['pre-B?', 'pro-B?'],
        'higher_in_pre_b_than_all_hspcs': ['pre-B?', 'pro-B?', 'all HSPCs'],
        'higher_in_b_than_pre_b': ['B', 'pre-B?', 'pro-B?'],
        'plasmablast_somewhat_specific': ['B', 'all HSPCs'],
        'higher_in_bemp_than_clp_m_and_nktdp': ['BEMP', 'CLP-M', 'CLP-L', 'NKTDP'],
        'higher_in_nktdp_than_clp': ['CLP-M', 'CLP-L', 'NKTDP'],
        'higher_in_nktdp_than_all_myeloid_hspcs': ['NKTDP', 'HSC', 'MPP', 'all myeloid HSPCs'],
        'higher_in_clp_m_than_all_myeloid_hspcs': ['CLP-M', 'CLP-L', 'HSC', 'MPP', 'all myeloid HSPCs'],
        'higher_in_plasmablast_than_clp': ['B', 'CLP-M', 'CLP-L'],
        'higher_in_bemp_than_gmp_l': ['BEMP', 'GMP-L'],
        'higher_in_mkp_than_mebemp_l_and_eryp': ['ERYP', 'MEBEMP-L'],
        'higher_in_bemp_than_mebemp_l_and_eryp': ['ERYP', 'MEBEMP-L', 'BEMP'],
        'higher_in_mpp_than_ep': ['HSC', 'MPP', 'MEBEMP-L', 'MEBEMP-E', 'ERYP'],
        'higher_in_ep_than_mpp': ['HSC', 'MPP', 'MEBEMP-L', 'MEBEMP-E', 'ERYP'],
        'higher_in_mebemp_l_than_mpp': ['HSC', 'MPP', 'MEBEMP-L', 'MEBEMP-E'],
    },
    'col_name_to_name_in_paper': {
        'blood_sample_id': 'Blood sample ID',
        'nimrod_compatible_donor_id': 'Individual ID',
        'donor_sex': 'sex',
        'donor_age': 'age',
        'bleeding_date_as_date': 'scRNA sample blood draw date',
        'facs_bm_date_repr': 'BM examination date',
        'scrna_exp_id': '10x library ID',
        'diagnosis_class': 'Diagnosis',
        'in_orig_ref_model': 'scRNA sample also used in original reference model',
        'used_to_fit_hspc_compos_knn': 'Used to calculate composition abnormality score',
        'in_fig4_ref_model': 'In Figure 4 reference model',
        'in_mds_cyto_groups_analysis': 'In Figure 4B-E',
        'in_train_set': 'In MDS classification cohort 1',
        'in_test_set': 'In MDS classification cohort 2',
        'in_blast_vs_clp_e_scatter': 'In Figure 4H',
        'cbc_date': 'CBC date',
        'gene': 'Gene',
        'CHR': 'Chromosome',
        'POS': 'Position',
        'REF': 'Reference allele',
        'ALT': 'Alternate allele',
        'mean_VAF': 'Average VAF',

        'first_gene_end_pos_in_chrom': 'Approximate start position',
        'last_gene_end_pos_in_chrom': 'Approximate end position',

        'log_c_BEMP': 'log2(BEMP frequency)',
        'log_c_CLP': 'log2(CLP frequency)',
        'log_c_ERYP': 'log2(ERYP frequency)',
        'log_c_GMP-L': 'log2(GMP-L frequency)',
        'log_c_HSC_MPP': 'log2(HSC_MPP frequency)',
        'log_c_MEBEMP-L': 'log2(MEBEMP-L frequency)',
        'log_c_MKP': 'log2(MKP frequency)',
        'log_c_NKTDP': 'log2(NKTDP frequency)',
        'log_c_high_MPP_sig_ERYP': 'log2(high_MPP_sig_ERYP frequency)',
        'log_c_pre-B?': 'log2(pre-B? frequency)',
        'log_c_pro-B?': 'log2(pro-B? frequency)',

        'num_BEMP': 'BEMP cell count (by metacell annotation)',
        'num_EP': 'ERYP cell count (by metacell annotation)',
        'num_MEBEMP.L': 'MEBEMP-L cell count (by metacell annotation)',
        'num_MEBEMP.E': 'MEBEMP-E cell count (by metacell annotation)',
        'num_GMP.E': 'GMP-E cell count (by metacell annotation)',
        'num_MPP': 'MPP cell count (by metacell annotation)',
        'num_HSC': 'HSC cell count (by metacell annotation)',
        'num_CLP.E': 'CLP-E cell count (by metacell annotation)',
        'num_CLP.M': 'CLP-M cell count (by metacell annotation)',
        'num_CLP.L': 'CLP-L cell count (by metacell annotation)',
        'num_NKTDP': 'NKTDP cell count (by metacell annotation)',

        'LMNA_mebemp_score': 'LMNA signature in MEBEMPs',
        'LMNA_clp_score': 'LMNA signature in CLPs',
        'sync_score': 'sync-score',

        'max_mean_VAF': 'maximal CH VAF',
        'mean_near_neighbor_dist': 'composition abnormality score',
        'cna_count': '#CNAs',
        'under_manual_exp_min_umi_count_count': '#cells excluded due to low original UMI count',
        'too_few_non_excluded_umis_count': '#cells excluded due to low UMI count after gene exclusion',
        'too_high_mito_and_malat1_count': '#cells excluded due to high mitochonrial and MALAT1 expression',
        'too_low_ribo_prot_count': '#cells excluded due to low ribosomal expression',
        'hba_hbb_contaminated_count': '#cells excluded due to RBC contamination',
        'neut_or_neut_contaminated_count': '#cells excluded due to neutrophil contamination',
        'platelet_or_platelet_contaminated_count': '#cells excluded due to platelet contamination',
        'too_low_norm_nuclear_expr_count': '#cells excluded due to low nuclear expression',
        'mebemp_monocyte_doublet_count': '#cells excluded as they appear to be MEBEMP/monocyte doublets',
        'enriched_genotype_doublet_mc_count': '#cells excluded due to appearing in a doublet-enriched metacell',
        'high_log_norm_umi_count_count': '#cells excluded due to high original UMI count',
        'final_c_hspc_count': '#HSPCs',
        
        'seq_platform': 'Sequencing platform',

        'manual_comment_for_blood_aging_paper': 'Notes',
    },
    'sig_name_to_name_in_paper': {
        'higher_in_dc_than_all_hspcs': 'DC signature',
        'nkt_somewhat_specific': 'NKT signature',
        'higher_in_b_than_all_other_hspcs': 'B signature',
        'endothel_somewhat_specific': 'Endothel signature',
        'higher_in_pre_b_than_pro_b': 'pre-B pro-B DEGs',
        'higher_in_pre_b_than_all_hspcs': 'pre-B signature',
        'higher_in_clp_m_than_all_myeloid_hspcs': 'CLP signature',
        'higher_in_nktdp_than_all_myeloid_hspcs': 'NKTDP signature',
        'higher_in_nktdp_than_clp': 'NKTDP CLP DEGs',
        'higher_in_bemp_than_mebemp_l_and_eryp': 'BEMP MEBEMP-L/ERYP DEGs',
        'higher_in_gmp_l_than_all_other_hspcs': 'GMP-L signature',
        'higher_in_mpp_than_ep': 'MPP ERYP DEGs',
        'higher_in_mebemp_l_than_mpp': 'MEBEMP-L MPP DEGs',
        'higher_in_ep_than_mpp': 'ERYP MPP DEGs',
        'higher_in_mkp_than_mebemp_l_and_eryp': 'MKP MEBEMP-L/ERYP DEGs',

        'gmp_l_ahnak_sig': 'GMP-L AHNAK signature',
        'gmp_l_elane_sig': 'GMP-L ELANE signature',
        'nktdp_lmna_sig': 'NKTDP LMNA signature',
        'nktdp_maf_sig': 'NKTDP MAF signature',
        'nktdp_samhd1_sig': 'NKTDP SAMHD1 signature',
        'clp_id2_runx2_sig': 'CLP ID2/RUNX2 signature',
        'clp_dntt_sig': 'CLP DNTT signature',
        'prss2_sig': 'PRSS2 signature',
        'clp_vim_sig': 'CLP VIM signature',
        'mpp_mebemp_vim_sig': 'MPP/MEBEMP VIM signature',
        'eryp_apoc1_sig': 'ERYP APOC1 signature',
        'eryp_plek_sig': 'ERYP PLEK signature',
        'bemp_cd74_sig': 'BEMP early signature',
        'bemp_ms4a2_sig': 'BEMP MS4A2 signature',
        's_phase_sig': 'S-phase signature',
        'mpp_slc40a1_sig': 'MPP SLC40A1 signature',
        'mebemp_l_hla_sig': 'MEBEMP-L MHC-II signature',
    },
    'pb_sig_name_to_info': {
        # c_state_sets_to_calc_statistics defaults to [set(relevant_states)]
        'gmp_l_ahnak_sig': {
            'relevant_states': ['GMP-L'],
        },
        'gmp_l_elane_sig': {
            'relevant_states': ['GMP-L'],
        },
        'nktdp_lmna_sig': {
            'relevant_states': ['NKTDP'],
        },
        'nktdp_maf_sig': {
            'relevant_states': ['NKTDP'],
        },
        'nktdp_samhd1_sig': {
            'relevant_states': ['NKTDP'],
        },
        # 'nktdp_fam30a_sig': { # commented out, but it looks like NKTDP MCs high in nktdp_fam30a_sig are contaminated with non-NKTDP cells...
        #     'relevant_states': ['NKTDP'],
        # },
        'clp_id2_runx2_sig': {
            'relevant_states': ['CLP'],
            'nimrod_148_relevant_states': ['CLP-M', 'CLP-L'],
        },
        'clp_dntt_sig': {
            'relevant_states': ['CLP'],
            'nimrod_148_relevant_states': ['CLP-M', 'CLP-L'],
        },
        'prss2_sig': {
            'relevant_states': ['CLP'],
            'nimrod_148_relevant_states': ['CLP-M', 'CLP-L'],
        },
        'clp_vim_sig': {
            'relevant_states': ['CLP'],
            'nimrod_148_relevant_states': ['CLP-M', 'CLP-L'],
        },
        'mpp_mebemp_vim_sig': {
            'relevant_states': ['MPP', 'MEBEMP-E', 'MEBEMP-L'],
            'c_state_sets_to_calc_statistics': [
                {'HSC_MPP'},
                {'MEBEMP-L'},
            ],
        },
        'eryp_apoc1_sig': {
            'relevant_states': ['ERYP'],
        },
        'eryp_plek_sig': {
            
            'relevant_states': ['ERYP'],
        },
        'bemp_cd74_sig': {
            'relevant_states': ['BEMP'],
        },
        'bemp_ms4a2_sig': {
            'relevant_states': ['BEMP'],
        },
        's_phase_sig': {
            'relevant_states': [
                ['MPP', 'MEBEMP-E', 'MEBEMP-L'],
                ['MEBEMP-L', 'ERYP'],
                ['BEMP'],
            ],
            'c_state_sets_to_calc_statistics': [
                # {'HSC_MPP', 'MEBEMP-L'},
                {'HSC_MPP'},
                {'MEBEMP-L'},
                {'ERYP'},
                {'BEMP'}, # better without it??
            ],
        },
        # 'higher_in_mpp_than_mebemp_l': {
        #     'relevant_states': ['MPP', 'MEBEMP-E', 'MEBEMP-L'],
        #     'c_state_sets_to_calc_statistics': [
        #         {'HSC_MPP', 'MEBEMP-L'},
        #     ],
        # },
        'mpp_slc40a1_sig': {
            'relevant_states': ['MPP'],
            'c_state_sets_to_calc_statistics': [
                {'HSC_MPP'},
            ],
        },
        'mebemp_l_hla_sig': {
            'relevant_states': ['MEBEMP-L'],
        },
    },
    
    'constant_random_seed_for_reproducibility': 1234,
    'donor_id_to_comments': {
        'N227': [
            'Higher S100A4 overall?',
            'High HBA1 and HBG2 in some MKPs?',
        ],
    },
    'hg_genes_gtf_file_path': os.path.join(mds_in_out_dir_paths.HUMAN_GENOME_DIR_PATH, 'genes.gtf'),
    'hg_fasta_file_path': os.path.join(mds_in_out_dir_paths.HUMAN_GENOME_DIR_PATH, 'genome.fa'),
    
    'arch_mutations': {
        'minimal_mutation_df_csv_file_path': os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240801_minimal_mutation_df.csv'),
        'arch4_mutation_file_paths': [
            (os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240621_ARCH_SNPs_June2024.tsv'), 'hg38'),
            (os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240623_ARCH_Indels_June2024.csv'), 'hg38'),
            (os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240625_ARCH_CALR_SNP_MPN.csv'), 'irrelevant_because_only_gene_names'),
        ],
        'mostly_arch3_mutation_file_paths': [
            (os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240727_hg37_nature_med_individual_mutations_from_nimrod.csv'), 'hg37'),
            (os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240727_hg38_nature_med_individual_mutations_from_nimrod.csv'), 'hg38'),
        ],
        'other_arch_mutation_file_paths': [
            (os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240621_CALR_by_nili.xlsx'), 'hg37'),
            (os.path.join(mds_in_out_dir_paths.ARCH_MUTATION_DIR_PATH, '240621_SRSF2.xlsx'), 'hg37'),
        ],
        
        'mds_irrelevant_shlush_donor_ids': [
            'NS10', 'NS17', 'NS21', 'NS23', 'NS24', 'NS6', 'NS13', 'NS4', 'NS9',
        ],
        'healthy_atlas_irrelevant_shlush_donor_ids': [
            'N55', 
        ],
    },

    
    'gene_exclusion': {
        # nimrodra forbidden genes (sent to me on 221117): mitochondrial (^MT-, ^MTRNR), MALAT1, NEAT1, XIST
        'min_num_of_total_umis_per_gene': 1,
        'excluded_gene_names': [
            'MALAT1', # highest normalized UMIs (around 3e-2). localized to the nucleus. when comparing batches, looks somewhat batchy

            # 'NEAT1', # normalized UMIs around 6e-4. localized to the nucleus. doesn't look particularly batchy when compraing batches.

            *MITOCHONDRIALLY_ENCODED_GENE_NAMES,
        ],
        'excluded_gene_name_patterns': [
            # Why do I see a high level of mitochondrial gene expression? (https://kb.10xgenomics.com/hc/en-us/articles/360001086611):
            # "High expression levels of mitochondrial genes could be an indicator of: (1) Poor sample quality,
            # leading to a high fraction of apoptotic or lysing cells. (2) Biology of the particular sample, for example, tumor biopsies,
            # may have increased mitochondrial gene expression due to metabolic activity and/or necrosis."
            '^MT-.*',
            # '^MTRNR.*', # MT-RNR2 Like xxx (Pseudogene - according to gene cards) - I guess all these are pseudogenes of MT-RNR2, which on 220712 I didn't see at all in multiple batches (I guess because 10x ignores rRNA, including mitochondrial rRNA). also, almost all of genes seem to me to be lowly expressed, so removing them from the excluded list.
        ],
        # 'examining_bursty_lonely_genes': { # 230706: seems like it is now automatically done by mc.pl.exclude_genes, or was it always automatic??
        #     'downsample_min_samples': 750, # the default on metacells v0.9.0
        #     'downsample_min_cell_quantile': 0.05, # the default on metacells v0.9.0
        #     'downsample_max_cell_quantile': 0.5, # the default on metacells v0.9.0
        #     # 221207: found no noisy lonely genes in some of the batches of the MDS data.
        # },


        'old': {
            'top_20_genes_used_to_compare_with_potential_genes_to_exclude': [
                'RPS8', 'RPL32', 'RPL34', 'RPS3A', 'RPS23', 'EEF1A1', 'RPS12', 'RPL39', 'RPL10', 'RPS6',
                'RPS24', 'RPL41', 'TPT1', 'RPLP1', 'RPL13',
                'MALAT1', 'MT-CO1', 'MT-ATP6', 'MT-CO3', 'MT-ND4',
            ],
            'gene_modules_of_most_prevalent_genes_containing_excluded_genes': {
                'mod_MT-CO2': [
                    'GNB1', # https://www.proteinatlas.org/ENSG00000078369-GNB1: "G protein subunit beta 1"
                    'PPP3CA', # https://www.proteinatlas.org/ENSG00000138814-PPP3CA: "Protein phosphatase 3 catalytic subunit alpha" and "Calcium-dependent, calmodulin-stimulated protein phosphatase which plays an essential role in the transduction of intracellular Ca(2+)-mediated signals"
                    'NIPBL', # https://www.proteinatlas.org/ENSG00000164190-NIPBL: "cohesin loading factor" and "Plays an important role in the loading of the cohesin complex on to DNA. Forms a heterodimeric complex (also known as cohesin loading complex) with MAU2/SCC4 which mediates the loading of the cohesin complex onto chromatin 1, 2. Plays a role in cohesin loading at sites of DNA damage."
                    'AP3B1', # https://www.proteinatlas.org/ENSG00000132842-AP3B1: "Adaptor related protein complex 3 subunit beta 1"
                    'MEF2C',
                    'NR3C1',
                    'RIPOR2',
                    'IKZF1',
                    'SRPK2',
                    'MKLN1',
                    'DOCK11',
                    'KDM2A',
                    'UVRAG',
                    'APBB1IP',
                    'ATF7IP',
                    'WDFY2',
                    'MYCBP2',
                    'TBC1D22A', # https://www.proteinatlas.org/ENSG00000054611-TBC1D22A: "TBC1 domain family member 22A" and "May act as a GTPase-activating protein for Rab family protein(s)." and "Molecular function (UniProt)i	GTPase activation"
                    'MT-CO2',
                ],
                'mod_MALAT1': [
                    'CD52', # https://www.proteinatlas.org/ENSG00000169442-CD52: "May play a role in carrying and orienting carbohydrate, as well as having a more specific role." - seems somewhat lymphoid specific.
                    'MALAT1',
                    # looked at its gene_umi_normalized_count_dist, and it neither looks particularly batchy, nor extremely highly expressed (say, around 4e-4).
                    'DDX5', # https://www.proteinatlas.org/ENSG00000108654-DDX5: "DEAD-box helicase 5" and "Involved in the alternative regulation of pre-mRNA splicing". https://www.proteinatlas.org/ENSG00000108654-DDX5/subcellular: "Main locationi Localized to the Nucleoplasm (supported), Nucleoli (supported)"
                    'CD37', # https://www.proteinatlas.org/ENSG00000104894-CD37. https://www.proteinatlas.org/ENSG00000104894-CD37/single+cell+type: "Single cell type expression clusteri	B-cells - Immune response (mainly)"
                ],
                'mod_MT-ATP6': [
                    'MT-ATP6',
                    'MT-CO3',
                    'MT-ND3',
                    'MT-ND4',
                    'MT-ND5',
                    'MT-CYB',
                ],
                'mod_MTRNR2L12': [
                    'OST4', # https://www.proteinatlas.org/ENSG00000228474-OST4: "Oligosaccharyltransferase complex subunit 4" and "Subunit of the oligosaccharyl transferase (OST) complex that catalyzes the initial transfer of a defined glycan (Glc(3)Man(9)GlcNAc(2) in eukaryotes) from the lipid carrier dolichol-pyrophosphate to an asparagine residue within an Asn-X-Ser/Thr consensus motif in nascent polypeptide chains, the first step in protein N-glycosylation."
                    'ARL6IP5', # https://www.proteinatlas.org/ENSG00000144746-ARL6IP5: "ADP ribosylation factor like GTPase 6 interacting protein 5" and "Regulates intracellular concentrations of taurine and glutamate. Negatively modulates SLC1A1/EAAC1 glutamate transport activity by decreasing its affinity for glutamate in a PKC activity-dependent manner. Plays a role in the retention of SLC1A1/EAAC1 in the endoplasmic reticulum."
                    'MTRNR2L12', # https://www.genecards.org/cgi-bin/carddisp.pl?gene=MTRNR2L12: "MT-RNR2 Like 12 (Pseudogene)"
                    'MTPN', # https://www.proteinatlas.org/ENSG00000105887-MTPN: "Myotrophin" and "The transcript produced from this gene is bi-cistronic and can encode both myotrophin and leucine zipper protein 6. The myotrophin protein is associated with cardiac hypertrophy, where it is involved in the conversion of NFkappa B p50-p65 heterodimers to p50-p50 and p65-p65 homodimers. This protein also has a potential function in cerebellar morphogenesis, and it may be involved in the differentiation of cerebellar neurons, particularly of granule cells."
                    'MT-ND1',
                    'MT-ND2',
                    'MT-CO1',
                    'MT-ATP8',
                ],
                'mod_NEAT1': [
                    'LAPTM5',
                    'S100A10',
                    'S100A11',
                    'S100A6',
                    'S100A4',
                    'LST1',
                    'IFITM3',
                    'NEAT1',
                    'NPC2',
                    'EVI2B',
                    'CST3',
                    'FTL',
                    'TSPO',
                ],
            },

        },

    },
    'downsampling_after_excluding_genes': {
        'downsample_min_samples': 750, # the default on metacells v0.9.0
        'downsample_min_cell_quantile': 0.05, # the default on metacells v0.9.0
        'downsample_max_cell_quantile': 0.5, # the default on metacells v0.9.0
    },
    'cell_exclusion': {
        # 'max_mitochondrial_to_ribosomal_protein_gene_umi_log_ratio': np.inf,
        # 'max_mitochondrial_to_ribosomal_protein_gene_umi_log_ratio': 1,

        'min_num_of_non_excluded_umis': 500, # NOTE (IMPORTANT): if you increase this, you will have to update which exps cannot be assumed to have unbiased composition.
        # why ribo_and_malat1 rather than all excluded genes? because what if we have IGLL1 as excluded? or maybe some ribosomals excluded because they had high ult-ill fold change?? we shouldn't penalize cells for that is if it is the same as mito. i am not even sure treating malat1 and mito the same is ok...
        'max_mito_and_malat1_umi_frac': 0.2,
        
        'min_ribo_prot_umi_frac': 0.1,
    },
    'bursty_lonely_genes_allowed_to_be_excluded': [
        
    ],
    'bursty_lonely_genes_not_allowed_to_be_excluded': [
        'HBB', # i assume it would be useful for ambient noise removal. and anyway i mark it as lateral and noisy, and it is only rarely extremely high (and i am probably going to remove these cells anyway).
    ],
    'lateral_gene_names_that_might_be_missing_from_cells_ad': lateral_and_noisy_genes.LATERAL_GENES_THAT_MIGHT_BE_MISSING_FROM_C_AD,
    
    # old metacell version - target is num of UMIs. in the new one it is num of cells.
    # 'target_num_of_umis_in_metacell_in_all_donor_model': int(320e3),
    'target_num_of_umis_in_metacell_in_all_donor_model': int(400e3), 
    # 'target_num_of_umis_in_metacell_in_all_donor_model': int(700e3), # 230629: already 5k metacells when i ask for 450k...
    'target_num_of_umis_in_metacell_in_single_donor_model': int(160e3),
    # 'target_num_of_umis_in_metacell_in_single_donor_model': int(100e3),
    'target_num_of_umis_in_metacell_in_single_donor_small_mc_model': int(80e3),
    # 'target_num_of_cells_in_metacell_in_all_donor_model': 100,
    'target_num_of_cells_in_metacell_in_all_donor_model': 140,
    'target_num_of_cells_in_metacell_in_single_donor_model': 60,
    'target_num_of_cells_in_metacell_in_single_donor_small_mc_model': 30,

    'max_num_of_cells_for_smaller_mcs': int(5e3),
    'min_cell_count_to_compute_mcs': 200,

    # 'forbidden_genes_from_nimrod_csv_file_path': '/dummy/dummy/dummy/tanay_group/mds/221117_nimrod_forbidden_blood_aging_genes.txt',
    'gene_batch_kruskal_pvals_from_nimrod_csv_file_path': os.path.join(
        mds_in_out_dir_paths.MISC_INPUT_DIR_PATH, '230309_10x_batch_kruskal_pvals_from_nimrod.csv'),
    
    'extra_lateral_gene_names_when_skipping_mcnoise': [
        # TODO: maybe should always be noisy and lateral, if indeed it might be such droplets rather than ambient noise...
    ],
    'lateral_gene_names_file_paths': [
        
        # '/dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/intermediate_output/230128_batchy_genes_to_mark_as_lateral.txt', # before excluding ribosomal protein genes

        # AC_AND_AL_CLONE_BASED_GENE_NAMES_FILE_PATH, # NOTE: decided against marking these as lateral. only if they come up as potentially lateral, i think we should mark them as such...
    ],
    'noisy_gene_names_file_paths': [
    ],
    'lateral_gene_name_patterns': [
        '^IGKV.*',
        '^IGHV.*',
        '^IGHJ.*',
        '^IGLJ.*',
        '^IGLV.*',
        '^IGKC.*',
        '^IGHA.*',
        '^IGHD.*',
        '^IGHE.*',
        '^IGHG.*',
        '^IGHM.*',
        '^IGLC.*',
        '^HLA-.*',
    ],
    'lateral_gene_names': lateral_and_noisy_genes.LATERAL_GENE_DICT,
    'noisy_gene_names': lateral_and_noisy_genes.NOISY_GENE_DICT,
    'suspect_lateral_gene_name_patterns': [
        '^ISG.*',

        # from https://github.com/tanaylab/metacells/blob/master/vignettes/Metacells_Vignette.ipynb
        '^IFI.*',
        'MCM[0-9]',
        'SMC[0-9]',
    ],

    'analyzing_feature_genes': {
        'gene_log_normalized_umis_similarity': {
            'epsilon_to_add_to_fraction_before_log': 1e-5,
            'linkage_method': 'ward',
            'linkage_metric': 'euclidean',
        },
    },
    'projection': {
        'project_to_healthy_blood_aging_atlas': True,
        'use_correction': True,
    },
    'epsilon_for_mitochondrial_to_ribosomal_protein_gene_umi_log_ratio': 1e-5,
    'epsilon_to_add_to_mc_gene_norm_expression': 1e-5,
    'epsilon_to_add_to_mc_gene_mean_expression': 1e-5,
    'epsilon_for_donor_log_c_state_freq': 0.02,
    'epsilon_for_facs_blast_log': 0.005,
    'max_bm_facs_blast_10x_sample_days_dist': 365,
    'min_hspc_count_for_composition_analysis': 500,
    'nearest_neighbor_count_for_composition_analysis': 4,
    # 'ugly_mcview_config_csv_file_path': '/dummy/dummy/dummy/tanay_group/mds/mcview_proj_config.csv',

    'single_donor_mc_projection': {
        'max_ultima_illumina_log_ratio_range': 0.5,
        'epsilon_for_log_projected_fold': 1e-5,
        'num_of_gc_bins': 10,
        'low_log_norm_expression_thresh': -15.5,
        'my_projected_correlation_thresh': 0.85,
    },

    'all_cna_attrs_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'all_cna_attrs_df.csv'),
    'donor_agg_cna_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'donor_agg_cna_df.csv'),

    'all_numbered_donor_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'all_numbered_donor_df.csv'),
    'all_ext_donor_feature_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'all_ext_donor_feature_df.csv'),
    'final_feature_names_and_types_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'final_feature_names_and_types_df.csv'),
    'with_vaf_classifier_final_feature_names_and_types_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'with_vaf_classifier_final_feature_names_and_types_df.csv'),
    'fig4_donor_composition_class_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'fig4_donor_composition_class_df.csv'),
    'usage_in_analysis_masks_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240803_usage_in_analysis_masks_df.csv'),

    'filtering_cells_before_computing_metacells': {
        'cells_to_discard_paths_for_final_models': [
        ],
        'for_removing_doublet_metacells': {
            'donors_to_keep_attributes': {
                'assigned_and_doublets', # i.e., discard unassigned
            },
        },
        'for_ambient_noise_estimation': {
            # NOTE: dont add specific_exps_to_exclude here that you just want to not use MCNoise on but do want them later (simply without using MCNoise to clean them). because if you add exps here in specific_exps_to_exclude, they will be excluded from all later analysis. instead, specify them as exps_to_exclude_from_mcnoise_and_continue_with_them_as_if_they_have_no_ambient_noise.
            'donors_to_keep_attributes': {
                'assigned', # i.e., discard unassigned and doublets
            },
        },
        # 'for_ambient_noise_estimation_without_accidents': {
        #     'donors_to_keep_attributes': {
        #         'assigned', # i.e., discard unassigned and doublets
        #     },
        #     'specific_exps_to_exclude': [
        #         'demux_28_11_21_1',
        #         'demux_07_02_22_1',
        #     ],
        # },
        'mk_only': {
            'donors_to_keep_attributes': {
                'assigned', # i.e., discard unassigned and doublets
            },
        },
        'final_only_specific_donors_template': {
        },
        'final_all_unhealthy': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'unhealthy', # this implicitly discards unassigned and doublets
                # 'normal', # this implicitly discards unassigned and doublets
            },
            'specific_donors_to_exclude': [
            ],
            'max_num_of_cells_to_keep_as_a_factor_of_median': 4,
        },
        'final_all_unhealthy_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'unhealthy', # this implicitly discards unassigned and doublets
                # 'normal', # this implicitly discards unassigned and doublets
            },
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],
            'max_num_of_cells_to_keep_as_a_factor_of_median': 4,
        },
        'final_all_normal_pb_except_problematics_and_outliers': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'normal', # this implicitly discards unassigned and doublets
            },
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
                *SPECIFIC_OUTLIER_EXP_DONORS_NAME_PAIRS,
            ],
            'only_pb': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_bm': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_exps_to_exclude': [x for x in SPECIFIC_EXPS_TO_EXCLUDE_NAMES if '_bm_' in x],
            'specific_exp_donors_to_exclude': [
                *[x for x in SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS if '_bm_' in x[0]],
            ],
            'only_bm': True,

            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_healthy_bm': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', # not "healthy", because N365 and N196 don't have a "normal" diagnosis...
            },
            'specific_donors_to_include': [
                'N365', 
                'N196',
            ],
            'specific_exps_to_exclude': [x for x in SPECIFIC_EXPS_TO_EXCLUDE_NAMES if '_bm_' in x],
            'specific_exp_donors_to_exclude': [
                *[x for x in SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS if '_bm_' in x[0]],
            ],
            'only_bm': True,

            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'pb_non_aml_donors_with_meg3_cells': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_donors_to_include': [
                'N212',
                'N198',
                'N251',
                'N265',
                'N262',
                'N245',
            ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],
            'only_pb': True,
            # 'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_pb_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_donors_to_exclude': [
            ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],
            'only_pb': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_pb_except_problematics_and_outliers': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_donors_to_exclude': [
            ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
                *SPECIFIC_OUTLIER_EXP_DONORS_NAME_PAIRS,
            ],
            'only_pb': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_cytopenic_pb_except_problematics_and_outliers': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'cytopenic', 
            },
            'specific_donors_to_exclude': [
                'N199', # TODO: remove when all files are updated so N199 is PV everywhere.
                'N334', # TODO: remove when all files are updated so N334 is not MDS everywhere.
            ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
                *SPECIFIC_OUTLIER_EXP_DONORS_NAME_PAIRS,
            ],
            'only_pb': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_mpn_and_associated_normal_pb_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'mpn_and_normal', 
            },
            # 'specific_donors_to_exclude': [
            # ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
                # *SPECIFIC_OUTLIER_EXP_DONORS_NAME_PAIRS,
            ],
            'only_pb': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_ultima_pb_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_donors_to_exclude': [
            ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],
            'only_pb': True,
            'only_ultima': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_illumina_pb_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_donors_to_exclude': [
            ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],
            'only_pb': True,
            'only_non_ultima': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_nimrod_atlas_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_donors_to_exclude': [
                '130', '151', '48', # not healthy. (N130, N151, N48)
                '154', # not healthy B cells? (N154)
            ],
            'only_pb': True,
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_and_doublets_pb_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned_and_doublets', # i.e., discard unassigned
            },
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],
            'only_pb': True,
            # 'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_and_doublets_bm_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned_and_doublets', # i.e., discard unassigned
            },
            'specific_exps_to_exclude': [x for x in SPECIFIC_EXPS_TO_EXCLUDE_NAMES if '_bm_' in x],
            'specific_exp_donors_to_exclude': [
                *[x for x in SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS if '_bm_' in x[0]],
            ],
            'only_bm': True,
            # 'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_and_doublets_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned_and_doublets', # i.e., discard unassigned
            },
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],
            # 'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned_except_problematics': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', 
            },
            'specific_donors_to_exclude': [
                # 'N279', # want to see if i see novel monocytes states also without it. yes. saw them. good.
                # 'N328', # JMML. and does mix with others (i.e., not only metacells dominated by her)
                
            ],
            'specific_exps_to_exclude': SPECIFIC_EXPS_TO_EXCLUDE_NAMES,
            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS,
            ],

            'specific_exp_donors_to_exclude': [
                *SPECIFIC_EXP_DONORS_TO_EXCLUDE_NAME_PAIRS
            ],
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_assigned': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', # i.e., discard unassigned and doublets
            },
            'specific_exps_to_exclude': [
                # TODO: remove this. shouldnt start with a cell matrix without BM...
                # 'N257_bm_06_10_22_1_illumina',
                # 'N279_bm_26_10_22_1_illumina',
            ],
            
            # 'specific_exp_donors_to_exclude': [ # TODO: remove this. not sure though. because of num_of_non_excluded_umis_per_barcode.png. choosing here to discard accidents with lowest num of non excluded UMIs.
            #     ('demux_09_06_22_accident', 'N251'),
                
            #     ('demux_15_12_22_accident', 'N328'),
                
            #     # ('demux_14_11_22_accident', 'N286'),
            #     # ('demux_14_11_22_accident', 'N288'),
            #     # ('demux_14_11_22_accident', 'N290'),
            # ],

            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_all_except_other_exp_condition': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', # i.e., discard unassigned and doublets
                'except_other_exp_condition',
            },
            'specific_exps_to_exclude': [

            ],
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_normal_pb_atlas': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', # i.e., discard unassigned and doublets
            },
            'specific_exps_to_exclude': [

            ],
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
        'final_mds_cyto_normal_excluding_atlas': {
            # 'cells_to_discard_paths': [

            # ],
            'donors_to_keep_attributes': {
                'assigned', # i.e., discard unassigned and doublets
            },
            'specific_exps_to_exclude': [

            ],
            'max_num_of_cells_to_keep_as_a_factor_of_median': 3,
        },
    },

    'hca_bm_cd34_filtered': {
        
        'cells_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_akhiad/HCA_immun/scrna_db/mc1_converted_to_h5ad/HCA_org_hsc_mat.h5ad',
        
        

    },
    'hca_bm': {
        # copied from ~obk/data/HCA.BM.ATLAS/few.320k/
        # 'cells_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_orenBK/HCA.BM.ATLAS/few.320k/cells.h5ad',
        'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_orenBK/HCA.BM.ATLAS/few.320k/converted_c.h5ad', # created by using 231213_convert_0.8_to_0.9.py on /dummy/dummy/dummy/raid/mds/from_orenBK/HCA.BM.ATLAS/few.320k/cells.h5ad and /dummy/dummy/dummy/raid/mds/from_nimrod/bm_data/metacells.h5ad
        # 'metacells_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_orenBK/HCA.BM.ATLAS/few.320k/metacells.h5ad',
        
        
        # 'metacells_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/bm_data/metacells.h5ad',
        'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/bm_data/converted_mc.h5ad', # created by using 231213_convert_0.8_to_0.9.py on /dummy/dummy/dummy/raid/mds/from_orenBK/HCA.BM.ATLAS/few.320k/cells.h5ad and /dummy/dummy/dummy/raid/mds/from_nimrod/bm_data/metacells.h5ad
    },
    'onek1k_pbmc': {
        # data from Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease | Science (https://www.science.org/doi/10.1126/science.abf3041, 2022)
        # "We present the OneK1K cohort, which consists of single-cell RNA sequencing (scRNA-seq) data from 1.27 million peripheral blood mononuclear cells (PMBCs) collected from 982 donors"
        # from supplementary materials: "Following library preparation using the published 10x Chromium Single Cell 3’ V2 Solution protocol". so i guess this is V2.

        # orig data obs had a column 'assay' with a single val: "10x 3' v2"
        # orig data obs had a column 'self_reported_ethnicity' with a single val: "European"
        # orig data obs had a column 'observation_joinid' which seemed like a unique identifier, so i removed it.

        # orig data uns['citation']: 'Publication: https://doi.org/10.1126/science.abf3041 Dataset Version: https://datasets.cellxgene.cziscience.com/e1239e77-732c-4dc6-9a4c-c9bb72607c61.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1'

        'c_ad_file_path': '/dummy/dummy/dummy/raid/OneK1K_pbmc/final_c_240207.h5ad', 
        'mc_ad_file_path': '/dummy/dummy/dummy/raid/OneK1K_pbmc/final_mc_240207.h5ad', 

        'atlas_cell_type_colors_csv_file_path': DEFAULT_CELL_STATE_COLORS_CSV_FILE_PATH,

        'cell_state_and_info_list': [
            ('NK', {
                'rules': [
                    [('expr', 'CD3D', -np.inf, k1k_pb_mc_score_threshs.CD3D_THRESH)],
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.NK_NKG7_THRESH, np.inf)],
                    [('expr', 'CD79A', -np.inf, k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH)],
                ],
            }),
            ('naive-CD8-T', {
                'rules': [
                    [('expr', 'CD8A', k1k_pb_mc_score_threshs.NAIVE_CD8_T_CD8A_THRESH, np.inf)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.NAIVE_CD8_T_NKG7_THRESH)],
                    [('expr', 'CD79A', -np.inf, k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH)],
                ],
            }),
            ('effector-memory-CD8-T', {
                'rules': [
                    [('expr', 'CD3D', k1k_pb_mc_score_threshs.CD3D_THRESH, np.inf)],
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.NKT_NKG7_THRESH, np.inf)],
                    [('expr', 'CD79A', -np.inf, k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH)],
                ],
            }),
            ('GZMK-high-NKG7-mid-T', {
                # Granzyme K+ CD8 T cells form a core population in inflamed human tissue (https://www.science.org/doi/10.1126/scitranslmed.abo0686, 2022): "GzmK is typically considered to be transiently expressed in memory CD8 T cell populations and subsequently down-regulated as CD8 T cells differentiate toward the GzmB+ effector cell phenotype seen in cytotoxic T lymphocytes (CTLs) (11, 12)"
                'rules': [
                    [('expr', 'CD3D', k1k_pb_mc_score_threshs.CD3D_THRESH, np.inf)],
                    [('expr', 'GZMK', k1k_pb_mc_score_threshs.GZMK_HIGH_NKG7_MID_T_GZMK_THRESH, np.inf)],
                    [('expr', 'CCL5', k1k_pb_mc_score_threshs.GZMK_HIGH_NKG7_MID_T_CCL5_THRESH, np.inf)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.NKT_NKG7_THRESH)],
                    [('expr', 'CD79A', -np.inf, k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH)],
                ],
            }),
            ('naive-CD4-T', {
                'rules': [
                    [('expr', 'CD3D', k1k_pb_mc_score_threshs.CD3D_THRESH, np.inf)],
                    [('expr', 'CD8A', -np.inf, k1k_pb_mc_score_threshs.NAIVE_CD8_T_CD8A_THRESH)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.NAIVE_MEMORY_CD4_NKG7_THRESH)],
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4693724/: "We report here that S100A4 is specifically expressed by memory T cells"
                    [('expr', 'S100A4', -np.inf, k1k_pb_mc_score_threshs.MEMORY_CD4_S100A4_THRESH)],
                    [('expr', 'HBB', -np.inf, k1k_pb_mc_score_threshs.HBA_HBB_CONTAMINATED_HBB_THRESH)],
                    # [('expr', 'CD14', -np.inf, pb_mc_score_threshs.CD14_NKG7_DOUBLET_CD14_THRESH)],
                    [('expr', 'CD79A', -np.inf, k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH)],
                    [('expr', 'GZMK', -np.inf, k1k_pb_mc_score_threshs.CD4_T_GZMK_THRESH)],
                    [('expr', 'GZMA', -np.inf, k1k_pb_mc_score_threshs.CD4_T_GZMA_THRESH)],
                    [('expr', 'CCL5', -np.inf, k1k_pb_mc_score_threshs.CD4_T_CCL5_THRESH)],
                ],
            }),
            ('memory-CD4-T', {
                'rules': [
                    [('expr', 'CD3D', k1k_pb_mc_score_threshs.CD3D_THRESH, np.inf)],
                    [('expr', 'CD8A', -np.inf, k1k_pb_mc_score_threshs.NAIVE_CD8_T_CD8A_THRESH)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.NAIVE_MEMORY_CD4_NKG7_THRESH)],
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4693724/: "We report here that S100A4 is specifically expressed by memory T cells"
                    [('expr', 'S100A4', k1k_pb_mc_score_threshs.MEMORY_CD4_S100A4_THRESH, np.inf)],
                    [('expr', 'CD79A', -np.inf, k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH)],
                    [('expr', 'FCGR3A', -np.inf, k1k_pb_mc_score_threshs.CD4_T_FCGR3A_THRESH)],
                    [('expr', 'GZMK', -np.inf, k1k_pb_mc_score_threshs.CD4_T_GZMK_THRESH)],
                    [('expr', 'GZMA', -np.inf, k1k_pb_mc_score_threshs.CD4_T_GZMA_THRESH)],
                    [('expr', 'CCL5', -np.inf, k1k_pb_mc_score_threshs.CD4_T_CCL5_THRESH)],
                ],
            }),
            ('naive-B', {
                'rules': [
                    # [('expr', 'LTB', pb_mc_score_threshs.B_LTB_THRESH, np.inf)],
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.B_CD79A_THRESH, np.inf)],
                    [('expr', 'TCL1A', k1k_pb_mc_score_threshs.NAIVE_B_TCL1A_THRESH, np.inf)],
                    [('expr', 'CD3D', -np.inf, k1k_pb_mc_score_threshs.B_CD3D_THRESH)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.B_NKG7_THRESH)],
                    [('expr', 'CCDC88A', -np.inf, k1k_pb_mc_score_threshs.CCDC88A_THRESH)],
                ],
            }),
            ('memory-B', {
                'rules': [
                    # [('expr', 'LTB', pb_mc_score_threshs.B_LTB_THRESH, np.inf)],
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.B_CD79A_THRESH, np.inf)],
                    [('expr', 'TCL1A', -np.inf, k1k_pb_mc_score_threshs.NAIVE_B_TCL1A_THRESH)],
                    [('expr', 'CD3D', -np.inf, k1k_pb_mc_score_threshs.B_CD3D_THRESH)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.B_NKG7_THRESH)],
                    [('expr', 'ZEB2', -np.inf, k1k_pb_mc_score_threshs.AGED_B_ZEB2_THRESH)],
                    [('expr', 'MZB1', -np.inf, k1k_pb_mc_score_threshs.PLASMABLAST_MZB1_THRESH)],
                    [('expr', 'HBB', -np.inf, k1k_pb_mc_score_threshs.HBA_HBB_CONTAMINATED_HBB_THRESH)],
                ],
            }),
            ('CCDC88A-high-B', {
                'rules': [
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.B_CD79A_THRESH, np.inf)],
                    [('expr', 'CCDC88A', k1k_pb_mc_score_threshs.CCDC88A_THRESH, np.inf)],
                ],
            }),
            ('aged-B?', {
                'rules': [
                    # [('expr', 'LTB', pb_mc_score_threshs.B_LTB_THRESH, np.inf)],
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.B_CD79A_THRESH, np.inf)],
                    [('expr', 'ZEB2', k1k_pb_mc_score_threshs.AGED_B_ZEB2_THRESH, np.inf)],
                    # [('expr', 'VPREB3', -np.inf, pb_mc_score_threshs.AGED_B_VPREB3_THRESH)],
                    [('expr', 'CD3D', -np.inf, k1k_pb_mc_score_threshs.B_CD3D_THRESH)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.B_NKG7_THRESH)],
                ],
            }),
            ('plasmablast', {
                'rules': [
                    [('expr', 'LTB', -np.inf, k1k_pb_mc_score_threshs.PLASMABLAST_LTB_THRESH)],
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.B_CD79A_THRESH, np.inf)],
                    [('expr', 'MZB1', k1k_pb_mc_score_threshs.PLASMABLAST_MZB1_THRESH, np.inf)],
                    [('expr', 'CD3D', -np.inf, k1k_pb_mc_score_threshs.B_CD3D_THRESH)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.B_NKG7_THRESH)],
                ],
            }),
            ('pDC', {
                # NOTE: split to high and low PTGDS? i guess not. https://www.ncbi.nlm.nih.gov/gene?Term=5730 says LCNL1 is a potential readthrough! so maybe PTGDS counts as a bursty-lonely?
                'rules': [
                    # GZMB is nice, but also strong in some T cells...
                    [('expr', 'IRF8', k1k_pb_mc_score_threshs.P_DC_IRF8_THRESH, np.inf)],
                    [('expr', 'TCF4', k1k_pb_mc_score_threshs.P_DC_TCF4_THRESH, np.inf)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.IRF8_NKG7_DOUBLET_NKG7_THRESH)],
                    # [('expr', 'CD79A', pb_mc_score_threshs.CD79A_THRESH, np.inf)],
                    # [('expr', 'JCHAIN', pb_mc_score_threshs.PLASMABLAST_JCHAIN_THRESH, np.inf)],
                ],
            }),
            ('cDC?', {
                # NOTE: Classical DC2 subsets and monocyte-derived DC: Delineating the developmental and functional relationship (https://onlinelibrary.wiley.com/doi/10.1002/eji.202149548, 2023) "the functional relationship between cDC and MoDC is not fully understood, as the overlapping phenotypes of certain type 2 cDC (cDC2) subsets and MoDC do not allow satisfactory distinction of these cells in the tissue, particularly during inflammation" though "This review will revise murine cDC2 and MoDC biology".
                # NOTE: split to high and low PTGDS? i guess not. https://www.ncbi.nlm.nih.gov/gene?Term=5730 says LCNL1 is a potential readthrough! so maybe PTGDS counts as a bursty-lonely?
                'rules': [
                    # GZMB is nice, but also strong in some T cells...
                    [('expr', 'CLEC10A', k1k_pb_mc_score_threshs.C_DC_CLEC10A_THRESH, np.inf)],
                    [('expr', 'FCER1A', k1k_pb_mc_score_threshs.C_DC_FCER1A_THRESH, np.inf)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.CLEC10A_NKG7_DOUBLET_NKG7_THRESH)],
                    # [('expr', 'CD79A', pb_mc_score_threshs.CD79A_THRESH, np.inf)],
                    # [('expr', 'JCHAIN', pb_mc_score_threshs.PLASMABLAST_JCHAIN_THRESH, np.inf)],
                ],
            }),
            ('AS-DC', {
                # Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors | Science (https://www.science.org/doi/full/10.1126/science.aah4573, 2017): "A new DC subset, representing 2 to 3% of the DC populations across all 10 donors tested, characterized by the expression of AXL, SIGLEC1, and SIGLEC6 antigens, named AS DCs."

                'rules': [
                    # GZMB is nice, but also strong in some T cells...
                    [('expr', 'AXL', k1k_pb_mc_score_threshs.AS_DC_AXL_THRESH, np.inf)],
                    [('expr', 'SIGLEC6', k1k_pb_mc_score_threshs.AS_DC_SIGLEC6_THRESH, np.inf)],
                    # [('expr', 'CD79A', pb_mc_score_threshs.CD79A_THRESH, np.inf)],
                    # [('expr', 'JCHAIN', pb_mc_score_threshs.PLASMABLAST_JCHAIN_THRESH, np.inf)],
                ],
            }),
            
            ('ncMonocyte', {
                # Yazar et al annotated these as 'CD14-low, CD-16-positive monocyte'
                'rules': [
                    [('expr', 'MS4A7', k1k_pb_mc_score_threshs.NC_MONOCYTE_MS4A7_THRESH, np.inf)],
                    [('expr', 'FCGR3A', k1k_pb_mc_score_threshs.NC_MONOCYTE_FCGR3A_THRESH, np.inf)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.MONOCYTE_NKG7_THRESH)],
                ],
            }),
            ('cMonocyte', {
                'rules': [
                    [('expr', 'CD14', k1k_pb_mc_score_threshs.CD14_THRESH, np.inf)],
                    [('expr', 'NKG7', -np.inf, k1k_pb_mc_score_threshs.CD14_NKG7_DOUBLET_NKG7_THRESH)],
                ],
            }),
            ('ERYP', {
                'rules': [
                    [('expr', 'CD34', k1k_pb_mc_score_threshs.CD34_THRESH, np.inf)],
                    [('expr', 'BLVRB', k1k_pb_mc_score_threshs.EP_BLVRB_THRESH, np.inf)],
                ],
            }),
            ('HSC', {
                'rules': [
                    [('expr', 'AVP', k1k_pb_mc_score_threshs.HSC_AVP_THRESH, np.inf)],
                ],
            }),
            ('CLP-M', {
                'rules': [
                    [('expr', 'HOXA9', k1k_pb_mc_score_threshs.CLP_HOXA9_THRESH, np.inf)],
                ],
            }),
            ('NKTDP', {
                'rules': [
                    [('expr', 'ACY3', k1k_pb_mc_score_threshs.NKTDP_ACY3_THRESH, np.inf)],
                ],
            }),
            ('GMP-L', {
                'rules': [
                    [('expr', 'MPO', k1k_pb_mc_score_threshs.GMP_L_MPO_THRESH, np.inf)],
                ],
            }),
            ('BEMP', {
                'rules': [
                    [('expr', 'HDC', k1k_pb_mc_score_threshs.BEMP_HDC_THRESH, np.inf)],
                ],
            }),
            ('MEBEMP', {
                'rules': [
                    [('expr', 'GATA2', k1k_pb_mc_score_threshs.MEBEMP_GATA2_THRESH, np.inf)],
                    [('expr', 'HDC', -np.inf, k1k_pb_mc_score_threshs.BEMP_HDC_THRESH)],
                    [('expr', 'BLVRB', -np.inf, k1k_pb_mc_score_threshs.EP_BLVRB_THRESH)],
                ],
            }),
            ('MPP', {
                'rules': [
                    [('expr', 'CD34', k1k_pb_mc_score_threshs.CD34_THRESH, np.inf)],
                    [('expr', 'MPO', -np.inf, k1k_pb_mc_score_threshs.GMP_L_MPO_THRESH)],
                    [('expr', 'AVP', -np.inf, k1k_pb_mc_score_threshs.HSC_AVP_THRESH)],
                    [('expr', 'CSF3R', k1k_pb_mc_score_threshs.MPP_CSF3R_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD14', k1k_pb_mc_score_threshs.CD14_NKG7_DOUBLET_CD14_THRESH, np.inf)],
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.CD14_NKG7_DOUBLET_NKG7_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'IRF8', k1k_pb_mc_score_threshs.IRF8_NKG7_DOUBLET_IRF8_THRESH, np.inf)],
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.IRF8_NKG7_DOUBLET_NKG7_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CLEC10A', k1k_pb_mc_score_threshs.CLEC10A_NKG7_DOUBLET_CLEC10A_THRESH, np.inf)],
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.CLEC10A_NKG7_DOUBLET_NKG7_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH, np.inf)],
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.CD79A_NKG7_DOUBLET_NKG7_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.CD79A_NKG7_DOUBLET_CD79A_THRESH, np.inf)],
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.B_NKG7_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.NK_T_CD79A_THRESH, np.inf)],
                    [('expr', 'CD3D', k1k_pb_mc_score_threshs.CD79A_CD3D_DOUBLET_CD3D_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD79A', k1k_pb_mc_score_threshs.CD79A_CD3D_DOUBLET_CD79A_THRESH, np.inf)],
                    [('expr', 'CD3D', k1k_pb_mc_score_threshs.B_CD3D_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'NKG7', k1k_pb_mc_score_threshs.MONOCYTE_NKG7_THRESH, -7.6)],
                    [('expr', 'FCGR3A', k1k_pb_mc_score_threshs.NC_MONOCYTE_FCGR3A_THRESH, np.inf)],
                ],
            }),
            ('HBA_HBB_contaminated', {
                'rules': [
                    [('expr', 'HBB', k1k_pb_mc_score_threshs.HBA_HBB_CONTAMINATED_HBB_THRESH, np.inf)],
                ],
            }),
        ],

        'state_column_name': 'state',
    },
    'palantir_bm': {
        # 240130: ugh. looking at metadata, seems like it is CD34+ CD38- BM cells!
        # IIUC, this is data from Characterization of cell fate probabilities in single-cell data with Palantir (https://www.nature.com/articles/s41587-019-0068-4, 2019). In Furer & Rappoport, referred to as Setty et al. (e.g., in sup fig 2).
        # "Cryopreserved bone marrow stem/progenitor CD34+ cells from healthy donors were purchased from AllCells, LLC. (catalog no. ABM022F) and stored in vapor phase nitrogen until use."

        # created by using 231213_convert_0.8_to_0.9.py on /dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/c.h5ad and /dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/mc.h5ad, which were copied from /home/nimrodra/proj/blood_aging/metacells_workdir/palantir_all_500000_cdata_output.h5ad and /home/nimrodra/proj/blood_aging/metacells_workdir/palantir_all_500000_mdata_output.h5ad.
        # 'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/orenmil_processing/converted_c.h5ad', 
        # 'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/orenmil_processing/converted_mc.h5ad', 
        # 'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/orenmil_processing/converted_c_with_states.h5ad', 
        # 'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/orenmil_processing/converted_mc_with_states.h5ad', 
        'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/orenmil_processing/c_240210.h5ad', 
        'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/palantir_bm_data/orenmil_processing/mc_240210.h5ad', 

        'atlas_cell_type_colors_csv_file_path': DEFAULT_CELL_STATE_COLORS_CSV_FILE_PATH,

        'cell_state_and_info_list': [
            

            ('MKP', {
                'rules': [
                    [('expr', 'ITGA2B', -11.5, np.inf)],
                    [('expr', 'PPBP', -np.inf, -12.5)], # not post-MKP
                    [('expr', 'AHSP', -np.inf, -15)], # chose MKP over EP
                ],
                # 'override_others': True,
            }),
            ('post-MKP', {
                'rules': [
                    [('expr', 'PPBP', -12.5, np.inf)],
                ],
                # 'override_others': True,
            }),
            ('ERYP', {
                'rules': [
                    [('expr', 'PLEK', -np.inf, -15.4)],
                    [('expr', 'BLVRB', -12.3, np.inf)],
                    [('expr', 'AHSP', -np.inf, -12.2)],
                    
                    # [('expr', 'ITGA2B', -14.2, np.inf)], # 240107: makes less sense to decrease PLEK and then increase it later. maybe there are still some transcription factors that increase ITGA2B, and maybe it is even necessary for something??? 231231: seems to me like it makes no sense for EP to increase ITGA2B, as later it would go down (it makes sense to increase it before giving up on becoming MKP (and maybe also BEMPs don't lower it later?). though of course maybe it should be transiently high for EP themselves??
                ],
                # 'override_others': True,
            }),
            ('post-EP', {
                'rules': [
                    [('expr', 'AHSP', -12.2, np.inf)],
                    [('expr', 'ITGA2B', -np.inf, -14.2)],
                ],
                # 'override_others': True,
            }),
            
            ('HSC', {
                'rules': [
                    [('expr', 'AVP', -12, np.inf)],
                    [('expr', 'HLF', -16, np.inf)],
                ],
                # 'override_others': True,
            }),
            
            ('GMP-E', {
                'rules': [
                    [('expr', 'MPO', -14, -10)],
                    [('expr', 'S100A8', -np.inf, -11)], # not monocytes
                    [('expr', 'PPBP', -np.inf, -12.5)], # not post-MKP
                    [('expr', 'JCHAIN', -np.inf, -15.2)], # not DC-prog
                    [('expr', 'LYZ', -np.inf, -11)], # not post-GMP2
                    [('expr', 'EPX', -np.inf, -15)], # not eosinophil-prog
                    [('expr', 'NKTR', -np.inf, -11.7)], # not pro-B-NKTR-outlier?
                    [('expr', 'FAHD2A', -np.inf, -14)], # not B-FAHD2A-outlier
                ],
                # 'override_others': True,
            }),
            # GMP-L and post-GMP1 have large variance in ['ELANE','AZU1','PRTN3']. not sure whether this is lateral.
            ('GMP-L', {
                'rules': [
                    [('expr', 'MPO', -10, np.inf)],
                    [('expr', 'ANXA2', -np.inf, -14)],
                    [('expr', 'NKTR', -np.inf, -11.7)], # not pro-B-NKTR-outlier?
                ],
                # 'override_others': True,
            }),
            ('Monocyte-prog?', {
                'rules': [
                    [('expr', 'MPO', -10, np.inf)],
                    [('expr', 'ANXA2', -14, np.inf)],
                    [('expr', 'S100A8', -np.inf, -11)], # not Monocytes or Monocyte-prog
                    [('expr', 'NKTR', -np.inf, -11.7)], # not pro-B-NKTR-outlier?
                ],
            }),
            ('Monocyte-prog?', {
                'rules': [
                    [('expr', 'ANXA2', -13, np.inf)],
                    
                    [('expr', 'IRF7', -np.inf, -12)], # not DC
                    [('expr', 'JCHAIN', -np.inf, -15.2)], # not DC-prog
                    [('expr', 'S100A8', -np.inf, -11)], # not Monocytes or Monocyte-prog
                    [('expr', 'NKTR', -np.inf, -11.7)], # not pro-B-NKTR-outlier?
                ],
            }),
            ('Monocyte-prog-L?', {
                'rules': [
                    # [('expr', 'CD14', -15, np.inf)], # works for one of the 2 MCs.
                    [('expr', 'S100A8', -11.3, np.inf)],
                ],
            }),
            ('DC', {
                'rules': [
                    [('expr', 'IRF7', -12, np.inf)],
                    # [('expr', 'IRF8', -11, np.inf)], # not needed
                ],
                # 'override_others': True,
            }),
            ('DC-prog?', {
                'rules': [
                    [('expr', 'JCHAIN', -15.2, np.inf)], 
                    # [('expr', 'IRF8', -11.5, np.inf)], 
                    # [('expr', 'SPINK2', -15, np.inf)], 
                    [('expr', 'DNTT', -np.inf, -11)], # not pro-B
                    [('expr', 'IRF7', -np.inf, -12)], # not DC
                    [('expr', ['CD79A', 'CD79B'], -np.inf, -13)], # not B
                    [('expr', 'KIAA0087', -np.inf, -13.5)], # not CLP-M
                    # [('expr', 'LYZ', -np.inf, -12)], # not post-GMP
                ],
                # 'override_others': True,
            }),
            
            ('eosinophil-prog?', {
                'rules': [
                    # NOTE: MPO is a bit high...
                    [('expr', 'EPX', -15, np.inf)],
                    # [('expr', 'CEBPE', -14, np.inf)], # not needed
                    # [('expr', 'PRG2', -14.5, np.inf)], # not needed
                    # [('expr', 'MPO', -14, np.inf)],
                    
                ],
                # 'override_others': True,
            }),
            ('BEMP', {
                'rules': [
                    [('expr', 'LMO4', -11.5, np.inf)],
                ],
                # 'override_others': True,
            }),
            ('MEBEMP-L', {
                'rules': [
                    [('expr', 'GATA1', -15.7, np.inf)],
                    [('expr', 'HOPX', -np.inf, -14.7)],
                    [('expr', 'LMO4', -np.inf, -11.5)], # not BEMP
                    [
                        ('expr', 'BLVRB', -np.inf, -12.3),
                        ('expr', 'PLEK', -15.4, np.inf),
                    ], # not EP
                    [('expr', 'ITGA2B', -np.inf, -11.5)], # not MKP
                    [('expr', 'EPX', -np.inf, -15)], # not eosinophil-prog
                ],
                # 'override_others': True,
            }),
            ('MPP', {
                'rules': [
                    [('expr', 'HOPX', -14.7, np.inf)],
                    [
                        # not HSC
                        ('expr', 'AVP', -np.inf, -12),
                        ('expr', 'HLF', -np.inf, -16),
                    ],
                    [('expr', 'MPO', -np.inf, -14)], # not GMP-E
                    [('expr', 'JCHAIN', -np.inf, -15.2)], # not DC-prog
                    [('expr', 'DNTT', -np.inf, -10.5)], # not pro-B
                    [('expr', 'DNTT', -np.inf, -13.9)], # not MPP-pro-or-pre-B-MC_Doublet?
                ],
                # 'override_others': True,
            }),

            ('VPREB1-B', {
                # NOTE: as far as i can tell, VPREB1 is not expressed in mature B cells.
                'rules': [
                    [('expr', ['CD79A', 'CD79B'], -12, np.inf)],
                    [('expr', 'DNTT', -np.inf, -11)],
                    [('expr', 'MPO', -np.inf, PALANTIR_B_MAX_MPO_EXPR)],
                ],
                # 'override_others': True,
            }),
            
            # Regulation of B-cell proliferation and differentiation by pre-B-cell receptor signalling (https://www.nature.com/articles/nri2491, 2009): "recent data show that the expression of cyclin D3 during B-cell development normally becomes induced in large pre-B cells and is directly downregulated in small pre-B cells, suggesting that cyclin D3 has an important role in the receptor-mediated proliferation of pre-B cells104. Indeed, deletion of CCND3 (which encodes cyclin D3) in mice results in a block of B-cell development at the pro-B-cell to pre-B-cell transition owing to impaired cell cycle progression of the large pre-B-cell population"
            # RNAi Screening Identifies A Novel Role for A-Kinase Anchoring Protein 12 (AKAP12) in B Cell Development and Function (https://ashpublications.org/blood/article/120/21/855/87490/RNAi-Screening-Identifies-A-Novel-Role-for-A, 2012)
            # https://www.immunology.org/public-information/bitesized-immunology/cells/b-cells says: "Signalling through the pre-BCR drives intense proliferation and differentiation into the small pre-B cell stage. Quiescent small pre-B cells then undergo V-J rearrangement of the Ig light chain [...]".
            # Unraveling B cell trajectories at single cell resolution (https://www.cell.com/trends/immunology/fulltext/S1471-4906%2822%2900003-5, 2022): fig 1 is nice.
            # Single-Cell Trajectory Detection Uncovers Progression and Regulatory Coordination in Human B Cell Development (https://www.sciencedirect.com/science/article/pii/S0092867414004711, 2014): fig2 shows that CD34 is pretty silent in pro-B, but then slightly goes up again in pre-B.
            ('small-pre-B?', {
                # highest CD9 here
                'rules': [
                    [('expr', ['CD79A', 'CD79B'], -8.2, np.inf)],
                    [('expr', 'PCLAF', -np.inf, -13.7)],
                    [('expr', 'MPO', -np.inf, PALANTIR_B_MAX_MPO_EXPR)],
                ],
                # 'override_others': True,
            }),
            ('large-pre-B?', {
                'rules': [
                    # [('expr', ['CD79A', 'CD79B'], -11.5, np.inf)], no need
                    [('expr', 'DNTT', -11, np.inf)],
                    [('expr', 'VPREB3', -14, np.inf)],
                    [('expr', 'PCLAF', -13.7, np.inf)],
                    [('expr', 'MPO', -np.inf, PALANTIR_B_MAX_MPO_EXPR)],
                    [('expr', 'AVP', -np.inf, -15.7)],
                ],
                # 'override_others': True,
            }),
            ('pro-B?', {
                'rules': [
                    
                    
                    # NOTE: we have two MCs here with KIAA0087 (and FAM30A and ACY3) a bit higher than others. maybe they recently differentiated from CLP-M to pro-B?? anyway, they don't seem different enough from other pro-B to call them pre-pro-B or split pro-B to early and late according to KIAA0087 etc...
                    [('expr', 'DNTT', -10.5, np.inf)],
                    [('expr', 'VPREB3', -np.inf, -14)],
                    [('expr', 'MPO', -np.inf, PALANTIR_B_MAX_MPO_EXPR)],
                ],
                # 'override_others': True,
            }),
            
            ('CLP', {
                'rules': [
                    [('expr', 'KIAA0087', -13.5, np.inf)],
                ],
                # 'override_others': True,
            }),

            
            
            ('MPP-pro-or-pre-B-MC_Doublet?', {
                # NOTE: the cells might be completely fine, just clustered to a metacell even though they are of different states. DNTT might be the best separator between cells in these MCs (can do diff exp between cells in these MCs with high and low DNTT...)
                'rules': [
                    [('expr', 'HOPX', -14.7, np.inf)],
                    [('expr', 'DNTT', -13.9, -10.5)], # not pro-B (<-10.5)
                    [('expr', 'MPO', -np.inf, -14)], # not GMP-E
                    [('expr', 'JCHAIN', -np.inf, -15.2)], # not DC-prog
                ],
                # 'override_others': True,
            }),
            ('SELL-outlier', {
                # NOTE: why can't it be CLP-E? because of MPO=-10. also high CD38.
                'rules': [
                    [('expr', 'SELL', -10, np.inf)],
                    # also outlier of: 'HSPA5', 'SLC3A2', 'ICAM3', 'ANKRD28', 'APLP2', 'CTSD', 'SPNS3', 'PTPRC'
                ],
                # 'override_others': True,
            }),
            ('pre-B-FAHD2A-outlier', {
                'rules': [
                    # also highest AKNA, MOB4
                    [('expr', 'FAHD2A', -14, np.inf)],
                    # [('expr', ['CD79A', 'CD79B'], -11, np.inf)],
                ],
                # 'override_others': True,
            }),
            ('pro-B-NKTR-outlier?', {
                # also high JADE1 and a bit high MPO.
                'rules': [
                    [('expr', 'DNTT', -10.5, np.inf)],
                    [('expr', 'NKTR', -11.7, np.inf)],
                ],
                # 'override_others': True,
            }),
            ('small-pre-B-dividing_MPP-Doublet?', {
                # also, very high in JADE1.
                'rules': [
                    [('expr', 'VPREB3', -9.5, np.inf)],
                    [('expr', 'AVP', -15.7, np.inf)],
                ],
                # 'override_others': True,
            }),
            ('large-pre-B-MPP-Doublet?', {
                'rules': [
                    [('expr', 'VPREB3', -12, np.inf)],
                    [('expr', 'AVP', -15.1, np.inf)],
                ],
                # 'override_others': True,
            }),
        ],


        'state_column_name': 'state',
    },

    'max_batch_kruskal_pval_to_ignore_gene_in_mc_calculation': 1e-6,

    'nimrod_oren_atlas': {
        # when people ask me to send:
        # cells: /net/dummy/dummy/dummy/export/tgdata/users/orenmil/mds/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_c.h5ad
        # metacells: /net/dummy/dummy/dummy/export/tgdata/users/orenmil/mds/240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_mc.h5ad
        
        'dataset_name': 'illu_mds',
        'mc_model_dir_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas'),
        'mc_ad_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_mc.h5ad'), 
        'c_ad_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240618_pb_illu_mds_cytopenia_normal/intermediate_output/mc_models/final_normal_pb_atlas/240626_cHSPC_79_normal_illu_atlas_c.h5ad'), 
        'atlas_cell_type_colors_csv_file_path': DEFAULT_CELL_STATE_COLORS_CSV_FILE_PATH,
        'state_column_name': 'state',


        'cell_state_and_info_list': [
            ('NKTDP', {
                'rules': [
                    [
                        ('expr', 'ID2', pb_cd34_mc_score_threshs.NKTDP_ID2_THRESH, np.inf),
                        ('expr', 'ACY3', pb_cd34_mc_score_threshs.NKTDP_ACY3_THRESH, np.inf),
                    ], 

                    [('expr', 'LAMP3', -np.inf, pb_cd34_mc_score_threshs.OUTLIER_LAMP3_THRESH)], 
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.NK_NKG7_THRESH)], 
                    [('expr', 'LYZ', -np.inf, pb_cd34_mc_score_threshs.NKTDP_LYZ_THRESH)], 
                ],
            }),
            ('CLP', {
                'rules': [
                    [('expr', 'HOXA9', pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH, np.inf)],
                    
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.NKTDP_ID2_THRESH)],
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.NKTDP_ACY3_THRESH)],
                ],
            }),
            ('CLP_NKTDP_intermediate?', {
                'rules': [
                    [('expr', 'HOXA9', pb_cd34_mc_score_threshs.CLP_NKTDP_INTERMEDIATE_HOXA9_THRESH, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'AVP', -np.inf, pb_cd34_mc_score_threshs.CLP_NKTDP_INTERMEDIATE_AVP_THRESH)], 
                    
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.NKTDP_ID2_THRESH)],
                    [('expr', 'ACY3', pb_cd34_mc_score_threshs.CLP_NKTDP_INTERMEDIATE_ACY3_THRESH, pb_cd34_mc_score_threshs.NKTDP_ACY3_THRESH)],
                ],
            }),

            ('CLP-E', {
                'rules': [
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'LTB', pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH, np.inf)],
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH, pb_cd34_mc_score_threshs.HSC_AVP_THRESH)], 
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.NKTDP_ID2_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GNLY', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_GNLY_THRESH)], 
                    [('expr', 'CD79A', -np.inf, pb_cd34_mc_score_threshs.B_CD79A_THRESH)], 
                    [('expr', 'APOC1', -np.inf, pb_cd34_mc_score_threshs.B_NKTDP_APOC1_THRESH)], 
                    # [('expr', 'MS4A1', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_MS4A1_THRESH)], 
                ],
            }),
            ('pre-B?', {
                # NOTE: written before 240208: we have MCs here with KIAA0087 (and FAM30A and ACY3) higher than we might expect. yet i guess this is one of the problems of negative (rather than positive) markers - they are much more dependent on mRNA decay and cell division rates. Maybe the ones with higher than expected CLP-M genes are those that performed relatively few cell divisions, if any, since differentiating from CLP-M???
                'rules': [
                    [('expr', 'DNTT', pb_cd34_mc_score_threshs.PRO_B_DNTT_THRESH, np.inf)],
                    [('expr', 'VPREB3', pb_cd34_mc_score_threshs.PRE_B_VPREB3_THRESH, np.inf)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.PRO_B_IL7R_THRESH)], 
                ],
            }),
            ('pro-B?', {
                # NOTE: 240115: this seems later than what we annotate as "pro-B?" in palantir BM. Specifically see IL7R.
                'rules': [
                    [('expr', 'DNTT', pb_cd34_mc_score_threshs.PRO_B_DNTT_THRESH, np.inf)],
                    [('expr', 'VPREB3', pb_cd34_mc_score_threshs.CLP_VPREB3_THRESH, pb_cd34_mc_score_threshs.PRE_B_VPREB3_THRESH)], 
                    [('expr', 'IL7R', pb_cd34_mc_score_threshs.PRO_B_IL7R_THRESH, np.inf)],
                ],
            }),
            ('pDC', {
                'rules': [
                    [('expr', 'IRF8', pb_cd34_mc_score_threshs.P_DC_IRF8_THRESH, np.inf)], 
                    [('expr', 'AXL', -np.inf, pb_cd34_mc_score_threshs.AS_DC_AXL_THRESH)], 
                ],
            }),
            ('cDC?', {
                'rules': [
                    [('expr', 'CLEC10A', k1k_pb_mc_score_threshs.C_DC_CLEC10A_THRESH, np.inf)],
                    [('expr', 'AXL', -np.inf, pb_cd34_mc_score_threshs.AS_DC_AXL_THRESH)], 
                ],
            }),
            ('AS-DC', {
                'rules': [
                    [('expr', 'AXL', pb_cd34_mc_score_threshs.AS_DC_AXL_THRESH, np.inf)], 
                    [('expr', 'SIGLEC6', pb_cd34_mc_score_threshs.AS_DC_SIGLEC6_THRESH, np.inf)], 
                ],
            }),
            ('NK', {
                'rules': [
                    [('expr', 'NKG7', pb_cd34_mc_score_threshs.NK_NKG7_THRESH, np.inf)], 
                    [('expr', 'CD3D', -np.inf, pb_cd34_mc_score_threshs.T_CD3D_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.NKT_HOXA9_THRESH)], 
                ],
            }),
            ('T', {
                # 240623: there is a T metacell with extremely high CCL4 and also CD34=-14.7, so maybe a doublet? anyway, not important so i leave it as it is.
                'rules': [
                    [('expr', 'CD3D', pb_cd34_mc_score_threshs.T_CD3D_THRESH, np.inf)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.NKT_HOXA9_THRESH)], 
                ],
            }),
            ('HSC', {
                'rules': [
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.HSC_AVP_THRESH, np.inf)], 
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH)],
                ],
            }),
            ('MPP', {
                'rules': [
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.MPP_AVP_THRESH, pb_cd34_mc_score_threshs.HSC_AVP_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH)],
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_NKG7_THRESH)],
                    [('expr', 'IGHA1', -np.inf, pb_cd34_mc_score_threshs.MYELOID_HSPC_IGHA1_THRESH)], 
                ],
            }),
            ('MEBEMP-E', {
                'rules': [
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH, pb_cd34_mc_score_threshs.MPP_AVP_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'BLVRB', -np.inf, pb_cd34_mc_score_threshs.EP_BLVRB_THRESH)], 
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH)],
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_NKG7_THRESH)],
                    [('expr', 'LMO4', -np.inf, pb_cd34_mc_score_threshs.NIMROD_OREN_BEMP_LMO4_THRESH)], 
                    [('expr', 'MS4A1', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_MS4A1_THRESH)], 
                    [('expr', 'IGHA1', -np.inf, pb_cd34_mc_score_threshs.MYELOID_HSPC_IGHA1_THRESH)], 
                    [('expr', 'SPARCL1', -np.inf, pb_cd34_mc_score_threshs.HSPC_SPARCL1_THRESH)], 
                ],
            }),
            ('MEBEMP-L', {
                'rules': [
                    [('expr', 'SLC40A1', pb_cd34_mc_score_threshs.MEBEMP_L_SLC40A1_THRESH, np.inf)], 
                    [('expr', 'AVP', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'BLVRB', -np.inf, pb_cd34_mc_score_threshs.EP_BLVRB_THRESH)], 
                    [('expr', 'LMO4', -np.inf, pb_cd34_mc_score_threshs.BEMP_LMO4_THRESH)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_ACY3_THRESH)], 
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_NKG7_THRESH)],
                    [('expr', 'MS4A1', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_MS4A1_THRESH)], 
                    [('expr', 'IGHA1', -np.inf, pb_cd34_mc_score_threshs.MYELOID_HSPC_IGHA1_THRESH)], 
                    [('expr', 'SPARCL1', -np.inf, pb_cd34_mc_score_threshs.HSPC_SPARCL1_THRESH)], 
                    # [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_LTB_THRESH)],
                ],
            }),
            ('ERYP', {
                'rules': [
                    [('expr', 'BLVRB', pb_cd34_mc_score_threshs.EP_BLVRB_THRESH, np.inf)], 
                    [('expr', 'LMO4', -np.inf, pb_cd34_mc_score_threshs.BEMP_LMO4_THRESH)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.CD14_THRESH)],
                    # [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_LTB_THRESH)],
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_ACY3_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.HSPC_CD14_THRESH)],
                    [('expr', 'MS4A1', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_MS4A1_THRESH)], 
                    [('expr', 'IGHA1', -np.inf, pb_cd34_mc_score_threshs.MYELOID_HSPC_IGHA1_THRESH)], 
                ],
            }),
            ('BEMP', {
                'rules': [
                    # [('expr', 'AVP', pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH, pb_cd34_mc_score_threshs.MPP_AVP_THRESH)], 
                    [('expr', 'LMO4', pb_cd34_mc_score_threshs.NIMROD_OREN_BEMP_LMO4_THRESH, np.inf)], 
                    [('expr', 'HDC', pb_cd34_mc_score_threshs.BEMP_HDC_THRESH, np.inf)], 
                    [('expr', 'KIAA0087', -np.inf, pb_cd34_mc_score_threshs.BEMP_KIAA0087_THRESH)], 
                    # [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_LTB_THRESH)],
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_ACY3_THRESH)], 
                    [('expr', 'HBB', -np.inf, pb_cd34_mc_score_threshs.BEMP_HBB_THRESH)], 
                    [('expr', 'MS4A1', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_MS4A1_THRESH)], 
                    [('expr', 'IGHA1', -np.inf, pb_cd34_mc_score_threshs.MYELOID_HSPC_IGHA1_THRESH)], 
                    [('expr', 'FAM30A', -np.inf, pb_cd34_mc_score_threshs.BEMP_FAM30A_THRESH)], 
                    [('expr', 'MAF', -np.inf, pb_cd34_mc_score_threshs.BEMP_GMP_L_MAF_THRESH)], 
                ],
            }),
            ('Monocyte-prog?', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_L_MPO_THRESH, np.inf)], 
                    [('expr', 'ANXA2', pb_cd34_mc_score_threshs.MONOCYTE_PROG_ANXA2_THRESH, np.inf)], 
                ],
            }),

            ('GMP-L', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_L_MPO_THRESH, np.inf)], 
                    [('expr', 'ANXA2', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_PROG_ANXA2_THRESH)], 
                    [('expr', 'HDC', -np.inf, pb_cd34_mc_score_threshs.GMP_HDC_THRESH)], 
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.GMP_ACY3_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.HSPC_CD14_THRESH)], 
                ],
            }),
            ('GMP-E', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH, pb_cd34_mc_score_threshs.GMP_L_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'CD79A', -np.inf, pb_cd34_mc_score_threshs.B_CD79A_THRESH)], 
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.GMP_E_DNTT_THRESH)], 
                    [('expr', 'HDC', -np.inf, pb_cd34_mc_score_threshs.GMP_HDC_THRESH)], 
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.GMP_ACY3_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.HSPC_CD14_THRESH)], 
                ],
            }),
            ('CCDC88A-high-B', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'CCDC88A', pb_cd34_mc_score_threshs.CCDC88A_THRESH, np.inf)], 
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'SLC40A1', -np.inf, pb_cd34_mc_score_threshs.B_SLC40A1_THRESH)], 
                ],
            }),
            ('plasmablast', {
                'rules': [
                    [('expr', 'MZB1', pb_cd34_mc_score_threshs.PLASMABLAST_MAZB1_THRESH, np.inf)],
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'SLC40A1', -np.inf, pb_cd34_mc_score_threshs.B_SLC40A1_THRESH)], 
                ],
            }),
            ('naive-B', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'TCL1A', pb_cd34_mc_score_threshs.NAIVE_B_TCL1A_THRESH, np.inf)],
                    [('expr', 'CCDC88A', -np.inf, pb_cd34_mc_score_threshs.CCDC88A_THRESH)], 
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'SLC40A1', -np.inf, pb_cd34_mc_score_threshs.B_SLC40A1_THRESH)], 
                ],
            }),
            ('memory-B', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'TCL1A', -np.inf, pb_cd34_mc_score_threshs.NAIVE_B_TCL1A_THRESH)],
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'CCDC88A', -np.inf, pb_cd34_mc_score_threshs.CCDC88A_THRESH)], 
                    [
                        ('expr', 'ZEB2', -np.inf, pb_cd34_mc_score_threshs.AGED_B_ZEB2_THRESH),
                        ('expr', 'VPREB3', pb_cd34_mc_score_threshs.AGED_B_VPREB3_THRESH, np.inf),
                    ],
                    [('expr', 'MZB1', -np.inf, pb_cd34_mc_score_threshs.PLASMABLAST_MAZB1_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'LAMP3', -np.inf, pb_cd34_mc_score_threshs.OUTLIER_LAMP3_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'SLC40A1', -np.inf, pb_cd34_mc_score_threshs.B_SLC40A1_THRESH)], 
                ],
            }),
            ('aged-B?', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'CCDC88A', -np.inf, pb_cd34_mc_score_threshs.CCDC88A_THRESH)], 
                    [('expr', 'ZEB2', pb_cd34_mc_score_threshs.AGED_B_ZEB2_THRESH, np.inf)],
                    [('expr', 'VPREB3', -np.inf, pb_cd34_mc_score_threshs.AGED_B_VPREB3_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'SLC40A1', -np.inf, pb_cd34_mc_score_threshs.B_SLC40A1_THRESH)], 
                ],
            }),
            ('ncMonocyte', {
                'rules': [
                    [('expr', 'MS4A7', pb_cd34_mc_score_threshs.NC_MONOCYTE_MS4A7_THRESH, np.inf)],
                    [('expr', 'FCGR3A', pb_cd34_mc_score_threshs.NC_MONOCYTE_FCGR3A_THRESH, np.inf)],
                    [('expr', 'DTL', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_DTL_THRESH)], 
                ],
            }),
            ('cMonocyte', {
                'rules': [
                    [('expr', 'CD14', pb_cd34_mc_score_threshs.CD14_THRESH, np.inf)],
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_PRSS57_THRESH)], 
                    [('expr', 'DTL', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_DTL_THRESH)], 
                    [('expr', 'FCGR3A', -np.inf, pb_cd34_mc_score_threshs.C_MONOCYTE_FCGR3A_THRESH)],
                ],
            }),
            ('intermMonocyte', {
                'rules': [
                    [('expr', 'CD14', pb_cd34_mc_score_threshs.CD14_THRESH, np.inf)],
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_PRSS57_THRESH)], 
                    [('expr', 'DTL', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_DTL_THRESH)], 
                    [('expr', 'FCGR3A', pb_cd34_mc_score_threshs.C_MONOCYTE_FCGR3A_THRESH, np.inf)],
                ],
            }),
            ('Endothel', {
                'rules': [
                    [('expr', 'VWF', pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH, np.inf)], 
                ],
            }),

            ('mregDC', {
                'rules': [
                    [('expr', 'LAMP3', pb_cd34_mc_score_threshs.M_REG_DC_LAMP3_THRESH, np.inf)], 
                ],
            }),
            
            ('Doublet', {
                # mregDC-NKTDP doublets?
                'rules': [
                    [('expr', 'LAMP3', pb_cd34_mc_score_threshs.NOT_M_REG_DC_LAMP3_THRESH, pb_cd34_mc_score_threshs.M_REG_DC_LAMP3_THRESH)], 
                    [('expr', 'ACY3', pb_cd34_mc_score_threshs.M_REG_DC_ACY3_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # Monocyte/DC-NKTDP doublets?
                'rules': [
                    [('expr', 'LYZ', pb_cd34_mc_score_threshs.NKTDP_LYZ_THRESH, np.inf)], 
                    [('expr', 'ACY3', pb_cd34_mc_score_threshs.MONO_DC_ACY3_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # NKT-GMP-L doublets?
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.NKT_MPO_THRESH, np.inf)], 
                    [('expr', 'IL32', pb_cd34_mc_score_threshs.GMP_L_IL32_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # ERYP-BEMP doublets?
                'rules': [
                    [('expr', 'HBD', pb_cd34_mc_score_threshs.BEMP_HBD_THRESH, np.inf)], 
                    [('expr', 'LMO4', pb_cd34_mc_score_threshs.ERYP_LMO4_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # MPP/CLP-BEMP doublets?
                'rules': [
                    [('expr', 'FAM30A', pb_cd34_mc_score_threshs.BEMP_FAM30A_THRESH, np.inf)], 
                    [('expr', 'MS4A2', pb_cd34_mc_score_threshs.MPP_CLP_MS4A2_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # NKTDP-GMP-L-BEMP multiplets?
                'rules': [
                    [('expr', 'MAF', pb_cd34_mc_score_threshs.BEMP_GMP_L_MAF_THRESH, np.inf)], 
                    [('expr', 'CLEC11A', pb_cd34_mc_score_threshs.NKTDP_CLEC11A_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # MEBEMP-B doublets?
                'rules': [
                    [('expr', 'MS4A1', pb_cd34_mc_score_threshs.MEBEMP_MS4A1_THRESH, np.inf)], 
                    [
                        ('expr', 'APOC1', pb_cd34_mc_score_threshs.B_NKTDP_APOC1_THRESH, np.inf),
                        ('expr', 'SLC40A1', pb_cd34_mc_score_threshs.MEBEMP_L_SLC40A1_THRESH, np.inf),
                    ], 
                ],
            }),
            ('Doublet', {
                # MEBEMP-B doublets?
                'rules': [
                    [('expr', 'MS4A1', pb_cd34_mc_score_threshs.NKTDP_MS4A1_THRESH, np.inf)], 
                    [('expr', 'SLC40A1', pb_cd34_mc_score_threshs.B_SLC40A1_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # MEBEMP-CLP doublets?
                'rules': [
                    [('expr', 'HOXA9', pb_cd34_mc_score_threshs.MEBEMP_HOXA9_THRESH, np.inf)], 
                    [('expr', 'APOC1', pb_cd34_mc_score_threshs.B_NKTDP_APOC1_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # MEBEMP-NKT doublets?
                'rules': [
                    [('expr', 'GNLY', pb_cd34_mc_score_threshs.MEBEMP_GNLY_THRESH, np.inf)], 
                    [('expr', 'SLC40A1', pb_cd34_mc_score_threshs.MEBEMP_L_SLC40A1_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # B-NKT doublets?
                'rules': [
                    [('expr', 'GNLY', pb_cd34_mc_score_threshs.B_GNLY_THRESH, np.inf)], 
                    [('expr', 'BASP1', pb_cd34_mc_score_threshs.NKT_BASP1_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # NKTDP-MPP-GMP-E-CLP multiplets?
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.NKTDP_MPO_THRESH, np.inf)], 
                    [('expr', 'SAMHD1', pb_cd34_mc_score_threshs.GMP_SAMHD1_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # NKTDP-B doublets?
                'rules': [
                    [('expr', 'ID2', pb_cd34_mc_score_threshs.B_ID2_THRESH, np.inf)], 
                    [('expr', 'IGHA1', pb_cd34_mc_score_threshs.CLP_NKTDP_IGHA1_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # MPP-B doublets?
                'rules': [
                    [('expr', 'MS4A1', pb_cd34_mc_score_threshs.NKTDP_MS4A1_THRESH, np.inf)], 
                    [('expr', 'PRSS57', pb_cd34_mc_score_threshs.B_PRSS57_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # DC-CLP/NKTDP doublets?
                'rules': [
                    [('expr', 'IRF8', pb_cd34_mc_score_threshs.NKTDP_IRF8_THRESH, np.inf)], 
                    [('expr', 'SPINK2', pb_cd34_mc_score_threshs.DC_SPINK2_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                # Endothel-MPP/MEBEMP doublets?
                'rules': [
                    [('expr', 'SPARCL1', pb_cd34_mc_score_threshs.HSPC_SPARCL1_THRESH, np.inf)], 
                    [('expr', 'PRSS57', pb_cd34_mc_score_threshs.ENDOTHEL_PRSS57_THRESH, np.inf)], 
                ],
            }),
        ],
    },
    'nimrod_148_atlas': {
        'mc_ad_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'from_nimrod/230918_frozen_model/orenmil_processing/240117_mc_with_fixed_type.h5ad'),
        'c_ad_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'from_nimrod/230918_frozen_model/all_cells_fil3_cdata_output.h5ad'),
        'atlas_cell_type_colors_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240620_cell_type_colors_for_furer_rapp_publication.csv'),
        'state_column_name': 'type',
        'orig_state_column_name': 'orig_type',
        'cell_state_and_info_list': [
            ('pro-B?', {
                'rules': [
                    [('expr', 'DNTT', pb_cd34_mc_score_threshs.NIMROD_148_PRO_B_DNTT_THRESH, np.inf)], 
                    [('expr', 'VPREB3', pb_cd34_mc_score_threshs.NIMROD_148_PRO_B_VPREB3_THRESH, np.inf)], 
                    [('expr', 'IL7R', pb_cd34_mc_score_threshs.NIMROD_148_PRO_B_IL7R_THRESH, np.inf)], 
                ],
            }),
            ('pre-B?', {
                'rules': [
                    [('expr', 'DNTT', pb_cd34_mc_score_threshs.NIMROD_148_PRO_B_DNTT_THRESH, np.inf)], 
                    [('expr', 'VPREB3', pb_cd34_mc_score_threshs.NIMROD_148_PRE_B_VPREB3_THRESH, np.inf)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.NIMROD_148_PRO_B_IL7R_THRESH)], 
                ],
            }),
            ('mregDC', {
                'rules': [
                    [('expr', 'LAMP3', pb_cd34_mc_score_threshs.NIMROD_148_LAMP3_MREG_DC_THRESH, np.inf)], 
                ],
            }),
        ],
    },
    'nimrod_atlas': {
        # 'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/231004_all_cells_fil2_500000_cdata_output_with_cell_type_exp_indiv.h5ad',
        # 'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/231004_all_cells_fil2_500000_mdata_output_with_type.h5ad',
        
        # 'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/c_240102.h5ad', # NOTE: created by mds/mds_one_offs/231231_fix_220411_nimrod_model.ipynb
        # 'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/mc_240102.h5ad', # NOTE: created by mds/mds_one_offs/231231_fix_220411_nimrod_model.ipynb
        
        # 'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/c_no_dominating_sick_and_clp_e_240115.h5ad', # NOTE: created by mds/mds_one_offs/231231_fix_220411_nimrod_model.ipynb
        # 'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/mc_no_dominating_sick_and_clp_e_240115.h5ad', # NOTE: created by mds/mds_one_offs/231231_fix_220411_nimrod_model.ipynb
        
        'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/c_no_sick_merged_myel_p_with_lymph_p_etc_240214.h5ad', # NOTE: created by mds/mds_one_offs/231231_fix_220411_nimrod_model.ipynb
        'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/220411_frozen_model/orenmil_processing/mc_no_sick_merged_myel_p_with_lymph_p_etc_240214.h5ad', # NOTE: created by mds/mds_one_offs/231231_fix_220411_nimrod_model.ipynb
        'atlas_cell_type_colors_csv_file_path': DEFAULT_CELL_STATE_COLORS_CSV_FILE_PATH,
        
        # 'state_column_name': 'type', # for first run of 231231_fix_220411_nimrod_model.ipynb
        'state_column_name': 'state',
        'orig_state_column_name': 'orig_type',


        
        # copied from /home/nimrodra/proj/blood_aging/metacells_workdir/all_cells2_500000_cdata_output.h5ad
        'atlas_cells_before_platelet_filtering_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/all_cells2_500000_cdata_output.h5ad',


        'cell_state_and_info_list': [
            # NOTE: 240626: watch out. some values in pb_cd34_mc_score_threshs were changed. i guess i should have a different such file for each atlas, but nevermind now. just take values from an earlier version of pb_cd34_mc_score_threshs to make the annotations here work (e.g., from 240620) 

            # TODO: M2.88 (endothel) probably contains two different cell states (can't remember why exactly i concluded that. maybe hist of endothel markers?) - so could clean or split to two MCs?
            # diff_exp:
            # n_c_ad.obs['metacell_name'].isin(['M2.88']) & (n_c_ad.obs['endothel_somewhat_specific'] > -6.8),
            # n_c_ad.obs['metacell_name'].isin(['M2.88']) & (n_c_ad.obs['endothel_somewhat_specific'] < -6.8),
            
            # TODO: M1401.79 (DC) - hist of higher_in_b_than_all_other_hspcs. in other DCs, we don't see so many high B-sig cells...
            # diff exp:
            # n_c_ad.obs['metacell_name'].isin(['M1401.79']) & (n_c_ad.obs['higher_in_b_than_all_other_hspcs'] > -7.7),
            # n_c_ad.obs['metacell_name'].isin(['M1401.79']) & (n_c_ad.obs['higher_in_b_than_all_other_hspcs'] < -7.7),
            #
            # diff exp (specifically see AXL and SIGLEC6 - so it seems like the non-B in M1401.79 are AS-DCs?? though there is some diff-exp. not very strong though. still TCF4, IRF8 lower in "fixed" M1401.79...):
            # n_c_ad.obs['metacell_name'].isin(['M1401.79']) & (n_c_ad.obs['higher_in_b_than_all_other_hspcs'] < -7.7),
            # n_c_ad.obs['metacell_name'].isin(['M1207.86']),
            
            # TODO: M1199.40 (T) - hist of ['NKG7','GNLY','CCL5']. 
            # diff exp:
            # n_c_ad.obs['metacell_name'].isin(['M1199.40']) & (mc_utils.get_genes_expr(n_c_ad, ['NKG7','GNLY','CCL5']) < -8.6),
            # n_c_ad.obs['metacell_name'].isin(['M1199.40']) & (mc_utils.get_genes_expr(n_c_ad, ['NKG7','GNLY','CCL5']) > -8.6),
            # similarly for M1198.53 (T) (though maybe another threshold than -8.6 is better)
            # n_c_ad.obs['metacell_name'].isin(['M1198.53']) & (mc_utils.get_genes_expr(n_c_ad, ['NKG7','GNLY','CCL5']) < -8.6),
            # n_c_ad.obs['metacell_name'].isin(['M1198.53']) & (mc_utils.get_genes_expr(n_c_ad, ['NKG7','GNLY','CCL5']) > -8.6),

            ('NKTDP', {
                'rules': [
                    # why SAMHD1 and not ACY3 as the main marker? look at HOXA9/KIAA0087 vs SAMHD1/NFKB1
                    # 240212: maybe the transition from CLP to NKTDP is best shown in a gene gene plot of RORC vs IRF8, colored by HOXA9, showing only CLP and NKTDP?
                    [('expr', 'SAMHD1', pb_cd34_mc_score_threshs.NKTDP_SAMHD1_THRESH, np.inf)], 
                    [('expr', 'ACY3', pb_cd34_mc_score_threshs.NKTDP_ACY3_THRESH, np.inf)],
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'LAMP3', -np.inf, pb_cd34_mc_score_threshs.OUTLIER_LAMP3_THRESH)], 
                    [('expr', 'HBD', -np.inf, pb_cd34_mc_score_threshs.CLP_NKTDP_HBD_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'CD79A', -np.inf, pb_cd34_mc_score_threshs.NKTDP_CD79A_THRESH)], 
                    [('expr', 'KLRB1', -np.inf, pb_cd34_mc_score_threshs.NKT_PROG_KLRB1_THRESH)], 
                    [('expr', 'TOX2', -np.inf, pb_cd34_mc_score_threshs.NKT_PROG_TOX2_THRESH)], 
                ],
            }),
            ('CLP', {
                'rules': [
                    [('expr', 'HOXA9', pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH, np.inf)],
                    [('expr', 'AVP', -np.inf, pb_cd34_mc_score_threshs.CLP_AVP_THRESH)], 
                    [('expr', 'VPREB3', -np.inf, pb_cd34_mc_score_threshs.CLP_VPREB3_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    # [('expr', 'SAMHD1', -np.inf, pb_cd34_mc_score_threshs.NKTDP_SAMHD1_THRESH)],
                    [('expr', 'CD79A', -np.inf, pb_cd34_mc_score_threshs.B_CD79A_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'HBD', -np.inf, pb_cd34_mc_score_threshs.CLP_NKTDP_HBD_THRESH)], 
                    [('expr', 'CD3D', -np.inf, pb_cd34_mc_score_threshs.CLP_CD3D_THRESH)], 
                    [('expr', 'IGHA1', -np.inf, pb_cd34_mc_score_threshs.CLP_IGHA1_THRESH)], 
                    [('expr', 'KLRB1', -np.inf, pb_cd34_mc_score_threshs.NKT_PROG_KLRB1_THRESH)], 
                    [('expr', 'TOX2', -np.inf, pb_cd34_mc_score_threshs.NKT_PROG_TOX2_THRESH)], 
                ],
            }),
            ('pre-B?', {
                # NOTE: written before 240208: we have MCs here with KIAA0087 (and FAM30A and ACY3) higher than we might expect. yet i guess this is one of the problems of negative (rather than positive) markers - they are much more dependent on mRNA decay and cell division rates. Maybe the ones with higher than expected CLP-M genes are those that performed relatively few cell divisions, if any, since differentiating from CLP-M???
                'rules': [
                    [('expr', 'DNTT', pb_cd34_mc_score_threshs.PRO_B_DNTT_THRESH, np.inf)],
                    [('expr', 'VPREB3', pb_cd34_mc_score_threshs.PRE_B_VPREB3_THRESH, np.inf)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.PRO_B_IL7R_THRESH)], 
                ],
            }),
            ('pro-B?', {
                # NOTE: 240115: this seems later than what we annotate as "pro-B?" in palantir BM. Specifically see IL7R.
                'rules': [
                    [('expr', 'DNTT', pb_cd34_mc_score_threshs.PRO_B_DNTT_THRESH, np.inf)],
                    [('expr', 'VPREB3', pb_cd34_mc_score_threshs.CLP_VPREB3_THRESH, pb_cd34_mc_score_threshs.PRE_B_VPREB3_THRESH)], 
                    [('expr', 'IL7R', pb_cd34_mc_score_threshs.PRO_B_IL7R_THRESH, np.inf)],
                ],
            }),
            ('pDC', {
                'rules': [
                    [('expr', 'IRF8', pb_cd34_mc_score_threshs.P_DC_IRF8_THRESH, np.inf)], 
                    [('expr', 'AXL', -np.inf, pb_cd34_mc_score_threshs.AS_DC_AXL_THRESH)], 
                ],
            }),
            ('AS-DC', {
                'rules': [
                    [('expr', 'AXL', pb_cd34_mc_score_threshs.AS_DC_AXL_THRESH, np.inf)], 
                    [('expr', 'SIGLEC6', pb_cd34_mc_score_threshs.AS_DC_SIGLEC6_THRESH, np.inf)], 
                ],
            }),
            ('NK', {
                'rules': [
                    [('expr', 'NKG7', pb_cd34_mc_score_threshs.NK_NKG7_THRESH, np.inf)], 
                    [('expr', 'CD3D', -np.inf, pb_cd34_mc_score_threshs.T_CD3D_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.NKT_HOXA9_THRESH)], 
                ],
            }),
            ('T', {
                'rules': [
                    [('expr', 'CD3D', pb_cd34_mc_score_threshs.T_CD3D_THRESH, np.inf)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.NKT_HOXA9_THRESH)], 
                ],
            }),
            ('NKT-prog?', {
                'rules': [
                    [('expr', 'SOX4', pb_cd34_mc_score_threshs.NKT_PROG_SOX4_THRESH, np.inf)], 
                    [
                        ('expr', 'KLRB1', pb_cd34_mc_score_threshs.NKT_PROG_KLRB1_THRESH, np.inf),
                        ('expr', 'TOX2', pb_cd34_mc_score_threshs.NKT_PROG_TOX2_THRESH, np.inf),
                    ], 
                    [('expr', 'CD3D', -np.inf, pb_cd34_mc_score_threshs.T_CD3D_THRESH)], 
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.NK_NKG7_THRESH)], 
                ],
            }),
            ('HSC', {
                'rules': [
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.HSC_AVP_THRESH, np.inf)], 
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH)],
                ],
            }),
            ('MPP', {
                'rules': [
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.MPP_AVP_THRESH, pb_cd34_mc_score_threshs.HSC_AVP_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH)],
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_NKG7_THRESH)],
                ],
            }),
            ('MEBEMP-E', {
                'rules': [
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH, pb_cd34_mc_score_threshs.MPP_AVP_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'BLVRB', -np.inf, pb_cd34_mc_score_threshs.EP_BLVRB_THRESH)], 
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH)],
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_NKG7_THRESH)],
                    # [('expr', 'LMO4', -np.inf, pb_cd34_mc_score_threshs.BEMP_LMO4_THRESH)], 
                ],
            }),
            ('MEBEMP-L', {
                'rules': [
                    [('expr', 'SLC40A1', pb_cd34_mc_score_threshs.MEBEMP_L_SLC40A1_THRESH, np.inf)], 
                    [('expr', 'AVP', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'BLVRB', -np.inf, pb_cd34_mc_score_threshs.EP_BLVRB_THRESH)], 
                    [('expr', 'LMO4', -np.inf, pb_cd34_mc_score_threshs.BEMP_LMO4_THRESH)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_LTB_THRESH)],
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_ACY3_THRESH)], 
                    [('expr', 'NKG7', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_NKG7_THRESH)],
                ],
            }),
            ('ERYP', {
                'rules': [
                    [('expr', 'BLVRB', pb_cd34_mc_score_threshs.EP_BLVRB_THRESH, np.inf)], 
                    [('expr', 'LMO4', -np.inf, pb_cd34_mc_score_threshs.BEMP_LMO4_THRESH)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.CD14_THRESH)],
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_LTB_THRESH)],
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_ACY3_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.HSPC_CD14_THRESH)],
                ],
            }),
            ('BEMP', {
                'rules': [
                    [('expr', 'LMO4', pb_cd34_mc_score_threshs.BEMP_LMO4_THRESH, np.inf)], 
                    [('expr', 'KIAA0087', -np.inf, pb_cd34_mc_score_threshs.BEMP_KIAA0087_THRESH)], 
                    [('expr', 'LTB', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_LTB_THRESH)],
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_ACY3_THRESH)], 
                    [('expr', 'HBB', -np.inf, pb_cd34_mc_score_threshs.BEMP_HBB_THRESH)], 
                ],
            }),
            ('Monocyte-prog?', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_L_MPO_THRESH, np.inf)], 
                    [('expr', 'ANXA2', pb_cd34_mc_score_threshs.MONOCYTE_PROG_ANXA2_THRESH, np.inf)], 
                ],
            }),
            ('GMP-L', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_L_MPO_THRESH, np.inf)], 
                    [('expr', 'ANXA2', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_PROG_ANXA2_THRESH)], 
                    [('expr', 'HDC', -np.inf, pb_cd34_mc_score_threshs.GMP_HDC_THRESH)], 
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.GMP_ACY3_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.HSPC_CD14_THRESH)], 
                ],
            }),
            ('GMP-E', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH, pb_cd34_mc_score_threshs.GMP_L_MPO_THRESH)], 
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'CD79A', -np.inf, pb_cd34_mc_score_threshs.B_CD79A_THRESH)], 
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.GMP_E_DNTT_THRESH)], 
                    [('expr', 'HDC', -np.inf, pb_cd34_mc_score_threshs.GMP_HDC_THRESH)], 
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.GMP_ACY3_THRESH)], 
                    [('expr', 'CD14', -np.inf, pb_cd34_mc_score_threshs.HSPC_CD14_THRESH)], 
                ],
            }),
            ('CCDC88A-high-B', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'CCDC88A', pb_cd34_mc_score_threshs.CCDC88A_THRESH, np.inf)], 

                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.B_PRSS57_THRESH)], 
                ],
            }),
            ('plasmablast', {
                'rules': [
                    [('expr', 'MZB1', pb_cd34_mc_score_threshs.PLASMABLAST_MAZB1_THRESH, np.inf)],
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.B_PRSS57_THRESH)], 
                ],
            }),
            ('naive-B', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'TCL1A', pb_cd34_mc_score_threshs.NAIVE_B_TCL1A_THRESH, np.inf)],
                    [('expr', 'CCDC88A', -np.inf, pb_cd34_mc_score_threshs.CCDC88A_THRESH)], 
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.B_PRSS57_THRESH)], 
                ],
            }),
            ('memory-B', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'TCL1A', -np.inf, pb_cd34_mc_score_threshs.NAIVE_B_TCL1A_THRESH)],
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'CCDC88A', -np.inf, pb_cd34_mc_score_threshs.CCDC88A_THRESH)], 
                    [
                        ('expr', 'ZEB2', -np.inf, pb_cd34_mc_score_threshs.AGED_B_ZEB2_THRESH),
                        ('expr', 'VPREB3', pb_cd34_mc_score_threshs.AGED_B_VPREB3_THRESH, np.inf),
                    ],
                    [('expr', 'MZB1', -np.inf, pb_cd34_mc_score_threshs.PLASMABLAST_MAZB1_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'LAMP3', -np.inf, pb_cd34_mc_score_threshs.OUTLIER_LAMP3_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'VPREB1', -np.inf, pb_cd34_mc_score_threshs.B_VPREB1_THRESH)], 
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.B_PRSS57_THRESH)], 
                ],
            }),
            ('aged-B?', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.B_DNTT_THRESH)],
                    [('expr', 'CCDC88A', -np.inf, pb_cd34_mc_score_threshs.CCDC88A_THRESH)], 
                    
                    [('expr', 'ZEB2', pb_cd34_mc_score_threshs.AGED_B_ZEB2_THRESH, np.inf)],
                    [('expr', 'VPREB3', -np.inf, pb_cd34_mc_score_threshs.AGED_B_VPREB3_THRESH)],
                    [('expr', 'ID2', -np.inf, pb_cd34_mc_score_threshs.B_ID2_THRESH)],
                    [('expr', 'MTSS1', -np.inf, pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH)], 
                    [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.B_IL7R_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'GATA1', -np.inf, pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH)], 
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.B_PRSS57_THRESH)], 
                ],
            }),
            ('ncMonocyte', {
                'rules': [
                    [('expr', 'MS4A7', pb_cd34_mc_score_threshs.NC_MONOCYTE_MS4A7_THRESH, np.inf)],
                    [('expr', 'FCGR3A', pb_cd34_mc_score_threshs.NC_MONOCYTE_FCGR3A_THRESH, np.inf)],
                    [('expr', 'DTL', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_DTL_THRESH)], 
                ],
            }),
            ('cMonocyte', {
                'rules': [
                    [('expr', 'CD14', pb_cd34_mc_score_threshs.CD14_THRESH, np.inf)],
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_PRSS57_THRESH)], 
                    [('expr', 'DTL', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_DTL_THRESH)], 
                    [('expr', 'FCGR3A', -np.inf, pb_cd34_mc_score_threshs.C_MONOCYTE_FCGR3A_THRESH)],
                ],
            }),
            ('intermMonocyte', {
                'rules': [
                    [('expr', 'CD14', pb_cd34_mc_score_threshs.CD14_THRESH, np.inf)],
                    [('expr', 'PRSS57', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_PRSS57_THRESH)], 
                    [('expr', 'DTL', -np.inf, pb_cd34_mc_score_threshs.MONOCYTE_DTL_THRESH)], 
                    [('expr', 'FCGR3A', pb_cd34_mc_score_threshs.C_MONOCYTE_FCGR3A_THRESH, np.inf)],
                ],
            }),
            ('Endothel', {
                'rules': [
                    [('expr', 'VWF', pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH, np.inf)], 
                ],
            }),
                    

            # ('Doublet', {
            #     'rules': [
            #         [('expr', 'GATA1', pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH, np.inf)], 
            #         [
            #             ('expr', 'HOXA9', pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH, np.inf),
            #             ('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf),
            #         ],
            #     ],
            # }),
            # ('Doublet', {
            #     'rules': [
            #         [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)],
            #         [('expr', 'DNTT', pb_cd34_mc_score_threshs.B_DNTT_THRESH, pb_cd34_mc_score_threshs.PRO_B_DNTT_THRESH)],
            #     ],
            # }),
            # ('Doublet', {
            #     'rules': [
            #         [('expr', 'DNTT', pb_cd34_mc_score_threshs.B_DNTT_THRESH, pb_cd34_mc_score_threshs.PRO_B_DNTT_THRESH)],
            #         [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)],
            #     ],
            # }),
            # ('Doublet', {
            #     'rules': [
            #         [('expr', 'HBD', pb_cd34_mc_score_threshs.CLP_NKTDP_HBD_THRESH, np.inf)], 
            #         [('expr', 'ACY3', pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_ACY3_THRESH, np.inf)], 
            #     ],
            # }),
            # ('Doublet', {
            #     'rules': [
            #         [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
            #         [('expr', 'ID2', pb_cd34_mc_score_threshs.B_ID2_THRESH, np.inf)],
            #         [('expr', 'VPREB3', -np.inf, pb_cd34_mc_score_threshs.CLP_VPREB3_THRESH)], 
            #         [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
            #         [('expr', 'IL7R', -np.inf, pb_cd34_mc_score_threshs.PRO_B_IL7R_THRESH)],
            #     ],
            # }),
            # ('Doublet', {
            #     'rules': [
            #         [('expr', 'ACY3', pb_cd34_mc_score_threshs.ACY3_IGHA1_DOUBLET_ACY3_THRESH, np.inf)],
            #         [('expr', 'IGHA1', pb_cd34_mc_score_threshs.ACY3_IGHA1_DOUBLET_IGHA1_THRESH, np.inf)],
            #     ],
            # }),

            ('mregDC', {
                'rules': [
                    [('expr', 'LAMP3', pb_cd34_mc_score_threshs.OUTLIER_LAMP3_THRESH, np.inf)], 
                ],
            }),
            # ('MTSS1_CLEC4A_B_outlier', {
            #     'rules': [
            #         [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)], 
            #         [('expr', 'MTSS1', pb_cd34_mc_score_threshs.B_OUTLIER_MTSS1_THRESH, np.inf)], 
            #     ],
            # }),

            ('Doublet', {
                'rules': [
                    [('expr', 'GATA1', pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH, np.inf)], 
                    [('expr', 'LTB', pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH, np.inf)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'GATA1', pb_cd34_mc_score_threshs.CLP_AND_T_AND_B_GATA1_THRESH, np.inf)], 
                    [('expr', 'LTB', pb_cd34_mc_score_threshs.MEBEMP_L_EP_BEMP_LTB_THRESH, np.inf)], 
                    [('expr', 'AVP', -np.inf, pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf)],
                    [('expr', 'PRSS57', pb_cd34_mc_score_threshs.B_PRSS57_THRESH, np.inf)],
                    [('expr', 'DNTT', -np.inf, pb_cd34_mc_score_threshs.PRO_B_DNTT_THRESH)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'LMO4', pb_cd34_mc_score_threshs.BEMP_LMO4_THRESH, np.inf)], 
                    [('expr', 'HBB', pb_cd34_mc_score_threshs.BEMP_HBB_THRESH, np.inf)], 
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH, np.inf)], 
                    [
                        ('expr', 'HDC', pb_cd34_mc_score_threshs.GMP_HDC_THRESH, np.inf),
                        ('expr', 'ACY3', pb_cd34_mc_score_threshs.GMP_ACY3_THRESH, np.inf),
                    ], 
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'AVP', pb_cd34_mc_score_threshs.MEBEMP_E_AVP_THRESH, np.inf)], 
                    [('expr', 'LTB', pb_cd34_mc_score_threshs.MEBEMP_LTB_THRESH, np.inf)], 
                    [('expr', 'CD79A', -np.inf, pb_cd34_mc_score_threshs.B_CD79A_THRESH)],
                    [('expr', 'HOXA9', -np.inf, pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH)],
                    [('expr', 'MPO', -np.inf, pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH)], 
                    [('expr', 'ACY3', -np.inf, pb_cd34_mc_score_threshs.NKTDP_ACY3_THRESH)], 
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'PRSS57', pb_cd34_mc_score_threshs.MONOCYTE_PRSS57_THRESH, np.inf)], 
                    [('expr', 'CD14', pb_cd34_mc_score_threshs.HSPC_CD14_THRESH, np.inf)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                ],
            }),
            ('Doublet', {
                'rules': [
                    [
                        ('expr', 'IGHA1', pb_cd34_mc_score_threshs.CLP_IGHA1_THRESH, np.inf),
                        ('expr', 'CD3D', pb_cd34_mc_score_threshs.CLP_CD3D_THRESH, np.inf),
                    ], 
                    [('expr', 'HOXA9', pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD3D', pb_cd34_mc_score_threshs.CLP_CD3D_THRESH, np.inf)], 
                    [('expr', 'HOXA9', pb_cd34_mc_score_threshs.NKT_HOXA9_THRESH, np.inf)],
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'CD14', pb_cd34_mc_score_threshs.HSPC_CD14_THRESH, np.inf)],
                    [('expr', 'DTL', pb_cd34_mc_score_threshs.MONOCYTE_DTL_THRESH, np.inf)], 
                    [('expr', 'VWF', -np.inf, pb_cd34_mc_score_threshs.ENDOTHEL_VWF_THRESH)], 
                ],
            }),
            ('Doublet', {
                'rules': [
                    [('expr', 'MPO', pb_cd34_mc_score_threshs.GMP_E_MPO_THRESH, np.inf)], 
                    [
                        ('expr', 'HOXA9', pb_cd34_mc_score_threshs.CLP_HOXA9_THRESH, np.inf),
                        ('expr', 'CD79A', pb_cd34_mc_score_threshs.B_CD79A_THRESH, np.inf), 
                        ('expr', 'DNTT', pb_cd34_mc_score_threshs.GMP_E_DNTT_THRESH, np.inf), 
                    ],
                ],
            }),
            
            
        ],
    },
    'dirty_nimrod_atlas': {
        
        # copied from /home/nimrodra/proj/blood_aging/metacells_workdir/all_cells_fil2_500000_cdata_output.h5ad
        # 'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/all_cells_fil2_500000_cdata_output.h5ad',
        
        # 'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/some_processing_by_orenmil/all_cells_fil2_500000_cdata_output_with_cell_type_exp_indiv.h5ad',
        'c_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/230918_frozen_model/orenmil_processing/final_c_231224.h5ad', # NOTE: created by mds/mds_one_offs/231122_fix_230918_nimrod_model.ipynb

        # from my first days in the lab, but should be fine to just take the donor id from it...
        
        # f47b44e7c9d3d8999c7324ba630bdbd3  /dummy/dummy/dummy/raid/symlink_to_nimrod_blood_aging_dir/metacell/mat.all_our_data_good_libs_fil2.Rda
        # 02947255bdd8dd722ccea01ae7b0d3cf  /dummy/dummy/dummy/raid/blood_aging_from_raw_data/MC1_raw_data_dir/mat.all_our_data_good_libs_fil2.Rda
        # 'atlas_cells_with_donor_id_ad_file_path': '/dummy/dummy/dummy/raid/blood_aging_from_raw_data/MC1_raw_data_dir/raw_converted_to_h5ad/all_our_data_good_libs_fil2.h5ad', # NOTE: 230703: i think we want atlas_cells_ad_file_path instead


        
        'mc_ad_file_path': '/dummy/dummy/dummy/raid/mds/from_nimrod/230918_frozen_model/orenmil_processing/final_mc_231224.h5ad', # NOTE: created by mds/mds_one_offs/231122_fix_230918_nimrod_model.ipynb
        
        # 'atlas_cell_types_csv_file_path': '/dummy/dummy/dummy/tanay_group/misc_standalones/221208_nimrod_metacell_types.tsv', # not needed anymore because now the atlas_ad file already contains the type.
        # 'cell_state_column_name': 'type',
        'state_column_name': 'state',
        'orig_state_column_name': 'orig_type',
        # 'atlas_cell_type_colors_csv_file_path': '/dummy/dummy/dummy/tanay_group/mds/cell_type_colors-2022-12-26.csv',
        'atlas_cell_type_colors_csv_file_path': DEFAULT_CELL_STATE_COLORS_CSV_FILE_PATH,

        'max_batch_kruskal_pval_to_ignore_gene_in_projection_to_atlas': 1e-6,

        # for the atlas, but probably not including cells before platelet filtering etc.
        # Individual's age and sex:
        # /home/nimrodra/proj/blood_aging/supp_tables/s1_individual_metadata.tsv
        # individual's cbcs:
        # /home/nimrodra/proj/blood_aging/supp_tables/s2_individual_cbc.tsv
        'project_while_ignoring_ultima_illumina_problematic_genes': True,

        'cell_state_and_info_list': [
            # NOTE: this is dirty_nimrod_atlas.
            # TODO: if needed, update this according to nimrod_atlas
            ('ERYP', {
                'rules': [
                    [('expr', 'AHSP', -15.2, np.inf)],
                    [('expr', 'PLEK', -np.inf, -13.3)],
                ],
                # 'override_others': True,
            }),
            ('MEBEMP-L', {
                'rules': [
                    [
                        ('expr', 'PLEK', -13.3, np.inf),
                        # ('expr', 'CA1', -np.inf, -13.7),
                        ('expr', 'AHSP', -np.inf, -15.2),
                    ],
                ],
                'obs_attr_val': [
                    ('orig_type', ['ERYP']),
                ],
                # 'override_others': True,
            }),
            # ('MEBEMP-L', {
            #     'rules': [
            #         [('AHSP', -np.inf, -15)],
            #     ],
            #     'obs_attr_val': [
            #         ('orig_type', ['ERYP']),
            #     ],
            #     # 'override_others': True,
            # }),
            ('pre-pro-B?', {
                'rules': [
                    [('DNTT', -10.5, np.inf)], # seems like VPREB1 starts rising after DNTT
                    [('VPREB3', -np.inf, -14)], 
                ],
                # 'obs_attr_val': [
                #     ('orig_type', ['CLP-L']),
                # ],
                # 'override_others': True,
            }),
            ('pre-B?', {
                # NOTE: we have multiple MCs here with KIAA0087 (and FAM30A and ACY3) higher than we might expect. yet i guess this is one of the problems of negative (rather than positive) markers - they are much more dependent on mRNA decay and cell division rates. Maybe the ones with higher than expected CLP-M genes are those that performed relatively few cell divisions, if any, since differentiating from CLP-M???
                'rules': [
                    [('DNTT', -10.5, np.inf)],
                    [('VPREB3', -14, np.inf)], 
                ],
                # 'obs_attr_val': [
                #     ('orig_type', ['CLP-L']),
                # ],
                # 'override_others': True,
            }),
            ('CLP-M', {
                'rules': [
                    [('DNTT', -np.inf, -10.5)],
                ],
                'obs_attr_val': [
                    ('orig_type', ['CLP-L']),
                ],
                # 'override_others': True,
            }),
            ('Endothel-SELE-outlier?', {
                'rules': [
                    [('SELE', -11, np.inf)],
                ],
                # 'override_others': True,
            }),
            ('Doublet', {
                'rules': [
                    [('CCL4', -13.3, np.inf)],
                ],
                'obs_attr_val': [
                    ('orig_type', ['MPP']),
                ],
                # 'override_others': True,
            }),
        ],
    },
    'mds_healthy_controls_atlas': {
        'atlas_mc_model_name': 'final_all_healthy',
        'cell_state_column_name': 'state',
        'file_path_attr_prefix': 'mds_healthy_controls_',
    },

    'cell_type_colors_csv_file_path': DEFAULT_CELL_STATE_COLORS_CSV_FILE_PATH,

    'min_dominant_donor_fraction': 0.9,
    
    'pb_cd34_enriched_mask_and_c_info_list': [
        # ('not_dc_nkt_endothel_b', {
        #     # NOTE: IMPORTANT: this is also excluding a lot of Monocytes!!! didnt use this in the end because of that...
        #     'rules': [
        #         [('metadata', 'higher_in_dc_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH)],
        #         [('metadata', 'higher_in_b_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_ALL_OTHER_HSPCS_THRESH)],
        #         [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
        #         [('metadata', 'nkt_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.NKT_SOMEWHAT_SPECIFIC_THRESH)],
        #     ],
        # }),
        ('B', {
            'rules': [
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'higher_in_dc_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH)],
                [('metadata', 'c_intermediate_nc_monocyte_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.C_INTERMEDIATE_NC_MONOCYTE_SOMEWHAT_SPECIFIC_THRESH)],
                [('metadata', 'strong_c_monocyte_specific', -np.inf, pb_cd34_c_score_threshs.STRONG_C_MONOCYTE_SPECIFIC_THRESH)],
                [('metadata', 'ighm_ighg', -np.inf, pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHM_IGHG_THRESH)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH)],

                [('metadata', 'higher_in_b_than_all_other_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_ALL_OTHER_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_b_than_nkt', pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_NKT_THRESH, np.inf)],
                [('metadata', 'higher_in_b_than_pre_b', pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_PRE_B_THRESH, np.inf)],
                [('metadata', 'higher_in_b_than_plasmablast', pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_PLASMABLAST_THRESH, np.inf)],
                
                [('metadata', 'higher_in_monocyte_than_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTE_THAN_B_THRESH)],
                [('metadata', 'higher_in_b_than_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH, np.inf)],
            ],
        }),
        ('monocyte', {
            'rules': [
                [('metadata', 'B', [False])],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                [('metadata', 'ighm_ighg', -np.inf, pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHM_IGHG_THRESH)],
                [('metadata', 'higher_in_c_dc1_than_as_dc_and_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC1_THAN_AS_DC_AND_C_DC2_THRESH)],
                [('metadata', 'higher_in_interm_monocyte_than_c_dc2', pb_cd34_c_score_threshs.HIGHER_IN_INTERM_MONOCYTE_THAN_C_DC2_THRESH, np.inf)],
                [('metadata', 'higher_in_c_dc2_than_interm_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC2_THAN_INTERM_MONOCYTE_THRESH)],
                [('metadata', 'c_intermediate_nc_monocyte_somewhat_specific', pb_cd34_c_score_threshs.C_INTERMEDIATE_NC_MONOCYTE_SOMEWHAT_SPECIFIC_THRESH, np.inf)],
                # [('metadata', 'strong_c_monocyte_specific', pb_cd34_c_score_threshs.STRONG_C_MONOCYTE_SPECIFIC_THRESH, np.inf)],
                
                [('metadata', 'higher_in_monocyte_than_b', pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTE_THAN_B_THRESH, np.inf)],
                [('metadata', 'higher_in_b_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH)],
            ],
        }),
        ('cDC1', {
            'rules': [
                [('metadata', 'B', [False])],
                [('metadata', 'higher_in_b_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                [('metadata', 'higher_in_nc_monocyte_than_c_and_interm_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NC_MONOCYTE_THAN_C_MONOCYTE_THRESH)],
                [('metadata', 'strong_c_monocyte_specific', -np.inf, pb_cd34_c_score_threshs.STRONG_C_MONOCYTE_SPECIFIC_THRESH)],
                [('metadata', 'higher_in_interm_monocyte_than_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_INTERM_MONOCYTE_THAN_C_DC2_THRESH)],
                [('metadata', 'higher_in_c_dc1_than_as_dc_and_c_dc2', pb_cd34_c_score_threshs.HIGHER_IN_C_DC1_THAN_AS_DC_AND_C_DC2_THRESH, np.inf)],
                [('metadata', 'higher_in_dc_than_all_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_dc_than_nkt', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_NKT_THRESH, np.inf)],

                [('metadata', 'higher_in_p_dc_than_as_and_c_dc', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_P_DC_THAN_AS_AND_C_DC_THRESH)],
                [('metadata', 'higher_in_monocyte_than_c_dc1', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTE_THAN_C_DC1_THRESH)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_DC_PROG_THRESH)],
            ],
        }),
        ('cDC2-prog?', { # 240403: left this one without testing much. i guess it is not very specific (i.e., cells of other states might get c_state=cDC2)...
            'rules': [
                # # Human Dendritic Cell Subsets, Ontogeny, and Impact on HIV Infection (https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01088/full, 2019): "In peripheral blood, what originally comprised of 3 key DC subsets (CD141+ cDC1, CD1c+ cDC2, and CD123+ pDC) has now been expanded to 6 putative subsets (cDC1, cDC2-A/B, CD16+ DC, Axl+ DC, and pDC) which can be distinguished by expression of CD11c, CD16, Clec9a/CADM1, CD1c, CD32b, CD163, Axl, Siglec6, and CD123."
                # # "Conventional DCs (cDCs), also known as myeloid DCs, can be defined as CD11c+ CD123−"
                # # "cDC2s represent the major subset of myeloid DC in blood"
                # 'THBD', # AKA CD141
                # 'CD1C', 
                # 'IL3RA', # AKA CD123

                # 'FCGR3A', # AKA CD16
                # 'AXL',
                # 'ITGAX', # AKA CD11c
                # 'CLEC9A',
                # 'CADM1',
                # 'FCGR2B', # AKA CD32b
                # 'CD163',
                # 'SIGLEC6',

                # # my additions
                # 'ENHO', # https://www.proteinatlas.org/ENSG00000168913-ENHO: "Immune cell enriched (myeloid DC)"
                # 'MRC1', # https://www.proteinatlas.org/ENSG00000260314-MRC1: "Immune cell enriched (myeloid DC)"
                # 'NCKAP5', # https://www.proteinatlas.org/ENSG00000176771-NCKAP5/immune+cell: in the monaco dataset myeloid DCs have highest NCKAP5, though still pretty low.
                # 'CLEC10A', # https://www.proteinatlas.org/ENSG00000132514-CLEC10A: "Immune cell enriched (myeloid DC)"
                # ###'XCR1', # ? not higher in cDCs in 1k1k
                # ###'IDO1', # ? not higher in cDCs in 1k1k
                # ###'CD207', # ? not higher in cDCs in 1k1k

                [('metadata', 'B', [False])],
                [('metadata', 'higher_in_b_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                
                [('metadata', 'c_intermediate_nc_monocyte_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.C_INTERMEDIATE_NC_MONOCYTE_SOMEWHAT_SPECIFIC_C_DC2_THRESH)],
                [('metadata', 'strong_c_monocyte_specific', -np.inf, pb_cd34_c_score_threshs.STRONG_C_MONOCYTE_SPECIFIC_C_DC2_THRESH)],
                [('metadata', 'higher_in_interm_monocyte_than_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_INTERM_MONOCYTE_THAN_C_DC2_THRESH)],
                
                [('metadata', 'higher_in_c_dc2_than_all_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_C_DC2_THAN_ALL_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_dc_than_nkt', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_NKT_THRESH, np.inf)],

                [('metadata', 'higher_in_p_dc_than_as_and_c_dc', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_P_DC_THAN_AS_AND_C_DC_THRESH)],
                [('metadata', 'higher_in_c_dc1_than_as_dc_and_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC1_THAN_AS_DC_AND_C_DC2_THRESH)],
                [('metadata', 'higher_in_as_dc_than_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_AS_DC_THAN_C_DC2_THRESH)],
                [('metadata', 'higher_in_c_dc2_than_interm_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_C_DC2_THAN_INTERM_MONOCYTE_THRESH, np.inf)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_DC_PROG_THRESH)],
            ],
        }),
        ('high_cMonocyte_sig_cDC2-prog?', {
            'rules': [
                [('metadata', 'B', [False])],
                [('metadata', 'higher_in_b_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH)],
                [('metadata', 'higher_in_b_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                [('metadata', 'strong_c_monocyte_specific', pb_cd34_c_score_threshs.STRONG_C_MONOCYTE_SPECIFIC_C_DC2_THRESH, np.inf)],
                
                [('metadata', 'higher_in_c_dc2_than_all_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_C_DC2_THAN_ALL_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_dc_than_nkt', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_NKT_THRESH, np.inf)],

                [('metadata', 'higher_in_p_dc_than_as_and_c_dc', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_P_DC_THAN_AS_AND_C_DC_THRESH)],
                [('metadata', 'higher_in_c_dc1_than_as_dc_and_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC1_THAN_AS_DC_AND_C_DC2_THRESH)],
                [('metadata', 'higher_in_as_dc_than_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_AS_DC_THAN_C_DC2_THRESH)],
                [('metadata', 'higher_in_c_dc2_than_interm_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_C_DC2_THAN_INTERM_MONOCYTE_THRESH, np.inf)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_DC_PROG_THRESH)],
            ],
        }),
        ('pDC_AS-DC', {
            'rules': [
                [('metadata', 'B', [False])],
                [('metadata', 'higher_in_b_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH)],
                [('metadata', 'c_intermediate_nc_monocyte_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.C_INTERMEDIATE_NC_MONOCYTE_SOMEWHAT_SPECIFIC_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                [('metadata', 'higher_in_nc_monocyte_than_c_and_interm_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NC_MONOCYTE_THAN_C_MONOCYTE_THRESH)],
                [('metadata', 'strong_c_monocyte_specific', -np.inf, pb_cd34_c_score_threshs.STRONG_C_MONOCYTE_SPECIFIC_THRESH)],
                [('metadata', 'higher_in_interm_monocyte_than_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_INTERM_MONOCYTE_THAN_C_DC2_THRESH)],

                [('metadata', 'higher_in_dc_than_all_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_dc_than_nkt', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_NKT_THRESH, np.inf)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_DC_PROG_THRESH)],
            ],
        }),
        ('mregDC', {
            'rules': [
                # https://www.tandfonline.com/doi/full/10.1080/2162402X.2023.2294564: "population known as mregDCs (also called LAMP3+DCs).Citation32,Citation33 Remarkably, this population shows the highest expression of maturation markers, such as, LAMP3, CCR7, CD83, BIRC3 and MARCKSL1"
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9937888/ - specifically table 2.
                [('metadata', 'B', [False])],
                [('metadata', 'higher_in_b_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH)],
                [('metadata', 'monocyte', [False])],
                [('metadata', 'lamp3_outlier_specific', pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH, np.inf)],
            ],
        }),
        ('DC_except_mregDC', {
            'rules': [
                [
                    ('metadata', 'cDC1', [True]),
                    ('metadata', 'cDC2-prog?', [True]),
                    ('metadata', 'pDC_AS-DC', [True]),
                    ('metadata', 'high_cMonocyte_sig_cDC2-prog?', [True]),
                ],
            ],
        }),
        ('DC', {
            'rules': [
                [
                    ('metadata', 'DC_except_mregDC', [True]),
                    ('metadata', 'mregDC', [True]),
                ],
            ],
        }),
        ('not_dc_nkt_monocyte_endothel_b_lamp3', {
            'rules': [
                [('metadata', 'c_intermediate_nc_monocyte_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.C_INTERMEDIATE_NC_MONOCYTE_SOMEWHAT_SPECIFIC_THRESH)],
                [('metadata', 'strong_c_monocyte_specific', -np.inf, pb_cd34_c_score_threshs.STRONG_C_MONOCYTE_SPECIFIC_THRESH)],
                [('metadata', 'higher_in_dc_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH)],
                [('metadata', 'higher_in_c_dc2_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC2_THAN_ALL_HSPCS_THRESH)],
                [('metadata', 'higher_in_b_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_ALL_OTHER_HSPCS_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'nkt_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.NKT_SOMEWHAT_SPECIFIC_THRESH)],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                [('metadata', 'plasmablast_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.PLASMABLAST_SOMEWHAT_SPECIFIC_THRESH)],
                [('metadata', 'ighm_ighg', -np.inf, pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHM_IGHG_THRESH)],
            ],
        }),
        ('myeloid_hspc', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_nktdp_than_all_myeloid_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_ALL_MYELOID_HSPCS_THRESH)],
                [('metadata', 'higher_in_pre_b_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH)],
                [('metadata', 'higher_in_clp_m_than_all_myeloid_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_THRESH)],
            ],
        }),
        ('HSC_CLP_traj', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_nktdp_than_clp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_CLP_THRESH)],
                [('metadata', 'higher_in_pre_b_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH)],
                [('metadata', 'higher_in_mkp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_THRESH)],
                [('metadata', 'higher_in_mebemp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MEBEMP_L_THAN_MPP_THRESH)],
                [('metadata', 'higher_in_bemp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_THRESH)],
                [('metadata', 'higher_in_ep_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH)],
                [('metadata', 'higher_in_gmp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_MPP_THRESH)],
                # [('metadata', 'higher_in_gmp_l_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_MEBEMP_L_THRESH)],
            ],
        }),
        ('CLP_NKTDP_traj', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_clp_m_and_nktdp_than_mpp', pb_cd34_c_score_threshs.HIGHER_IN_CLP_AND_NKTDP_THAN_MPP_PRESERVING_CLP_NKTDP_TRAJ_THRESH, np.inf)],
                [('metadata', 'higher_in_bemp_than_clp_m_and_nktdp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_CLP_AND_NKTDP_PRESERVING_CLP_NKTDP_TRAJ_THRESH)],
                [('metadata', 'higher_in_pre_b_than_clp_m_and_nktdp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_CLP_AND_NKTDP_PRESERVING_CLP_NKTDP_TRAJ_THRESH)],
            ],
        }),
        ('CLP_pro_B_traj', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_clp_m_and_pro_b_than_mpp', pb_cd34_c_score_threshs.HIGHER_IN_CLP_AND_PRO_B_THAN_MPP_PRESERVING_CLP_PRO_B_TRAJ_THRESH, np.inf)],
                [('metadata', 'higher_in_nktdp_than_clp_m_and_pro_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_CLP_AND_PRO_B_PRESERVING_CLP_PRO_B_TRAJ_THRESH)],
                [('metadata', 'higher_in_pre_b_than_pro_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_PRO_B_THRESH)],
            ],
        }),
        ('pro_B_pre_B_traj', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_pre_b_than_all_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_b_than_pre_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_PRE_B_THRESH)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH)],
            ],
        }),
        ('HSC_GMP_traj', {
            'rules': [
                [('metadata', 'myeloid_hspc', [True])],
                [('metadata', 'higher_in_bemp_than_gmp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_GMP_L_THRESH)],
                [('metadata', 'higher_in_mkp_than_mebemp_l_and_eryp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_AND_ERYP_THRESH)],
                [('metadata', 'higher_in_eryp_than_gmp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ERYP_THAN_GMP_L_THRESH)],
                [('metadata', 'higher_in_mebemp_l_than_gmp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MEBEMP_L_THAN_GMP_L_PRESERVING_HSC_GMP_TRAJ_THRESH)],
            ],
        }),
        ('HSC_MEBEMP_and_later_traj', {
            'rules': [
                [('metadata', 'myeloid_hspc', [True])],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH)],
            ],
        }),
        ('MEBEMP_BEMP_traj', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_and_later_traj', [True])],
                [('metadata', 'higher_in_mkp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_THRESH)],
                [('metadata', 'higher_in_ep_than_bemp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ERYP_THAN_BEMP_PRESERVING_MEBEMP_BEMP_TRAJ_THRESH)],
                [('metadata', 'higher_in_mpp_than_bemp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_BEMP_PRESERVING_MEBEMP_BEMP_TRAJ_THRESH)],
            ],
        }),
        ('MEBEMP_MKP_traj', { # NOTE: this only takes MEBEMP-L with somewhat high higher_in_mkp_than_mebemp_l_and_eryp...
            'rules': [
                [('metadata', 'HSC_MEBEMP_and_later_traj', [True])],
                [('metadata', 'higher_in_bemp_than_mebemp_l_and_eryp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_AND_ERYP_THRESH)],
                [('metadata', 'higher_in_mkp_than_mebemp_l_and_eryp', pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_AND_ERYP_PRESERVING_MEBEMP_MKP_TRAJ_THRESH, np.inf)],
                # [('metadata', 'higher_in_eryp_than_mkp_and_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ERYP_THAN_MKP_AND_MEBEMP_L_PRESERVING_MEBEMP_MKP_TRAJ_THRESH)],
            ],
        }),
        ('HSC_MEBEMP_ERYP_traj', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_and_later_traj', [True])],
                [('metadata', 'higher_in_bemp_than_mebemp_l_and_eryp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_AND_ERYP_THRESH)],
                [('metadata', 'higher_in_mkp_than_mebemp_l_and_eryp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_AND_ERYP_THRESH)],
            ],
        }),
        # ('HSC_MEBEMP_traj', { # 240625: seems like i am not really using it...
        #     'rules': [
        #         [('metadata', 'HSC_MEBEMP_ERYP_traj', [True])],
        #         [('metadata', 'higher_in_ep_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH)],
        #     ],
        # }),
        ('MEBEMP_ERYP_traj', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_ERYP_traj', [True])],
                [('metadata', 'higher_in_mpp_than_ep', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_ERYP_PRESERVING_MEBEMP_ERYP_TRAJ_THRESH)],
            ],
        }),
        ('CLP-E', {
            'rules': [
                [('metadata', 'HSC_CLP_traj', [True])],
                [('metadata', 'higher_in_clp_m_than_all_myeloid_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_CLP_E_LOW_THRESH, pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_CLP_E_HIGH_THRESH)],
            ],
        }),
        ('CLP-E2', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_nktdp_than_clp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_CLP_THRESH)],
                [('metadata', 'higher_in_pre_b_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH)],
                [('metadata', 'higher_in_plasmablast_than_clp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PLASMABLAST_THAN_CLP_THRESH)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH)],

                [('metadata', 'higher_in_mebemp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MEBEMP_L_THAN_MPP_THRESH)],
                # [('metadata', 'higher_in_bemp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_THRESH)],
                # [('metadata', 'higher_in_ep_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH)],
                # [('metadata', 'higher_in_gmp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_MPP_THRESH)],
                
                [('metadata', 'higher_in_clp_m_than_all_myeloid_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_CLP_E_LOW_THRESH, pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_CLP_E_HIGH_THRESH)],
            ],
        }),
    ],
    'disease_pb_cd34_enriched_cell_state_and_c_info_list': [
        ('N257_MEBEMP', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_and_later_traj', [False])],
                [('metadata', 'higher_in_mebemp_l_than_mpp', pb_cd34_c_score_threshs.HIGHER_IN_MEBEMP_L_THAN_MPP_THRESH, np.inf)],
                # [('metadata', 'higher_in_dc_than_nkt', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_NKT_THRESH, np.inf)],
            ],
        }),
    ],
    'pb_cd34_enriched_cell_state_and_c_info_list': [
        ('monocyte_B_doublet', {
            'rules': [
                [('metadata', 'higher_in_monocyte_than_b', pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTE_THAN_B_THRESH, np.inf)],
                [('metadata', 'higher_in_b_than_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_MONOCYTE_THRESH, np.inf)],
            ],
        }),
        ('pDC', {
            'rules': [
                [('metadata', 'pDC_AS-DC', [True])],
                [('metadata', 'higher_in_p_and_as_dc_than_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_P_AND_AS_DC_THAN_MONOCYTE_THRESH, np.inf)],
                [('metadata', 'higher_in_p_dc_than_as_and_c_dc', pb_cd34_c_score_threshs.HIGHER_IN_P_DC_THAN_AS_AND_C_DC_THRESH, np.inf)],
                [('metadata', 'higher_in_as_dc_than_p_dc', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_AS_DC_THAN_P_DC_THRESH)],
            ],
        }),
        ('AS-DC', {
            'rules': [
                [('metadata', 'pDC_AS-DC', [True])],
                [('metadata', 'higher_in_p_and_as_dc_than_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_P_AND_AS_DC_THAN_MONOCYTE_THRESH, np.inf)],
                [('metadata', 'higher_in_p_dc_than_as_and_c_dc', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_P_DC_THAN_AS_AND_C_DC_THRESH)],
                [('metadata', 'higher_in_c_dc_and_c_monocyte_than_as_dc', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC_THAN_AS_DC_THRESH)],
                [('metadata', 'higher_in_c_dc1_than_as_dc_and_c_dc2', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC1_THAN_AS_DC_AND_C_DC2_THRESH)],
                [('metadata', 'higher_in_as_dc_than_p_dc', pb_cd34_c_score_threshs.HIGHER_IN_AS_DC_THAN_P_DC_THRESH, np.inf)],
                [('metadata', 'higher_in_as_dc_than_c_dc2', pb_cd34_c_score_threshs.HIGHER_IN_AS_DC_THAN_C_DC2_THRESH, np.inf)],
            ],
        }),
        ('cDC1', {
            'rules': [
                [('metadata', 'cDC1', [True])],
            ],
        }),
        ('cDC2-prog?', {
            'rules': [
                [('metadata', 'cDC2-prog?', [True])],
            ],
        }),
        ('DC-prog?', {
            'rules': [
                [('metadata', 'monocyte', [False])],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'higher_in_dc_than_all_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_DC_PROG_THRESH, np.inf)],
            ],
        }),
        ('high_cMonocyte_sig_cDC2-prog?', {
            'rules': [
                [('metadata', 'high_cMonocyte_sig_cDC2-prog?', [True])],
            ],
        }),
        ('Endothel', {
            'rules': [
                [('metadata', 'B', [False])],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],
                [('metadata', 'endothel_somewhat_specific', pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH, np.inf)],
            ],
        }),
        ('mregDC', {
            'rules': [
                [('metadata', 'mregDC', [True])],
            ],
        }),
        ('cMonocyte', {
            'rules': [
                [('metadata', 'monocyte', [True])],
                [('metadata', 'higher_in_nc_monocyte_than_c_and_interm_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NC_MONOCYTE_THAN_C_MONOCYTE_THRESH)],
                [('metadata', 'higher_in_interm_monocyte_than_c_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_INTERM_MONOCYTE_THAN_C_MONOCYTE_THRESH)],
                # [('metadata', 'higher_in_dc_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_MONOCYTE_THRESH)],
                # [('metadata', 'higher_in_monocytes_than_c_dc', pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTES_THAN_C_DC_THRESH, np.inf)],
            ],
        }),
        ('intermMonocyte', {
            'rules': [
                [('metadata', 'monocyte', [True])],
                [('metadata', 'higher_in_nc_monocyte_than_c_and_interm_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NC_MONOCYTE_THAN_C_MONOCYTE_THRESH)],
                [('metadata', 'higher_in_interm_monocyte_than_c_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_INTERM_MONOCYTE_THAN_C_MONOCYTE_THRESH, np.inf)],
                # [('metadata', 'higher_in_dc_than_monocyte', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_MONOCYTE_THRESH)],
                # [('metadata', 'higher_in_monocytes_than_c_dc', pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTES_THAN_C_DC_THRESH, np.inf)],
            ],
        }),
        ('ncMonocyte', {
            'rules': [
                [('metadata', 'monocyte', [True])],
                [('metadata', 'higher_in_nc_monocyte_than_c_and_interm_monocyte', pb_cd34_c_score_threshs.HIGHER_IN_NC_MONOCYTE_THAN_C_MONOCYTE_THRESH, np.inf)],
                # [('metadata', 'higher_in_nc_monocyte_than_nkt', pb_cd34_c_score_threshs.HIGHER_IN_NC_MONOCYTE_THAN_NKT_THRESH, np.inf)],
                # [('metadata', 'higher_in_c_dc_than_nc_monocytes', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC_THAN_NC_MONOCYTES_THRESH)],
                
            ],
        }),
        ('NKT', {
            'rules': [
                [('metadata', 'B', [False])],
                [('metadata', 'monocyte', [False])],
                [('metadata', 'higher_in_dc_than_nkt', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_NKT_THRESH)],
                [('metadata', 'higher_in_b_than_nkt', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_NKT_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'lamp3_outlier_specific', -np.inf, pb_cd34_c_score_threshs.LAMP3_OUTLIER_SPECIFIC_THRESH)],

                [('metadata', 'nkt_somewhat_specific', pb_cd34_c_score_threshs.NKT_SOMEWHAT_SPECIFIC_THRESH, np.inf)],
                [('metadata', 'higher_in_nkt_than_nktdp', pb_cd34_c_score_threshs.HIGHER_IN_NKT_THAN_NKTDP_THRESH, np.inf)],
            ],
        }),
        ('naive-B', {
            'rules': [
                [('metadata', 'B', [True])],
                [('metadata', 'higher_in_naive_b_than_memory_b', pb_cd34_c_score_threshs.HIGHER_IN_NAIVE_B_THAN_MEMORY_B_THRESH, np.inf)],
                [('metadata', 'higher_in_ccdc88a_high_b_than_other_mature_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_CCDC88A_HIGH_B_THAN_OTHER_MATURE_B_THRESH)],
            ],
        }),
        ('memory-B', {
            'rules': [
                [('metadata', 'B', [True])],
                [('metadata', 'higher_in_naive_b_than_memory_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NAIVE_B_THAN_MEMORY_B_THRESH)],
                [('metadata', 'higher_in_aged_b_than_memory_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_AGED_B_THAN_MEMORY_B_THRESH)],
                [('metadata', 'higher_in_ccdc88a_high_b_than_other_mature_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_CCDC88A_HIGH_B_THAN_OTHER_MATURE_B_THRESH)],
            ],
        }),
        ('aged-B?', {
            'rules': [
                [('metadata', 'B', [True])],
                [('metadata', 'higher_in_naive_b_than_memory_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NAIVE_B_THAN_MEMORY_B_THRESH)],
                [('metadata', 'higher_in_ccdc88a_high_b_than_other_mature_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_CCDC88A_HIGH_B_THAN_OTHER_MATURE_B_THRESH)],
                [('metadata', 'higher_in_aged_b_than_memory_b', pb_cd34_c_score_threshs.HIGHER_IN_AGED_B_THAN_MEMORY_B_THRESH, np.inf)],
            ],
        }),
        ('CCDC88A-high-B', {
            'rules': [
                [('metadata', 'B', [True])],
                [('metadata', 'higher_in_ccdc88a_high_b_than_other_mature_b', pb_cd34_c_score_threshs.HIGHER_IN_CCDC88A_HIGH_B_THAN_OTHER_MATURE_B_THRESH, np.inf)],
            ],
        }),
        ('plasmablast', {
            'rules': [
                [('metadata', 'ighm_ighg', -np.inf, pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHM_IGHG_THRESH)],
                [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
                [('metadata', 'plasmablast_somewhat_specific', pb_cd34_c_score_threshs.PLASMABLAST_SOMEWHAT_SPECIFIC_THRESH, np.inf)],
                [('metadata', 'higher_in_b_than_plasmablast', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_PLASMABLAST_THRESH)],
                [('metadata', 'higher_in_monocyte_than_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTE_THAN_B_THRESH)],
            ],
        }),
        ('plasmablast', {
            'rules': [
                [('metadata', 'ighm_ighg', pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHM_IGHG_THRESH, np.inf)],
                [('metadata', 'igha_iglc', pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHA_IGLC_THRESH, np.inf)],
            ],
        }),
        ('plasmablast_ighm_ighg', {
            'rules': [
                [('metadata', 'ighm_ighg', pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHM_IGHG_THRESH, np.inf)],
                [('metadata', 'igha_iglc', -np.inf, pb_cd34_c_score_threshs.IGHM_IGHG_PLASMABLAST_IGHA_IGLC_THRESH)],
            ],
        }),
        ('pro-B?', {
            'rules': [
                [('metadata', 'pro_B_pre_B_traj', [True])],
                [('metadata', 'higher_in_pre_b_than_pro_b', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_PRO_B_THRESH)],
            ],
        }),
        ('pre-B?', {
            'rules': [
                [('metadata', 'pro_B_pre_B_traj', [True])],
                [('metadata', 'higher_in_pre_b_than_pro_b', pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_PRO_B_THRESH, np.inf)],
            ],
        }),
        ('NKTDP', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_bemp_than_clp_m_and_nktdp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_CLP_M_AND_NKTDP_THRESH)],
                [('metadata', 'higher_in_pre_b_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH)],

                [('metadata', 'higher_in_nktdp_than_all_myeloid_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_ALL_MYELOID_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_nktdp_than_clp', pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_CLP_THRESH, np.inf)],
                # [('metadata', 'higher_in_nktdp_than_gmp_l', pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_GMP_L_THRESH, np.inf)], # seems to change pretty much nothing (even in MDS ult), so better without it, i think.
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH)],
            ],
        }),
        ('CLP', {
            'rules': [
                [('metadata', 'not_dc_nkt_monocyte_endothel_b_lamp3', [True])],
                [('metadata', 'higher_in_nktdp_than_clp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_CLP_THRESH)],
                [('metadata', 'higher_in_pre_b_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH)],
                [('metadata', 'higher_in_plasmablast_than_clp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PLASMABLAST_THAN_CLP_THRESH)],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH)],
                
                [('metadata', 'higher_in_clp_m_than_all_myeloid_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_THRESH, np.inf)],
            ],
        }),
        ('GMP-L', {
            'rules': [
                [('metadata', 'myeloid_hspc', [True])],
                [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH, np.inf)],
                [('metadata', 'higher_in_bemp_than_gmp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_GMP_L_THRESH)],
                # [('metadata', 'higher_in_cfd_tryptase_gmp_l_than_gmp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_CFD_TRYPTASE_GMP_L_THAN_GMP_L_THRESH)],
            ],
        }),
        # ('CFD_tryptase_GMP-L', {
        #     'rules': [
        #         [('metadata', 'myeloid_hspc', [True])],
        #         [('metadata', 'higher_in_gmp_l_than_all_other_hspcs', pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_ALL_OTHER_HSPCS_THRESH, np.inf)],
        #         [('metadata', 'higher_in_monocyte_than_gmp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MONOCYTE_THAN_GMP_L_THRESH)],
        #         [('metadata', 'higher_in_dc_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH)],
        #         [('metadata', 'higher_in_c_dc2_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_C_DC2_THAN_ALL_HSPCS_THRESH)],
        #         [('metadata', 'higher_in_cfd_tryptase_gmp_l_than_gmp_l', pb_cd34_c_score_threshs.HIGHER_IN_CFD_TRYPTASE_GMP_L_THAN_GMP_L_THRESH, np.inf)],
        #         # [('metadata', 'higher_in_bemp_than_gmp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_GMP_L_THRESH)],
        #     ],
        # }),
        ('BEMP', {
            # N151 has high FKBP5, which i guess leads to his BEMPs dominating some metacells. there doesn't seem to be very strong diff exp between these BEMPs and normal ones...
            'rules': [
                [('metadata', 'HSC_MEBEMP_and_later_traj', [True])],
                [('metadata', 'higher_in_mkp_than_mebemp_l_and_eryp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_AND_ERYP_THRESH)],
                [('metadata', 'higher_in_bemp_than_mebemp_l_and_eryp', pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_AND_ERYP_THRESH, np.inf)],
                [('metadata', 'higher_in_bemp_than_clp_m_and_nktdp', pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_CLP_M_AND_NKTDP_THRESH, np.inf)],
            ],
        }),
        ('MKP', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_and_later_traj', [True])],
                [('metadata', 'higher_in_bemp_than_mebemp_l_and_eryp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_AND_ERYP_THRESH)],
                [('metadata', 'higher_in_mkp_than_mebemp_l_and_eryp', pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_AND_ERYP_THRESH, np.inf)],
                # [('metadata', 'higher_in_mpp_than_mkp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_MKP_THRESH)],
            ],
        }),
        # ('high_MPP_sig_MKP', {
        #     'rules': [
        #         [('metadata', 'HSC_MEBEMP_and_later_traj', [True])],
        #         [('metadata', 'higher_in_bemp_than_mebemp_l_and_eryp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_AND_ERYP_THRESH)],
        #         [('metadata', 'higher_in_mkp_than_mebemp_l_and_eryp', pb_cd34_c_score_threshs.HIGHER_IN_MKP_THAN_MEBEMP_L_AND_ERYP_THRESH, np.inf)],
        #         [('metadata', 'higher_in_mpp_than_mkp', pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_MKP_THRESH, np.inf)],
        #     ],
        # }),
        ('ERYP', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_ERYP_traj', [True])],
              
                [('metadata', 'higher_in_mpp_than_ep', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_EP_THRESH)],
                [('metadata', 'higher_in_ep_than_mpp', pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH, np.inf)],
            ],
        }),
        ('high_MPP_sig_ERYP', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_ERYP_traj', [True])],
                
                [('metadata', 'higher_in_mpp_than_ep', pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_EP_THRESH, np.inf)],
                [('metadata', 'higher_in_ep_than_mpp', pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH, np.inf)],
            ],
        }),
        ('MEBEMP-L', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_ERYP_traj', [True])],
                # [('metadata', 'higher_in_mpp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_MEBEMP_L_THRESH)],
                [('metadata', 'higher_in_ep_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH)],
                [('metadata', 'higher_in_mebemp_l_than_mpp', pb_cd34_c_score_threshs.HIGHER_IN_MEBEMP_L_THAN_MPP_THRESH, np.inf)],
            ],
        }),
        ('HSC_MPP', {
            'rules': [
                [('metadata', 'HSC_MEBEMP_ERYP_traj', [True])],
                [('metadata', 'higher_in_ep_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH)],
                [('metadata', 'higher_in_mebemp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MEBEMP_L_THAN_MPP_THRESH)],
                
                # [('metadata', 'higher_in_mpp_than_mebemp_l', pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_MEBEMP_L_THRESH, np.inf)],
            ],
        }),
        # ('GMP-E', {
        #     'rules': [
        #         [('metadata', 'higher_in_mpp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MPP_THAN_MEBEMP_L_THRESH)],
        #         [('metadata', 'higher_in_mebemp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_MEBEMP_L_THAN_MPP_THRESH)],
                
        #         [('metadata', 'higher_in_bemp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_THRESH)],
        #         [('metadata', 'higher_in_ep_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH)],
        #         [('metadata', 'higher_in_gmp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_MPP_THRESH)],
        #         [('metadata', 'higher_in_gmp_l_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_MEBEMP_L_THRESH)],
        #         [('metadata', 'higher_in_nktdp_than_all_myeloid_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_ALL_HSPCS_THRESH)],
        #         [('metadata', 'higher_in_clp_m_than_all_myeloid_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_THRESH)],
        #         [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
        #         [('metadata', 'nkt_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.NKT_SOMEWHAT_SPECIFIC_THRESH)],
        #         [('metadata', 'higher_in_b_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_ALL_OTHER_HSPCS_THRESH)],
        #         [('metadata', 'higher_in_dc_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH)],
        #         [('metadata', 'strong_c_monocyte_specific', -np.inf, pb_cd34_c_score_threshs.STRONG_MONOCYTE_SPECIFIC_THRESH)],
        #         [('metadata', 'higher_in_pre_b_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH)],
        #     ],
        # }),
        # ('HSC_MPP_MEBEMP', {
        #     'rules': [
        #         # [('metadata', 'higher_in_gmp_l_than_ep', pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_EP_THRESH, np.inf)],
                
        #         [('metadata', 'higher_in_bemp_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_BEMP_THAN_MEBEMP_L_THRESH)],
        #         [('metadata', 'higher_in_ep_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_EP_THAN_MPP_THRESH)],
        #         [('metadata', 'higher_in_gmp_l_than_mpp', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_MPP_THRESH)],
        #         [('metadata', 'higher_in_gmp_l_than_mebemp_l', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_GMP_L_THAN_MEBEMP_L_THRESH)],
        #         [('metadata', 'higher_in_nktdp_than_all_myeloid_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_NKTDP_THAN_ALL_HSPCS_THRESH)],
        #         [('metadata', 'higher_in_clp_m_than_all_myeloid_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_THRESH)],
        #         [('metadata', 'endothel_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_ENDOTHEL_THAN_ALL_EXCEPT_UNKNOWN_THRESH)],
        #         [('metadata', 'nkt_somewhat_specific', -np.inf, pb_cd34_c_score_threshs.NKT_SOMEWHAT_SPECIFIC_THRESH)],
        #         [('metadata', 'higher_in_b_than_all_other_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_B_THAN_ALL_OTHER_HSPCS_THRESH)],
        #         [('metadata', 'higher_in_dc_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_DC_THAN_ALL_HSPCS_THRESH)],
        #         [('metadata', 'strong_c_monocyte_specific', -np.inf, pb_cd34_c_score_threshs.STRONG_MONOCYTE_SPECIFIC_THRESH)],
        #         [('metadata', 'higher_in_pre_b_than_all_hspcs', -np.inf, pb_cd34_c_score_threshs.HIGHER_IN_PRE_B_THAN_ALL_HSPCS_THRESH)],
        #     ],
        # }),
    ],
    'mds_and_mpn_pb_cd34_enriched_cell_state_and_info_list': [
        ('N205-MDS', {
            'rules': [
                [('expr', 'GALNTL6', -14, np.inf)],
                
                # ('at_least', 1, [
                #     # ('expr', 'GALNTL6', -12.5, np.inf),
                #     # ('expr', 'PRDX5', -11.5, np.inf),
                    
                #     # by the way, both GALNTL6 and CACNB4 are in too_often_pretty_silent_in_both_query_and_atlas_gene_mask.
                #     ('expr', 'GALNTL6', -14.3, np.inf),
                #     # ('expr', 'CACNB4', -11.7, np.inf),
                # ]),
            ],
            # 'metadata_attr_range': [
            #     ('donor_id_N205_fraction', 0.7, 1),
            # ],
            'override_others': True,
        }),
        ('N186-MDS', {
            'rules': [
                [
                    # ('expr', 'DLK1', -7.5, np.inf), # the smaller clone?
                    # ('expr', 'ATP5F1A', -9.5, np.inf),
                    # ('expr', 'XRN2', -9.5, np.inf),
                    ('expr', 'TPGS2', -11.5, np.inf),
                    ('expr', 'PDRG1', -12.5, np.inf),
                ],
                [('expr', 'HBD', -np.inf, -12.5)], # 230704: some Erythroid progenitors seem to have somewhat high TPGS2, and N186 does not seem to have such erythroid progenitors...
            ],
            # 'metadata_attr_range': [
            #     ('donor_id_N186_fraction', 0.7, 1),
            # ],
            'override_others': True,
        }),
        ('MKP', {
            'rules': [
                
                # https://www.proteinatlas.org/ENSG00000085721-RRN3: "RRN3 homolog, RNA polymerase I transcription factor" and "RRN3 (DKFZp566E104, TIF-IA)" and "Required for efficient transcription initiation by RNA polymerase I".
                # https://www.proteinatlas.org/ENSG00000105849-POLR1F: "RNA polymerase I subunit F" and "POLR1F (A43, RPA43, TWISTNB)" and "Component of RNA polymerase I which synthesizes ribosomal RNA precursors. Through its association with RRN3/TIF-IA may be involved in recruitment of Pol I to rDNA promoters".
                # tried to look at all RNA pol I genes but they didn't seem higher in our MKPs....
                # 'POLR1A',
                # 'POLR1B',
                # 'POLR1C',
                # 'POLR1D',
                # 'POLR1E',
                # 'TWISTNB', # POLR1F
                # 'CD3EAP', # POLR1G
                # 'ZNRD1', # POLR1H
                
                # old metacell version:
                # [('expr', ['RRN3', 'GSPT1'], -11.5, np.inf)],
                # [('expr', 'PF4', -13, np.inf)],
                
                # a later old version. we don't use GSPT1 anymore because its correlation with ITGA2B acorss MKPs is not very high.
                # [('expr', ['RRN3', 'GSPT1'], -12.5, np.inf)],

                # from a later old draft
                # [('expr', GENES_AT_LEAST_64_FOLD_LOWER_IN_MKP_THAN_IN_PLATELETS_AND_LOW_IN_MKP, -np.inf, -15.5)], # not contaminated with platelets.
                # [('expr', GENES_LOWER_IN_MKP_COMPARED_TO_OTHER_MEBEMPS, -np.inf, -11)], # not other mebemps
                
                # the correlation between genes low in platelets (most of which make no sense in platelets, e.g. cell cycle) and ITGA2B is very reassuring.
                [('expr', GENES_HIGHER_IN_MKP_COMPARED_TO_EP_AND_MEBEMP_L_AND_VERY_POS_CORRELATED_WITH_ITGA2B_ACROSS_MKPS_AND_LOW_IN_PLATELETS, -11.5, np.inf)],
                [('expr', 'ITGA2B', -13, np.inf)], # AKA CD41
                [('expr', 'CD79A', -np.inf, -13.5)], # not B
                [('expr', 'HBB', -np.inf, -10)], # not post-EP
                # [('expr', 'CD34', -15, np.inf)],
            ],
        }),
    ],


    'donor_sex_palette': {'male': 'blue', 'female': 'red'},

    'metacell_metadata_statistic_type_to_column_names_in_cells_ad_obj': {
        'bool_fraction': [
            'is_ultima',
            'is_bm',
        ],
        'category_enrichment_and_mode_and_more': [
        ],
        'str_mode_mode2_and_fracs_and_nunique': [
            'exp_name',
            'diagnosis',
            'diagnosis_class',
            'donor_sex',
            'donor_id',
            'tenx_technology',
            'donor_id_and_date',
            # *(['genotype_doublet'] if clean_ad.obs['genotype_doublet'].any() else [])
        ],
        # 'category_fraction_diff': [('exp_name', 'indiv_id')],
        'nan_mean': [
            'soup_log_ratio_doublet_singlet_prob',
            'vireo_doublet_prob',
            'vireo_doublet_log_likelihood_ratio',
            'log_cell_ranger_umi_count',
            'log_num_of_non_excluded_umis',
            'log_norm_cell_ranger_umi_count',
            'excluded_umi_frac',
            'mito_umi_frac',
            'malat1_umi_frac',
            'ribo_prot_umi_frac',
            'mito_to_ribo_prot_umi_log_ratio',
            'donor_age',
        ],
        # 'nan_mean_and_median_and_std': [
        # ],
    },

    'forced_cell_types_by_metadata': {
        # 'ERYP': [
    },

    'doublet_and_multiplet_identification': {
        'mc_model_with_genotype_doublets': {
            'target_metacell_umis': int(170e3),
            'target_metacell_size': 80,
        },
        'mc_genotype_doublet_enrichment': {
            'dataset_out_dir_path_to_max_mc_genotype_doublet_freq': {
                # os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240222_pb_exps_potentially_containing_mpn_illu'): 0.24,
                # os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240222_pb_exps_potentially_containing_mpn_ult'): 0.45,
                # os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240222_pb_exps_potentially_containing_mds_or_cytopenia_illu'): 0.33,
                os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240618_pb_illu_mds_cytopenia_normal'): 0.33,
                os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240623_pb_ult_mds_cytopenia_normal'): 0.45,
                os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, '240222_pb_exps_potentially_containing_mds_or_cytopenia_ult'): 0.45,
            },
            'min_num_of_filtered_cells_in_mc': 10,
            'min_num_of_donors_identified_in_exp': 2,
        },
        
        
        # 'max_soup_or_vireo_doublet_fraction_per_metacell': 0.75,
        'excluding_barcodes_with_high_norm_cell_ranger_umi_count': {
            'num_of_clusters': 10,
            'min_exp_cluster_cell_count': 40,
            'max_log_norm_cell_ranger_umi_count': 1.3,
        },
    },
    'traj_analyses': {
        'analysis_name_to_params': {
            # NOTE: maybe better to instead just use union of cell states instead of trajectories??? i guess often not. this might work when c_state is somewhat hierarchical (e.g. for ERYP and high_MPP_sig_ERYP)..

            # good ones
            'CLP_sig_across_hsc_clp_traj': dict(
                traj_pos_obs_col_name='higher_in_clp_m_than_all_myeloid_hspcs',
                c_mask_col='HSC_CLP_traj',
                
                x_is_cell_freq=True,
                # x_is_cell_freq_or_cum_freq_out_of_background=True,
                # background_c_mask_for_freq_col='not_dc_nkt_monocyte_endothel_b_lamp3',
                
                traj_pos_min_truncate_val=-10,
                traj_pos_max_truncate_val=-5.5,
                bin_count=9,
                min_donor_cell_count=100,
            ),
            'higher_in_gmp_l_than_mpp_across_hsc_gmp_traj': dict(
                traj_pos_obs_col_name='higher_in_gmp_l_than_mpp',
                c_mask_col='HSC_GMP_traj',
                traj_pos_min_truncate_val=-9.5,
                traj_pos_max_truncate_val=-5.5,
                x_is_cell_freq=True,
                bin_count=9,
                min_donor_cell_count=100,
            ),
            'higher_in_mpp_than_mebemp_l_across_hsc_mpp_and_mebemp_l': dict(
                traj_pos_obs_col_name='higher_in_mpp_than_mebemp_l',
                # c_mask_col='HSC_MEBEMP_traj',
                c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-6,
                
                # x_is_cell_freq=True,
                # bin_count=6,
                # min_donor_cell_count=50,
                
                x_is_cell_cum_freq=True,
                bin_count=100,
                min_donor_cell_count=10,
            ),
            'higher_in_mpp_than_mebemp_l_stratified_across_hsc_mpp_mebemp_l': dict(
                traj_pos_obs_col_name='higher_in_mebemp_l_than_mpp',
                c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
                traj_pos_min_truncate_val=-8.6,
                traj_pos_max_truncate_val=-5.8,
                
                # x_is_cell_freq=True,
                
                x_obs_col_name='higher_in_mpp_than_mebemp_l',
                norm_each_traj_bin=True,
                
                bin_count=5,
                min_bin_cell_count=5,
                min_donor_cell_count=1,
            ),
            's_phase_sig_stratified_across_hsc_mpp_mebemp_l': dict(
                traj_pos_obs_col_name='higher_in_mpp_than_mebemp_l',
                # c_mask_col='HSC_MEBEMP_ERYP_traj',
                c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-6,
                
                # x_is_cell_freq=True,
                
                x_obs_col_name='s_phase_sig',
                norm_each_traj_bin=True,
                
                bin_count=4,
                min_bin_cell_count=10,
                min_donor_cell_count=1,
            ),
            's_phase_sig_across_hsc_mpp_mebemp_l': dict(
                traj_pos_obs_col_name='s_phase_sig',
                c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
                # c_mask_col_and_vals=('c_state', ['ERYP']), # not bad.
                # c_mask_col_and_vals=('c_state', ['BEMP']), # does not look good.
                # c_mask_col='HSC_MEBEMP_ERYP_traj',
                # c_mask_col='HSC_MEBEMP_traj',
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-6.5,
                # x_is_cell_freq=True,
                x_is_cell_cum_freq=True,
                # bin_count=6,
                bin_count=100,
                min_donor_cell_count=10,
            ),
            'vim_sig_stratified_across_hsc_mpp_mebemp_l': dict(
                traj_pos_obs_col_name='higher_in_mpp_than_mebemp_l',
                # c_mask_col='HSC_MEBEMP_ERYP_traj',
                c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-6,
                
                # x_is_cell_freq=True,
                
                x_obs_col_name='mpp_mebemp_vim_sig',
                # norm_each_traj_bin=True,
                norm_each_donor=True,
                
                bin_count=4,
                min_bin_cell_count=10,
                min_donor_cell_count=1,
            ),
            'vim_sig_across_hsc_mpp': dict(
                traj_pos_obs_col_name='mpp_mebemp_vim_sig',
                # c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
                c_mask_col_and_vals=('c_state', ['HSC_MPP']),
                # c_mask_col_and_vals=('c_state', ['ERYP']), # not bad.
                # c_mask_col_and_vals=('c_state', ['BEMP']), # does not look good.
                # c_mask_col='HSC_MEBEMP_ERYP_traj',
                # c_mask_col='HSC_MEBEMP_traj',
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-6.3,
                
                x_is_cell_freq=True,
                bin_count=6,
                min_donor_cell_count=50,
                
                # x_is_cell_cum_freq=True,
                # bin_count=100,
                # min_donor_cell_count=10,
            ),
            'vim_sig_across_clp': dict(
                traj_pos_obs_col_name='mpp_mebemp_vim_sig',
                c_mask_col_and_vals=('c_state', ['CLP']),
                # c_mask_col='HSC_MEBEMP_ERYP_traj',
                # c_mask_col='HSC_MEBEMP_traj',
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-6.3,
                
                # x_is_cell_freq=True,
                # bin_count=6,
                # min_donor_cell_count=50,
                
                x_is_cell_cum_freq=True,
                bin_count=100,
                min_donor_cell_count=10,
            ),
            'id2_runx2_sig_across_clp': dict(
                traj_pos_obs_col_name='clp_id2_runx2_sig',
                c_mask_col_and_vals=('c_state', ['CLP']),
                # c_mask_col='HSC_MEBEMP_ERYP_traj',
                # c_mask_col='HSC_MEBEMP_traj',
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-6,
                
                # x_is_cell_freq=True,
                # bin_count=6,
                # min_donor_cell_count=50,
                
                x_is_cell_cum_freq=True,
                bin_count=100,
                min_donor_cell_count=10,
            ),
            'higher_in_bemp_than_mebemp_l_across_mebemp_bemp_traj': dict(
                traj_pos_obs_col_name='higher_in_bemp_than_mebemp_l',
                # c_mask_col='MEBEMP_BEMP_traj',
                c_mask_col_and_vals=('c_state', ['MEBEMP-L', 'BEMP']),
                traj_pos_min_truncate_val=-10,
                traj_pos_max_truncate_val=-5,
                
                # x_is_cell_freq=True,
                # bin_count=7,
                # min_donor_cell_count=100,
                
                x_is_cell_cum_freq=True,
                bin_count=100,
                min_donor_cell_count=10,
            ),
            'higher_in_ms4a2_high_bemp_than_other_bemp_across_bemp': dict(
                traj_pos_obs_col_name='bemp_ms4a2_sig',
                c_mask_col_and_vals=('c_state', ['BEMP']),
                # c_mask_col_and_vals=('c_state', ['MEBEMP-L', 'BEMP']),
                # c_mask_col_and_vals=('c_state', ['MEBEMP-L']),
                # c_mask_col='MEBEMP_BEMP_traj',
                # c_mask_col='mebemp_with_high_bemp_sig',
                traj_pos_min_truncate_val=-8.5, # for higher_in_ms4a2_high_bemp_than_other_bemp
                traj_pos_max_truncate_val=-5.5, # for higher_in_ms4a2_high_bemp_than_other_bemp
                # traj_pos_min_truncate_val=-8, # for higher_in_bemp_than_mebemp_l_and_eryp
                # traj_pos_max_truncate_val=-5, # for higher_in_bemp_than_mebemp_l_and_eryp
                # x_is_cell_freq=True,
                x_is_cell_cum_freq=True,
                bin_count=100,
                # x_obs_col_name='higher_in_ms4a2_low_bemp_than_other_bemp',
                # bin_count=5,
                min_donor_cell_count=10,
                # ignore_lower_than_truncate_val=True,
            ),
            'higher_in_pro_b_than_clp_across_clp_pro_b_traj': dict(
                traj_pos_obs_col_name='higher_in_pro_b_than_clp_m',
                c_mask_col='CLP_pro_B_traj',
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-5.5,
                bin_count=9,
                ignore_lower_than_truncate_val=True,
                min_donor_cell_count=50,
                
                x_is_cell_freq=True,
                # x_is_cell_freq_or_cum_freq_out_of_background=True,
                # background_c_mask_for_freq_col='not_dc_nkt_monocyte_endothel_b_lamp3',
            ),
            
            'prss2_sig_across_clp': dict(
                # traj_pos_obs_col_name='higher_in_nktdp_than_clp_m_and_pro_b',
                traj_pos_obs_col_name='prss2_sig',
                # c_mask_col='CLP_NKTDP_traj',
                c_mask_col_and_vals=('c_state', ['CLP']),

                traj_pos_min_truncate_val=-9.5,
                traj_pos_max_truncate_val=-5.5,
                
                # x_is_cell_freq=True,
                # bin_count=6,
                # min_donor_cell_count=40,
                
                x_is_cell_cum_freq=True,
                bin_count=100,
                min_donor_cell_count=10,
            ),
            # interesting only in ultima currently?
            'higher_in_nktdp_than_clp_m_and_pro_b_across_clp_nktdp_traj': dict(
                traj_pos_obs_col_name='higher_in_nktdp_than_clp_m_and_pro_b',
                # c_mask_col='CLP_NKTDP_traj',
                c_mask_col_and_vals=('c_state', ['CLP']),

                traj_pos_min_truncate_val=-9.5,
                traj_pos_max_truncate_val=-5.5,
                
                x_is_cell_freq=True,
                bin_count=6,
                min_donor_cell_count=50,
            ),
            
            # maybe good?
            'higher_in_mpp_than_mkp_across_hsc_mebemp_mkp_traj': dict(
                traj_pos_obs_col_name='higher_in_mpp_than_mkp',
                c_mask_col='MEBEMP_MKP_traj',
                traj_pos_min_truncate_val=-9.5,
                traj_pos_max_truncate_val=-5.5,
                x_is_cell_freq=True,
                bin_count=7,
                min_donor_cell_count=100,
            ),
            'higher_in_mkp_than_mebemp_l_across_mebemp_mkp_traj': dict(
                traj_pos_obs_col_name='higher_in_mkp_than_mebemp_l',
                c_mask_col='MEBEMP_MKP_traj',
                traj_pos_min_truncate_val=-8,
                traj_pos_max_truncate_val=-5,
                x_is_cell_freq=True,
                bin_count=6,
                min_donor_cell_count=100,
            ),
            'higher_in_mpp_than_ep_across_eryp': dict(
                traj_pos_obs_col_name='higher_in_mpp_than_ep',
                c_mask_col_and_vals=('c_state', ['ERYP', 'high_MPP_sig_ERYP']),
                traj_pos_min_truncate_val=-9,
                traj_pos_max_truncate_val=-5.5,
                x_is_cell_freq=True,
                bin_count=8,
                min_donor_cell_count=100,
            ),
            # is higher_in_ep_than_mebemp_l_across_mebemp_eryp_traj batchy??
            'higher_in_ep_than_mebemp_l_across_mebemp_eryp_traj': dict(
                traj_pos_obs_col_name='higher_in_ep_than_mebemp_l',
                c_mask_col='MEBEMP_ERYP_traj',
                traj_pos_min_truncate_val=-10,
                traj_pos_max_truncate_val=-6,
                x_is_cell_freq=True,
                bin_count=6,
                min_donor_cell_count=100,
            ),
            
            # nothing interesting?
            # nothing here...

            # not enough donors with enough cells
            'higher_in_pre_b_than_pro_b_across_pro_b_pre_b_traj': dict(
                traj_pos_obs_col_name='higher_in_pre_b_than_pro_b',
                c_mask_col='pro_B_pre_B_traj',
                traj_pos_min_truncate_val=-9.5,
                traj_pos_max_truncate_val=-4.5,
                x_is_cell_freq=True,
                bin_count=6,
                min_donor_cell_count=20,
            ),

            # similar to one of the above, but seems worse
            # 'higher_in_mpp_than_ep_across_s_phase_sig_across_mebemp_l': dict(
            #     # traj_pos_obs_col_name='s_phase_sig',
            #     traj_pos_obs_col_name='s_phase_sig',
            #     c_mask_col_and_vals=('c_state', ['MEBEMP-L']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-6.5,
            #     # x_is_cell_freq=True,
            #     x_obs_col_name='higher_in_mpp_than_ep',
            #     # norm_each_traj_bin=True,
            #     # min_bin_cell_count=1,
            #     bin_count=6,
            #     min_donor_cell_count=50,
            # ),
            # 'higher_in_mpp_than_ep_across_mebemp_l': dict(
            #     traj_pos_obs_col_name='higher_in_mpp_than_ep',
            #     c_mask_col_and_vals=('c_state', ['MEBEMP-L']),
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-5.5,
            #     x_is_cell_freq=True,
            #     bin_count=6,
            #     min_donor_cell_count=50,
            # ),
            # 'higher_in_mebemp_l_than_mpp_across_hsc_mpp_and_mebemp_l': dict( # meh. higher_in_mpp_than_mebemp_l_across_hsc_mpp_and_mebemp_l is better, i think.
            #     traj_pos_obs_col_name='higher_in_mebemp_l_than_mpp',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
            #     traj_pos_min_truncate_val=-8.8,
            #     traj_pos_max_truncate_val=-5.5,
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),
            # 'higher_in_mpp_than_ep_across_hsc_mebemp_eryp_traj': dict(
            #     traj_pos_obs_col_name='higher_in_mpp_than_ep',
            #     c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-5.5,
            #     x_is_cell_freq=True,
            #     bin_count=8,
            #     min_donor_cell_count=100,
            # ),
            # 'higher_in_mpp_than_ep_across_mebemp_eryp': dict(
            #     traj_pos_obs_col_name='higher_in_mpp_than_ep',
            #     c_mask_col_and_vals=('c_state', ['MEBEMP-L', 'ERYP', 'high_MPP_sig_ERYP']),
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-5.5,
            #     x_is_cell_freq=True,
            #     bin_count=8,
            #     min_donor_cell_count=100,
            # ),
            # 'pcna_sig_across_mebemp_bemp_traj': dict(
            #     traj_pos_obs_col_name='higher_in_bemp_than_mebemp_l',
            #     c_mask_col='MEBEMP_BEMP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-5,
            #     x_obs_col_name='pcna_sig',
            #     bin_count=5,
            #     norm_each_traj_bin=True,
            #     min_bin_cell_count=1,
            #     min_donor_cell_count=50,
            # ),
            # 'pcna_sig_across_mebemp_eryp_traj': dict(
            #     traj_pos_obs_col_name='higher_in_ep_than_mebemp_l',
            #     c_mask_col='MEBEMP_ERYP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-5.5,
            #     x_obs_col_name='pcna_sig',
            #     bin_count=5,
            #     norm_each_traj_bin=True,
            #     min_bin_cell_count=1,
            #     min_donor_cell_count=50,
            # ),
            # 's_phase_sig_across_mebemp_l': dict( # s_phase_sig_across_hsc_mpp_mebemp_l is better, i think.
            #     traj_pos_obs_col_name='s_phase_sig',
            #     c_mask_col_and_vals=('c_state', ['MEBEMP-L']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-6.5,
            #     # x_is_cell_freq=True,
            #     x_is_cell_cum_freq=True,
            #     # bin_count=6,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),
            # 's_phase_sig_across_hsc_mpp': dict( # s_phase_sig_across_hsc_mpp_mebemp_l is better, i think.
            #     traj_pos_obs_col_name='s_phase_sig',
            #     c_mask_col_and_vals=('c_state', ['HSC_MPP']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-6.5,
            #     # x_is_cell_freq=True,
            #     x_is_cell_cum_freq=True,
            #     # bin_count=6,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),
            # 'pcna_sig_across_hsc_gmp_traj': dict( # meh, i think
            #     traj_pos_obs_col_name='higher_in_gmp_l_than_mpp',
            #     c_mask_col='HSC_GMP_traj',
            #     traj_pos_min_truncate_val=-9.5,
            #     traj_pos_max_truncate_val=-5.5,
            #     x_obs_col_name='pcna_sig',
            #     bin_count=5,
            #     norm_each_traj_bin=True,
            #     min_bin_cell_count=1,
            #     min_donor_cell_count=50,
            # ),
            # 'pcna_sig_across_hsc_clp_traj': dict( # meh, i think
            #     traj_pos_obs_col_name='higher_in_clp_m_than_all_myeloid_hspcs',
            #     c_mask_col='HSC_CLP_traj',
            #     traj_pos_min_truncate_val=-9.5,
            #     traj_pos_max_truncate_val=-5.5,
            #     x_obs_col_name='pcna_sig',
            #     bin_count=6,
            #     norm_each_traj_bin=True,
            #     min_bin_cell_count=2,
            #     min_donor_cell_count=50,
            # ),
            # 'pcna_sig_across_eryp': dict(
            #     traj_pos_obs_col_name='pcna_sig',
            #     c_mask_col_and_vals=('c_state', ['ERYP']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-6,
            #     x_is_cell_freq=True,
            #     bin_count=6,
            #     min_donor_cell_count=50,
            # ),
            # 'higher_in_mpp_than_mebemp_l_across_mebemp_l': dict( # 240415: looked worse than higher_in_mpp_than_mebemp_l_across_hsc_mpp_and_mebemp_l
            #     traj_pos_obs_col_name='higher_in_mpp_than_mebemp_l',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     c_mask_col_and_vals=('c_state', ['MEBEMP-L']),
            #     traj_pos_min_truncate_val=-8.5,
            #     traj_pos_max_truncate_val=-6,
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),
            # 'higher_in_mebemp_l_than_mpp_across_hsc_mpp': dict( # seems like higher_in_mebemp_l_than_mpp_across_hsc_mpp_and_mebemp_l is better
            #     traj_pos_obs_col_name='higher_in_mebemp_l_than_mpp',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     c_mask_col_and_vals=('c_state', ['HSC_MPP']),
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-6.5, # if only looking at HSC_MPP, can't reach lower than thresh that defines HSC_MPP (-6.8 on 240415)...
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),
            # 'higher_in_mebemp_l_than_mpp_across_mebemp_l': dict( # seems like higher_in_mebemp_l_than_mpp_across_hsc_mpp_and_mebemp_l is better
            #     traj_pos_obs_col_name='higher_in_mebemp_l_than_mpp',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     c_mask_col_and_vals=('c_state', ['MEBEMP-L']),
            #     traj_pos_min_truncate_val=-6.8,
            #     traj_pos_max_truncate_val=-5.5,
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),
            # 'higher_in_mpp_than_mebemp_l_across_hsc_mpp': dict( # 240415: looked worse than higher_in_mpp_than_mebemp_l_across_mebemp_l
            #     traj_pos_obs_col_name='higher_in_mpp_than_mebemp_l',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     c_mask_col_and_vals=('c_state', ['HSC_MPP', 'MEBEMP-L']),
            #     traj_pos_min_truncate_val=-8.5,
            #     traj_pos_max_truncate_val=-6,
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),
            # 'higher_in_ms4a2_high_bemp_than_other_bemp_across_bemp': dict(
            #     traj_pos_obs_col_name='bemp_ms4a2_sig',
            #     # traj_pos_obs_col_name='higher_in_bemp_than_mebemp_l_and_eryp',
            #     c_mask_col_and_vals=('c_state', ['BEMP']),
            #     # c_mask_col='MEBEMP_BEMP_traj',
            #     # c_mask_col='mebemp_with_high_bemp_sig',
            #     traj_pos_min_truncate_val=-8.5, # for higher_in_ms4a2_high_bemp_than_other_bemp
            #     traj_pos_max_truncate_val=-5.5, # for higher_in_ms4a2_high_bemp_than_other_bemp
            #     # traj_pos_min_truncate_val=-8, # for higher_in_bemp_than_mebemp_l_and_eryp
            #     # traj_pos_max_truncate_val=-5, # for higher_in_bemp_than_mebemp_l_and_eryp
            #     # x_is_cell_freq=True,
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     # x_obs_col_name='higher_in_ms4a2_low_bemp_than_other_bemp',
            #     # bin_count=5,
            #     min_donor_cell_count=10,
            #     # ignore_lower_than_truncate_val=True,
            # ),
            # 'higher_in_ep_than_mpp_across_hsc_mebemp_eryp_traj': dict(
            #     traj_pos_obs_col_name='higher_in_mpp_than_ep',
            #     c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-5.5,
            #     x_obs_col_name='higher_in_ep_than_mpp',
            #     bin_count=6,
            #     norm_each_traj_bin=True,
            #     min_bin_cell_count=3,
            #     min_donor_cell_count=50,
            # ),

            
            # meh
            'nktdp_maf_sig_across_nktdp': dict(
                traj_pos_obs_col_name='nktdp_maf_sig',
                c_mask_col_and_vals=('c_state', ['NKTDP']),
                # c_mask_col='HSC_MEBEMP_ERYP_traj',
                # c_mask_col='HSC_MEBEMP_traj',
                traj_pos_min_truncate_val=-10,
                traj_pos_max_truncate_val=-7.2,
                # x_obs_col_agg_func=lambda x: np.quantile(x, 0.9),
                
                # x_is_cell_freq=True,
                # bin_count=6,
                # min_donor_cell_count=50,
                
                x_is_cell_cum_freq=True,
                bin_count=100,
                min_donor_cell_count=10,
            ),            
            # 'nktdp_fam30a_sig_across_nktdp': dict(
            #     traj_pos_obs_col_name='nktdp_fam30a_sig',
            #     c_mask_col_and_vals=('c_state', ['NKTDP']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-11,
            #     traj_pos_max_truncate_val=-7,
            #     # x_obs_col_agg_func=lambda x: np.quantile(x, 0.9),
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),            
            # 'nktdp_lmna_sig_across_nktdp': dict(
            #     traj_pos_obs_col_name='nktdp_lmna_sig',
            #     c_mask_col_and_vals=('c_state', ['NKTDP']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-11,
            #     traj_pos_max_truncate_val=-7,
            #     # x_obs_col_agg_func=lambda x: np.quantile(x, 0.9),
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),            
            # 'nktdp_samhd1_sig_across_nktdp': dict(
            #     traj_pos_obs_col_name='nktdp_samhd1_sig',
            #     c_mask_col_and_vals=('c_state', ['NKTDP']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-9.2,
            #     traj_pos_max_truncate_val=-5.2,
            #     # x_obs_col_agg_func=lambda x: np.quantile(x, 0.9),
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),            
            # 'nktdp_nfil3_sig_across_nktdp': dict(
            #     traj_pos_obs_col_name='nktdp_nfil3_sig',
            #     c_mask_col_and_vals=('c_state', ['NKTDP']),
            #     # c_mask_col='HSC_MEBEMP_ERYP_traj',
            #     # c_mask_col='HSC_MEBEMP_traj',
            #     traj_pos_min_truncate_val=-9,
            #     traj_pos_max_truncate_val=-6.5,
            #     # x_obs_col_agg_func=lambda x: np.quantile(x, 0.9),
                
            #     # x_is_cell_freq=True,
            #     # bin_count=6,
            #     # min_donor_cell_count=50,
                
            #     x_is_cell_cum_freq=True,
            #     bin_count=100,
            #     min_donor_cell_count=10,
            # ),            
            
            
            # exploring
        },
    },
    'ambient_noise_removal': {
        'mc_model': {
            'target_metacell_umis': int(170e3),
            'target_metacell_size': 80,
        },
        'exps_to_exclude_from_mcnoise_and_continue_with_them_as_if_they_have_no_ambient_noise': [
            'N328_15_12_22_1', # 230531: MCNoise estimated ~70% noise for each of the 3 depth bins, so I guess the assumption that the native expression is the same did not hold.
        ],
        
        'projected_and_assigned_types_to_ignore': [
            'N186-MDS',
            'N205-MDS',
            'N191-AML',
            # 'N224-high-ribo',

            'Doublet',
            'Mixture',
            'Dissimilar',
        ],
        'min_max_donor_id_fraction_to_ignore_metacell': 0.95,
        'min_max_donor_id_log_enrichment_to_ignore_metacell': 5,
        'min_max_exp_name_fraction_to_ignore_metacell': 0.95,
        'min_max_exp_name_log_enrichment_to_ignore_metacell': 3.5,
        'extra_genes_to_forbid': [ # unused if we use manually chosen clusters
            'ANXA2',
            'LMO4', # potentially a nice BEMP. but it isn't low in any MC (after filtering AML etc)...
            'ALOX5AP', # potentially a nice BEMP. but also B and monocytes, and not as strong as MS4A2 and HDC, so i guess we can give up on it?
            'MZB1',
            'FHL1', 'HEMGN', 'PCDH9',
            'HPGD',
            'LTB', 'IGHM',
            'HOPX', 'SPINK2',
            'FAM30A',
            'CSF3R',
            'TYMP', 'ITGB2',
            'TYROBP',
            'CST3', 'SAMHD1',
            'GRN',
            'CCR7',
            'ITGAL',
            'TNFAIP2', 'IGHD', 'DMXL2',
            'ZFP36L1',
            'LYZ', # so sad.
            'JAML', 'PSAP', 'PTPRE',
            'CD36',
            'S100A10',
            'S100A11',
            'S100A6',
            'FCER1G',
            'FYB1',
            'F13A1',
            'CPVL',
            'AOAH',

            'LGALS1',

            'CD52',
            'CD48',
            'CCDC50',
            'BANK1',
            'TMEM154',
            'MARCKS', 'STK17A',
            'CTSS',
            'RALGPS2', 'LBH', 'AFF3', 'COBLL1', 'TTN',
            'HLA-DQA1', 'HLA-DQB1', 'HLA-DMB',
            'BIRC3',
            'ETS1', 'FAM107B',
            'COTL1',
            'P2RX5',
            'EVI2B',
            'CD79B', 'CTSZ', 'KLF2',
            'POU2F2',

            'GALNTL6',
        ],
        'use_anyway_gene_names': [ # unused if we use manually chosen clusters
            # GMP-L
            'CFD', 'MPO',

            # B
            'BLK', 'PAX5',
            'CD79A', 'CD79B',
            'MS4A1', 'IGKC', 'BANK1',

            # Monocytes
            'S100A9', 'S100A8',
            'LYZ', # seems to not be silent for most types

            'MARCH1',
            'SCIMP', 'CTSH', 'CYBB',

            # EP
            'HBB',
            'HBD',

            # HSPC
            'CD34', 'EGFL7', 'NPR3',
        ],
        'gene_clusters': [
            [ # CLP (doesnt look very good for MCNoise, but whatever)
                'SLC2A5', 'ID2', 'MME', 'CLNK', 'PALLD', 'KCNK17', 'KIAA0087', 
                
                'ACY3', 'DDIT4', 'DNTT', 'ITGAL', 'CYTH4',
                'MAF', 'JCHAIN',
            ],

            [ # EP
                'HBB',
                'HBD',
                'APOC1',
            ],

            [ # GMP
                'CLEC12A', 
                'TNFSF13B', 
                'MGST1',
                'RNASE2', 'ELANE',
                'CFD', # maybe too noisy?
                'MPO',
                'CTSG',
            ],

            [ # monocytes and B
                'TBC1D9', 'MARCH1',
                'LY86',
                'NCF1',
                'IL10RA',
                'CTSH',
                'SCIMP',
            # ],
            # [ # B
                'ID3',
                'LINC01781',
                'FCRL5', 'FCRL3', 'FCRL2',
                'AIM2',
                'FCRLA',
                'LINC01857',
                'SP140', 'SH3BP5', 'OSBPL10',
                'TLR10', 'ARHGAP24',
                'CPNE5', 'CD24',
                'PAX5',
                'MS4A1',
                'CD27', 'CLECL1',
                'HVCN1', 'GPR183',
                'LINC00926',
                'IRF8', # though maybe too weak?
                'TNFRSF13B', 'IKZF3',
                'FAM129C', 'CD22', 'CD79A', 'SPIB',
                'VPREB3', 'TNFRSF13C',
            # ],
            # [ # monocytes
                'TNFRSF1B',
                'FGR',
                'S100A9', 'S100A12', 'S100A8',
                'MNDA',
                'FCGR2A',
                'NCF2',
                'SLC11A1',
                'CSTA', 'P2RY13',
                'VCAN', 'CD14',
                'FGL2',
                'GIMAP4', 'CYBB', 'CFP', 'CEBPD',
                'FCN1', 'MPEG1', 'MS4A6A', 'MS4A7',
                'CLEC7A', 'LRRK2',
                'KCTD12',
                'SLC7A7', 'LGALS3',
                'SERPINA1',
                'IGSF6',
                'CD68',
                'CD300E', 'RAB31', 'CD93',
                'FPR1', 'LILRB3', 'LILRB2', 'LILRB1', 'IL17RA', 'ADA2',
                'LGALS2',
                'RGS2', 'TGFBI',
                'MAFB',
            ],
            [ # BEMP
                'MS4A3', 'MS4A2', 'HDC',
            ],
            [ # not B
                'ANXA1',
                'ATP8B4',
                'GIMAP7', 'DUSP6',
            ],
            [ # not B and not monocytes (i.e., HSPC??)
                'COL24A1', 'TSPAN2', 'CD34', 'GUCY1A1', 'NPR3', 'LAPTM4B', 'BAALC', 'ANGPT1',
                'MDK', 'ZBTB16',
                'SMIM24', 'LINC02573',
                'CRHBP', 'ADGRG6',
                'MEG3', 'AC011139.1',
                'CPA3',
                'TFPI',
                # 'ANKRD28', # maybe not low enough in B and monocytes...

                # not very clean, but might be good
                'FCER1A', 'CNRIP1', 'HPGDS', 'GATA1', 'AL157895.1', 'MINPP1', 'KLF1',
                # 'CSF2RB',
                'PBX1',
                # 'IL1B',
                # 'SLC40A1',
                'ZNF385D', 'GATA2', 'CYTL1', 'KIT', 'EREG',
                #  'SOX4',
                'MYB',
                #  'CDK6',
                'AC002454.1',
                # 'NUCB2',
                # 'AC084033.3',
            ],
        ],
        'umi_depth_number_of_bins': 3,
        # 'number_of_metacells_clusters': 10,
        # 'number_of_genes_clusters': 10,
        'max_relative_expression_for_metacells_genes_to_use': -3,
        'max_relative_expression_for_metacells_genes_to_use_step': 4,
    },
    'gene_group_name_to_color': {
        'platelet_related': 'purple',
        'cell_cycle_related': 'black',
        'misc': 'grey',
        'presumably_stress_or_apoptosis_related': 'orange',
        'relatively_noisy': 'blue',
        'highly_polymorphic_genes': 'red',
        'ribosomal': 'green',
        'uncharacterized': 'yellow',
    },
    'metacell_attrs_to_copy_to_cells_ad': [
        'projected_type',
        'state',
        'state_color',
        # 'similar',
    ],
    'nkt_hba_contamination': {
        'strong_nkt_specific_genes': ['IL32', 'GZMA', 'CD3D', 'CD3G', 'KLRD1', 'GZMK'],
        'min_nkt_specific_expr': -8,
        'min_hba_expr': -10,
    },
    'platelet_and_mk': {
        'characterizing_platelets': {
            'axiomatic_platelet_marker_gene_names': ['PF4', 'PPBP'], # must stay only PF4 and PPBP here to avoid circularity. after we characterize platelets, we can use PROBABLY_PLATELET_SPECIFIC_GENES everywhere. 
            'platelet_min_total_axiomatic_platelet_marker_downsampled_umis': 15, # TODO: add num of UMIs in downsample. 790?
        },
        
        'strong_platelet_specific_marker_gene_names': ['PPBP'], 
        'mkp_or_platelet_contaminated_min_platelet_specific_norm_expression': 2e-3, # 231109: this is based on the observation that in http://tanaywiz.wisdom2k3.weizmann.ac.il/hg_hematopoiesis/hsc_stimulation/combined_1_2_3_4_5/ and http://tanaywiz.wisdom2k3.weizmann.ac.il/hg_hematopoiesis//GSE192519_hsc_expansion_SolPlus/ PPBP is virtually zero in mostly everything.
        
        # 'min_num_of_non_excluded_umis_for_exemption_according_to_mkp_markers': 1000,
        # 'mkp_marker_gene_names': GENES_HIGHER_IN_MKP_COMPARED_TO_PLATELTS_AND_ALSO_EP_AND_MEBEMP_L,
        # 'min_mkp_downsampled_mkp_specific_umi_count': 3,
        # 'min_mkp_specific_norm_expression': 3 / 1000,

        # 'platelet_min_downsampled_platelet_specific_umi_count': 20, # the expression
        'platelet_min_platelet_specific_norm_expression': 2**(-6),
    },
    'min_num_of_gene_reads_to_consider_in_read2_gc_calc': 10,
    'ultima_illumina_technical_repeats': {
        'gene_ratio_df_csv_file_path': '/dummy/dummy/dummy/raid/ultima_vs_illumina_in_the_6_technical_repeat_experiments/gene_ratio_df.csv',
        'gene_read_gc_df_csv_file_path': os.path.join(mds_in_out_dir_paths.GC_BIAS_ANALYSIS_DIR_PATH, 'gene_read_gc_df.csv'),
        'gene_read_per_umi_df_csv_file_path': '/dummy/dummy/dummy/raid/ultima_vs_illumina_in_the_6_technical_repeat_experiments/gc_bias_analysis/all_gene_read_per_umi_df.csv',
        'log_pseudocount_for_log': -18,
        'min_log_ratio_range_to_forbid_from_features': 1,
        'min_abs_median_log_ratio_to_forbid_from_features': 1,
    },
    'reads_per_umi': {
        'min_num_of_reads_per_gene': 1000,
    },
    'batch_effects': {
        'max_exp_mc_binom_pval_to_consider': 1e-2,
        'min_exp_mc_min_donor_binom_pval': 0.05,
        'min_num_of_cells_of_exp_in_metacell_for_suspicious_exp_mc_pairs': 30,


        'max_exp_ks_pval_to_also_perform_ks_ignoring_each_donor': 1e-20,
        'min_num_of_metacells_with_any_cell_of_exp_or_donor_to_perform_ks': 10,
        'min_num_of_cells_per_group_to_perform_ks': 150,
        'max_gene_exp_cell_type_max_ks_pval_to_consider': 1e-35,
        'max_gene_exp_cell_type_max_ks_pval_to_mark_as_lateral': 1e-50,
        'max_gene_exp_cell_type_min_median_abs_diff_to_mark_as_lateral': 0.5,
        'max_exp_ks_pval_to_also_perform_ks_ignoring_each_donor': 1e-20,

        # 'min_num_of_metacells_with'
    },

    # NOTE: use mds_analysis_params.NON_STEM_CELL_STATE_NAMES instead.
    # 'differentiated_cell_types': [
    #     'B',
    #     'cMonocyte',
    #     'DC',
    #     'NKT',
    # ],

    'lmna_sig': {
        # 
        'gene_names': ['LMNA', 'AHNAK', 'MYADM', 'TSPAN2', 'ANXA1', 'ANXA2'],
        'cell_types': [
            'MPP', 'MEBEMP-M', 'MEBEMP-L',
            # 'GMP-E',
        ],
        'metacell_max_donor_id_fraction': 0.5,
        'num_of_bins': 10,
        'genes_to_bin_by': 'AVP',
        'min_num_of_cells_to_calculate_score': 50,
        'geomean_across_pooled_cells': True, 
        
    },


    'karyotype_estimation': {
        'ignore_problematic_ultima_illumina_genes': True,
        'ordered_chrom_names_for_karyotype_estimation': [
            *(str(x) for x in range(1,23)), 
            'X', 
            # 'Y', # 230402: nothing left here if not considering sex diff expressed genes
        ],

        'seemingly_batchy_cna_infos': [
            {
                'name': 'region_containing_some_HIST_and_HLA',
                'chrom': '6',
                'first_gene': 'HIST1H1C',
                'last_gene': 'RING1',
            },
        ],
        
        'low_log_norm_expression_thresh': -15.3, # 230705: -15 was the value for a long time. now seemed like we are left with a bit small number of genes, so i tried a lower number.
        
        # 'min_total_num_of_donors': 5,
        'min_total_num_of_donors_for_extra_ignored_genes': 20,
        
        # 'min_total_num_of_donors_for_ignoring_genes_by_low_projected_correlation': 50,
        # 'min_gene_projected_correlation': 0.1, # seems like there are very few genes with low projected correlation that aren't ignored by other criteria.

        ### strict
        # 'max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression': 0.01, # 230330: compared to 0.9, using 0.01 here seemed much better.
        # 'max_max_abs_log_norm_expression_ratio_20_80_percentiles': 0.5,

        ### relaxed
        # 'max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression': 0.98, # i think 0.98 was too high for wide scanning, e.g. N199, N186
        'max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression': 0.5, # 0.5 seemed the best for N200 del_5_except_5p_start
        'max_max_abs_log_norm_expression_ratio_20_80_percentiles': 1,
        
        # anyway we take stratified to bins, i think, because this inhibits the bias (e.g., due to GC).
        'projected_fold_type_for_cna_estimation': 'geomean_of_projected_folds', # geomean implicitly by taking the mean of logs
        # 'projected_fold_type_for_cna_estimation': 'projected_fold_of_mean', # i think here zeros would have a lower weight (as usual in arithmetic mean).

        'use_median_for_gc_bin_correction': True, # alternative is to use mean, which seems to work worse
        # 'max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression': 0.01,
        # 'max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression': 0.98, 

        'epsilon_for_log_projected_fold': 1e-5,
        # 'epsilon_for_log_projected_fold': 5e-5,
        'num_of_gc_bins': 10, # we don't want too many bins, because there might be a lot of CNAs in some donors. 10 is reasonable, i guess?
        

        # 'num_of_gene_stratification_bins': {
        #     # TODO: go back to same params for cells and metacells, and also no bins for expr level...
        #     'cells': {
        #         # TODO: completely ignore cells with any empty bin? (for which we set np.nan instead of projected fold)
        #         'mean_illumina_mean_read2_gc_fraction': 3,
        #         'median_log_norm_expression': 3, 
        #     },
        #     'metacells': {
        #         # 'mean_illumina_mean_read2_gc_fraction': 10, # seems like a value higher than 10 here is not making a lot of difference
        #         'mean_illumina_mean_read2_gc_fraction': 5, # seems like there isn't much difference between 5 and 10, though maybe 10 is a bit better. still, for uniformity when looking at cells and metacells, i think i prefer to use 5.
        #         'median_log_norm_expression': 5, # this helps also with N223, for example.
        #     }
        # },
        
        # 'wanted_num_of_genes_in_bin': 10,
        'wanted_num_of_genes_in_bin': 20,
        
        'moving_window_size': 21,

        
        'default_min_num_of_cells_in_cluster_to_show_in_cluster_moving_median_plot': 100,
        'default_min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 4,
        
        'default_min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist': 4,


        # NOTE: normalize_each_mc_projected_fold_row_by_median does not affect the clustering by pearson, because pearson does not care about adding a constant to each vector.
        'normalize_each_mc_projected_fold_row_by_median': True,
        # 'normalize_each_mc_projected_fold_row_by_median': False,
        
        'normalize_each_cell_projected_fold_row_by_mean': True,
        # 'normalize_each_cell_projected_fold_row_by_mean': False,
        
        'clustering_cells_or_metacells': {
            'linkage_method': 'average', # is this appropriate? tools_by_others/MCNoise-main/MCNoise/AmbientNoiseFinder.py uses it. and some page on the internet said 'average' and 'complete' are the two most popular ones for hierarchical clustering.
            
            'linkage_metric': 'correlation',
            # 'linkage_metric': 'euclidean',


            'default_num_of_mcs_per_cluster_for_identifying_clones': 10,
            'default_num_of_mc_cna_clusters': 10,
            'min_num_of_clusters_for_identifying_clones': 2,
            # 'normalized_mc_projected_fold_truncation_range_for_clustering': (-0, 0), # tried different values hoping for N200 to become sane, until i got an error for 0.4, i guess due to an all-zero row (though the error was: "ValueError: The condensed distance matrix must contain only finite values.").

            'default_mc_cna_clusters_threshold_dist_quantile': 0.8,
            'default_mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.8,
        },
        
        # 'very_low_log_norm_expression_thresh': -16,
        # 'max_abs_log_norm_expression_ratio': 1,
        # 'atlas_to_use': 'mds_healthy_controls_atlas',
        'atlas_to_use': 'nimrod_oren_atlas',

        'donor_id_to_clone_and_cna_info': {
            # useful: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38
            # useful: https://www.nature.com/articles/s41591-021-01411-9/figures/9


            # 95's CNAs were the strongest (and only) CNAs i found in the healthy donors in our original atlas (if we exclude N48 as unhealthy). So why don't we see more CNAs, like Ogawa (Combined landscape of single-nucleotide variants and copy number alterations in clonal hematopoiesis | Nature Medicine (https://www.nature.com/articles/s41591-021-01411-9, 2021))? Ogawa identified CNAs in 2032 out of 10562 (EDF 1), so ~19%. The cohort age distribution was 45.7% 60-69yo, 40.6% 70-79yo, 12.8% 80-89yo, 0.9% 90-101yo (supp table 1), so older than our original atlas (but not much older, i think). it seems that the CNA VAFs he reported are mostly very low (supp fig S1C, most common VAF for CNA in their data is ~0.006), and looking at EDF 3, it seems that most of the CNAs Ogawa identified are copy number neutral LOH events (which we don't identify currently), as well as many CNAs of "unclassifiable" type. Even if the "unclassifiable" are not part of the 19%, i guess that given the mostly low VAFs that Ogawa identified, our results are not that surprising. If we include N48, then we have identified CNAs in both N48 and 95, i.e., 2%. There were also multiple donors with a single metacell that seemed to have abnormal expression, but i just ignore these (and i wouldn't guess these are CNAs).
            # Ogawa sampled blood without CD34+ enrichment ("Our analyses of CNAs are based on results from a previous publication21" (citing Chromosomal alterations among age-related haematopoietic clones in Japan | Nature (https://www.nature.com/articles/s41586-020-2426-2, 2020), which says: "DNA was obtained from blood samples")
                # In each of N48 and 95, we identified CNAs only in B cells, i think. So maybe many of the CNAs that Ogawa identified are in B cells?
                # I guess less likely, maybe HSPCs with CNAs have an expanded progeny relative to wild-type HSPCs?



            'N212': {
                'num_of_mc_clusters_for_identifying_clones': 7,
                # 'consider_small_clusters_combined': True,

                'cna_clustermap_width': 3,
                'cna_clustermap_height': 8,
                'cna_clustermap_colors_ratio': 0.12,


                'cna_cluster_linkage_metric': 'euclidean',
                'num_of_mc_cna_clusters': 3,

                'cna_infos': [
                    {
                        # feels borderline, as the metacell hist has a peak around -0.65 rather than -1.
                        'name': 'N212_del_5q_middle',
                        'chrom': '5',
                        # 'first_gene': 'WDR36',
                        # 'first_gene': 'SRP19',
                        # 'first_gene': 'MAN2A1',
                        'first_gene': 'C5orf15',
                        # 'last_gene': 'DIAPH1',
                        # 'last_gene': 'NDFIP1',
                        'last_gene': 'CTNNA1',
                        # 'last_gene': 'RNF14',
                        # 'SRP19', 'DCP2', 'PGGT1B', 'ATG12', 'AP3S1', 'DMXL1', 'TNFAIP8', 'HSD17B4', 'SRFBP1', 'SNX2', 'CSNK1G3', 'PHAX', 'PRRC1', 'HINT1', 'CDC42SE2', 'FNIP1', 'IRF1', 'UQCRQ', 'AFF4', 'ZCCHC10', 'HSPA4', 'VDAC1', 'SKP1', 'PPP2CA', 'UBE2B', 'SAR1B', 'DDX46', 'C5orf24', 'TXNDC15', 'SMAD5', 'HNRNPA0', 'FAM13B', 'BRD8', 'KDM3B', 'ETF1', 'HSPA9', 'PAIP2', 'UBE2D2', 'CXXC5', 'PURA', 'PFDN1', 'NDUFA2', 'IK', 'HARS', 'ZMAT2', 'TAF7', 'DIAPH1', 'DELE1', 'NDFIP1'

                        'max_median_cluster_median_projected_fold_threshold': -0.48, 
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.3, 

                        'mc_median_normalized_projected_fold_hist_bins': 15,
                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.2,
                        # 'metacell_median_projected_fold_threshold': -0.3,
                        # 'cell_mean_projected_fold_diff_threshold': -1, # awful separation.
                    },
                    {
                        'name': 'N212_del_20q_middle',
                        'chrom': '20',
                        # 'first_gene': 'CHD6',
                        'first_gene': 'SRSF6',
                        # 'last_gene': 'DPM1',
                        'last_gene': 'NFATC2',
                        
                        'max_median_cluster_median_projected_fold_threshold': -0.4, 
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.25, 
                        'mc_median_normalized_projected_fold_hist_bins': 15,

                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.3,
                        # 'metacell_median_projected_fold_threshold': -0.4,
                        # 'cell_mean_projected_fold_diff_threshold': -1.25, # awful separation.
                    },
                    # {
                    #     # NOTE: very borderline: pretty small region, only a few metacells, and doesn't seem like a deletion hotspot in https://www.nature.com/articles/s41591-021-01411-9/figures/9. so commenting out the threshold (and not assigning a clone for it). also, IIRC, it seemed batchy in this experiment or another experiment.
                    #     'skip_cna': True,

                    #     'name': 'N212_del_17q_near_end',
                    #     'chrom': '17',
                    #     'first_gene': 'NPLOC4',
                    #     'last_gene': 'DUS1L',
                        
                    #     # 'max_median_cluster_median_projected_fold_threshold': -0.5, 

                    #     # 'expected_projected_fold': -1,
                    #     # 'metacell_median_projected_fold_diff_threshold': -0.3,
                    #     # 'metacell_median_projected_fold_threshold': -0.4,
                    #     # 'cell_mean_projected_fold_diff_threshold': -1.25, # awful separation.
                    # },
                ],
                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N212_1': [
                #         ('N212_del_5q_middle', 'partially'),
                #         ('N212_del_20q_middle', 'partially'), 
                #     ],
                #     'N212_2': [
                #         ('N212_del_5q_middle', 'yes'),
                #         ('N212_del_20q_middle', 'yes'), 
                #     ],
                # },
            },
            'N211': {
                'num_of_mc_clusters_for_identifying_clones': 5,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 3,

                # 'consider_small_clusters_combined': True,

                'single_cna_median_projected_fold_bins': (-np.inf, -0.45, np.inf),
                'cna_clustermap_width': 3,
                'cna_clustermap_height': 10,
                'cna_clustermap_colors_ratio': 0.15,

                'cna_infos': [
                    {
                        'name': 'N211_del_5q_middle',
                        'chrom': '5',
                        'first_gene': 'MCTP1',
                        'last_gene': 'UBLCP1',

                        'max_median_cluster_median_projected_fold_threshold': -0.45, 

                        'exp_names': ['demux_16_01_22_1'],
                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.5,
                        # 'metacell_median_projected_fold_threshold': -0.5,
                        # 'cell_mean_projected_fold_diff_threshold': -0.8, # awful separation.
                    },
                ],
                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N211_0': [
                #         ('N211_del_5q_middle', 'no'),
                #     ],
                #     'N211_1': [
                #         ('N211_del_5q_middle', 'yes'),
                #     ],
                # },
            },
            'N187': {
                'num_of_mc_clusters_for_identifying_clones': 2,
                
                'cna_clustermap_width': 3,
                'cna_clustermap_height': 8,
                'cna_clustermap_colors_ratio': 0.12,

                'cna_cluster_linkage_metric': 'euclidean',
                'num_of_mc_cna_clusters': 2,

                'cna_infos': [
                    {
                        'name': 'N187_del_5q_middle',
                        'chrom': '5',
                        # TODO: make sure the first and last genes are good. maybe should not do this manually but implement an HMM for all...
                        'first_gene': 'EPB41L4A-AS1',
                        'last_gene': 'MRPL22',
                        # NOTE: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4653765/ (with regard to lenalidomide mechanism of action): "In del(5q) MDS, the ubiquitination and degradation of casein kinase 1A1 causes the death of del(5q) cells, because they express this protein at haploinsufficient levels". assuming this is the only pathway by which lenalidomide helps del(5q) patients, then we can check whether CSNK1A1 is deleted or not. here, for example, it is deleted.
                        # 'EPB41L4A-AS1', 'APC', 'SRP19', 'DCP2', 'PGGT1B', 'ATG12', 'AP3S1', 'DMXL1', 'TNFAIP8', 'HSD17B4', 'SRFBP1', 'SNX2', 'CSNK1G3', 'PHAX', 'PRRC1', 'HINT1', 'CDC42SE2', 'FNIP1', 'IRF1', 'UQCRQ', 'AFF4', 'ZCCHC10', 'HSPA4', 'VDAC1', 'SKP1', 'PPP2CA', 'UBE2B', 'SAR1B', 'DDX46', 'C5orf24', 'TXNDC15', 'SMAD5', 'HNRNPA0', 'FAM13B', 'BRD8', 'KDM3B', 'ETF1', 'HSPA9', 'PAIP2', 'UBE2D2', 'CXXC5', 'PURA', 'PFDN1', 'NDUFA2', 'IK', 'HARS', 'ZMAT2', 'TAF7', 'DIAPH1', 'DELE1', 'NDFIP1', 'NR3C1', 'YIPF5', 'LARS', 'TCERG1', 'CSNK1A1', 'HMGXB3', 'TCOF1', 'RPS14', 'RBM22', 'TNIP1', 'CCDC69', 'ATOX1', 'G3BP1', 'LARP1', 'CNOT8', 'MRPL22'
                        
                        'max_median_cluster_median_projected_fold_threshold': -0.4, 

                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.3,
                        # 'metacell_median_projected_fold_threshold': -0.6,
                        # 'cell_mean_projected_fold_diff_threshold': -0.7, # awful separation.
                    },
                    {
                        'name': 'N187_del_11q_near_end',
                        'chrom': '11',
                        # TODO: make sure the first and last genes are good. maybe should not do this manually but implement an HMM for all...
                        'first_gene': 'NKAPD1',
                        # 'last_gene': 'TRAPPC4',
                        'last_gene': 'TBRG1',
                        # 'last_gene': 'RPS25',
                        # 'NKAPD1', 'TIMM8B', 'SDHD', 'REXO2', 'PAFAH1B2', 'PCSK7', 'UBE4A', 'ATP5MG', 'KMT2A', 'DDX6', 'RPS25', 'TBRG1'
                        
                        'max_median_cluster_median_projected_fold_threshold': -0.4, 

                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.3,
                        # 'metacell_median_projected_fold_threshold': -0.4,
                        # 'cell_mean_projected_fold_diff_threshold': -1.2, # awful separation.
                    },
                ],

                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N187_0': [
                #         ('N187_del_5q_middle', 'no'),
                #         ('N187_del_11q_near_end', 'no'),
                #     ],
                #     'N187_1': [
                #         ('N187_del_5q_middle', 'yes'),
                #         ('N187_del_11q_near_end', 'yes'),
                #     ],
                # },
            },
            'N227': {
                # 'num_of_mc_clusters_for_identifying_clones': 23,
                'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.8,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 10,
                # 'consider_small_clusters_combined': True,

                'cna_clustermap_width': 3,
                'cna_clustermap_height': 8,
                'cna_clustermap_colors_ratio': 0.12,

                'cna_cluster_linkage_metric': 'euclidean',
                'num_of_mc_cna_clusters': 4,
                'num_of_cell_cna_clusters': 8,

                'donor_general_manual_comment': 'clinical data says (among others) "t(1q;21p)". we dont see the t(21p). assuming the t(21p) is real, i guess we miss it because of the low number of genes there',

                'cna_infos': [
                    {
                        'name': 'N227_dup_1q_except_start',
                        'chrom': '1',
                        # 'first_gene': 'SCNM1',
                        'first_gene': 'PLEKHO1',
                        
                        'min_median_cluster_median_projected_fold_threshold': 0.26, 

                        # 'expected_projected_fold': 0.58,
                        # 'metacell_median_projected_fold_diff_threshold': 0.3,
                        # 'metacell_median_projected_fold_threshold': 0.3,
                        # 'cell_mean_projected_fold_diff_threshold': 0.25, # awful separation. though probably very few or none at all healthy cells?

                    },
                    {
                        'name': 'N227_del_7',
                        'chrom': '7',

                        'max_median_cluster_median_projected_fold_threshold': -0.43,
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.25,
                        
                        'max_mean_cell_cluster_mean_projected_fold_threshold': -0.23,

                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.6,
                        # 'metacell_median_projected_fold_threshold': -0.7,
                        # 'cell_mean_projected_fold_diff_threshold': -0.5, 

                    },

                    # {
                    #     'name': 'fake_neg_control_N227_dup_8',
                    #     'chrom': '8',
                        
                    #     'min_median_cluster_median_projected_fold_threshold': 0.46, 

                    #     # 'expected_projected_fold': 0.58,
                    #     # 'metacell_median_projected_fold_diff_threshold': 0.25,
                    #     # 'metacell_median_projected_fold_threshold': 0.35,
                    #     # 'cell_mean_projected_fold_diff_threshold': 0.75, 
                    # },

                    # {
                    #     'name': 'fake_neg_control_N227_del_5_except_5p_start',
                    #     'chrom': '5',
                    #     # 'first_gene': 'GPBP1',
                    #     'first_gene': 'SREK1IP1',

                    #     'max_median_cluster_median_projected_fold_threshold': -0.5, 

                    #     # 'expected_projected_fold': -1,
                    #     # 'metacell_median_projected_fold_diff_threshold': -0.5,
                    #     # 'metacell_median_projected_fold_threshold': -0.5,
                    #     # 'cell_mean_projected_fold_diff_threshold': -0.7, 
                    # },
                ],
                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N227_0': [
                #         ('N227_dup_1q_except_start', 'no'),
                #         ('N227_del_7', 'no'),
                #     ],
                #     'N227_1': [
                #         ('N227_dup_1q_except_start', 'yes'),
                #         ('N227_del_7', 'no'),
                #     ],
                #     'N227_2': [
                #         ('N227_dup_1q_except_start', 'yes'),
                #         ('N227_del_7', 'partially'),
                #     ],
                #     'N227_3': [
                #         ('N227_dup_1q_except_start', 'yes'),
                #         ('N227_del_7', 'yes'),
                #     ],
                # },
            },
            'N200': {
                'min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist': 8,
                # 'mc_cna_clusters_threshold_dist_quantile': 0.99,
                'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.3,

                'cna_clustermap_width': 4,
                'cna_clustermap_height': 8,
                'cna_clustermap_colors_ratio': 0.11,

                'num_of_cell_cna_clusters': 2,
                'cna_cluster_linkage_metric': 'euclidean',

                'cna_infos': [
                    {
                        'name': 'N200_del_3p_end_and_3q_start',
                        'chrom': '3',
                        # 'first_gene': 'ZNF445',
                        'first_gene': 'KIAA1143',
                        # 'last_gene': 'NEPRO',
                        'last_gene': 'SPICE1',
                        # 'ZNF445', 'KIAA1143', 'SACM1L', 'CCDC12', 'SETD2', 'SMARCC1', 'DHX30', 'MAP4', 'NME6', 'TMA7', 'SHISA5', 'IP6K2', 'PRKAR2A', 'ARIH2', 'WDR6', 'NDUFAF3', 'IMPDH2', 'QARS', 'USP4', 'RHOA', 'RBM6', 'RBM5', 'GNAI2', 'RASSF1', 'TEX264', 'ABHD14B', 'RPL29', 'WDR82', 'PBRM1', 'GNL3', 'SPCS1', 'TKT', 'DCP1A', 'SELENOK', 'FAM208A', 'APPL1', 'SLMAP', 'PDHB', 'THOC7', 'PSMD6', 'UBA3', 'ARL6IP5', 'FOXP1', 'RYBP', 'PPP4R2', 'CHMP2B', 'ZNF654', 'TFG', 'SENP7', 'TRMT10C', 'RPL24', 'BBX', 'CD47', 'ATG3', 'NEPRO'

                        'max_median_cluster_median_projected_fold_threshold': -0.5, 
                        
                        'max_mean_cell_cluster_mean_projected_fold_threshold': -0.5, 

                        'exp_names': ['demux_11_04_21_2_illumina', 'demux_11_04_21_2', 'demux_10_01_22_1', 'demux_26_03_23_1'],
                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.5,
                        # 'metacell_median_projected_fold_threshold': -0.5,
                        # 'cell_mean_projected_fold_diff_threshold': -0.85, 
                    },
                    {
                        'name': 'N200_del_5q',
                        'chrom': '5',
                        # 'first_gene': 'GPBP1',
                        'first_gene': 'SREK1IP1',

                        'max_median_cluster_median_projected_fold_threshold': -0.5, 
                        
                        'max_mean_cell_cluster_mean_projected_fold_threshold': -0.38, 
                        'exp_names': ['demux_11_04_21_2_illumina', 'demux_11_04_21_2', 'demux_10_01_22_1', 'demux_26_03_23_1'],

                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.5,
                        # 'metacell_median_projected_fold_threshold': -0.5,
                        # 'cell_mean_projected_fold_diff_threshold': -0.7, 
                    },
                    {
                        'name': 'N200_dup_8',
                        'chrom': '8',
                        
                        'min_median_cluster_median_projected_fold_threshold': 0.46, 

                        'min_mean_cell_cluster_mean_projected_fold_threshold': 0.45, 
                        'exp_names': ['demux_11_04_21_2_illumina', 'demux_11_04_21_2', 'demux_10_01_22_1', 'demux_26_03_23_1'],

                        # 'expected_projected_fold': 0.58,
                        # 'metacell_median_projected_fold_diff_threshold': 0.25,
                        # 'metacell_median_projected_fold_threshold': 0.35,
                        # 'cell_mean_projected_fold_diff_threshold': 0.75, 
                    },
                    {
                        'name': 'N200_del_10q_middle',
                        'chrom': '10',
                        'first_gene': 'ANAPC16',
                        'last_gene': 'VDAC2',

                        'max_median_cluster_median_projected_fold_threshold': -0.48, 
                        
                        'skip_cna': True, # it doesn't satisfy the requirements of the protocol (too small a region and not extremely strong). it jumps out as it seems to be a part of the clear clone, but the cna_clustermap suggests that it is not exactly the same clone...
                    },
                ],
            },
            'N186': {
                'num_of_mc_clusters_for_identifying_clones': 4,
                'consider_small_clusters_combined': True,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 1,

                'cna_clustermap_width': 7,
                'cna_clustermap_height': 10,
                'cna_clustermap_colors_ratio': 0.05,
                # 'cna_clustermap_abs_vmin_vmax': 3,
                'num_of_mc_cna_clusters': 7,

                'cna_infos': [
                    {
                        'name': 'N186_del_5q_middle',
                        'chrom': '5',
                        'first_gene': 'SMAD5',
                        # 'last_gene': 'SIL1',
                        'last_gene': 'HARS',

                        # the following mc thresholds are because of a single B-cell metacell (see N186_del_5q_middle_mc_vs_cell_mean_normalized_projected_fold_scatter.png). i guess it contains some cells without a del(5q).
                        'max_median_cluster_median_projected_fold_threshold': -0.5, 
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.3, 
                    },
                    {
                        'name': 'N186_dup_dup_9p',
                        'chrom': '9',
                        # 'last_gene': 'PTAR1',
                        'last_gene': 'DCAF10',

                        'min_mean_cell_cluster_mean_projected_fold_threshold': 0.5,

                        'min_median_cluster_median_projected_fold_threshold': 0.75, 
                        'min_median_cluster_median_projected_fold_partial_threshold': 0.4, 
                    },
                    # {
                    #     'name': 'N186_dup_9q_middle',
                    #     'chrom': '9',
                    #     # 'first_gene': 'SMC5', # I guess not.
                    #     'first_gene': 'TRIM14',
                    #     # 'last_gene': 'PTBP3', 
                    #     'last_gene': 'INIP', 
                    #     # 
                    #     'min_median_cluster_median_projected_fold_threshold': 0.55, 
                    #     'min_median_cluster_median_projected_fold_partial_threshold': 0.33, 

                    # },
                    {
                        'name': 'N186_dup_9q',
                        'chrom': '9',
                        # 'first_gene': 'SMC5', # I guess not.
                        'first_gene': 'FAM122A',
                        # 'last_gene': 'PTBP3', 
                        # 'last_gene': 'INIP', 
                        # 
                        'min_median_cluster_median_projected_fold_threshold': 0.2, 
                        # 'min_median_cluster_median_projected_fold_partial_threshold': 0.33, 
                    },
                    {
                        # feels a bit borderline, as the metacell hist has a peak around 0.5 rather than 0.58.
                        'name': 'N186_dup_10',
                        'chrom': '10',
                        # 

                        'min_median_cluster_median_projected_fold_threshold': 0.25, 
                    },
                    {
                        'name': 'N186_del_17p_start',
                        'chrom': '17',
                        'last_gene': 'STX8',

                        'max_median_cluster_median_projected_fold_threshold': -0.4, 
                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.4, 
                    },
                    # {
                    #     'name': 'N186_dup_18q_middle',
                    #     'chrom': '18',
                    #     'first_gene': 'KIAA1328',
                    #     # 'first_gene': 'GALNT1',
                    #     'last_gene': 'SMAD2', 
                    #     # 

                    #     'min_median_cluster_median_projected_fold_threshold': 0.7, 
                    #     'min_median_cluster_median_projected_fold_partial_threshold': 0.27, 
                    # },
                    {
                        'name': 'N186_dup_dup_18q_middle',
                        'chrom': '18',
                        # 'first_gene': 'INO80C',
                        'first_gene': 'C18orf21',
                        # 'first_gene': 'GALNT1',
                        'last_gene': 'SMAD2', 
                        # 
                        
                        'min_mean_cell_cluster_mean_projected_fold_threshold': 1.5,

                        'min_median_cluster_median_projected_fold_threshold': 2, 
                        'min_median_cluster_median_projected_fold_partial_threshold': 0.8, 


                        'manual_comment': 'the number of genes we consider is very small, so maybe it is form of regional transcriptional regulation?',
                    },
                    
                    {
                        'name': 'N186_del_20p_start',
                        'chrom': '20',
                        # 'first_gene': 'GALNT1',
                        'last_gene': 'RALGAPA2', 
                        
                        'max_median_cluster_median_projected_fold_threshold': -0.45, 
                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.25, # 230427: not relevant to any clone currently

                    },
                    {
                        'name': 'N186_del_20q_middle',
                        'chrom': '20',
                        'first_gene': 'COMMD7',
                        # 'first_gene': 'NXT1', # NOTE: i think this is a great example for a case in which we are not sure about the exact edge of the deletion, and seeing genes strongly upregulated seems like (strong?) evidence against these genes being part of the deletion (HM13, PDRG1, TM9SF4). I.e., strongly upregulated despite a deletion seems to me more surprising than strongly downregulated despite a duplication. thus i think NXT1 being the first gene in the deletion is a mistake.
                        'last_gene': 'PIGT', 
                        # 
                        
                        'max_median_cluster_median_projected_fold_threshold': -0.45, 
                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.35, 

                    },
                    {
                        # feels borderline, but the metacell hist has a peak very close to 0.58, so i guess it is legit.
                        'name': 'N186_dup_22',
                        'chrom': '22',
                        # 

                        
                        'min_median_cluster_median_projected_fold_threshold': 0.3, 
                        # 'min_median_cluster_median_projected_fold_partial_threshold': 0.25, 
                        
                    },
                    
                    # {
                    #     # maybe true in some subclusters, but not when considering the two big ones
                    #     'name': 'N186_del_4q_middle',
                    #     'chrom': '4',
                    #     # 
                    #     'first_gene': 'HADH',
                    #     # 'last_gene': 'PAPSS1', 
                    #     # 'last_gene': 'JADE1', 
                    #     'last_gene': 'SCLT1', 

                        
                    #     # 'min_median_cluster_median_projected_fold_threshold': 0.4, 
                        
                    # },
                ],
                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N186_2': [
                #         ('N186_del_5q_middle', 'yes'),
                #         ('N186_dup_dup_9p_and_9q_start', 'yes'),
                #         ('N186_dup_9q_middle', 'yes'),
                #         ('N186_dup_10', 'yes'),
                #         ('N186_del_17p_start', 'yes'),
                #         ('N186_dup_dup_18q_middle', 'yes'),
                #         ('N186_del_20p_start', 'yes'),
                #         ('N186_del_20q_middle', 'yes'),
                #         ('N186_dup_22', 'yes'),
                #     ],
                #     'N186_1_1': [
                #         ('N186_del_5q_middle', 'yes'),
                #         ('N186_dup_dup_9p_and_9q_start', 'no'),
                #         ('N186_dup_9q_middle', 'partially'),
                #         ('N186_dup_10', 'no'),
                #         ('N186_del_17p_start', 'yes'),
                #         ('N186_dup_dup_18q_middle', 'no'), # the only (major) difference between N186_1_1 and N186_1_2
                #         ('N186_del_20p_start', 'no'),
                #         ('N186_del_20q_middle', 'yes'),
                #         ('N186_dup_22', 'no'),
                #     ],
                #     'N186_1_2': [
                #         ('N186_del_5q_middle', 'yes'),
                #         ('N186_dup_dup_9p_and_9q_start', 'no'),
                #         ('N186_dup_9q_middle', 'partially'),
                #         ('N186_dup_10', 'no'),
                #         ('N186_del_17p_start', 'yes'),
                #         ('N186_dup_dup_18q_middle', 'partially'),
                #         ('N186_del_20p_start', 'no'),
                #         ('N186_del_20q_middle', 'yes'),
                #         ('N186_dup_22', 'no'),
                #     ],
                # },
            },
            'N48': {
                # NOTE: the signal seems way too weak...
                
                'num_of_mc_clusters_for_identifying_clones': 7, # to have the metacells with a seeming weak signal of dup_12 be in a separate cluster
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 9,
                # 'consider_small_clusters_combined': True,

                'single_cna_median_projected_fold_bins': (-np.inf, 0.4, np.inf), # dummy so we don't crash.

                'cna_clustermap_width': 3,
                'cna_clustermap_height': 10,
                'cna_clustermap_colors_ratio': 0.15,                

                'min_num_of_cells_in_cluster_to_show_in_cluster_moving_median_plot': 30,
                'cna_cluster_linkage_metric': 'euclidean',

                'cna_infos': [
                    
                    # https://en.wikipedia.org/wiki/Chronic_lymphocytic_leukemia: "Chronic lymphocytic leukemia (CLL) is a type of cancer in which the bone marrow makes too many lymphocytes (a type of white blood cell)" and "CLLs are, in virtually all cases, preceded by a particular subtype of monoclonal B-cell lymphocytosis (MBL). This subtype, termed chronic lymphocytic leukemia-type MBL (CLL-type MBL) is an asymptomatic, indolent, and chronic disorder in which people exhibit a mild increase in the number of circulating B-cell lymphocytes. These B-cells are abnormal: they are monoclonal, i.e. produced by a single ancestral B-cell".
                    # https://ashpublications.org/blood/article/123/26/4101/32798/Trisomy-12-chronic-lymphocytic-leukemia-cells: "Chronic lymphocytic leukemia (CLL) is a disease of considerable clinical and genetic heterogeneity. Trisomy 12 is the third most common cytogenetic abnormality and has several distinguishing features including abnormal morphology and increased prevalence of NOTCH1 mutations.1,2  Although trisomy 12 is present in approximately 16% of cases of CLL, the prevalence of this cytogenetic abnormality is significantly higher in small lymphocytic lymphoma (SLL) where it is present in 28% of cases.3  Furthermore, acquisition of trisomy 12 also has been recently implicated in a third of cases of Richter’s transformation"
                    {
                        # seems that it is only in B cells
                        'name': 'N48_dup_12',
                        'chrom': '12',

                        # 'min_median_cluster_median_projected_fold_partial_threshold': 0.15,
                        'min_median_cluster_median_projected_fold_partial_threshold': 0.24,
                        # 'min_median_cluster_median_projected_fold_threshold': 0.3, # dummy
                        # 'min_mean_cell_cluster_mean_projected_fold_threshold': 0.35,
                        'exp_names': ['demux_03_11_21_1'],
                        'manual_comment': 'borderline. pretty clear in demux_03_11_21_1 (in B cells). maybe also in demux_22_02_21_1, but at least when looking at metacells (in demux_22_02_21_1) the evidence seems too weak',
                        'manual_comment_for_blood_aging_paper': 'seems to appear only in B cells',
                    },
                ],
            },
            'N243': {
                'num_of_mc_clusters_for_identifying_clones': 3,
                
                # 'single_cna_median_projected_fold_bins': (-np.inf, 0.4, np.inf), # dummy so we don't crash.
                'cna_clustermap_width': 3,
                'cna_clustermap_height': 10,
                'cna_clustermap_colors_ratio': 0.15,

                'min_num_of_cells_in_cluster_to_show_in_cluster_moving_median_plot': 30,
                
                'cna_infos': [
                    {
                        # too many genes have negative projected fold values, so commenting out the threshold (which was too low anyway).
                        'name': 'N243_dup_3',
                        'chrom': '3',
                        
                        # 'min_median_cluster_median_projected_fold_threshold': 0.4, # dummy
                        'min_median_cluster_median_projected_fold_partial_threshold': 0.2, 

                        'min_mean_cell_cluster_mean_projected_fold_threshold': 0.3, # looks real? but why so few cells?

                        # 'manual_comment': 'very borderline.',
                        'skip_cna': True, # 231004: on second thought, not even strong enough to be considered borderline.
                    },
                    {
                        # too many genes have negative projected fold values, so commenting out the threshold (which was too low anyway).
                        'name': 'N243_del_2q_middle',
                        'chrom': '2',
                        'first_gene': 'GYPC',
                        'last_gene': 'ZEB2',
                        
                        'max_median_cluster_median_projected_fold_threshold': -0.4, # dummy
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.37, 

                        'skip_cna': True, # the number of genes is too small, and the effect is not extremely strong.

                        # 'min_mean_cell_cluster_mean_projected_fold_threshold': 0.3, # looks real? but why so few cells?
                    },

                    # {
                    #     # NOTE: added it just because https://www.uptodate.com/contents/chronic-myelomonocytic-leukemia-clinical-features-evaluation-and-diagnosis/print says "The most common cytogenetic abnormalities in CMML are trisomy 8 and various abnormalities of chromosome 7"
                    #     'name': 'N243_dup_8',
                    #     'chrom': '8',
                        
                    #     'min_median_cluster_median_projected_fold_threshold': 0.4, # dummy
                    # },
                    
                    # { too weak, and too small a region, to consider.
                    #     'name': 'N243_del_4p_middle',
                    #     'chrom': '4',
                        
                    #     # 'min_median_cluster_median_projected_fold_threshold': 0.4, 
                    # },
                ],
            },
            'N201': {
                'num_of_mc_clusters_for_identifying_clones': 6,
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.44,
                # 'consider_small_clusters_combined': True,

                'cna_clustermap_width': 3,
                'cna_clustermap_height': 8,
                'cna_clustermap_colors_ratio': 0.15,
                # 'cna_clustermap_abs_vmin_vmax': 1,
                
                # 'single_cna_median_projected_fold_bins': (-np.inf, -0.44, -0.27, np.inf),
                'single_cna_median_projected_fold_bins': (-np.inf, -0.44, np.inf),

                'cna_infos': [
                    {
                        # the signal is weak, and the region is very small, but it seems to me like exactly the common region of almost all del(5q) CNAs that we see, so I am not commenting it out.
                        'name': 'N201_del_5q_middle',
                        'chrom': '5',
                        'first_gene': 'PHAX',
                        'last_gene': 'DUSP1',
                        
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.25, 
                        # 'max_median_cluster_median_projected_fold_threshold': -0.7, # not really relevant to our data.

                        # 'mc_median_normalized_projected_fold_hist_bins': 15,
                        # 'mc_median_normalized_projected_fold_hist_bins': 30,

                        # 'num_of_rpt_non_edge_bps': 2,

                        # 'skip_cna': True,
                        
                        # NOTE: after looking at cell hists, seems like we don't have enough evidence for this one. maybe even my best guess would be that there is no CNA at all here. yet being exactly in the common region is suspicious.
                        'manual_comment': 'somewhat weak evidence, yet the region seems to contain the common region for almost all del(5q) CNAs that we see. N201_del_5q_middle_uncorrected_X_cna_region_expr_log_ratio_hist.png seems like the strongest evidence',
                    },
                    # {
                    #     # doesnt look real
                    #     'name': 'N201_dup_18_middle',
                    #     'chrom': '18',
                    #     'first_gene': 'ROCK1',
                    #     'last_gene': 'LMAN1',
                        
                    #     # 'min_median_cluster_median_projected_fold_threshold': 0.4, 

                    #     # 'mc_median_normalized_projected_fold_hist_bins': 15,
                    #     # 'mc_median_normalized_projected_fold_hist_bins': 30,
                    # },
                ],
                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N201_1': [
                #         ('N201_del_5q_middle', 'yes'), 
                #     ],
                # },

            },
            'N198': {
                # 'num_of_mc_clusters_for_identifying_clones': 1,
                'num_of_mc_clusters_for_identifying_clones': 20,
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.78,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 20, 

                # TODO: demux_02_01_22_1: something in chrom 6 in a subset of MKPs?
            },
            'N199': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N220': {
                'num_of_mc_clusters_for_identifying_clones': 2,
            },
            'N238': {
                'num_of_mc_clusters_for_identifying_clones': 10,
                # 'consider_small_clusters_combined': True,
            },
            'N239': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N240': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N241': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N75': {
                # 'num_of_mc_clusters_for_identifying_clones': 65,
                'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.2,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 3,
                
                # 'cna_infos': [
                #     {
                #         'name': 'N75_dup_19_i_guess_RNA_only',
                #         'chrom': '19',
                        

                #         # 'min_mean_cell_cluster_mean_projected_fold_threshold': 1.2,
                        
                #         'min_median_cluster_median_projected_fold_partial_threshold': 0.25, 
                #         # 'min_median_cluster_median_projected_fold_partial_threshold': 0.4, 

                #         'manual_comment': 'I guess this is only in RNA and not DNA. also much stronger in demux_23_11_20_1_illumina than in demux_23_11_20_1. looks like a different state?',
                        
                #         # 'skip_cna': True, # not convincing enough.
                #     },
                # ],
            },
            'N295': {
                # 'num_of_mc_clusters_for_identifying_clones': 65,
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.2,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 3,
                
                # 'cna_infos': [ # NOTE: that was due to some weird batchy cells with low nuclear mRNA levels. maybe their nuclei failed to lyse?
                #     {
                #         'name': 'N75_dup_19_i_guess_RNA_only',
                #         'chrom': '19',
                        

                #         # 'min_mean_cell_cluster_mean_projected_fold_threshold': 1.2,
                        
                #         'min_median_cluster_median_projected_fold_partial_threshold': 0.25, 
                #         # 'min_median_cluster_median_projected_fold_partial_threshold': 0.4, 

                #         'manual_comment': 'I guess this is only in RNA and not DNA. also much stronger in demux_23_11_20_1_illumina than in demux_23_11_20_1. looks like a different state?',
                        
                #         # 'skip_cna': True, # not convincing enough.
                #     },
                # ],
            },
            'N242': {
                # 'num_of_mc_clusters_for_identifying_clones': 1,
                'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.65,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 15,
                
                # TODO: need to analyze again with a proper atlas for BM... i suspect the stuff i see is not CNAs are only differences in RNA and not in DNA...
            },
            'N244': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N12': {
                'num_of_mc_clusters_for_identifying_clones': 2,
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.78,
                'cna_infos': [
                    {
                        'name': 'N12_del_17q_near_start',
                        'chrom': '17',
                        
                        'first_gene': 'POLDIP2',
                        'last_gene': 'TAOK1',

                        # 'min_mean_cell_cluster_mean_projected_fold_threshold': 1.2,
                        
                        # 'min_median_cluster_median_projected_fold_partial_threshold': 0.24, 
                        # 'min_median_cluster_median_projected_fold_partial_threshold': 0.4, 
                        
                        'skip_cna': True, # not convincing enough.
                    },

                ],
            },
            'N60': {
                'num_of_mc_clusters_for_identifying_clones': 3,
            },
            'N209': {
                'num_of_mc_clusters_for_identifying_clones': 3,
            },
            'N224': {
                'num_of_mc_clusters_for_identifying_clones': 7,
                # 'consider_small_clusters_combined': True,
            },
            'N225': {
                'num_of_mc_clusters_for_identifying_clones': 4,
                # 'consider_small_clusters_combined': True,
            },
            'N226': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N185': {
                'num_of_mc_clusters_for_identifying_clones': 3,
            },
            'N188': {
                'num_of_mc_clusters_for_identifying_clones': 3,
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.86,

                'min_num_of_cells_in_cluster_to_show_in_cluster_moving_median_plot': 10,

                'cna_infos': [
                    {
                        'name': 'N188_dup_14',
                        'chrom': '14',
                        
                        # 'min_mean_cell_cluster_mean_projected_fold_threshold': 1.2,
                        
                        'min_median_cluster_median_projected_fold_partial_threshold': 0.1, 
                        
                        'skip_cna': True, # definitely not convincing enough.
                    },
                    {
                        'name': 'N188_dup_xp_near_start',
                        'chrom': 'X',
                        
                        'first_gene': 'TRAPPC2',
                        'last_gene': 'ZFX',

                        'min_mean_cell_cluster_mean_projected_fold_threshold': 1.2,
                        
                        'min_median_cluster_median_projected_fold_partial_threshold': 0.24, 
                        # 'min_median_cluster_median_projected_fold_partial_threshold': 0.4, 
                        
                        'skip_cna': True, # not convincing enough.
                    },

                ],
            },
            'N203': {
                'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.74,
                # 'num_of_mc_clusters_for_identifying_clones': 6, # to have the cluster that is composed of monocytes and a single CLP-M, for which a weak dup_8 signal seems to exist.

                'min_num_of_cells_in_cluster_to_show_in_cluster_moving_median_plot': 15,

                'cna_infos': [
                    {
                        # borderline, but N203_dup_8_cell_total_umis_min_quantile_0.8_mean_normalized_projected_fold_hist.png is pretty convincing.
                        'name': 'N203_dup_8',
                        'chrom': '8',

                        'mc_median_normalized_projected_fold_hist_bins': 15,

                        'min_mean_cell_cluster_mean_projected_fold_threshold': 0.55,
                        
                        'min_median_cluster_median_projected_fold_partial_threshold': 0.13, 

                        'manual_comment': 'very borderline.',
                        'manual_comment_for_blood_aging_paper': 'somewhat borderline',
                    },
                ],
            },
            'N204': {
                'num_of_mc_clusters_for_identifying_clones': 3,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 3, 

                'cna_infos': [
                    {
                        'name': 'N204_del_xq_near_start',
                        'chrom': 'X',
                        'first_gene': 'IL2RG',
                        'last_gene': 'ABCB7',

                        'skip_cna': True, # not convincing enough. probably shouldn't even consider it if i work according to the protocol.
                        # 'max_median_cluster_median_projected_fold_threshold': -0.8, 
                    },
                ],
            },
            'N206': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N205': {
                'num_of_mc_clusters_for_identifying_clones': 4,
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.98,
                'mc_cna_clusters_threshold_dist_quantile': 0.99,
                # 'consider_small_clusters_combined': True,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 2, # to show the presumably "healthy" clone in the end.
                
                'cna_clustermap_width': 5,
                'cna_clustermap_height': 10,
                'cna_clustermap_colors_ratio': 0.07,
                
                'cna_infos': [
                    {
                        # a very small number of genes
                        'name': 'N205_del_3p_middle',
                        'chrom': '3',
                        'first_gene': 'TRANK1',
                        'last_gene': 'WDR48',

                        'max_median_cluster_median_projected_fold_threshold': -0.2, 

                        'skip_cna': True, # the number of genes is too small, and the effect is not extremely strong.
                    },
                    {
                        # feels borderline
                        'name': 'N205_dup_3p_middle',
                        'chrom': '3',
                        'first_gene': 'RPSA',
                        'last_gene': 'ARF4',
                        
                        # 'first_gene': 'SMARCC1',
                        # 'last_gene': 'MAPKAPK3',

                        'min_median_cluster_median_projected_fold_threshold': 0.12, 

                        'expected_projected_fold': 0.58,

                        
                        'manual_comment': 'the fluctuating pattern in chromosome 3 is very suspicious. could this be some form of regional transcriptional regulation? (keeping this CNA because the number of genes in it is not too small). hypothesis (not mine): this might be inversion 3 it causes a unique gene expression pattern',
                    },
                    {
                        # very small number of genes
                        'name': 'N205_del_3q_middle',
                        'chrom': '3',
                        'first_gene': 'MFSD1',
                        'last_gene': 'MECOM',

                        'max_median_cluster_median_projected_fold_threshold': -0.3,
                        
                        'skip_cna': True, # the number of genes is too small, and the effect is not extremely strong.
                    },

                    {
                        # borderline - the MC histogram has a peak around 0.45 instead of 0.58
                        'name': 'N205_dup_3q_end',
                        'chrom': '3',
                        # 'first_gene': 'TNFSF10',
                        # 'first_gene': 'PIK3CA',
                        'first_gene': 'MYNN',
                        # 'first_gene': 'ALG3',

                        'min_median_cluster_median_projected_fold_threshold': 0.12, 
                        
                        'expected_projected_fold': 0.58,

                        
                        'manual_comment': 'the fluctuating pattern in chromosome 3 is very suspicious. could this be some form of regional transcriptional regulation? (keeping this CNA because the number of genes in it is not too small). hypothesis (not mine): this might be inversion 3 it causes a unique gene expression pattern',
                    },
                    
                    {
                        # borderline - the MC histogram has a peak around 0.46 instead of 0.58
                        'name': 'N205_dup_8',
                        'chrom': '8',
                        # 'first_gene': 'TNFSF10',
                        
                        'min_median_cluster_median_projected_fold_threshold': 0.36, 

                    },

                    {
                        'name': 'N205_dup_9p_and_9q_start',
                        'chrom': '9',
                        # 'first_gene': 'BICD2',
                        'last_gene': 'NOL8',
                       
                        # 'max_median_cluster_median_projected_fold_threshold': -0.5, 

                        'skip_cna': True, # the signal is too weak.

                        # 'manual_comment': 'very borderline, as the number of genes is small, and the fluctuating pattern in chromosome 9 is very suspicious. could this be some form of regional transcriptional regulation?',
                    },
                    {
                        'name': 'N205_del_9q_middle',
                        'chrom': '9',
                        'first_gene': 'BICD2',
                        'last_gene': 'TRIM14',
                       
                        'max_median_cluster_median_projected_fold_threshold': -0.5, 

                        'manual_comment': 'very borderline, as the number of genes is small, and the fluctuating pattern in chromosome 9 is very suspicious. could this be some form of regional transcriptional regulation?',
                    },
                    {
                        'name': 'N205_del_9q_near_end',
                        'chrom': '9',
                        'first_gene': 'GTF3C4',
                        'last_gene': 'BRD3',
                       
                        'max_median_cluster_median_projected_fold_threshold': -0.8, # dummy
                        
                        'skip_cna': True, # the number of genes is too small, and the effect is not extremely strong.
                    },
                    
                    {
                        # borderline - the MC histogram has a peak around 0.44 instead of 0.58
                        'name': 'N205_dup_11',
                        'chrom': '11',

                        'min_median_cluster_median_projected_fold_threshold': 0.35, 
                        
                    },
                    
                    {
                        'name': 'N205_del_12p_middle',
                        'chrom': '12',
                        'first_gene': 'KLRG1',
                        'last_gene': 'EMP1',
                       
                        'max_median_cluster_median_projected_fold_threshold': -0.8, # dummy
                        'skip_cna': True, # the number of genes is too small, and the effect is not extremely strong.
                    },
                    {
                        # borderline - the MC histogram has a peak around 0.45 instead of 0.58
                        'name': 'N205_dup_13',
                        'chrom': '13',
                        
                        'min_median_cluster_median_projected_fold_threshold': 0.35, 

                    },
                    
                    {
                        'name': 'N205_dup_14',
                        'chrom': '14',

                        'min_median_cluster_median_projected_fold_threshold': 0.37, 

                    },
                
                ],
                
                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N205_0': [
                #         ('N205_del_3p_middle', 'no'),
                #         ('N205_dup_3p_middle', 'no'),
                #         ('N205_del_3q_middle', 'no'),
                #         ('N205_dup_3q_end', 'no'),
                #         ('N205_dup_8', 'no'),
                #         ('N205_del_9q_middle', 'no'),
                #         ('N205_del_9q_near_end', 'no'),
                #         ('N205_dup_11', 'no'),
                #         ('N205_del_12p_middle', 'no'),
                #         ('N205_dup_13', 'no'),
                #         ('N205_dup_14', 'no'),
                #     ],
                #     'N205_1': [
                #         ('N205_del_3p_middle', 'yes'),
                #         ('N205_dup_3p_middle', 'yes'),
                #         ('N205_del_3q_middle', 'yes'),
                #         ('N205_dup_3q_end', 'yes'),
                #         ('N205_dup_8', 'yes'),
                #         ('N205_del_9q_middle', 'yes'),
                #         ('N205_del_9q_near_end', 'yes'),
                #         ('N205_dup_11', 'yes'),
                #         ('N205_del_12p_middle', 'yes'),
                #         ('N205_dup_13', 'yes'),
                #         ('N205_dup_14', 'yes'),
                #     ],
                # },
    
            },
            'N190': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N192': {
                'num_of_mc_clusters_for_identifying_clones': 2,
                # 'cna_infos': [
                #     {
                #         'name': 'N192_dup_3',
                #         'chrom': '3',
                #         'skip_cna': True,
                #     },
                # ],
            },
            'N193': {
                # 'num_of_mc_clusters_for_identifying_clones': 1,
                'num_of_mc_clusters_for_identifying_clones': 6,

                'donor_general_manual_comment': 'expression differs from atlas in many short regions, none very clear, so cant identify CNAs (or their absence) confidently',
            },
            'N191': {
                'num_of_mc_clusters_for_identifying_clones': {
                    ('final_only_N191_only_demux_28_11_21_1', 'demux_28_11_21_1'): 20, # ugh. 10 was also ugh.
                    ('final_only_N191_only_demux_n_04_12_23_1', 'demux_n_04_12_23_1'): 5,
                },
                

                'cna_clustermap_width': 5,
                'cna_clustermap_height': 10,
                'cna_clustermap_colors_ratio': 0.07,
                'mc_cna_clusters_threshold_dist_quantile': 0.87,
                
                'cna_infos': [
                    # {
                    #     # NOTE: appears in all donors in demux_28_11_21_1 (i.e., batchy), so commenting out the threshold. maybe upregulated as part of the (viral) stress response?
                    #     # NOTE: contains 8 HLA genes
                    #     'name': 'N191_dup_6p_middle',
                    #     'chrom': '6',
                    #     # 'first_gene': 'HLA-DRA',
                    #     'first_gene': 'LTB',
                    #     'last_gene': 'HLA-DPB1',

                    #     # 'min_median_cluster_median_projected_fold_threshold': 0.45, 
                    # },
                    # {
                    #     # NOTE: appears in all donors in demux_28_11_21_1 (i.e., batchy), so commenting out the threshold. maybe upregulated as part of the (viral) stress response?
                    #     # NOTE: contains 8 HLA genes
                    #     'name': 'N191_dup_6p_middle_shorter',
                    #     'chrom': '6',
                    #     # 'first_gene': 'HLA-DRA',
                    #     'first_gene': 'NELFE',
                    #     'last_gene': 'HLA-DPB1',
                        
                    #     # 'min_median_cluster_median_projected_fold_threshold': 0.45, 
                    # },
                    
                    
                    {
                        # N191_del_5q_middle_cell_total_umis_min_quantile_0.9_mean_normalized_projected_fold_hist.png is quite convincing that this is real.
                        'name': 'N191_del_5q_middle',
                        'chrom': '5',
                        # 'first_gene': 'CHD1',
                        # 'last_gene': 'SPARC', # 'last_gene': 'MED7',
                        
                        # looks better in demux_n_04_12_23_1 (chosen by the rpt lib)
                        'first_gene': 'COX7C',
                        'last_gene': 'G3BP1',

                        # 'num_of_rpt_non_edge_bps': 2,

                        'manual_comment': 'observed in demux_28_11_21_1, and much more clearly in demux_n_04_12_23_1',

                        'exp_names': ['demux_28_11_21_1', 'demux_n_04_12_23_1'],
                        

                        'max_median_cluster_median_projected_fold_partial_threshold': -0.2, 
                    },

                    # ugh chrom 4 and 7. i guess i just shouldn't trust the colors.
                    { # NOTE: in the original mc clustermap, it seems obvious that there is something here, and also in cells with high umi count, it kind of seems like this is real. but it also seems like it is not the whole chromosome.
                        'name': 'N191_del_7',
                        'chrom': '7',
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.25, 

                        'exp_names': ['demux_28_11_21_1'],
                        'manual_comment': 'observed in demux_28_11_21_1, and only maybe in demux_n_04_12_23_1',
                    },

                    {
                        'name': 'N191_dup_13',
                        'chrom': '13',

                        'min_median_cluster_median_projected_fold_threshold': 0.48, 

                        'exp_names': ['demux_28_11_21_1'],
                        'manual_comment': 'observed in demux_28_11_21_1, but not in demux_n_04_12_23_1',
                    },

                    {
                        'name': 'N191_del_17p_start',
                        'chrom': '17',
                        # 'last_gene': 'MAP2K4',
                        'last_gene': 'ELAC2',

                        'max_median_cluster_median_projected_fold_threshold': -0.4, 
                        
                        'exp_names': ['demux_28_11_21_1'],
                        'manual_comment': 'observed in demux_28_11_21_1, but not in demux_n_04_12_23_1',
                    },
                    
                    {
                        'name': 'N191_del_21',
                        'chrom': '21',

                        'max_median_cluster_median_projected_fold_threshold': -0.4, 
                        
                        'exp_names': ['demux_28_11_21_1'],
                        'manual_comment': 'observed in demux_28_11_21_1, but not in demux_n_04_12_23_1',
                    },

                    {
                        'name': 'N191_del_xq_middle',
                        'chrom': 'X',
                        # 'first_gene': 'ITM2A',
                        'first_gene': 'NAP1L3',
                        'last_gene': 'TCEAL3',
                        # 'last_gene': 'TCEAL1',

                        'max_median_cluster_median_projected_fold_threshold': -0.8, 
                        
                        'exp_names': ['demux_28_11_21_1'],
                        'manual_comment': 'observed in demux_28_11_21_1, but not in demux_n_04_12_23_1. the number of genes we consider is very small, so maybe it is form of regional transcriptional regulation?',

                        # 'expected_projected_fold': -1,
                        # # NOTE: very weird pair of metacell and cell hists. the metacell one suggests a deletion in the active X chromosome (N191 is female), leading to very low expression of the region. yet the cells seem like some of them actually have pretty high expression level of this region. weird. perhaps this region is often somewhat skewed to the right? it is also this way for N218.
                        # 'metacell_median_projected_fold_diff_threshold': -0.5, 
                        # 'metacell_median_projected_fold_threshold': -0.6,
                        # 'cell_mean_projected_fold_diff_threshold': -1.5, # pretty bad separation. 
                    },
                ],

                # 'clone_name_to_cna_names_and_detection_levels': {
                #     'N191_1': [
                #         'N191_del_5q_middle',
                #     ],
                #     'N191_2': [
                #         'N191_del_7q_middle',
                #         'N191_dup_13',
                #         'N191_del_17p_start',
                #         'N191_del_21',
                #         'N191_del_xq_middle',
                #     ],
                # },
            },
            
            'N228': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N229': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N230': {
                'num_of_mc_clusters_for_identifying_clones': 2,
            },
            'N207': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N208': {
                'num_of_mc_clusters_for_identifying_clones': 2,
            },
            'N210': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N231': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N232': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N233': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N234': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N213': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N214': {
                'num_of_mc_clusters_for_identifying_clones': 5,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 8,
                # 'consider_small_clusters_combined': True,
            },
            'N189': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N194': {
                'num_of_mc_clusters_for_identifying_clones': 1,
                'cna_infos': [
                    {
                        'name': 'N194_del_5q_middle',
                        'chrom': '5',
                        
                        # 'num_of_rpt_non_edge_bps': 2,
                        'first_gene': 'IRF1',
                        'last_gene': 'PAIP2',
                        
                        # 'num_of_rpt_non_edge_bps': 1,
                        # 'first_gene': 'PHAX',

                        'manual_comment': 'weak evidence',
                        'skip_cna': True, # way too weak signal
                        # 'max_median_cluster_median_projected_fold_threshold': -0.15, 
                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.15, 
                    },
                ],
            },
            'N195': {
                'num_of_mc_clusters_for_identifying_clones': 2,
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.
                'cna_infos': [
                    {
                        'name': 'N195_del_13',
                        'chrom': '13',
                        # 'first_gene': 'MCTP1',
                        # 'last_gene': 'UBLCP1',


                        'skip_cna': True, # way too weak signal
                        # 'max_median_cluster_median_projected_fold_threshold': -0.15, 

                        # 'expected_projected_fold': -1,
                        # 'metacell_median_projected_fold_diff_threshold': -0.5,
                        # 'metacell_median_projected_fold_threshold': -0.5,
                        # 'cell_mean_projected_fold_diff_threshold': -0.8, # awful separation.
                    },
                ],
            },
            'N196': {
                'num_of_mc_clusters_for_identifying_clones': 3,
            },
            'N197': {
                'num_of_mc_clusters_for_identifying_clones': 3,
            },
            'N235': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N215': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N217': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N218': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N223': {
                'num_of_mc_clusters_for_identifying_clones': 1,
            },
            'N245': {
                # 'num_of_mc_clusters_for_identifying_clones': 1,
                'mc_clusters_for_identifying_clones_threshold_dist_quantile': 0.3,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 10,
                # 'consider_small_clusters_combined': True,

                # 230601: demux_02_05_22_1_ultima__demux_02_05_22_2_ultima, demux_02_05_22_2_ultima, demux_02_05_22_1_ultima: seems like no CNAs here, according to protocol. but it does look like there are some abnormally expressed genes.
            },
            'N248': {
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 1,
                
                # 230601: demux_02_05_22_2_ultima: seems like no CNAs here, according to protocol. maybe maybe a deletion in 6q end, and a deletion in 14q middle, but only in a single metacell...
            },
            'N249': {
                'num_of_mc_clusters_for_identifying_clones': 6,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 25,
                # 'consider_small_clusters_combined': True,
                
                # 230601: demux_06_06_22_2_ultima: seems like no CNAs here, according to protocol.
                # 230601: demux_06_06_22_accident: seems like no CNAs here, according to protocol.
            },
            'N250': {
                'num_of_mc_clusters_for_identifying_clones': 4,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 1,
                # 'consider_small_clusters_combined': True,
                
                'cna_infos': [
                    {
                        # very very borderline. 230914: i guess i noticed it only because i knew about his lenalidomide treatment. 
                        'name': 'N250_del_5q_middle',
                        'chrom': '5',
                        'first_gene': 'CETN3',
                        'last_gene': 'MAT2B',

                        'max_median_cluster_median_projected_fold_partial_threshold': -0.2,
                        
                        'skip_cna': True, # doesn't look convincing enough.

                    },
                ],
                
                # 230602: demux_09_06_22_1_ultima: seems like no CNAs here, according to protocol.
                # 230602: demux_09_06_22_2_ultima: seems like no CNAs here, according to protocol.
            },
            'N251': {
                'num_of_mc_clusters_for_identifying_clones': 5,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 25,
                # 'consider_small_clusters_combined': True,
                
                # 230602: demux_09_06_22_1_ultima: seems like no CNAs here, according to protocol.
                # 230602: demux_09_06_22_2_ultima: seems like no CNAs here, according to protocol.
            },
            'N253': {
                'num_of_mc_clusters_for_identifying_clones': 1,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 25,
                'consider_small_clusters_combined': True,
                
                # 230602: demux_06_06_22_2_ultima: seems like no CNAs here, according to protocol.
                # 230602: demux_06_06_22_accident: seems like no CNAs here, according to protocol.
            },
            'N254': {
                'num_of_mc_clusters_for_identifying_clones': 1,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 25,
                'consider_small_clusters_combined': True,
                
                # 230602: : seems like no CNAs here, according to protocol.
                # 230602: : seems like no CNAs here, according to protocol.
            },
            'NS20': {
                'num_of_mc_clusters_for_identifying_clones': 6,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 15,
                'consider_small_clusters_combined': True,
                
                # 230601: NS20_29_09_21_1: seems like no CNAs here, according to protocol.
            },

            'N257': {
                'num_of_mc_clusters_for_identifying_clones': 3,
                
                
                'cna_infos': [
                    {
                        'name': 'N257_del_4q',
                        'chrom': '4',
                        'first_gene': 'BMP2K',
                        # 'last_gene': '',
                        # 'num_of_rpt_non_edge_bps': 1,

                        'max_median_cluster_median_projected_fold_partial_threshold': -0.23, 

                        'skip_cna': True, # doesn't look convincing enough, also maybe appears in N285 48h samples, so i guess it is some regulated response.
                    },
                ],
            },

            
            '01_01_23_1_a': {
                'num_of_mc_clusters_for_identifying_clones': 1,
                # 'num_of_mc_clusters_for_identifying_clones': 6,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 15,
                
                # 230604: quick look: seems like no CNAs here
            },

            'N367': {
                'num_of_mc_clusters_for_identifying_clones': 6,
                # 'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 15,


                'cna_infos': [
                    {
                        'name': 'N367_del_13q_middle',
                        'chrom': '13',
                        
                        'first_gene': 'HMGB1',
                        'last_gene': 'PCDH9',
                        # 'num_of_rpt_non_edge_bps': 2,

                        'manual_comment': r'contains RB1 (mentioning because it is well-known. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2234042/ says "It is estimated that RB is dysfunctional/inactivated in up to 40% of human leukemias")',

                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.43, 

                        # 'skip_cna': True, # doesn't look convincing enough.
                    },
                ],

            },
            'N369': {
                'num_of_mc_clusters_for_identifying_clones': 6,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 15,


                'cna_infos': [
                    {
                        'name': 'N369_del_8q_end',
                        'chrom': '8',
                        'first_gene': 'GSDMD',

                        'max_median_cluster_median_projected_fold_partial_threshold': -0.43, 

                        'skip_cna': True, # doesn't look convincing enough.
                    },
                    {
                        'name': 'N369_del_20q_near_start',
                        'chrom': '20',
                        # 'num_of_rpt_non_edge_bps': 2,
                        # 'first_gene': 'HM13',
                        # 'last_gene': 'ADNP',
                        'first_gene': 'EIF2S2',
                        'last_gene': 'ZFAS1',

                        'max_median_cluster_median_projected_fold_threshold': -0.6, 
                    },
                ],

            },
            'N280': {
                # 'mc_clusters_for_identifying_clones_threshold_dist_quantile': {
                #     'final_only_N280_only_N280_bm_06_11_22_1': 0.95,
                #     'final_only_N280_only_demux_06_11_22_1': 0.3,
                #     'final_only_N280_only_demux_n_15_08_23_2': 0.3, # dummy
                #     'final_only_N280_only_demux_n_bm_15_08_23_1': 0.3, # dummy
                # },
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 15,

                'cna_infos': [
                    {
                        'name': 'N280_del_7',
                        'chrom': '7',

                        'max_median_cluster_median_projected_fold_partial_threshold': -0.45, 

                        'manual_comment': 'weirdly, i can see this del(7) in cells of N280 in demux_06_11_22_1 and in demux_n_15_08_23_2 (cell hists are much clearer in demux_n_15_08_23_2, and seems like less cells have del(7) in demux_n_15_08_23_2, though N280_del_7_uncorrected_X_cna_region_expr_log_ratio_hist is weird), but not in N280_bm_06_11_22_1. but maybe this is only because projecting a BM sample on the PB atlas is not good enough? TODO: demux_n_bm_15_08_23_1',

                        'exp_names': ['demux_06_11_22_1', 'demux_n_15_08_23_2'],

                        # 'skip_cna': True, # doesn't look convincing enough.

                    },
                ],
            },
            'N376': {
                'num_of_mc_clusters_for_identifying_clones': 4,
                'min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot': 10,

                'cna_infos': [
                    {
                        'name': 'N376_del_7',
                        'chrom': '7',

                        'max_median_cluster_median_projected_fold_partial_threshold': -0.25, 

                        # 'skip_cna': True, # doesn't look convincing enough.
                    },
                    {
                        'name': 'N376_del_20q',
                        'chrom': '20',
                        
                        # 'num_of_rpt_non_edge_bps': 2,
                        'first_gene': 'EIF2S2',
                        'last_gene': 'ZFAS1',
                        'max_median_cluster_median_projected_fold_partial_threshold': -0.4, 

                        # 'skip_cna': True, # doesn't look convincing enough.
                    },
                    {
                        'name': 'N376_dup_in_X',
                        'chrom': 'X',
                        
                        # 'num_of_rpt_non_edge_bps': 2,

                        'first_gene': 'ZMAT1',
                        'last_gene': 'BEX3',

                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.33, 

                        'skip_cna': True, # doesn't look convincing enough.
                    },
                ],
            },
            'N335': {
                # 'num_of_mc_clusters_for_identifying_clones': 3,

                'cna_infos': [
                    # { # NOTE: that was due to some weird batchy neut cells or batchy monocyte-mebemp doublets
                    #     'name': 'N335_dup_4_except_start',
                    #     'chrom': '4',
                    #     'first_gene': 'CD38',

                    #     'min_median_cluster_median_projected_fold_partial_threshold': 0.25, 

                    #     # 'skip_cna': True, # doesn't look convincing enough.
                    # },
                    # {
                    #     'name': 'N335_dup_9p_start',
                    #     'chrom': '9',
                    #     'last_gene': 'KLHL9',

                    #     'min_median_cluster_median_projected_fold_partial_threshold': 0.5, 

                    #     # 'skip_cna': True, # doesn't look convincing enough.
                    # },
                ],
            },
            'N403': {
                'num_of_mc_clusters_for_identifying_clones': 1,

                'cna_infos': [
                    {
                        'name': 'N403_del_5q_middle',
                        'chrom': '5',

                        # 'num_of_rpt_non_edge_bps': 2,

                        'first_gene': 'TTC37',
                        'last_gene': 'MAT2B',

                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.25, 

                        # 'skip_cna': True, # doesn't look convincing enough.
                    },
                    {
                        'name': 'N403_del_7',
                        'chrom': '7',

                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.25, 

                        # 'skip_cna': True, # doesn't look convincing enough.
                    },
                ],
            },
            'N405': {
                'num_of_mc_clusters_for_identifying_clones': 1,

                'cna_infos': [
                    {
                        'name': 'N405_dup_21',
                        'chrom': '21',

                        # 'num_of_rpt_non_edge_bps': 2,

                        # 'first_gene': 'TTC37',
                        # 'last_gene': 'MAT2B',
                        
                        'manual_comment': 'hard to be sure as pretty much all metacells (and cells) of N405 in demux_n_30_11_23_1 look the same. but chrom 21 expression in N405 seems higher than others in demux_n_30_11_23_1',

                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.25, 

                        # 'skip_cna': True, # doesn't look convincing enough.
                    },
                ],
            },
            'G8': {
                # 'num_of_mc_clusters_for_identifying_clones': 1,
                'cna_infos': [
                    {
                        'name': 'G8_del_20q_middle',
                        'chrom': '20',
                        
                        # 'num_of_rpt_non_edge_bps': 2,
                        'first_gene': 'EIF2S2',
                        'last_gene': 'ZFAS1',
                        
                        # 'num_of_rpt_non_edge_bps': 1,
                        # 'first_gene': 'PHAX',

                        # 'manual_comment': 'weak evidence',
                        # 'skip_cna': True, # way too weak signal
                        # 'max_median_cluster_median_projected_fold_threshold': -0.15, 
                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.15, 
                    },
                ],
            },
            'N95': { # 73 years old male. demux_21_12_20_2 ('N96', 'N98', 'N95', 'N92'). 
                # donor 95 from healthy atlas (73yo male). other donors in the same experiment had very few (if at all) B cells. can't tell if there really are CNAs in 12 and 7q, or this is a B cell state uncommon in our healthy atlas.

                'cna_infos': [
                    {
                        'name': '95_dup_12q',
                        'chrom': '12',
                        
                        # 'num_of_rpt_non_edge_bps': 2,
                        # 'first_gene': 'TUBA1A',
                        # 'last_gene': 'PAWR',
                        
                        # 'num_of_rpt_non_edge_bps': 1,
                        'first_gene': 'TUBA1A',

                        'manual_comment': 'weak signal. maybe a B cell state uncommon in our healthy atlas?',
                        'skip_cna': True, # weak signal
                        # 'max_median_cluster_median_projected_fold_threshold': -0.15, 
                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.15, 
                    },
                    {
                        'name': '95_del_7q',
                        'chrom': '7',
                        
                        # 'num_of_rpt_non_edge_bps': 1,
                        'first_gene': 'SMIM30',

                        'manual_comment': 'weak signal. maybe a B cell state uncommon in our healthy atlas?',
                        'skip_cna': True, # weak signal
                        # 'max_median_cluster_median_projected_fold_threshold': -0.15, 
                        # 'max_median_cluster_median_projected_fold_partial_threshold': -0.15, 
                    },
                ],
            },

        },

        # batchy CNAs??
        # NOTE: could it be something specific to chromosome ends??? maybe some telomere related stress response that inhibits transcription nearby???
        # del(8q_end) (demux_07_02_22_1)
        # del(11_middle) (demux_07_02_22_1)
        # del(16p_start) (demux_07_02_22_1) (very clear)
        # del(17q_end) (demux_07_02_22_1)
        # del(19p_start) (demux_07_02_22_1) (very clear)
        # somewhat: del(4p_start) (demux_07_02_22_1)
        # somewhat: del(5q_end) (demux_07_02_22_1)
        # somewhat: del(9q_end) (demux_07_02_22_1)
        # somewhat: del(6p_middle) (demux_02_03_22_1)
        
    },

    'hba_hbb_contaminated_cells': {
        
        'min_hba_and_hbb_expr': -5,
    },
    'neut_or_neut_contaminated_cells': {
        'neut_specific_genes': ['CAMP', 'LTF', 'LCN2'],
        'min_neut_specific_expr': -7.5,
        'min_s100a8_s100a9_expr': -8,
    },
    'batchy_monocyte_mebemp_doublet': {
        'higher_in_mebemp_e_genes': HIGHER_IN_MEBEMP_E_THAN_MONOCYTE_GENES,
        'higher_in_monocyte_genes': STRONG_MONOCYTE_SPECIFIC_GENE_NAMES,
        'min_higher_in_mebemp_e_expr': -8,
        'min_higher_in_monocyte_expr': -10, # see higher_in_mono_hist_only_rel_balanced_exps_and_high_higher_in_mebemp.png
    },
    'min_norm_nuclear_expr_score': -2,

    'gene_loc_df_csv_file_path': os.path.join(mds_in_out_dir_paths.OUTPUT_ROOT_DIR_PATH, 'gene_loc_df.csv'),
}

def get_mds_path_dict(dataset_name=None):
    main_paths = MDS_ANALYSIS_PARAMS['dataset_name_to_main_paths'][dataset_name]
    mds_out_dir_path = main_paths['out_dir_path']
    intermediate_out_dir_path = os.path.join(mds_out_dir_path, 'intermediate_output')
    path_dict = {
        'mds_out_dir': mds_out_dir_path,
        'intermediate_out_dir': intermediate_out_dir_path,
        'raw_ad_file_path': main_paths['raw_ad_file_path'],
        'preprocessing_out_dir_file_path': main_paths['preprocessing_out_dir_file_path'],
        'preprocessing_out_exp_df_csv_file_path': main_paths['preprocessing_out_exp_df_csv_file_path'],
        
        'chrom_pos_df_csv': os.path.join(intermediate_out_dir_path, 'chrom_pos_df.csv'),
    
        # 'raw_with_downsampled_cell_umis_before_excluding_genes_ad_file_path': os.path.join(intermediate_out_dir_path, 'raw_with_downsampled_cell_umis.h5ad'), # seems to me like we don't need this.
        'raw_with_clean_gene_mask_ad_file_path': os.path.join(intermediate_out_dir_path, 'raw_with_clean_gene_mask.h5ad'),
        'too_few_non_excluded_umis_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'too_few_non_excluded_umis_count_df.csv'),
        'too_high_mito_and_malat1_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'too_high_mito_and_malat1_count_df.csv'),
        'too_low_ribo_prot_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'too_low_ribo_prot_count_df.csv'),
        'hba_hbb_contaminated_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'hba_hbb_contaminated_count_df.csv'),
        'neut_or_neut_contaminated_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'neut_or_neut_contaminated_count_df.csv'),
        'platelet_or_platelet_contaminated_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'platelet_or_platelet_contaminated_count_df.csv'),
        'too_low_norm_nuclear_expr_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'too_low_norm_nuclear_expr_count_df.csv'),
        'mebemp_monocyte_doublet_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'mebemp_monocyte_doublet_count_df.csv'),
        'enriched_genotype_doublet_mc_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'enriched_genotype_doublet_mc_count_df.csv'),
        'high_log_norm_umi_count_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'high_log_norm_umi_count_count_df.csv'),
        'dc_nkt_monocyte_endothel_b_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'dc_nkt_monocyte_endothel_b_count_df.csv'),
        'presumably_hspc_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'presumably_hspc_count_df.csv'),
        'cell_filtering_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'cell_filtering_count_df.csv'),

        'cells_discarded_due_to_many_cells_per_donor_obs_names_file_path': os.path.join(intermediate_out_dir_path, 'cells_discarded_due_to_many_cells_per_donor_obs_names.txt'),
        'too_high_norm_cr_num_of_umis_obs_names_file_path': os.path.join(intermediate_out_dir_path, 'too_high_norm_cr_num_of_umis_obs_names.txt'),
        'cells_discarded_due_to_donor_attributes_obs_names_file_path': os.path.join(intermediate_out_dir_path, 'cells_discarded_due_to_donor_attributes_obs_names.txt'),
        'for_removing_doublet_metacells_cells_discarded_due_to_donor_attributes_obs_names_file_path': os.path.join(intermediate_out_dir_path, 'for_removing_doublet_metacells_cells_discarded_due_to_donor_attributes_obs_names.txt'),
        'for_ambient_noise_cells_discarded_due_to_donor_attributes_obs_names_file_path': os.path.join(intermediate_out_dir_path, 'for_ambient_noise_cells_discarded_due_to_donor_attributes_obs_names.txt'),
        'cells_ad_before_basic_cell_exclusion_file_path': os.path.join(intermediate_out_dir_path, 'cells_ad_before_basic_cell_exclusion.h5ad'),
        'cells_ad_before_discarding_mk_file_path': os.path.join(intermediate_out_dir_path, 'cells_ad_before_discarding_mk.h5ad'),
        'clean_ad_file_path': os.path.join(intermediate_out_dir_path, 'raw_without_excluded_genes_and_cells.h5ad'),
        'excluded_gene_umi_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'excluded_gene_umi_df.csv'),
        'with_doublets_with_mc_c_ad_file_path': os.path.join(intermediate_out_dir_path, 'with_doublets_with_mc_c.h5ad'),
        'with_doublets_mc_ad_file_path': os.path.join(intermediate_out_dir_path, 'with_doublets_mc.h5ad'),
        'assigned_and_no_doublet_enirhced_mc_c_ad_file_path': os.path.join(intermediate_out_dir_path, 'assigned_and_no_doublet_enirhced_mc_c.h5ad'),
        'for_mcnoise_with_mc_c_ad_file_path': os.path.join(intermediate_out_dir_path, 'for_mcnoise_with_mc_c.h5ad'),
        'for_mcnoise_mc_ad_file_path': os.path.join(intermediate_out_dir_path, 'for_mcnoise_mc.h5ad'),



        # 'with_doublets_outliers_ad_file_path': os.path.join(intermediate_out_dir_path, 'with_doublets_outliers.h5ad'),
        # 'with_doublets_metacells_metadata_csv_file_path': os.path.join(intermediate_out_dir_path, 'with_doublets_metacells_metadata.csv'),
        'cells_with_too_high_norm_cr_num_of_umis_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'cells_with_too_high_norm_cr_num_of_umis_df.csv'),
        'cells_in_metacells_marked_as_doublet_metacells_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'cells_in_metacells_marked_as_doublet_metacells_df.csv'),
        'potential_batch_degs_dir_path': os.path.join(intermediate_out_dir_path, 'potential_batch_DEGs'),

        'platelet_pooled_downsampled_umis_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'platelet_pooled_downsampled_umis_df.csv'),

        'mk_only_cells_for_metacells_ad_file_path': os.path.join(intermediate_out_dir_path, 'mk_only_cells_for_metacells.h5ad'),
        'mk_only_cells_with_metacells_ad_file_path': os.path.join(intermediate_out_dir_path, 'mk_only_cells_with_metacells.h5ad'),
        'mk_only_bare_metacells_ad_file_path': os.path.join(intermediate_out_dir_path, 'mk_only_bare_metacells.h5ad'),
        'mk_only_bare_outliers_ad_file_path': os.path.join(intermediate_out_dir_path, 'mk_only_bare_outliers.h5ad'),
        
        'cells_after_doublet_cleanup_ad_file_path': os.path.join(intermediate_out_dir_path, 'cells_after_doublet_cleanup.h5ad'),
        'for_ambient_noise_estimation_bare_metacells_ad_file_path': os.path.join(intermediate_out_dir_path, 'for_ambient_noise_estimation_bare_metacells.h5ad'),
        'for_ambient_noise_estimation_bare_outliers_ad_file_path': os.path.join(intermediate_out_dir_path, 'for_ambient_noise_estimation_outliers.h5ad'),
        'exp_names_excluded_from_mcnoise': os.path.join(intermediate_out_dir_path, 'exp_names_excluded_from_mcnoise.txt'),
        'batch_empty_droplet_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'batch_empty_droplet_df.csv'),
        'ambient_noise_levels_per_exp_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'ambient_noise_levels_per_exp_df.csv'),
        'ambient_noise_anf_pickle_file_path': os.path.join(intermediate_out_dir_path, 'ambient_noise_anf.pickle'),
        'ambient_noise_ane_pickle_file_path': os.path.join(intermediate_out_dir_path, 'ambient_noise_ane.pickle'),
        'ambient_noise_cells_obs_mcnoise_columns_csv': os.path.join(intermediate_out_dir_path, 'ambient_noise_cells_obs_mcnoise_columns.csv'),
        # 'ambient_noise_cells_with_noise_level_estimations_ad_file_path': os.path.join(
        #     intermediate_out_dir_path, 'ambient_noise_cells_with_noise_level_estimations.h5ad'),

        'metacells_for_ambient_noise_estimation_with_uncorrected_projection_ad_file_path': os.path.join(intermediate_out_dir_path, 'metacells_for_ambient_noise_estimation_with_uncorrected_projection.h5ad'),
        'metacells_for_ambient_noise_estimation_with_corrected_projection_ad_file_path': os.path.join(intermediate_out_dir_path, 'metacells_for_ambient_noise_estimation_with_corrected_projection.h5ad'),

        'cells_with_metacell_attrs_ad_file_path': os.path.join(intermediate_out_dir_path, 'cells_with_metacell_attrs.h5ad'),
        'metacells_metadata_csv_file_path': os.path.join(intermediate_out_dir_path, 'metacells_metadata.csv'),

        # 'metacells_with_uncorrected_projection_ad_file_path': os.path.join(intermediate_out_dir_path, 'metacells_with_uncorrected_projection.h5ad'),
        # 'metacells_with_corrected_projection_ad_file_path': os.path.join(intermediate_out_dir_path, 'metacells_with_corrected_projection.h5ad'),
        # 'projection_weights_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'projection_weights_df.csv'),
        # 'cell_metadata_df_for_mcview_csv_file_path': os.path.join(intermediate_out_dir_path, 'cell_metadata_df_for_mcview.csv'),
        # 'cell_to_metacell_df_for_mcview_csv_file_path': os.path.join(intermediate_out_dir_path, 'cell_to_metacell_df_for_mcview.csv'),

        'metacells_with_correlations_ad_file_path': os.path.join(intermediate_out_dir_path, 'metacells_with_correlations.h5ad'),


        'correlated_feature_and_lateral_genes_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'correlated_feature_and_lateral_genes_df.csv'),
        'max_norm_expression_ratio_df': os.path.join(intermediate_out_dir_path, 'max_norm_expression_ratio_df.csv'),
        'sex_diff_genes_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'sex_diff_genes_df.csv'),
        'same_donor_doublets_identified_by_doublet_enriched_metacells_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'same_donor_doublets_identified_by_poibin_doublet_metacells_count_df.csv'),
        'same_donor_doublets_identified_by_too_high_norm_cr_num_of_umis_count_df_csv_file_path': os.path.join(intermediate_out_dir_path, 'same_donor_doublets_identified_by_too_high_norm_cr_num_of_umis_count_df.csv'),
        'batch_diff_expression_out_dir_path': os.path.join(intermediate_out_dir_path, 'batch_diff_expression'),
        'batchy_metacells_out_dir_path': os.path.join(intermediate_out_dir_path, 'batchy_metacells'),
        'excluding_cells_out_dir_path': os.path.join(intermediate_out_dir_path, 'excluding_cells'),
        'karyotype_estimation_out_dir_path': os.path.join(intermediate_out_dir_path, 'karyotype_estimation'),
        'all_my_single_donor_mc_projected_correlation_df_csv': os.path.join(intermediate_out_dir_path, 'all_my_single_donor_mc_projected_correlation_df.csv'),

        'exp_df_csv_file_path': os.path.join(mds_out_dir_path, 'exp_df1.csv'),
        'exp_df_after_cleaning_cells_csv_file_path': os.path.join(mds_out_dir_path, 'exp_df2.csv'),
        'mk_only_metacell_log_file_path': os.path.join(mds_out_dir_path, 'mk_only_mc_log.txt'),
        'metacell_log_file_path': os.path.join(mds_out_dir_path, 'mc_log.txt'),
        'metacell_projection_log_file_path': os.path.join(mds_out_dir_path, 'mc_projection_log.txt'),
        'ambient_noise_removal_log_file_path': os.path.join(mds_out_dir_path, 'ambient_noise_removal_log.txt'),
        'single_donor_mcview_proj_dir_path': os.path.join(mds_out_dir_path, 'single_donor_mcview_proj'),
        'multiple_donor_mcview_proj_dir_path': os.path.join(mds_out_dir_path, 'multiple_donor_mcview_proj'),
        'ugly_mcview_config_csv_file_paths_file_path': os.path.join(mds_out_dir_path, 'ugly_mcview_config_csv_file_paths.txt'),
        'out_figs_dir_path': os.path.join(mds_out_dir_path, 'out_figs'),
        'out_tables_dir_path': os.path.join(mds_out_dir_path, 'out_tables'),
    }

    path_dict['gene_ratio_from_technical_repeats_df'] = os.path.join(path_dict['potential_batch_degs_dir_path'], 'gene_ratio_df.csv')

    path_dict['orig_num_of_umis_per_cell_hist'] = os.path.join(path_dict['out_figs_dir_path'], 'orig_num_of_umis_per_cell_hist.png')
    path_dict['ribosomal_protein_gene_umi_fraction_hist'] = os.path.join(path_dict['out_figs_dir_path'], 'ribosomal_protein_gene_umi_fraction_hist.png')
    path_dict['mitochondrial_gene_umi_fraction_hist'] = os.path.join(path_dict['out_figs_dir_path'], 'mitochondrial_gene_umi_fraction_hist.png')
    path_dict['malat1_umi_fraction_hist'] = os.path.join(path_dict['out_figs_dir_path'], 'malat1_umi_fraction_hist.png')
    path_dict['mitochondrial_to_ribosomal_protein_gene_umi_log_ratio_hist'] = os.path.join(path_dict['out_figs_dir_path'], 'mitochondrial_to_ribosomal_protein_gene_umi_log_ratio_hist.png')
    path_dict['num_of_non_excluded_umis_hist'] = os.path.join(path_dict['out_figs_dir_path'], 'num_of_non_excluded_umis_hist.png')
    path_dict['excluded_umi_fraction_hist'] = os.path.join(path_dict['out_figs_dir_path'], 'excluded_umi_fraction_hist.png')
    path_dict['concise_qc_summary_dir_path'] = os.path.join(path_dict['out_figs_dir_path'], 'concise_QC_summary')
    path_dict['threshold_choice_summary_dir_path'] = os.path.join(path_dict['out_figs_dir_path'], 'threshold_choice_summary')

    path_dict['ambient_noise_removal_figs_dir'] = os.path.join(path_dict['out_figs_dir_path'], 'ambient_noise_removal')
    path_dict['with_doublets_metacells_figs_dir'] = os.path.join(path_dict['out_figs_dir_path'], 'with_doublets_metacells')
    path_dict['final_metacell_model_figs_dir'] = os.path.join(path_dict['out_figs_dir_path'], 'final_metacell_model')

    all_out_paths = list(path_dict.values())
    assert len(all_out_paths) == len(set(all_out_paths))
    
    for x in [
        intermediate_out_dir_path,
        mds_out_dir_path,
        path_dict['single_donor_mcview_proj_dir_path'],
        path_dict['multiple_donor_mcview_proj_dir_path'],
        path_dict['out_figs_dir_path'],
        path_dict['out_tables_dir_path'],
        path_dict['with_doublets_metacells_figs_dir'],
        path_dict['ambient_noise_removal_figs_dir'],
        path_dict['final_metacell_model_figs_dir'],
        path_dict['potential_batch_degs_dir_path'],
        path_dict['batch_diff_expression_out_dir_path'],
        path_dict['excluding_cells_out_dir_path'],
        path_dict['batchy_metacells_out_dir_path'],
        path_dict['concise_qc_summary_dir_path'],
        path_dict['threshold_choice_summary_dir_path'],
        path_dict['karyotype_estimation_out_dir_path'],
    ]:
        pathlib.Path(x).mkdir(parents=True, exist_ok=True)
    
    return path_dict




# MDS_ANALYSIS_PARAMS['doublets']['same_donor_doublets_df_csv_file_paths'] = [
#     MDS_ANALYSIS_PARAMS['out_paths']['same_donor_doublets_identified_by_doublet_enriched_metacells_count_df_csv_file_path'],
#     MDS_ANALYSIS_PARAMS['out_paths']['same_donor_doublets_identified_by_too_high_norm_cr_num_of_umis_count_df_csv_file_path'],
# ]


MDS_ANALYSIS_PARAMS['pb_hspc_cell_state_and_c_info_list'] = [
    x for x in MDS_ANALYSIS_PARAMS['pb_cd34_enriched_cell_state_and_c_info_list']
    if x[0] in PB_HSPC_STATE_NAMES
]

def get_metadata_cols_in_info_rules(state_info):
    metadata_cols_in_state_info = set()
    list_of_match_all_rules = state_info['rules']
    for list_of_match_any_rules in list_of_match_all_rules:
        for match_rule in list_of_match_any_rules:
            if match_rule[0] == 'metadata':
                metadata_cols_in_state_info.add(match_rule[1])
    return metadata_cols_in_state_info

def get_state_to_used_root_metadata_cols(state_and_c_info_list, mask_and_c_info_list):
    mask_to_used_cols = {
        mask: get_metadata_cols_in_info_rules(info)
        for mask, info in mask_and_c_info_list
    }
    mask_cols = set(mask_to_used_cols)


    state_to_used_cols = {
        state: get_metadata_cols_in_info_rules(info)
        for state, info in state_and_c_info_list
    }
    while 1:
        prev_state_to_used_cols = dict(state_to_used_cols)
        
        for state, used_cols in prev_state_to_used_cols.items():
            used_mask_cols = used_cols & mask_cols
            new_used_cols = set(itertools.chain.from_iterable([mask_to_used_cols[x] for x in used_mask_cols]))
            state_to_used_cols[state] = (state_to_used_cols[state] | new_used_cols) - used_mask_cols # watchout. don't change the original...

        if state_to_used_cols == prev_state_to_used_cols:
            break
    return state_to_used_cols

def get_gene_names_from_lateral_or_noisy_gene_names_file(lateral_gene_file_path):
    return generic_utils.read_text_file(lateral_gene_file_path).strip().split()

def get_nimrod_batchy_genes_by_kruskal_pval(max_kruskal_pval):
    gene_batch_kruskal_pval_df = pd.read_csv(MDS_ANALYSIS_PARAMS['gene_batch_kruskal_pvals_from_nimrod_csv_file_path'], sep='\t')
    batchy_gene_names = sorted(gene_batch_kruskal_pval_df.loc[gene_batch_kruskal_pval_df['kruskal_pval'] <= max_kruskal_pval, 'gene'].unique())
    return batchy_gene_names


def get_lateral_and_noisy_gene_names(is_ult_and_non_ult_mixed_model, lateral_gene_types_to_skip=set(), noisy_gene_types_to_skip=set(), add_all_noisy_genes_to_lateral_genes=True, skip_malat_ribo_mito_assert=False):
    lateral_gene_names = []
    for lateral_gene_type, curr_gene_names in MDS_ANALYSIS_PARAMS['lateral_gene_names'].items():
        if lateral_gene_type not in lateral_gene_types_to_skip:
            lateral_gene_names.extend(curr_gene_names)
    # assert len(lateral_gene_names) == len(set(lateral_gene_names)) # TUBB1 is in two groups, but i don't see a problem with that
    for lateral_gene_file_path in MDS_ANALYSIS_PARAMS['lateral_gene_names_file_paths']:
        lateral_gene_names.extend(get_gene_names_from_lateral_or_noisy_gene_names_file(lateral_gene_file_path))

    lateral_gene_names.extend(get_nimrod_batchy_genes_by_kruskal_pval(MDS_ANALYSIS_PARAMS['max_batch_kruskal_pval_to_ignore_gene_in_mc_calculation']))

    noisy_gene_names = []
    for noisy_gene_type, curr_gene_names in MDS_ANALYSIS_PARAMS['noisy_gene_names'].items():
        # print(noisy_gene_type)
        if noisy_gene_type not in noisy_gene_types_to_skip:
            noisy_gene_names.extend(curr_gene_names)
        if is_ult_and_non_ult_mixed_model:
            noisy_gene_names.extend(
                lateral_and_noisy_genes.GENES_WITH_MORE_THAN_4_ULT_ILL_FOLD
                + lateral_and_noisy_genes.GENES_WITH_HIGH_ULT_ILL_LOG_RATIO_STD_AND_AT_LEAST_SOMETIME_ULT_ILL_DIFF_EXP
            )
    for noisy_gene_file_path in MDS_ANALYSIS_PARAMS['noisy_gene_names_file_paths']:
        noisy_gene_names.extend(get_gene_names_from_lateral_or_noisy_gene_names_file(noisy_gene_file_path))

    if not skip_malat_ribo_mito_assert:
        # just make sure MALAT1, ribosomal and mitochondrial are either excluded or marked as lateral or noisy
        assert (
            set(['MALAT1'] + RIBOSOMAL_PROTEIN_GENE_NAMES + MITOCHONDRIALLY_ENCODED_GENE_NAMES) <=
            set(MDS_ANALYSIS_PARAMS['gene_exclusion']['excluded_gene_names'] + lateral_gene_names + noisy_gene_names)
        )

    lateral_gene_names = sorted(set(lateral_gene_names)) # because it is ok if i list a gene more than once.
    noisy_gene_names = sorted(set(noisy_gene_names)) # because it is ok if i list a gene more than once.
    if add_all_noisy_genes_to_lateral_genes:
        lateral_gene_names = sorted(set(lateral_gene_names + noisy_gene_names))
    return (lateral_gene_names, noisy_gene_names)

def get_lateral_gene_name_to_group_name():
    print(f'NOTE: this function does not consider lateral genes specified in lateral_gene_names_file_paths.')
    lateral_gene_name_to_group_name = {}
    for group_name, curr_gene_names in MDS_ANALYSIS_PARAMS['lateral_gene_names'].items():
        for gene_name in curr_gene_names:
            assert gene_name not in lateral_gene_name_to_group_name
            lateral_gene_name_to_group_name[gene_name] = group_name
    return lateral_gene_name_to_group_name

def get_karyotype_estimation_genes_to_ignore_names(sex_diff_gene_names, ignore_problematic_ultima_illumina_genes):
    genes_to_ignore_names = (
        []
        + sex_diff_gene_names
        # + RIBOSOMAL_PROTEIN_GENE_NAMES
        # + STRONGEST_PLATELET_GENES
        # + NON_RIBOSOMAL_VERY_HIGHLY_EXPRESSED_GENES
        # + TECH_PROBLEMATIC_GENE_NAMES
        # + get_lateral_gene_names() # e.g., for N220,N238,N239,N240 a batchy del(6p) pretty much disappears after adding this. (note that this includes ULTIMA_VS_ILLUMINA_PROBLEMATIC_GENE_NAMES_FILE_PATH)
        # + batchy_gene_names # looks cleaner without excluding these it.
        # + REPLICATION_DEPENDENT_HISTONE_GENE_NAMES # in N220,N238,N239,N240 ['HIST1H1C', 'HIST1H4C', 'HIST1H2BC', 'HIST1H2AC', 'HIST1H1E', 'HIST1H2BD', 'HIST1H4E', 'HIST1H1D'], seem much lower than projected, so it seems batchy. NOTE: commented out on 230419 as it seems pretty arbitrary.
    )
    if ignore_problematic_ultima_illumina_genes:
        genes_to_ignore_names += lateral_and_noisy_genes.ULT_ILL_PROBLEMATIC_GENES # i think this helps also when projecting to ultima controls

    return sorted(set(genes_to_ignore_names))
    
def get_karyotype_donor_exp_out_dir_path(donor_id, exp_names_repr, corrected_by_gc, mds_path_dict):
    donor_exp_out_dir_path = os.path.join(mds_path_dict['karyotype_estimation_out_dir_path'], donor_id, exp_names_repr)
    if corrected_by_gc:
        donor_exp_out_dir_path = os.path.join(donor_exp_out_dir_path, 'corrected_by_gc')
    pathlib.Path(donor_exp_out_dir_path).mkdir(parents=True, exist_ok=True)
    return donor_exp_out_dir_path

def get_karyotype_donor_cna_out_dir_path(donor_id, exp_names_repr, cna_info, corrected_by_gc):
    cna_name = cna_info['name']
    if 'generic_scan' in cna_info:
        relative_path = os.path.join('generic_scan', cna_name)
    else:
        relative_path = cna_name
    donor_cna_out_dir_path = os.path.join(get_karyotype_donor_exp_out_dir_path(donor_id, exp_names_repr, corrected_by_gc=corrected_by_gc), relative_path)
    pathlib.Path(donor_cna_out_dir_path).mkdir(parents=True, exist_ok=True)
    return donor_cna_out_dir_path

def get_donor_id_to_cna_infos(skip_cnas_with_skip_cna_attr=True):
    donor_id_to_cna_infos = {}
    for donor_id, clone_and_cna_info in MDS_ANALYSIS_PARAMS['karyotype_estimation']['donor_id_to_clone_and_cna_info'].items():
        if 'cna_infos' in clone_and_cna_info:
            cna_infos = clone_and_cna_info['cna_infos']
            if skip_cnas_with_skip_cna_attr:
                cna_infos = [x for x in cna_infos if 'skip_cna' not in x]
            if cna_infos:
                donor_id_to_cna_infos[donor_id] = cna_infos
    return donor_id_to_cna_infos

def get_cna_name_to_cna_info(donor_id):
    cna_infos = MDS_ANALYSIS_PARAMS['karyotype_estimation']['donor_id_to_clone_and_cna_info'][donor_id]['cna_infos']
    return {x['name']: x for x in cna_infos}

def get_names_of_genes_used_to_assign_given_types(cell_type_names):
    gene_names = set()
    for cell_type_to_assign_attr_infos in MDS_ANALYSIS_PARAMS['cell_type_assignment'].values():
        for cell_type, assign_attr_infos in cell_type_to_assign_attr_infos.items():
            if cell_type in cell_type_names:
                if 'expression' in assign_attr_infos:
                    list_of_genes_and_log_norm_expression_ranges = assign_attr_infos['expression']
                    for genes_and_log_norm_expression_ranges in list_of_genes_and_log_norm_expression_ranges:
                        if isinstance(genes_and_log_norm_expression_ranges, tuple) and genes_and_log_norm_expression_ranges[0] == 'at_least':
                            genes_and_log_norm_expression_ranges = genes_and_log_norm_expression_ranges[2]
                        gene_names |= {x[0] for x in genes_and_log_norm_expression_ranges}
    return sorted(gene_names)

def get_lower_and_upper_bounds_for_expected_num_of_yet_unidentified_same_donor_doublets():
    # TODO: this can be further refined to consider the expected num_of_yet_unidentified_same_donor_doublets combining different cell types. but probably in another function.

    exp_df = pd.read_csv(MDS_ANALYSIS_PARAMS['preprocessing_out_exp_df_csv_file_path'])
    skipped_exp_names = set()
    flat_dicts = []

    for _, row in exp_df.iterrows():
        if not isinstance(row['donor_id_to_num_of_cells'], str):
            assert np.isnan(row['donor_id_to_num_of_cells'])
            print(f'skipping {row["exp_name"]} because donor_id_to_num_of_cells is nan')
            skipped_exp_names.add(row['exp_name'])
            continue

        curr_flat_dict = {
            'exp_name': row['exp_name'],
        }
        for bound_type in ('upper', 'lower'):
            for donor_id_and_bound in row[f'donor_id_to_{bound_type}_bound_for_expected_num_of_same_donor_doublets'].split(';'):
                # print(donor_id_and_bound)
                donor_id, _, num_of_cells = donor_id_and_bound.partition(':')
                curr_flat_dict['donor_id'] = donor_id
                num_of_cells = np.nan if (num_of_cells == 'nan') else int(num_of_cells)
                curr_flat_dict[f'{bound_type}_bound_num_of_unidentified_same_donor_doublets'] = num_of_cells
        flat_dicts.append(curr_flat_dict)

    exp_donor_unidentified_doublets_df = pd.DataFrame(flat_dicts)

    print('ATTENTION: this code assumes we didnt remove the same cell twice, so make sure each analysis removes only cells that werent already removed')
    for same_donor_doublets_df_csv_file_path in MDS_ANALYSIS_PARAMS['doublets']['same_donor_doublets_df_csv_file_paths']:
        curr_df = pd.read_csv(same_donor_doublets_df_csv_file_path)
        assert set(list(curr_df)) == {'exp_name', 'donor_id', 'num_of_identified_same_donor_doublet_cells'}
        exp_donor_unidentified_doublets_df = exp_donor_unidentified_doublets_df.merge(curr_df, how='left')

        has_data_mask = ~exp_donor_unidentified_doublets_df['num_of_identified_same_donor_doublet_cells'].isna()
        for bound_type in ('upper', 'lower'):
            exp_donor_unidentified_doublets_df.loc[has_data_mask, f'{bound_type}_bound_num_of_unidentified_same_donor_doublets'] -= exp_donor_unidentified_doublets_df.loc[
                has_data_mask, 'num_of_identified_same_donor_doublet_cells']
        exp_donor_unidentified_doublets_df.drop('num_of_identified_same_donor_doublet_cells', axis=1, inplace=True)
    exp_donor_unidentified_doublets_df.sort_values('upper_bound_num_of_unidentified_same_donor_doublets', ascending=False, na_position='first', inplace=True)

    return exp_donor_unidentified_doublets_df

def get_mc_model_paths(mc_model_name_or_abs_path, dataset_name):
    if os.path.isabs(mc_model_name_or_abs_path):
        raise RuntimeError('i think i shouldnt use this')
        mc_model_dir_path = mc_model_name_or_abs_path
        mc_model_name = os.path.basename(mc_model_name_or_abs_path)
    else:
        mc_model_name = mc_model_name_or_abs_path
        mds_path_dict = get_mds_path_dict(dataset_name=dataset_name)
        mc_model_dir_path = os.path.join(mds_path_dict['intermediate_out_dir'], 'mc_models', mc_model_name)
    figs_dir_path = os.path.join(mc_model_dir_path, 'figs')
    pathlib.Path(figs_dir_path).mkdir(parents=True, exist_ok=True)
    mcview_proj_dir_path = os.path.join(mc_model_dir_path, 'mcview_proj')
    pathlib.Path(mcview_proj_dir_path).mkdir(parents=True, exist_ok=True)

    return {
        'mc_model_dir': mc_model_dir_path,
        
        'exp_df_before_computing_metacells_csv': os.path.join(mc_model_dir_path, 'exp_df_before_computing_metacells.csv'),
        'donor_df_before_computing_metacells_csv': os.path.join(mc_model_dir_path, 'donor_df_before_computing_metacells.csv'),
        'misc_numbered_donor_info_pickle': os.path.join(mc_model_dir_path, 'misc_numbered_donor_info.pickle'),
        'bare_metacells_before_denoising_ad': os.path.join(mc_model_dir_path, 'bare_metacells_before_denoising.h5ad'),
        'clean_with_metacells_ad': os.path.join(mc_model_dir_path, 'clean_with_metacells.h5ad'),
        'outliers_ad': os.path.join(mc_model_dir_path, 'outliers.h5ad'),
        'bare_metacells_ad': os.path.join(mc_model_dir_path, 'bare_metacells.h5ad'),
        
        'num_of_outliers': os.path.join(mc_model_dir_path, 'num_of_outliers.txt'),
        'num_of_metacells': os.path.join(mc_model_dir_path, 'num_of_metacells.txt'),

        'projection_weights_df_csv': os.path.join(mc_model_dir_path, 'projection_weights_df.csv'),
        'metacells_with_projection_ad': os.path.join(mc_model_dir_path, 'metacells_with_projection.h5ad'),
        'metacells_with_projection_checkpoint_ad': os.path.join(mc_model_dir_path, 'metacells_with_projection_checkpoint.h5ad'),
        'metacells_with_projection_checkpoint2_ad': os.path.join(mc_model_dir_path, 'metacells_with_projection_checkpoint2.h5ad'),
        'cells_with_metacell_attrs_ad': os.path.join(mc_model_dir_path, 'cells_with_metacell_attrs.h5ad'),
        
        'metacells_with_cna_attrs_etc_ad': os.path.join(mc_model_dir_path, 'metacells_with_cna_attrs_etc_ad.h5ad'),
        'cells_with_cna_attrs_etc_ad': os.path.join(mc_model_dir_path, 'cells_with_cna_attrs_etc_ad.h5ad'),

        'mutation_donor_count_df_csv': os.path.join(mc_model_dir_path, 'mutation_donor_count_df.csv'),
        
        # only healthy HSPC states (i.e., HSPC and appears in the atlas)
        'my_mc_projected_correlation_df_csv': os.path.join(mc_model_dir_path, 'my_mc_projected_correlation_df.csv'),

        'donor_diff_exp': os.path.join(mc_model_dir_path, 'donor_diff_exp'),
        
        'umi_count_per_donor_exp_state': os.path.join(mc_model_dir_path, 'umi_count_per_donor_exp_state'),

        # 'mds_healthy_controls_projection_weights_df_csv': os.path.join(mc_model_dir_path, 'mds_healthy_controls_projection_weights_df.csv'),
        # 'mds_healthy_controls_metacells_with_projection_ad': os.path.join(mc_model_dir_path, 'metacells_with_mds_healthy_controls_projection.h5ad'),
        # 'mds_healthy_controls_cells_with_metacell_attrs_ad': os.path.join(mc_model_dir_path, 'cells_with_metacell_attrs_of_mds_healthy_controls_projection.h5ad'),
        
        'cell_metadata_df_for_mcview_csv': os.path.join(mc_model_dir_path, 'cell_metadata_df_for_mcview.csv'),
        'cell_to_metacell_df_for_mcview_csv': os.path.join(mc_model_dir_path, 'cell_to_metacell_df_for_mcview.csv'),
        'ugly_mcview_config_csv': os.path.join(mc_model_dir_path, 'ugly_mcview_config.csv'),
        'figs_dir': figs_dir_path,
        'mcview_proj_dir': mcview_proj_dir_path,
        'ipssm_df_csv': os.path.join(mc_model_dir_path, 'ipssm_df.csv'),
        'ipssm_res_df_csv': os.path.join(mc_model_dir_path, 'ipssm_res_df.csv'),
        'ext_donor_feature_df': os.path.join(mc_model_dir_path, 'ext_donor_feature_df.csv'),
        'feature_cols': os.path.join(mc_model_dir_path, 'feature_cols.txt'),

        'genes_excluded_from_cna_due_to_consistent_obs_exp_diff': os.path.join(mc_model_dir_path, 'genes_excluded_from_cna_due_to_consistent_obs_exp_diff.txt'),
        'genes_excluded_from_cna_due_to_low_expr': os.path.join(mc_model_dir_path, 'genes_excluded_from_cna_due_to_low_expr.txt'),
        'genes_excluded_from_cna_due_to_missing_position_or_on_y': os.path.join(mc_model_dir_path, 'genes_excluded_from_cna_due_to_missing_position_or_on_y.txt'),
    }

def get_mc_model_minimal_donor_exp_info(mc_model_name, dataset_name):
    mc_model_paths = get_mc_model_paths(mc_model_name, dataset_name)
    
    exp_df = pd.read_csv(mc_model_paths['exp_df_before_computing_metacells_csv'])
    exp_names = sorted(exp_df.loc[~exp_df['donor_id_to_num_of_cells_for_which_metacells_were_calculated'].isna(), 'exp_name'].unique())
    exp_names_repr = '__'.join(sorted(exp_names))
    if len(exp_names) == 1:
        exp_name = exp_names[0]
    else:
        exp_name = None

    donor_df = pd.read_csv(mc_model_paths['donor_df_before_computing_metacells_csv'])
    donor_ids = sorted(donor_df['donor_id'].astype(str).unique())
    donor_ids_repr = '__'.join(sorted(donor_ids))
    if len(donor_ids) == 1:
        donor_id = donor_ids[0]
    else:
        donor_id = None

    list_of_exp_name_and_donor_id = []
    for _, row in exp_df.loc[
        ~exp_df['donor_id_to_num_of_cells_for_which_metacells_were_calculated'].isna(), 
        ['exp_name', 'donor_id_to_num_of_cells_for_which_metacells_were_calculated']
    ].iterrows():
        exp_name = row['exp_name']
        
        for donor_id_and_num_of_cells_repr in row['donor_id_to_num_of_cells_for_which_metacells_were_calculated'].split(';'):
            list_of_exp_name_and_donor_id.append((exp_name, donor_id_and_num_of_cells_repr.partition(':')[0]))
    
    if os.path.isfile(mc_model_paths['num_of_outliers']):
        num_of_outliers = int(generic_utils.read_text_file(mc_model_paths['num_of_outliers']))
    else:
        num_of_outliers = None
    if os.path.isfile(mc_model_paths['num_of_metacells']):
        num_of_metacells = int(generic_utils.read_text_file(mc_model_paths['num_of_metacells']))
    else:
        num_of_metacells = None

    return {
        'exp_names': exp_names,
        'exp_names_repr': exp_names_repr,
        'exp_name': exp_name,
        'donor_ids': donor_ids,
        'donor_ids_repr': donor_ids_repr,
        'donor_id': donor_id,
        'list_of_exp_name_and_donor_id': list_of_exp_name_and_donor_id,
        'num_of_outliers': num_of_outliers,
        'num_of_metacells': num_of_metacells,
    }

def add_missing_atlas_params(atlas_params):
    if 'mc_ad_file_path' not in atlas_params:
        atlas_mc_model_paths = get_mc_model_paths(atlas_params['atlas_mc_model_name'])
        atlas_params['mc_ad_file_path'] = atlas_mc_model_paths['metacells_with_projection_ad']

    if 'atlas_cell_type_colors_csv_file_path' not in atlas_params:
        atlas_params['atlas_cell_type_colors_csv_file_path'] = MDS_ANALYSIS_PARAMS['cell_type_colors_csv_file_path']


def get_all_existing_mc_model_names(mds_path_dict):
    return os.listdir(os.path.join(mds_path_dict['intermediate_out_dir'], 'mc_models'))


def write_my_ugly_mcview_config_csv(mc_model_name, mc_ad_file_path, bundle_name_and_mcview_proj_title=None):
    if bundle_name_and_mcview_proj_title is None:
        bundle_name_and_mcview_proj_title = mc_model_name
    mc_model_paths = get_mc_model_paths(mc_model_name)
    if MDS_ANALYSIS_PARAMS['projection']['project_to_healthy_blood_aging_atlas']:
        # metacell_type_field = 'projected_type'
        metacell_type_field = 'state'
        # mc_ad_file_path = mc_model_paths['metacells_with_cna_attrs_etc_ad']
    else:
        metacell_type_field = np.nan
        # mc_ad_file_path = mc_model_paths['bare_metacells_ad']

    donor_ids = get_mc_model_minimal_donor_exp_info(mc_model_name)['donor_ids']
    single_donor_model = len(donor_ids) == 1

    bundle_path_prefix = '/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tanaywiz/apps/hg_hematopoiesis/'
    if single_donor_model:
        # /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/single_donor_mcview_proj
        # mcview_proj_dir_path = MDS_ANALYSIS_PARAMS['out_paths']['single_donor_mcview_proj_dir_path']
        bundle_path_suffix = 'mds_single_donor_models'
    else:
        # /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/multiple_donor_mcview_proj
        # mcview_proj_dir_path = MDS_ANALYSIS_PARAMS['out_paths']['multiple_donor_mcview_proj_dir_path']
        bundle_path_suffix = 'mds_multiple_donor_models'
    bundle_path = os.path.join(bundle_path_prefix, bundle_path_suffix) # 230629: commented out but changed my mind and didn't try commenting this out yet
    # bundle_path = os.path.join(bundle_path_prefix, bundle_path_suffix, mc_model_name) # 230629: maybe this is what i need

    pd.DataFrame({
        'MCVIEW_PROJ_DIR_PATH': [mc_model_paths['mcview_proj_dir']], # this would not allow same URL with different datasets for different models.
        # 'MCVIEW_PROJ_DIR_PATH': [mcview_proj_dir_path], # need to have the same proj dir for different models if i want them to be in different datasets of the same URL. NOTE: seems like we don't want that because MCView would have memory problems due to loading all datasets...
        'METACELLS_FILE_PATH': [mc_ad_file_path],
        # 'METACELLS_METADATA_FILE_PATH': [mc_model_paths['metacells_metadata_csv_file_path']],
        'PROJECTION_WEIGHTS_FILE_PATH': [mc_model_paths['projection_weights_df_csv']],
        'CELL_METADATA_CSV_FILE_PATH': [mc_model_paths['cell_metadata_df_for_mcview_csv']],
        'CELL_TO_METACELL_CSV_FILE_PATH': [mc_model_paths['cell_to_metacell_df_for_mcview_csv']],
        'MCVIEW_PROJ_TITLE': [bundle_name_and_mcview_proj_title],
        'MC_TYPE_COLORS_CSV_FILE_PATH': [MDS_ANALYSIS_PARAMS['cell_type_colors_csv_file_path']],
        'MC_TYPE_FIELD': [metacell_type_field],
        'BUNDLE_NAME': [bundle_name_and_mcview_proj_title],
        'BUNDLE_PATH': [bundle_path],
        'ATLAS_PROJECT': ['/net/dummy/dummy/dummy/export/tgdata/db/tgdb/tanaywiz/apps/hg_hematopoiesis/240626_cHSPC_79_normal_illu_atlas/project/'],
        'ATLAS_DATASET': ['240626_cHSPC_79_normal_illu_atlas'],
    # }).to_csv(MDS_ANALYSIS_PARAMS['ugly_mcview_config_csv_file_path'], sep=',', index=False)
    }).to_csv(mc_model_paths['ugly_mcview_config_csv'], sep=',', index=False)

# def get_fixed_gene_loc_df():

    # gene_loc_df = gtf_interface_and_utils.fix_gene_loc_df_to_have_a_single_row_for_each_gene(unfixed_gene_loc_df)
    # assert gene_loc_df['gene_name'].is_unique
    # return gene_loc_df

def get_donor_id_to_genes_with_manual_arch_mutation_info():
    donor_id_to_genes_with_manual_arch_mutation_info = collections.defaultdict(set)
    for donor_id, manual_mutation_infos in MDS_ANALYSIS_PARAMS['arch_mutations']['donor_id_to_manual_arch_mutation_infos'].items():
        for mutation_info in manual_mutation_infos:
            donor_id_to_genes_with_manual_arch_mutation_info[donor_id].add(mutation_info['gene_name'])
    donor_id_to_genes_with_manual_arch_mutation_info = dict(donor_id_to_genes_with_manual_arch_mutation_info) # I don't want a defaultdict moving around.
    return donor_id_to_genes_with_manual_arch_mutation_info

def get_chrom_name_to_cna_infos(cna_infos):
    chrom_name_to_cna_infos = collections.defaultdict(list)
    for cna_info in cna_infos:
        chrom_name_to_cna_infos[cna_info['chrom']].append(cna_info)
    chrom_name_to_cna_infos = dict(chrom_name_to_cna_infos) # I don't want a defaultdict moving around.
    return chrom_name_to_cna_infos

def get_manual_arch_mutation_df():
    flat_dicts = []
    for donor_id, manual_mutation_infos in MDS_ANALYSIS_PARAMS['arch_mutations']['donor_id_to_manual_arch_mutation_infos'].items():
        flat_dicts.extend([
            {
                'donor_id': donor_id,
                **x,
            } for x in manual_mutation_infos
        ])
    df = pd.DataFrame(flat_dicts)
    df[['seq_batch_name', 'CHR', 'POS', 'REF', 'ALT']] = 'dummy'
    return df

def get_disease_cell_state_assignment_info(state_name):
    state_infos = MDS_ANALYSIS_PARAMS['cell_type_assignment']['disease']
    state_infos = [x for x in state_infos if x[0] == state_name]
    assert len(state_infos) == 1
    return state_infos[0][1]

def get_donor_seq_batch_name_and_sample_date_df():
    return pd.DataFrame(
        MDS_ANALYSIS_PARAMS['arch_mutations']['donor_seq_batch_name_and_sample_date'], columns=['donor_id', 'seq_batch_name', 'sample_date'])

# def get_genes_to_ignore_in_projection(atlas_name, mc_ad):
#     # TODO: ignore batchy genes according to our data (e.g., those i identify after differential gene expression analysis controlling for composition)?
#     atlas_params = MDS_ANALYSIS_PARAMS[atlas_name]
#     ignored_gene_names = set(mc_ad.var_names[mc_ad.var[['lateral_gene', 'noisy_gene']].any(axis=1)])
#     if 'max_batch_kruskal_pval_to_ignore_gene_in_projection_to_atlas' in atlas_params:
#         ignored_gene_names |= set(get_nimrod_batchy_genes_by_kruskal_pval(atlas_params['max_batch_kruskal_pval_to_ignore_gene_in_projection_to_atlas']))

#     return sorted(ignored_gene_names)

def get_single_donor_single_batch_mc_model_types_and_names(list_of_exp_name_and_donor_id=None, c_ad=None, c_mask=None):
    only_specific_donors_mc_model_name_prefix = 'final_only_'

    if list_of_exp_name_and_donor_id is None:
        list_of_exp_name_and_donor_id = mc_utils.get_list_of_exp_name_and_donor_id(c_ad, c_mask)
    else:
        assert isinstance(list_of_exp_name_and_donor_id, list)

    return sorted([
        ('final_only_specific_donors_in_specific_exps', f'{only_specific_donors_mc_model_name_prefix}{donor_id}_only_{exp_name}')
        for exp_name, donor_id in list_of_exp_name_and_donor_id
    ])

def get_single_donor_mc_model_types_and_names(list_of_exp_name_and_donor_id=None, c_ad=None):
    only_specific_donors_mc_model_name_prefix = 'final_only_'

    if list_of_exp_name_and_donor_id is None:
        list_of_exp_name_and_donor_id = mc_utils.get_list_of_exp_name_and_donor_id(c_ad)
    
    all_donor_ids = sorted({x[1] for x in list_of_exp_name_and_donor_id})

    return sorted([
        *[
            ('final_only_specific_donors', f'{only_specific_donors_mc_model_name_prefix}{x}') for x in 
            all_donor_ids
        ],
    ])

def get_cell_state_to_color(add_outlier_black_color=False):
    color_df = pd.read_csv(MDS_ANALYSIS_PARAMS['cell_type_colors_csv_file_path'])
    cell_state_to_color = {x['cell_type']: x['color'] for _, x in color_df.iterrows()}
    if add_outlier_black_color:
        cell_state_to_color['Outliers'] = 'black'
    return cell_state_to_color


def get_bio_rep_donor_bleeding_date_pairs(df, donor_id_to_bleeding_dates=None):
    # TODO: better to move this to sc_rna_seq_preprocessing_params, i guess?
    if donor_id_to_bleeding_dates:
        raise NotImplementedError('implement this - so you could choose manually choose for a donor with >2 bleeding dates the dates you want.')
    donor_date_df = df[['donor_id', 'bleeding_date']].drop_duplicates()
    bio_replicate_donor_ids = sorted((donor_date_df['donor_id'].value_counts() > 1).loc[lambda x: x].index)
    donor_date_df = donor_date_df[donor_date_df['donor_id'].isin(bio_replicate_donor_ids)]

    sort_date_key = lambda vec: vec.apply(lambda x: ''.join(x.split('.')[::-1]))
    return generic_utils.merge_preserving_df1_index_and_row_order(
        donor_date_df.sort_values('bleeding_date', ascending=True, key=sort_date_key).drop_duplicates(subset='donor_id'),
        donor_date_df.sort_values('bleeding_date', ascending=False, key=sort_date_key).drop_duplicates(subset='donor_id'),
        on='donor_id', suffixes=('1','2'),
        verbose=False,
    )



def get_traj_mask_names(mask_and_info_list):
    mask_names = [x[0] for x in mask_and_info_list]
    return [x for x in mask_names if x.endswith('_traj')]



# NOTE: pipeline stages
# - choose_excluded_genes.ipynb
# - choose_excluded_cells.ipynb
# - calculate_metacells_and_prepare_for_mcview.ipynb
# - calculate_metacells_and_prepare_for_mcview.ipynb (this time with cells filtered to keep only unhealthy (which implicitly discards soup or vireo doublets and unassigned, but also explicitly discard doublets and other problematic cells we identified))
# - project_to_healthy_blood_aging_atlas.ipynb

def identify_unused_sigs(sig_names, py_file_path='/dummy/dummy/dummy/tanay_group/mds/mds_analysis_params.py'):
    lines = generic_utils.read_text_file(py_file_path).splitlines()
    lines = [x.partition('#')[0] for x in lines]
    all_text = '\n'.join(lines)
    df = pd.DataFrame(
        [(sig_name, all_text.count(f'"{sig_name}"') + all_text.count(f"'{sig_name}'")) for sig_name in sig_names],
        columns=['sig_name', 'count'],
    )
    df.sort_values('count', inplace=True)
    return df

def get_within_state_sig_df():
    final_feature_names_and_types_df = pd.read_csv(MDS_ANALYSIS_PARAMS['final_feature_names_and_types_df_csv_file_path'])
    gene_prog_sig_names = list(MDS_ANALYSIS_PARAMS['pb_sig_name_to_info'])
    within_state_sig_df = pd.DataFrame(
        [(x, y, f'{x}_{y}_0.5') for x, y in itertools.product(gene_prog_sig_names, PB_HSPC_STATE_NAMES)], columns=['sig_name', 'state', 'feature'])
    within_state_sig_df = within_state_sig_df[within_state_sig_df['feature'].isin(final_feature_names_and_types_df['name'])]
    return within_state_sig_df

def get_gene_sig_score_to_name_in_paper():
    return {
        row['feature']: MDS_ANALYSIS_PARAMS['sig_name_to_name_in_paper'][row['sig_name']] + f' across {row["state"]}' 
        for _, row in get_within_state_sig_df().iterrows()
    }



print('mds_analysis_params was loaded/reloaded')