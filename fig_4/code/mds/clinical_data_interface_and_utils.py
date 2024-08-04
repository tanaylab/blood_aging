import pandas as pd
import numpy as np
import os
import re

from sc_rna_seq_preprocessing import sc_rna_seq_preprocessing_params
from mds import mds_analysis_params
from mds import mip_interface_and_utils
from generic import generic_utils

FIRST_ARCH4_DATE = '03.11.21'
FIRST_ARCH4_DATE_AS_DATE = pd.to_datetime(FIRST_ARCH4_DATE, format='%d.%m.%y', errors='raise')

CBC_COL_TO_PUBLISH_COL = {
    'Basophils#': 'Basophils#',
    'Basophils%': 'Basophils%',
    'Eosinophils#': 'Eosinophils#',
    'Eosinophils%': 'Eosinophils%',
    'Hematocrit (%)': 'Hematocrit (%)',
    'Hemoglobin (g/dl)': 'Hemoglobin (g/dl)',
    'Lympho#': 'Lymphocytes#',
    'Lympho%': 'Lymphocytes%',
    'MCH (pg)': 'MCH (pg)',
    'MCHC (g/dl)': 'MCHC (g/dl)',
    'MCV (fL)': 'MCV (fL)',
    'MPV (fL)': 'MPV (fL)',
    'Mono#': 'Monocytes#',
    'Mono%': 'Monocytes%',
    'Neutro#': 'Neutrophils#',
    'Neutro%': 'Neutrophils%',
    'Platelets (10^3/microliter)': 'Platelets (10^3/microliter)',
    'RBC (10^6/microliter)': 'RBC (10^6/microliter)',
    'RDW (%)': 'RDW (%)',
    'WBC (10^3/microliter)': 'WBC (10^3/microliter)',
}

MF_DIAGNOSES = [
    'MF',
    'MDS/CMML MF',
    'post-PV MF',
    'post-ET MF',
    'PMF',
    'PMF-CLL',
]
MF_SUSPECT_DIAGNOSES = [
    'MF?', # TODO: should be here?
    'PMF?', # TODO: should be here?
    'ET or MF', # TODO: should be here?
    'ET/MF',
    'ET/MF?',
    'PV/MF',
    'PV +/- MF',
]

ET_DIAGNOSES = [
    'ET',
]

NON_MPN_CYTOSIS = [
    # NOTE: really not sure about this
    'thrombocytosis',
    'polycythemia',
    'extreme thrombocytosis triple negative',
]

CMML_DIAGNOSES = [
    'CMML',
    'MDS/CMML'
]

DIAGNOSIS_VAL_TO_CORRECTED_DIAGNOSIS_VAL = {
    ' neutropenia': 'neutropenia',
    'ITP ': 'ITP',
    'CMML ': 'CMML',
    'Monocytosis ': 'monocytosis',
    'Myelofibrosis': 'MF',
    'myelofibrosis': 'MF',
    'mf': 'MF',
    'thrombocytosis ': 'thrombocytosis',
    'Polycythemia': 'erythrocytosis', # https://en.wikipedia.org/wiki/Polycythemia
    'polycythemia JAK2 negative': 'erythrocytosis, JAK2-',
    'polycythemia, JAK2-': 'erythrocytosis, JAK2-',
    'MDS/CMML Myelofibrosis': 'MDS/CMML MF',
    'AML ': 'AML',
    'PV ': 'PV',
    'pv': 'PV',
    'et': 'ET',
    'mds': 'MDS',
    'cmml': 'CMML',
    'thrombocytopenia neutropenia': 'neutropenia thrombocytopenia',
    'et_vs_mf': 'ET/MF',
    'pmf': 'PMF',
    'pmf_?': 'PMF?',
    'mpn_leukocytosis': 'MPN leukocytosis',
    'post pv MF': 'post-PV MF',
    'post_pv_mf': 'post-PV MF',
    'post-PV myelofibrosis': 'post-PV MF',
    'Copper deficiency?': 'copper deficiency?',
    'Aplastic Anemia': 'aplastic anemia', 
    'Aplastic anemia': 'aplastic anemia', 
    'Aplastic Anemia no ARCH mutations': 'aplastic anemia no ARCH mutations',
    'MDS no ARCH mutations': 'MDS no ARCH mutations',
    'PV/myelofibrosis': 'PV/MF',
    'ET/myelofibrosis': 'ET/MF',
    'ET/myelofibrosis?': 'ET/MF?',
    'ET or myelofibrosis': 'ET/MF',
    'myelofibrosis?': 'MF?',
    'primary myelofibrosis': 'PMF',
    'primary myelofibrosis?': 'PMF?',
    
    'healthy': 'normal', 
    'MDS/MPN overlap': 'MDS/MPN',
    
    'Chronic Neutrophilic Leukemia? no ARCH mutations': 'Chronic Neutrophilic Leukemia?',

    'post ET-MF': 'post-ET MF', 
    'pmf-cll': 'PMF-CLL',
    'healthy 24 hours': 'normal_24h',
    'healthy 36 hours': 'normal_36h',
    'healthy 48 hours': 'normal_48h',
    'pv ': 'PV',
}

MDS_AND_CYTOPENIA_RELATED_DIAGNOSES = [
    'monocytosis',
]
MPN_RELATED_DIAGNOSIS_CLASSES = [
    'thrombo_or_erythro_cytosis',
]


DIAGNOSIS_TO_DIAGNOSIS_CLASS = {
    'MDS?': 'MDS', # TODO: should be here?
    'MDS with excess blasts': 'MDS',
    'MDS no ARCH mutations': 'MDS',
    'hypoplastic MDS': 'MDS',
    'MDS (SF3B1, DNMT3A)': 'MDS',
    'MDS? awaiting sequencing and BM biopsy results': 'MDS',
    'hypoplastic MDS?': 'MDS',
    'MDS/smoldering MM': 'MDS', 
    

    'CMML': 'MDS/MPN',
    'CMML?': 'MDS/MPN',
    'MDS/MPN': 'MDS/MPN',
    'MDS/CMML': 'MDS/MPN',
    'MDS/CMML MF': 'MDS/MPN',

    'anemia': 'cytopenia',
    'neutropenia': 'cytopenia',
    'thrombocytopenia': 'cytopenia',
    'anemia thrombocytopenia': 'cytopenia',
    'pancytopenia': 'cytopenia',
    'neutropenia thrombocytopenia': 'cytopenia',
    'ITP': 'cytopenia', # TODO: should be here?
    'anemia thrombocytopenia thalassemia': 'cytopenia', # TODO: should be here?
    'anemia monocytosis': 'cytopenia', # TODO: should be here?
    'ICUS': 'cytopenia',
    'cytopenia?': 'cytopenia',
    'cytopenia': 'cytopenia',
    'cytopenia (MF?)': 'cytopenia', # here because the patient appeared in 240223_MDS_clinical_data_new.xlsx.
    
    'MPN': 'MPN',
    'MF': 'MPN',
    'MF?': 'MPN',
    'PV': 'MPN',
    'ET': 'MPN',
    'ET-preMF': 'MPN',
    'PMF': 'MPN',
    'PMF?': 'MPN',
    'MPN leukocytosis': 'MPN',
    'post-PV MF': 'MPN',
    'post-ET MF': 'MPN',
    'JAK2+, post PV MF?': 'MPN',
    'PV/MF': 'MPN',
    'ET/MF': 'MPN',
    'ET/MF?': 'MPN',
    'PV +/- MF, JAK2': 'MPN',
    'PV +/- MF': 'MPN',
    'PMF-CLL': 'MPN',
    
    'extreme thrombocytosis triple negative': 'thrombo_or_erythro_cytosis',
    'thrombocytosis': 'thrombo_or_erythro_cytosis',
    'thrombocytosis no mutations': 'thrombo_or_erythro_cytosis',
    'erythrocytosis': 'thrombo_or_erythro_cytosis',
    'erythrocytosis, JAK2-': 'thrombo_or_erythro_cytosis',

    'MDS/AML': 'AML',    
    'AML on vidaza/venetoclax': 'AML',    

    'normal': 'normal',
    'normal?' : 'normal',
    'AML': 'AML',
    'MDS': 'MDS',

    't cell disease?': 'other_disease',
    'APLA': 'other_disease',
    'monocytosis': 'other_disease', # TODO: should be here?
    'leukocytosis': 'other_disease', # TODO: should be here?
    'eosinophilia': 'other_disease', # TODO: should be here?
    'hyper-eosinophilic syndrome': 'other_disease', # TODO: should be here?
    'mastocytosis': 'other_disease', # TODO: should be here?
    'CML stop therapy': 'other_disease', # TODO: should be here?
    'some abnormally high counts': 'other_disease', # TODO: should be here?
    'abnormally high counts': 'other_disease', # TODO: should be here?
    
    'copper deficiency?': 'other_disease',
    'AML-rem': 'other_disease',
    'JMML': 'other_disease',
    'aplastic anemia': 'other_disease', 
    'aplastic anemia no ARCH mutations': 'other_disease', 
    'MM': 'other_disease',
    'liver disease': 'other_disease',
    'ICF': 'other_disease',
    'not_MDS': 'other_disease',
    'Amyloidosis': 'other_disease',
    'Chronic Neutrophilic Leukemia': 'other_disease',
    'Chronic Neutrophilic Leukemia?': 'other_disease',
    'Hairy Cell Leukemia': 'other_disease',
    'mature T cell leukemia': 'other_disease',
    'Multiple Myeloma': 'other_disease',
    '?': 'other_disease',
    
    'normal_24h': 'other_exp_condition',
    'normal_36h': 'other_exp_condition',
    'normal_48h': 'other_exp_condition',
    'normal_other_exp_condition': 'other_exp_condition',
    'other_exp_condition': 'other_exp_condition',
    'other_exp_condition_excluded': 'other_exp_condition',
    
    'nan': 'nan',
}

MDS_AND_MDS_MPN_AND_CYTOPENIA_DIAGNOSIS_CLASSES = ['cytopenia', 'ICUS', 'CCUS', 'MDS', 'MDS/MPN']

DIAGNOSIS_CLASS_TO_COLOR = {
    'AML': 'black',
    'MDS': 'tab:red',
    'MDS/MPN': 'tab:olive',
    'MPN': 'tab:purple',
    'normal': 'tab:green',
    'cytopenia': 'tab:blue',
    'CCUS': 'tab:brown',
    'ICUS': 'tab:cyan',
    'thrombo_or_erythro_cytosis': 'tab:pink',
    
    'PV': 'tab:red',
    'ET': 'tab:blue',
    'MF': 'tab:orange',
    'MF?': 'tab:cyan',
    
    'other_disease': '#A0A0A0',
    
    '?': '#A0A0A0',
    'nan': '#E0E0E0',
    'other_exp_condition': '#A0A0A0',
}


def get_sc_rna_seq_preprocessing_params():
    return sc_rna_seq_preprocessing_params.SC_RNA_SEQ_PREPROCESSING_PARAMS

def get_mds_params():
    return mds_analysis_params.MDS_ANALYSIS_PARAMS


def get_technical_replicates_df():
    source_file_path = get_sc_rna_seq_preprocessing_params()['donor_table_paths']['technical_replicates_xlsx']
    df = pd.read_excel(source_file_path)
    df['source_table_file_name'] = os.path.basename(source_file_path)

    final_rows = []
    for _, row in df.iterrows():
        if '/' in row['date']:
            row_with_date1 = row.copy()
            date1 = row['date'].split('/')[0]
            row_with_date1['date'] = date1
            row_with_date2 = row.copy()
            date2 = row['date'].split('/')[1]
            row_with_date2['date'] = date2
            
            # year1 = int(date1.rpartition('_')[-1])
            # year2 = int(date2.rpartition('_')[-1])
            # if year1 > year2:
            # nah. we just use 230501_MDS_10x_age_and_gender.xlsx, which has the ages of the donors in later experiments...
            # final_rows.extend([row_with_date1, row_with_date2])
            final_rows.extend([row_with_date1]) # 230614: i guess the second date always appears in the clinical data xlsx? if we also add the second date, then we will have duplicates when concatenating to the clinical data xlsx...
        else:
            final_rows.append(row)

    df = pd.DataFrame(final_rows)
    for num in (1,2,3):
        suffix = f'_p{num}'
        px_mask = df['date'].str.endswith(suffix)
        if px_mask.any():
            df.loc[px_mask, 'exp_num_in_bleeding_date'] = str(num)
            df.loc[px_mask, 'date'] = df.loc[px_mask, 'date'].str.partition(suffix, expand=True)[0]

    assert df['num'].str.fullmatch('N[0-9]+', na=False).all()

    # df['donor_cohort'] = 'healthy' # 
    # df.loc[df['num'].isin(['N180', 'N78']), 'donor_cohort'] = np.nan
    df.rename(columns={
        'date': 'exp_date', 
        'cohort': 'diagnosis',
        'num': 'donor_id', 
        'gender': 'donor_sex', 
        'age': 'donor_age'
    }, inplace=True)
    df['diagnosis'].replace({'healthy': 'normal'}, inplace=True)
    assert set(df['diagnosis'].unique()) <= {'normal', 'MDS'}

    return df

def get_prev_donor_df():
    source_file_path = get_sc_rna_seq_preprocessing_params()['donor_table_paths']['prev_donor_id_xlsx']
    df = pd.read_excel(source_file_path)
    return df.rename(columns={'num': 'donor_id', 'previous_num': 'prev_donor_id'})[['donor_id', 'prev_donor_id']]

def fix_and_verify_nili_longitudinal_cbc_df(df):
    new_donor_id_col = []
    curr_donor_id = None
    for _, row in df.iterrows():
        donor_id = row['10x_number']
        if pd.isna(donor_id):
            assert curr_donor_id is not None
        else:
            assert curr_donor_id != donor_id, f'curr_donor_id == donor_id: {curr_donor_id}, {donor_id}'
            curr_donor_id = donor_id
        new_donor_id_col.append(curr_donor_id)
    df.insert(0, 'donor_id', pd.Series(new_donor_id_col).astype(str))
    df.drop(columns='10x_number', inplace=True)

    df['CBC dates'] = df['CBC dates'].astype(str)
    for donor_id in sorted(df['donor_id'].unique()):
        if donor_id in [
            
            # 
        ]:
            continue
        # print(f'verifying {donor_id} CBC dates are as valid and ordered')
        curr_df = df[df['donor_id'] == donor_id]
        invalid_date_mask = curr_df['CBC dates'].str.count(r'\.') != 2
        if invalid_date_mask.any():
            print(f'donor_id: {donor_id}, skipping invalid dates: {curr_df.loc[invalid_date_mask, "CBC dates"]}')
            curr_df = curr_df[~invalid_date_mask]
        if curr_df.empty:
            continue
        curr_dates_df = curr_df['CBC dates'].str.split('.', expand=True).astype(int)
        before_2000_mask = curr_dates_df[2] > 40
        curr_dates_df.loc[before_2000_mask, 2] = curr_dates_df.loc[before_2000_mask, 2] + 1900
        curr_dates_df.loc[~before_2000_mask, 2] = curr_dates_df.loc[~before_2000_mask, 2] + 2000

        # print(curr_dates_df)
        assert (curr_dates_df.index == curr_dates_df.sort_values([2, 1, 0], ascending=False).index).all()
    return df

def get_nili_cbc_df():
    source_file_path = get_sc_rna_seq_preprocessing_params()['donor_table_paths']['nili_clinical_data_xlsx']
    cbc_df = pd.read_excel(source_file_path, sheet_name='longitudinal CBCs')
    cbc_df = fix_and_verify_nili_longitudinal_cbc_df(cbc_df)
    cbc_df.rename(columns={
        '10x_number': 'donor_id', 'sampling date': 'exp_date', 'Diagnosis': 'diagnosis',
    }, inplace=True)
    cbc_df['donor_id'].replace(get_sc_rna_seq_preprocessing_params()['other_donor_id_to_donor_id_i_use'], inplace=True)
    cbc_df['cbc_date'] = pd.to_datetime(cbc_df['CBC dates'], format='%d.%m.%y', errors='raise')
    cbc_df.drop_duplicates(inplace=True)

    return cbc_df

NOT_TREATED_TREATMENT_VALUES = [
    'not treated',
    'not treated (previously treated with Aranesp - no response)',
]
UNKNOWN_TREATMENT_VALUES = [
    'unclear', 
    'Unknown',
]
DRUG_TREATMENT_VALUES = [
    'Aranesp (EPO)',
    'EPO',
    'Plaquenil',
    'Prednisone (hemolytic anemia)',
    'Hydrea',
    'Revlimid',
    'vidaza+venetoclax (from Jan. 2022), Aranesp',
]

TREATMENT_VALUES = [
    *UNKNOWN_TREATMENT_VALUES,
    *NOT_TREATED_TREATMENT_VALUES,
    *DRUG_TREATMENT_VALUES,
]
def get_nili_treatment_df():
    source_file_path = get_sc_rna_seq_preprocessing_params()['donor_table_paths']['nili_treatment_xlsx']
    treatment_df = pd.read_excel(source_file_path)
    treatment_df = treatment_df[['donor_id', 'bleeding_date', 'treatment status']]
    treatment_df.rename(columns={'treatment status': 'treatment'}, inplace=True)
    unexpected_treatment_values = set(treatment_df['treatment']) - set(TREATMENT_VALUES)
    assert generic_utils.set_is_empty_or_contains_only_np_nan(unexpected_treatment_values), f'unexpected_treatment_values: {unexpected_treatment_values}'
    treatment_df['known_treatment'] = treatment_df['treatment'].isin(DRUG_TREATMENT_VALUES)
    return treatment_df

def identify_potentially_abnormal_cbcs():
    clinical_df = get_clinical_data_df()
    cbc_df = pd.read_csv(get_sc_rna_seq_preprocessing_params()['donor_table_paths']['minimal_cbc_df_csv'])
    # cbc_df.to_csv('temp/cbc/cbc_for_nili.csv', index=False)

    cbc_df_with_sex = generic_utils.merge_preserving_df1_index_and_row_order(cbc_df, clinical_df[['donor_id', 'donor_sex']].drop_duplicates(), how='left')
    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag = cbc_df_with_sex[
        ~cbc_df_with_sex['donor_id'].isin(clinical_df.loc[~clinical_df['diagnosis'].astype(str).isin(['normal', 'nan']), 'donor_id'])]

    missing_from_clinical_df_donor_ids = sorted(cbc_df_with_sex[cbc_df_with_sex['donor_sex'].isna()]['donor_id'].unique())
    if missing_from_clinical_df_donor_ids:
        raise RuntimeError(f'missing_from_clinical_df_donor_ids: {missing_from_clinical_df_donor_ids}')

    cbc_df_with_sex.loc[(cbc_df_with_sex['donor_sex'] == 'male') & (cbc_df_with_sex['Hemoglobin (g/dl)'] < 13.5), ['donor_id', 'CBC dates']].to_csv(
        'all_male_hgb_under_13.5.csv', index=False)
    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
        (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['donor_sex'] == 'male') & (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Hemoglobin (g/dl)'] < 13.5),
        ['donor_id', 'CBC dates']
    ].to_csv('temp/cbc/filtered_male_hgb_under_13.5.csv', index=False)

    cbc_df_with_sex.loc[(cbc_df_with_sex['donor_sex'] == 'female') & (cbc_df_with_sex['Hemoglobin (g/dl)'] < 12), ['donor_id', 'CBC dates']].to_csv(
        'all_female_hgb_under_12.csv', index=False)
    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
        (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['donor_sex'] == 'female') & (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Hemoglobin (g/dl)'] < 12),
        ['donor_id', 'CBC dates']
    ].to_csv('temp/cbc/filtered_female_hgb_under_12.csv', index=False)


    cbc_df_with_sex.loc[cbc_df_with_sex['Platelets (10^3/microliter)'] < 150, ['donor_id', 'CBC dates']].to_csv('temp/cbc/all_plt_under_150.csv', index=False)
    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
        cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Platelets (10^3/microliter)'] < 150,
        ['donor_id', 'CBC dates']
    ].to_csv('temp/cbc/filtered_plt_under_150.csv', index=False)


    cbc_df_with_sex.loc[cbc_df_with_sex['Neutro#'] < 1.8, ['donor_id', 'CBC dates']].to_csv('temp/cbc/all_neut_under_1.8.csv', index=False)
    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
        cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Neutro#'] < 1.8,
        ['donor_id', 'CBC dates']
    ].to_csv('temp/cbc/filtered_neut_under_1.8.csv', index=False)

    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
        (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['donor_sex'] == 'male') & (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Hemoglobin (g/dl)'] > 17.5),
        ['donor_id', 'CBC dates']
    ].to_csv('temp/cbc/filtered_male_hgb_above_17.5.csv', index=False)
    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
        (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['donor_sex'] == 'female') & (cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Hemoglobin (g/dl)'] > 16),
        ['donor_id', 'CBC dates']
    ].to_csv('temp/cbc/filtered_female_hgb_above_16.csv', index=False)

    cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
        cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Platelets (10^3/microliter)'] > 400,
        ['donor_id', 'CBC dates']
    ].to_csv('temp/cbc/filtered_plt_above_400.csv', index=False)

    for thresh in (6, 6.5, 7, 7.5, 8):
        cbc_df_with_sex_only_without_any_non_normal_or_nan_diag.loc[
            cbc_df_with_sex_only_without_any_non_normal_or_nan_diag['Neutro#'] > thresh,
            ['donor_id', 'CBC dates']
        ].to_csv(f'temp/cbc/filtered_neut_above_{thresh}.csv', index=False)

def fix_nili_exp_date(df):
    # df['exp_date'] = df['exp_date'].astype(str).str.replace('/', '.', regex=False)
    df['exp_date'].replace({
        # NOTE: this is before removing p1/p2/p3 from the exp_date, so what's here should contain all p1/p2/p3 parts of the original date...
        '4.12.23': '04.12.23',
        '5.12.23': '05.12.23',
        '06.06.22 PB . 06.10.22 BM': '06.06.22', # ah?? in excel it looks like "06.06.22 PB / 06.10.22 BM". i notice i am confused. ugh.
        '1.1.24': '01.01.24', # TODO: remove after this is fixed.
        '21.10.20/07.12.20_p2': '07.12.20',
        '30.11.20_p2/07.02.22': '30.11.20', # N78
        '30.11.20_p1/07.12.20_p1/03.07.22': '07.12.20', 
        '30.11.20_p2/30.01.22': '30.11.20', # 'N80'
        '30.11.20_p1/11.04.21_p1': '11.04.21', # 'N83'
        '07.12.20_p2/16.01.22': '07.12.20', # 'N88'
        '01.02.21_p2/14.02.22': '01.02.21', # 'N91'
        '21.12.20_p2/30.01.22': '21.12.20', # 'N98'
        '28.12.20_p1/14.02.22': '28.12.20', # 'N99'
        '01.02.21_p1/07.02.22': '01.02.21', # 'N131'
        '21.02.21_p2/28.02.22': '21.02.21', # 'N150'
        '22.02.21_p1/30.01.22': '22.02.21', # 'N157'
        '22.02.21_p2/14.02.22': '22.02.21', # 'N159'
        '01.03.21/07.02.22': '01.03.21', # 'N164'
        '01.03.21/07.02.22': '01.03.21', # 'N165'
        '07.03.21_p2/14.02.22': '07.03.21', # 'N175'
        '11.04.21_p1/03.03.22': '11.04.21', # 'N179'
        '11.04.21_p2/10.01.22': '11.04.21', # 'N180'
        '04.01.21_p1/02.03.22': '04.01.21', # 'N109'
        '11.01.21_p2/06.06.22': '11.01.21', # 'N122'
    }, inplace=True)

    df['exp_date'] = df['exp_date'].astype(str).str.replace(r'_p1$', '', regex=True)
    df['exp_date'] = df['exp_date'].astype(str).str.replace(r'_p2$', '', regex=True)
    df['exp_date'] = df['exp_date'].astype(str).str.replace(r'_p3$', '', regex=True)

def read_nili_clinical_data_xlsx():
    source_file_path = get_sc_rna_seq_preprocessing_params()['donor_table_paths']['nili_clinical_data_xlsx']
    else_df = pd.read_excel(source_file_path, sheet_name='all else')
    

    if 0:
        import openpyxl
        wb = openpyxl.load_workbook(source_file_path) 
        print(wb.sheetnames)

    assert 'Sampling date' in else_df.columns
    else_df.rename(columns={'Sampling date': 'Exp. Date'}, inplace=True)


    cytopenic_df = pd.read_excel(source_file_path)
    

    else_df = else_df[~else_df['Sample num'].isin([
        'N151', 'N130', 'N181', 'N116',

        'N87', 'N107', 'N81', 'N149', 'N59', 'N85', 'N82', 'N86',
    ])]

    donor_ids_in_both_cytopenic_and_else = set(else_df['Sample num']) & set(cytopenic_df['Sample num'])
    assert not donor_ids_in_both_cytopenic_and_else, f'donor_ids_in_both_cytopenic_and_else: {donor_ids_in_both_cytopenic_and_else}' # TODO: uncomment

    df = pd.concat([cytopenic_df, else_df], ignore_index=True)
    # assert df['Sample num'].is_unique # why was this needed? N48 doesn't satisfy this...
    df.rename(columns={
        'Sample num': 'donor_id', 'Previous num': 'prev_donor_id', 'Exp. Date': 'exp_date', 'Age': 'donor_age', 'Gender': 'donor_sex', 'Diagnosis': 'diagnosis',
    }, inplace=True)
    

    fix_nili_exp_date(df)
    
    if 0:
        df = generic_utils.merge_preserving_df1_index_and_row_order(df, cbc_df, on=['donor_id', 'exp_date'], how='left', indicator=True)
        df['same_date_cbc'] = df['_merge'] == 'both'
        df.drop(columns='_merge', inplace=True)
    

    df['prev_donor_id'].replace({
        'N180, N200': 'N200, N180',
        'No cells': np.nan,
    }, inplace=True)

    df = generic_utils.merge_preserving_df1_index_and_row_order(
        df, get_prev_donor_df().rename(columns={'prev_donor_id': 'new_prev_donor_id'}), how='left')
    
    new_previous_num_mask = df['prev_donor_id'].isna() & (~df['new_prev_donor_id'].isna())
    df.loc[new_previous_num_mask, 'prev_donor_id'] = df.loc[new_previous_num_mask, 'new_prev_donor_id']

    new_different_previous_num_mask = (~df['prev_donor_id'].isna()) & (~df['new_prev_donor_id'].isna()) & (df['prev_donor_id'] != df['new_prev_donor_id'])
    problematic_df = df.loc[new_different_previous_num_mask, ['donor_id', 'prev_donor_id', 'new_prev_donor_id']]
    assert not new_different_previous_num_mask.any(), problematic_df
    df.drop(columns='new_prev_donor_id', inplace=True)

    df['source_table_file_path'] = source_file_path
    
    blast_df = pd.read_excel(source_file_path, sheet_name='FACS blasts', header=None, names=['donor_id', 'FACS BM blasts (%)'])
    donor_ids_with_blasts_in_blast_sheet = set(blast_df.loc[blast_df['FACS BM blasts (%)'].notna(), 'donor_id'])
    donor_ids_with_blasts_in_cytopenic_sheet = set(blast_df.loc[blast_df['FACS BM blasts (%)'].notna(), 'donor_id'])
    assert not (donor_ids_with_blasts_in_cytopenic_sheet - donor_ids_with_blasts_in_blast_sheet), str((donor_ids_with_blasts_in_cytopenic_sheet - donor_ids_with_blasts_in_blast_sheet))
    assert not (donor_ids_with_blasts_in_blast_sheet - donor_ids_with_blasts_in_cytopenic_sheet), str((donor_ids_with_blasts_in_blast_sheet - donor_ids_with_blasts_in_cytopenic_sheet))

    return df

def read_gal_dadi_clinical_data_xlsx():
    source_file_path = get_sc_rna_seq_preprocessing_params()['donor_table_paths']['gal_dadi_clinical_data_xlsx']
    df1 = pd.read_excel(source_file_path)
    # df2 = pd.read_excel(source_file_path, sheet_name='CBC')
    # assert df1['num'].is_unique # G8 cells were used again to see whether frozen cells can be used...
    # assert df2['num'].is_unique #

    # df1_donor_id_bleeding_dates = df1[['num', 'exp_date']].to_records(index=False).tolist()
    # df2_donor_id_bleeding_dates = df2[['num', 'date']].to_records(index=False).tolist()
    # assert len(set(df1_donor_id_bleeding_dates)) == len(df1_donor_id_bleeding_dates) # G8 cells from 13/07/2023 were used again to see whether frozen cells can be used...
    # assert len(set(df2_donor_id_bleeding_dates)) == len(df2_donor_id_bleeding_dates)
    # assert set(df1_donor_id_bleeding_dates) == set(df2_donor_id_bleeding_dates)
    # df = generic_utils.merge_preserving_df1_index_and_row_order(df1, df2.drop(columns='diagnosis'), on=['date', 'num'])
    df = df1
    df['source_table_file_path'] = source_file_path
    df.rename(columns={
        'num': 'donor_id', 
        'previous_num': 'prev_donor_id', 
        'gender': 'donor_sex', 
        'age': 'donor_age'
    }, inplace=True)

    df = generic_utils.merge_preserving_df1_index_and_row_order(df, get_prev_donor_df().rename(columns={
        'prev_donor_id': 'new_prev_donor_id'}), how='left')

    new_previous_num_mask = df['prev_donor_id'].isna() & (~df['new_prev_donor_id'].isna())
    df.loc[new_previous_num_mask, 'prev_donor_id'] = df.loc[new_previous_num_mask, 'new_prev_donor_id']

    new_different_previous_num_mask = (~df['prev_donor_id'].isna()) & (~df['new_prev_donor_id'].isna()) & (df['prev_donor_id'] != df['new_prev_donor_id'])
    assert not new_different_previous_num_mask.any()
    df.drop(columns='new_prev_donor_id', inplace=True)


    df['exp_date'] = df['exp_date'].astype(str).str.replace('/', '.', regex=False)
    df['exp_date'] = df['exp_date'].astype(str).str.replace(r'\.2023$', r'.23', regex=True)
    df['exp_date'] = df['exp_date'].astype(str).str.replace(r'\.2024$', r'.24', regex=True)
    
    return df

def validate_exp_date_col(df):
    date_pattern = '[0-3][0-9]\.[01][0-9]\.[0-9][0-9]'
    invalid_exp_date_mask = df['exp_date'].apply(lambda x: re.fullmatch(date_pattern, x) is None)
    assert not invalid_exp_date_mask.any(), f'invalid exp dates:\n{df.loc[invalid_exp_date_mask, ["donor_id", "exp_date"]]}'

def get_nili_shlush_donor_id_and_exp_date_df_for_merging_with_arch_mut_dfs():
    donor_id_and_exp_date_df = read_nili_clinical_data_xlsx()[['donor_id', 'exp_date']]
    
    
    len_before = len(donor_id_and_exp_date_df)
    donor_id_and_exp_date_df = donor_id_and_exp_date_df[~(
        (donor_id_and_exp_date_df['donor_id'] == 'N48') &
        (donor_id_and_exp_date_df['exp_date'] == '22.02.21')
    )]
    assert len(donor_id_and_exp_date_df) == (len_before - 1)
    
    return donor_id_and_exp_date_df

def get_gal_dadi_shlush_donor_id_and_exp_date_df_for_merging_with_arch_mut_dfs():
    donor_id_and_exp_date_df = read_gal_dadi_clinical_data_xlsx()[['donor_id', 'exp_date', 'Notes']]
    
    len_before = len(donor_id_and_exp_date_df)
    donor_id_and_exp_date_df = donor_id_and_exp_date_df[
        (donor_id_and_exp_date_df['Notes'] != 'BM - no data as qc1 failed')
    ]
    assert len(donor_id_and_exp_date_df) == (len_before - 2)
    
    return donor_id_and_exp_date_df[['donor_id', 'exp_date']]

def get_gal_dadi_and_nili_shlush_donor_id_and_exp_date_df_for_merging_with_arch_mut_dfs():
    return pd.concat([
        get_gal_dadi_shlush_donor_id_and_exp_date_df_for_merging_with_arch_mut_dfs(),
        get_nili_shlush_donor_id_and_exp_date_df_for_merging_with_arch_mut_dfs(),
    ], ignore_index=True)

def get_arch4_donor_id_and_exp_date_df():
    df = get_gal_dadi_and_nili_shlush_donor_id_and_exp_date_df_for_merging_with_arch_mut_dfs()
    orig_donor_ids = df['donor_id'].drop_duplicates()
    arch4_orig_donor_ids = [x for x in orig_donor_ids if x.startswith('G') or (x.startswith('N') and int(x[1:]) >= 185)]
    df = df[df['donor_id'].isin(arch4_orig_donor_ids)].drop_duplicates()
    
    special_case_df = pd.DataFrame([(x, FIRST_ARCH4_DATE) for x in ('N12', 'N48', 'N60')], columns=['donor_id', 'exp_date'])
    df = pd.concat([df, special_case_df], ignore_index=True)
    df['donor_id'].replace(get_sc_rna_seq_preprocessing_params()['other_donor_id_to_donor_id_i_use'], inplace=True)
    df.drop_duplicates(inplace=True)
    return df

def get_minimal_clinical_data_df():
    return pd.read_csv(get_sc_rna_seq_preprocessing_params()['donor_table_paths']['nili_minimal_clinical_df_csv'])

    
def get_clinical_data_df():
    df = read_nili_clinical_data_xlsx()

    # df.loc[new_previous_num_mask, ['num', 'prev_donor_id', 'new_prev_donor_id']]
    # df.loc[new_different_previous_num_mask, ['num', 'prev_donor_id', 'new_prev_donor_id']]

    # assert (df['age'] == df['age'].astype(int)).all() # this is not true, so we use floats. (one donor is 2.5 years old)
    # df['age'] = df['age'].astype(int)

    df = pd.concat([df, read_gal_dadi_clinical_data_xlsx()], ignore_index=True)

    validate_exp_date_col(df)

    sc_rna_seq_preprocessing_params.verify_donor_id_couples_in_other_donor_id_to_donor_id_i_use(df, 'donor_id', 'prev_donor_id')
    df['donor_id'].replace(get_sc_rna_seq_preprocessing_params()['other_donor_id_to_donor_id_i_use'], inplace=True)

    for (donor_id, exp_date), forced_diagnosis_and_sex in get_sc_rna_seq_preprocessing_params()['clinical_table_donor_id_and_exp_date_to_forced_diagnosis_and_sex'].items():
        mask = (df['donor_id'] == donor_id) & (df['exp_date'] == exp_date)
        num_of_rows = mask.sum()
        assert num_of_rows <= 1
        if num_of_rows == 1:
            df.loc[mask, 'diagnosis'] = forced_diagnosis_and_sex[0]
            df.loc[mask, 'donor_sex'] = forced_diagnosis_and_sex[1]
        else:
            assert num_of_rows == 0
            print(f'WARNING: did not find {donor_id} {exp_date} in clinical data xlsx. adding a row for them in clinical_data_df with attrs obtained from clinical_table_donor_id_and_exp_date_to_forced_diagnosis_and_sex. it would be really better to add them to the clinical data xlsx.')

            df = pd.concat([df, pd.DataFrame([{
                'donor_id': donor_id, 'exp_date': exp_date, 
                'diagnosis': forced_diagnosis_and_sex[0], 'donor_sex': forced_diagnosis_and_sex[1]}])], ignore_index=True)

    df['diagnosis'] = df['diagnosis'].astype(str)

    multiple_rows_for_same_donor_exp_date_mask = df[['donor_id', 'exp_date']].duplicated(keep=False)
    collapsed_dfs = []
    for donor_id_and_exp_date in df.loc[multiple_rows_for_same_donor_exp_date_mask, ['donor_id', 'exp_date']].drop_duplicates().to_records(index=False):
        # print(donor_id_and_exp_date)
        # print(df[df[['donor_id', 'exp_date']].to_records(index=False) == donor_id_and_exp_date])
        collapsed_dfs.append(generic_utils.collapse_df_rows(
            df[df[['donor_id', 'exp_date']].to_records(index=False) == donor_id_and_exp_date], 
            allow_conflicts_cols=[
                'source_table_file_path',
                'pool',
                'Notes',
                'prev_donor_id',
            ],
        ))
    # print(collapsed_dfs)

    df = pd.concat([
        df[~multiple_rows_for_same_donor_exp_date_mask],
        *collapsed_dfs,
    ], ignore_index=True)

    df['diagnosis'].replace(DIAGNOSIS_VAL_TO_CORRECTED_DIAGNOSIS_VAL, inplace=True)

    unexpected_diagnosis_vals = set(df['diagnosis'].unique()) - set(DIAGNOSIS_TO_DIAGNOSIS_CLASS)
    assert not unexpected_diagnosis_vals, f'unexpected_diagnosis_vals: {unexpected_diagnosis_vals}'

    # df['diagnosis_class'] = df['diagnosis'].replace(DIAGNOSIS_TO_DIAGNOSIS_CLASS) # NOTE: commenting out on purpose, because we might overwrite diagnosis by exp_name and donor_id, so diagnosis_class should be determined only after that.
    # df.loc[~df['diagnosis_class'].isin(['normal', 'MDS', 'cytopenia', 'MPN', 'AML']), 'diagnosis_class'] = 'other_unhealthy'

    # df['diagnosis'].value_counts()
    # df['diagnosis_class'].value_counts()
    return df

def add_donor_closest_cbc_cols(df):
    assert 'donor_id' in df.columns
    assert 'exp_name' in df.columns
    assert 'bleeding_date' in df.columns

    assert not df[['donor_id', 'exp_name']].duplicated().any()

    orig_df_len = len(df)
    cbc_df = pd.read_csv(get_sc_rna_seq_preprocessing_params()['donor_table_paths']['minimal_cbc_df_csv'])
    df = df.merge(cbc_df.drop(columns='exp_date'), on='donor_id', how='left')
    df['cbc_10x_sample_days_diff'] = (
        generic_utils.convert_two_dot_two_dot_two_date_col_to_days_since_ref_date(df, 'CBC dates', (2020, 1, 1), skip_na_vals=True) -
        generic_utils.convert_two_dot_two_dot_two_date_col_to_days_since_ref_date(df, 'bleeding_date', (2020, 1, 1), skip_na_vals=False)
    )
    df['cbc_10x_sample_days_dist'] = df['cbc_10x_sample_days_diff'].abs()
    df = df.sort_values('cbc_10x_sample_days_dist', ascending=True).drop_duplicates(subset=['donor_id', 'exp_name'], keep='first')
    assert len(df) == orig_df_len, f'{len(df)}, {orig_df_len}'
    return df

def get_bm_blast_list_of_date_and_frac(raw_blast_str):
    list_of_date_and_percent_repr = raw_blast_str.split(', ')
    list_of_date_and_frac = []
    for date_and_percent_repr in list_of_date_and_percent_repr:
        date_and_percent_repr_split = date_and_percent_repr.split(' - ')
        if len(date_and_percent_repr_split) != 2:
            raise ValueError(f'len(date_and_percent_repr_split) != 2. date_and_percent_repr_split: {date_and_percent_repr_split}')
        date_repr, percent_repr = date_and_percent_repr_split
        assert percent_repr.endswith('%')
        frac = float(percent_repr[:-1]) / 100
        list_of_date_and_frac.append((date_repr, frac))
    return list_of_date_and_frac

def convert_raw_bm_facs_blast_to_mean_frac(raw):
    if raw == 'nan':
        return np.nan
    list_of_date_and_frac = get_bm_blast_list_of_date_and_frac(raw)
    mean_frac = np.mean([x[1] for x in list_of_date_and_frac])
    return mean_frac

def get_blast_df(clinical_df, epsilon_for_facs_blast_log):
    curr_df = clinical_df.loc[clinical_df['FACS BM blasts (%)'].notna(), ['donor_id', 'FACS BM blasts (%)']].copy()
    flat_dicts = []
    for _, row in curr_df.iterrows():
        donor_id = row['donor_id']
        list_of_date_and_frac = get_bm_blast_list_of_date_and_frac(row['FACS BM blasts (%)'])
        flat_dicts.extend([
            {'donor_id': donor_id, 'FACS_BM_date_raw': date, 'bm_FACS_blast_frac': frac}
            for date, frac in list_of_date_and_frac
        ])
    blast_df = pd.DataFrame(flat_dicts)
    blast_df['bm_FACS_blast_frac'] = blast_df['bm_FACS_blast_frac'].astype(float)
    blast_df['log_bm_FACS_blast_frac'] = np.log2(blast_df['bm_FACS_blast_frac'] + epsilon_for_facs_blast_log)
    year_only_mask = ~blast_df['FACS_BM_date_raw'].str.contains('.', regex=False)
    blast_df.loc[year_only_mask, 'FACS_BM_year'] = blast_df.loc[year_only_mask, 'FACS_BM_date_raw'].astype(int)
    assert blast_df.loc[year_only_mask, 'FACS_BM_year'].between(2000, 2040).all()
    
    # fake middle of the year
    blast_df.loc[year_only_mask, 'FACS_BM_month'] = 7
    blast_df.loc[year_only_mask, 'FACS_BM_day'] = 1
    
    # print(blast_df[(blast_df['FACS_BM_date_raw'].str.count('.') != 1) & ~year_only_mask]) 
    # raise
    assert (blast_df.loc[~year_only_mask, 'FACS_BM_date_raw'].str.count('\.') == 1).all()
    blast_df.loc[~year_only_mask, ['FACS_BM_year', 'FACS_BM_month']] = blast_df.loc[~year_only_mask, 'FACS_BM_date_raw'].str.split('.', expand=True).astype(int).rename(columns={0: 'FACS_BM_month', 1: 'FACS_BM_year'})
    assert blast_df.loc[~year_only_mask, 'FACS_BM_year'].between(0, 40).all()
    blast_df.loc[~year_only_mask, 'FACS_BM_year'] += 2000
    
    # fake middle of the month
    blast_df.loc[~year_only_mask, 'FACS_BM_day'] = 15

    
    
    blast_df['fake_FACS_BM_day'] = True
    blast_df['fake_FACS_BM_month'] = year_only_mask

    assert blast_df['FACS_BM_year'].notna().all()
    assert blast_df['FACS_BM_month'].notna().all()
    assert blast_df['FACS_BM_day'].notna().all()
    blast_df[["FACS_BM_year", "FACS_BM_month", "FACS_BM_day"]] = blast_df[["FACS_BM_year", "FACS_BM_month", "FACS_BM_day"]].astype(int)

    # print(pd.to_datetime(blast_df[["FACS_BM_year", "FACS_BM_month", "FACS_BM_day"]].iloc[0:1].rename(columns={'FACS_BM_year': 'year', 'FACS_BM_month': 'month', 'FACS_BM_day': 'day'}), errors='raise'))

    blast_df['FACS_BM_date'] = pd.to_datetime(blast_df[["FACS_BM_year", "FACS_BM_month", "FACS_BM_day"]].rename(columns={'FACS_BM_year': 'year', 'FACS_BM_month': 'month', 'FACS_BM_day': 'day'}), errors='raise')
    blast_df.drop(columns=['FACS_BM_year', 'FACS_BM_month', 'FACS_BM_day', 'FACS_BM_date_raw'], inplace=True)
    return blast_df


def add_donor_closest_bm_facs_blast_cols(df, clinical_df):
    blast_df = get_blast_df(clinical_df, epsilon_for_facs_blast_log=get_mds_params()['epsilon_for_facs_blast_log'])
    
    assert 'donor_id' in df.columns
    assert 'exp_name' in df.columns
    assert 'bleeding_date' in df.columns

    assert not df[['donor_id', 'exp_name']].duplicated().any()
    
    df['bleeding_date_as_date'] = pd.to_datetime(df['bleeding_date'], format='%d.%m.%y', errors='raise')


    orig_df_len = len(df)
    unexpected_common_cols = (set(df.columns) & set(blast_df.columns)) - {'donor_id'}
    assert not unexpected_common_cols, f'unexpected_common_cols: {unexpected_common_cols}'
    df = df.merge(blast_df, on='donor_id', how='left')
    df['bm_facs_blast_10x_sample_days_diff'] = (df['FACS_BM_date'] - df['bleeding_date_as_date']).dt.days
    df['bm_facs_blast_10x_sample_days_dist'] = df['bm_facs_blast_10x_sample_days_diff'].abs()
    df = df.sort_values('bm_facs_blast_10x_sample_days_dist', ascending=True).drop_duplicates(subset=['donor_id', 'exp_name'], keep='first')
    assert len(df) == orig_df_len, f'{len(df)}, {orig_df_len}'
    return df


def add_facs_blast_columns(df, epsilon_for_facs_blast_log=0.005):
    
    curr_df = df.copy()
    curr_df['bm_FACS_blast_frac'] = curr_df['FACS BM blasts (%)'].astype(str).apply(convert_raw_bm_facs_blast_to_mean_frac)

    curr_df['bm_FACS_blast_frac'] = curr_df['bm_FACS_blast_frac'].astype(float)
    curr_df['log_bm_FACS_blast_frac'] = np.log2(curr_df['bm_FACS_blast_frac'] + epsilon_for_facs_blast_log)
    return curr_df

def update_diagnosis_column_and_add_diagnosis_class_column(df, allow_df_to_contain_donor_id_and_exp_name=False):
    # why we need this function? because of the danger of map converting to nan...

    if not allow_df_to_contain_donor_id_and_exp_name and ({'donor_id', 'exp_name'} <= set(df.columns)):
        raise RuntimeError('you probably want to use get_df_with_diagnosis_and_diagnosis_class_by_exp_name_and_donor_id_and_clinical_df instead. otherwise, set allow_df_to_contain_donor_id_and_exp_name=True.')

    assert 'diagnosis' in df.columns
    df['diagnosis'] = df['diagnosis'].replace(DIAGNOSIS_VAL_TO_CORRECTED_DIAGNOSIS_VAL).astype(str)
    unexpected_diagnosis_vals = set(df['diagnosis'].unique()) - set(DIAGNOSIS_TO_DIAGNOSIS_CLASS)
    assert not unexpected_diagnosis_vals, f'unexpected_diagnosis_vals: {unexpected_diagnosis_vals}'
    df['diagnosis_class'] = df['diagnosis'].map(DIAGNOSIS_TO_DIAGNOSIS_CLASS).astype(str)


def get_df_with_diagnosis_and_diagnosis_class_by_exp_name_and_donor_id_and_clinical_df(df, clinical_df, dont_add_exp_date=True):
    # NOTE: this function does not update all metadata columns (e.g., if donor age was initially missing, this won't fix it)!
    drop_exp_date = False
    if 'exp_date' not in df:
        drop_exp_date = True
        sc_rna_seq_preprocessing_params.add_exp_date_column_to_df(df)
    
    df.drop(columns='diagnosis', inplace=True, errors='ignore')
    df = generic_utils.merge_preserving_df1_index_and_row_order(
        df, clinical_df[['donor_id', 'exp_date', 'diagnosis']], on=['donor_id', 'exp_date'], how='left')
    df = sc_rna_seq_preprocessing_params.overwrite_according_to_exp_name_and_donor_id(df, column_names_to_set={'diagnosis'})
    if drop_exp_date and dont_add_exp_date:
        df.drop(columns='exp_date', inplace=True)
    
    update_diagnosis_column_and_add_diagnosis_class_column(df, allow_df_to_contain_donor_id_and_exp_name=True)
    return df

