import pandas as pd
import numpy as np
import os
import re





from sc_rna_seq_preprocessing import sc_rna_seq_preprocessing_params
from mds import mds_analysis_params
from generic import liftover_interface_and_utils
from mds import mip_interface_and_utils
from generic import generic_utils


def get_sc_rna_seq_preprocessing_params():
    return sc_rna_seq_preprocessing_params.SC_RNA_SEQ_PREPROCESSING_PARAMS

def get_mds_params():
    return mds_analysis_params.MDS_ANALYSIS_PARAMS

def get_hotspot_mutation_desc(row):
    expected_prefix = row['gene'] + '-'
    hotspot_str = row['Hotspot']
    assert hotspot_str.startswith(expected_prefix), str((expected_prefix, hotspot_str))
    return row['gene'] + '__' + hotspot_str.partition(expected_prefix)[-1]

def get_srsf2_mutation_desc(row):
    expected_prefix = 'p.'
    hotspot_str = row['SRSF2']
    assert hotspot_str.startswith(expected_prefix)
    return 'SRSF2__' + hotspot_str.partition(expected_prefix)[-1]

def validate_exp_date_col(df):
    date_pattern = '[0-3][0-9]\.[01][0-9]\.[0-9][0-9]'
    invalid_exp_date_mask = df['exp_date'].apply(lambda x: re.fullmatch(date_pattern, x) is None)
    assert not invalid_exp_date_mask.any(), f'invalid exp dates:\n{df.loc[invalid_exp_date_mask, ["donor_id", "exp_date"]]}'

# https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=1788
HG37_DNMT3A_POS_RANGE = (25450743, 25565459)
HG38_DNMT3A_POS_RANGE = (25227874, 25342590)

# https://www.ncbi.nlm.nih.gov/gene/54790
HG37_TET2_POS_RANGE = (106067032, 106200960)
HG38_TET2_POS_RANGE = (105145875, 105279803)

# https://www.ncbi.nlm.nih.gov/gene/171023
HG37_ASXL1_POS_RANGE = (30946134, 31027122)
HG38_ASXL1_POS_RANGE = (32358331, 32439319)

# https://www.ncbi.nlm.nih.gov/gene/3418
HG37_IDH2_POS_RANGE = (90626277, 90645700)
HG38_IDH2_POS_RANGE = (90083045, 90102468)

# https://www.ncbi.nlm.nih.gov/gene/3845
HG37_KRAS_POS_RANGE = (25358180, 25403863)
HG38_KRAS_POS_RANGE = (25205246, 25250929)


def add_probably_hg37_38(mut_df):
    mut_df = mut_df.copy()
    
    dnmt3a_mask = mut_df['gene'] == 'DNMT3A'
    dnmt3a_df = mut_df[dnmt3a_mask]
    unexpected_dnmt3a_pos_mask = dnmt3a_df['POS'].notna() & ~(dnmt3a_df['POS'].between(*HG37_DNMT3A_POS_RANGE) | dnmt3a_df['POS'].between(*HG38_DNMT3A_POS_RANGE))
    assert not unexpected_dnmt3a_pos_mask.any(), dnmt3a_df[unexpected_dnmt3a_pos_mask]
    
    tet2_mask = mut_df['gene'] == 'TET2'
    tet2_df = mut_df[tet2_mask]
    unexpected_tet2_pos_mask = tet2_df['POS'].notna() & ~(tet2_df['POS'].between(*HG37_TET2_POS_RANGE) | tet2_df['POS'].between(*HG38_TET2_POS_RANGE))
    assert not unexpected_tet2_pos_mask.any(), tet2_df[unexpected_tet2_pos_mask]
    
    asxl1_mask = mut_df['gene'] == 'ASXL1'
    asxl1_df = mut_df[asxl1_mask]
    unexpected_asxl1_pos_mask = asxl1_df['POS'].notna() & ~(asxl1_df['POS'].between(*HG37_ASXL1_POS_RANGE) | asxl1_df['POS'].between(*HG38_ASXL1_POS_RANGE))
    assert not unexpected_asxl1_pos_mask.any(), asxl1_df[unexpected_asxl1_pos_mask]
    
    idh2_mask = mut_df['gene'] == 'IDH2'
    idh2_df = mut_df[idh2_mask]
    unexpected_idh2_pos_mask = idh2_df['POS'].notna() & ~(idh2_df['POS'].between(*HG37_IDH2_POS_RANGE) | idh2_df['POS'].between(*HG38_IDH2_POS_RANGE))
    assert not unexpected_idh2_pos_mask.any(), idh2_df[unexpected_idh2_pos_mask]
    
    kras_mask = mut_df['gene'] == 'KRAS'
    kras_df = mut_df[kras_mask]
    unexpected_kras_pos_mask = kras_df['POS'].notna() & ~(kras_df['POS'].between(*HG37_KRAS_POS_RANGE) | kras_df['POS'].between(*HG38_KRAS_POS_RANGE))
    assert not unexpected_kras_pos_mask.any(), kras_df[unexpected_kras_pos_mask]
    
    mut_df['probably_hg37'] = (
        # https://www.ncbi.nlm.nih.gov/clinvar/RCV001293765.1/
        ((mut_df['gene'] == 'SRSF2') & (mut_df['POS'] == 74732959))
        # https://www.ncbi.nlm.nih.gov/clinvar/RCV000424424/
        | ((mut_df['gene'] == 'SF3B1') & (mut_df['POS'] == 198266834))
        # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA149708
        | ((mut_df['gene'] == 'CALR') & (mut_df['POS'] == 13054564))
        | (dnmt3a_mask & mut_df['POS'].between(*HG37_DNMT3A_POS_RANGE))
        | (tet2_mask & mut_df['POS'].between(*HG37_TET2_POS_RANGE))
        | (asxl1_mask & mut_df['POS'].between(*HG37_ASXL1_POS_RANGE))
        | (idh2_mask & mut_df['POS'].between(*HG37_IDH2_POS_RANGE))
        | (kras_mask & mut_df['POS'].between(*HG37_KRAS_POS_RANGE))
    )
    mut_df['probably_hg38'] = (
        # https://www.ncbi.nlm.nih.gov/clinvar/RCV001293765.1/
        ((mut_df['gene'] == 'SRSF2') & (mut_df['POS'] == 76736877))
        # https://www.ncbi.nlm.nih.gov/clinvar/RCV000424424/
        | ((mut_df['gene'] == 'SF3B1') & (mut_df['POS'] == 197402110))
        # https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA149708
        | ((mut_df['gene'] == 'CALR') & (mut_df['POS'] == 12943750))
        | (dnmt3a_mask & mut_df['POS'].between(*HG38_DNMT3A_POS_RANGE))
        | (tet2_mask & mut_df['POS'].between(*HG38_TET2_POS_RANGE))
        | (asxl1_mask & mut_df['POS'].between(*HG38_ASXL1_POS_RANGE))
        | (idh2_mask & mut_df['POS'].between(*HG38_IDH2_POS_RANGE))
        | (kras_mask & mut_df['POS'].between(*HG38_KRAS_POS_RANGE))
    )
    return mut_df

def get_mutation_df(arch_mutation_file_path, hg_version, shlush_donor_id_and_exp_date_df, irrelevant_shlush_donor_ids):
    source_file_name = os.path.basename(arch_mutation_file_path)
    # arch_mutation_params = get_mds_params()['arch_mutations']
    
    if source_file_name.endswith('.xlsx'):
        sheet_name_to_df = pd.read_excel(arch_mutation_file_path, sheet_name=None)
        assert len(sheet_name_to_df) == 1
        df = next(iter(sheet_name_to_df.values()))

    elif source_file_name.endswith('.csv'):
        df = pd.read_csv(arch_mutation_file_path)
    elif source_file_name.endswith('.tsv'):
        df = pd.read_csv(arch_mutation_file_path, sep='\t')
    else:
        assert False

    df.rename(columns={
        'SAMPLE_ID': 'donor_id', 
        'Sample_ID': 'donor_id', 
        'sample_ID': 'donor_id', 
        'Individual ID': 'donor_id', 
        'Individual.ID': 'donor_id', 
        'Patient_ID': 'donor_id', 
        'Gene.refGene': 'gene', 
        'Gene': 'gene', 
        'samping date': 'exp_date',
        'sampling date': 'exp_date',
        'sampling_date': 'exp_date',
        'Chr': 'CHR',
        'Chromosome': 'CHR',
        'Start': 'POS',
        'Position': 'POS',
        'Ref': 'REF',
        'Reference allele': 'REF',
        'Reference.allele': 'REF',
        'Alt': 'ALT',
        'Alternate allele': 'ALT',
        'Alternate.allele': 'ALT',
        'Average VAF': 'AVG_VAF',
        'Average.VAF': 'AVG_VAF',
    }, inplace=True)

    if source_file_name == '240727_hg38_nature_med_individual_mutations_from_nimrod.csv':
        df = df[df['donor_id'] != 'N55'] # this is wrong. it is just a copy of N295, which is a later sample of the same donor, and appears in the ARCH4 file. only the 23.11.22 sample of N55 (aka N295) appears in the 148 atlas.
    
    # print(set(df.columns))
    assert set(shlush_donor_id_and_exp_date_df.columns) == {'donor_id', 'exp_date'}
    if not ({'donor_id', 'exp_date'} <= set(df.columns)):
        df = generic_utils.merge_preserving_df1_index_and_row_order(
            df, shlush_donor_id_and_exp_date_df, 
            # how='left',
        )
    
    df['exp_date'].replace({
        '4.12.23': '04.12.23',
        '5.12.23': '05.12.23',
        '22.06.23_p1': '22.06.23',
    }, inplace=True)
    validate_exp_date_col(df)

    if 'AVG_VAF' not in df.columns:
        assert 'VAF' in df.columns
        assert 'VAF_DUP' not in df.columns
        if set(df.columns) == {'exp_date', 'CHR', 'Depth', 'End', 'VAF', 'Count', 'donor_id', 'POS'}:
            group_by_cols = ['exp_date', 'donor_id', 'CHR', 'End', 'POS']
        elif set(df.columns) == {'AAChange.refGene', 'exp_date', 'Depth', 'End', 'cosmic70', 'VAF', 'Exonic Function', 'CHR', 'ALT', 'REF', 'donor_id', 'Function', 'avsnp150', 'gene', 'POS'}:
            group_by_cols = ['AAChange.refGene', 'exp_date', 'End', 'cosmic70', 'Exonic Function', 'CHR', 'ALT', 'REF', 'donor_id', 'Function', 'avsnp150', 'gene', 'POS']
        elif set(df.columns) == {'VAF', 'gene', 'donor_id', 'exp_date'}:
            group_by_cols = ['gene', 'donor_id', 'exp_date']
        else:
            group_by_cols = None
            should_be_unique_cols = ['exp_date', 'donor_id', 'CHR', 'POS']
            assert not df[should_be_unique_cols].duplicated().any()
            df.rename(columns={'VAF': 'mean_VAF'}, inplace=True)
        if group_by_cols is not None:
            assert (
                ({'exp_date', 'donor_id', 'CHR', 'POS'} <= set(group_by_cols)) # good
                | (set(df.columns) == {'VAF', 'gene', 'donor_id', 'exp_date'}) # ugh
            )
            df = df.groupby(group_by_cols)['VAF'].mean().reset_index(name='AVG_VAF')

    df.rename(columns={
        'AVG_VAF': 'mean_VAF',
    }, inplace=True)
    
    print(f'dropping the following donor_ids because they are irrelevant: {irrelevant_shlush_donor_ids}')
    df = df[~df['donor_id'].isin(irrelevant_shlush_donor_ids)]


    df['donor_id'].replace(get_sc_rna_seq_preprocessing_params()['other_donor_id_to_donor_id_i_use'], inplace=True)
    
    
    if source_file_name in {'240621_CALR_by_nili.xlsx'}:
        df['gene'] = 'CALR'
    for col in ['CHR', 'POS', 'REF', 'ALT']:
        if col not in df.columns:
            df[col] = np.nan


    u2af1_u2af1l5_mask = df['gene'] == 'U2AF1;U2AF1L5'
    if u2af1_u2af1l5_mask.any():
        df.loc[u2af1_u2af1l5_mask, 'gene'] = 'U2AF1'

    df['gene'].replace({
        'U2AF1;U2AF1L5': 'U2AF1',
        'CALR .K385fs': 'CALR',
        
    }, inplace=True)

    df['mutation_desc'] = np.nan # this is needed to avoid an error in the iloc that uses get_loc below.

    if 'ExonicFunc.refGene' in df.columns:
        # NOTE: i think we should be able to distinguish between hotspot mutations, stopgains, synonymous and nonsynonymous, because otherwise we might have false positives...
        unexpected_exonic_funcs = set(df['ExonicFunc.refGene'].unique()) - {'nonsynonymous SNV', 'stopgain', 'frameshift insertion', 'frameshift deletion', 'frameshift_variant', 'splicing', 0, '0'}
        if unexpected_exonic_funcs:
            print("df[df['ExonicFunc.refGene'].isin(list(unexpected_exonic_funcs))]")
            print(df[df['ExonicFunc.refGene'].isin(list(unexpected_exonic_funcs))])
        assert not unexpected_exonic_funcs, unexpected_exonic_funcs

        stopgain_or_frameshift_mask = df['ExonicFunc.refGene'].isin(['stopgain', 'frameshift insertion', 'frameshift deletion', 'frameshift_variant'])
        df.loc[stopgain_or_frameshift_mask, 'mutation_desc'] = df.loc[stopgain_or_frameshift_mask, 'gene'] + '__stopgain_or_frameshift'

        nonsynonymous_mask = df['ExonicFunc.refGene'] == 'nonsynonymous SNV'
        df.loc[nonsynonymous_mask, 'mutation_desc'] = df.loc[nonsynonymous_mask, 'gene'] + '__nonsynonymous'

        splicing_mask = df['ExonicFunc.refGene'] == 'splicing'
        df.loc[splicing_mask, 'mutation_desc'] = df.loc[splicing_mask, 'gene'] + '__splicing'
        
        unknown_exonic_func_mask = df['ExonicFunc.refGene'] == 0
        df.loc[unknown_exonic_func_mask, 'mutation_desc'] = df.loc[unknown_exonic_func_mask, 'gene'] + '__unknown'


    

    assert not df.duplicated(subset=['donor_id', 'gene', 'CHR', 'POS', 'REF', 'ALT', 'exp_date']).any()
    df = df[['donor_id', 'gene', 'mutation_desc', 'mean_VAF', 'exp_date', 'CHR', 'POS', 'REF', 'ALT']].copy()


    if 0:
        # i don't see why this should be done...    
        dropped_mutation_row_count_df = df[['donor_id', 'mutation_desc', 'exp_date']].value_counts() - 1
        dropped_mutation_row_count_df = dropped_mutation_row_count_df[dropped_mutation_row_count_df > 0]
        print(f'the following number of mutation rows were dropped because a donor had in the same (sample) date multiple mutations with the same mutation_desc. the row of the mutation with highest VAF was kept:\n{dropped_mutation_row_count_df}\n')
        df = df.sort_values('mean_VAF', ascending=False).drop_duplicates(
            subset=['donor_id', 'mutation_desc', 'exp_date'], keep='first')
        assert not df.duplicated(subset=['donor_id', 'mutation_desc', 'exp_date']).any()

    df['source_file_name'] = source_file_name

    assert hg_version in {'hg37', 'hg38', 'both', 'irrelevant_because_only_gene_names'}
    df['source_file_hg_version'] = hg_version
    if hg_version == 'hg37':
        # orig_pos_with_hg38_df = liftover_interface_and_utils.replace_hg37_with_hg38_coordinates(
        #     df=df[['CHR', 'POS']].drop_duplicates(), chr_col='CHR', start_col='POS', end_col=None) 
        df = liftover_interface_and_utils.replace_hg37_with_hg38_coordinates(
            df=df, chr_col='CHR', start_col='POS', end_col=None, convert_to_vireo_chr_names=True) 
    elif hg_version == 'both':
        print(f'WARNING: {arch_mutation_file_path} contains both hg37 and hg38 positions!')

    ext_df = add_probably_hg37_38(df)
    if ext_df['probably_hg37'].any():
        print(f'WARNING: {arch_mutation_file_path} probably contains hg37 positions!')
        print('positions that raise the suspicion:')
        print(ext_df[ext_df["probably_hg37"]])

    return df

def get_all_mutation_df(shlush_donor_id_and_exp_date_df):
    print('TODO: make sure to distinguish between zero VAF and nan (not enough DNA data)')
    mut_dfs = []
    arch4_file_paths = get_mds_params()['arch_mutations']['arch4_mutation_file_paths']
    other_arch_file_paths = get_mds_params()['arch_mutations']['other_arch_mutation_file_paths']
    mostly_arch3_file_paths = get_mds_params()['arch_mutations']['mostly_arch3_mutation_file_paths']
    get_paths_as_set = lambda x: {y[0] for y in x}
    generic_utils.assert_sets_are_disjoint([get_paths_as_set(arch4_file_paths), get_paths_as_set(other_arch_file_paths), get_paths_as_set(mostly_arch3_file_paths)])
    for mut_file_path, hg_version in arch4_file_paths + other_arch_file_paths + mostly_arch3_file_paths:
        print(mut_file_path)
        mut_df = get_mutation_df(
            mut_file_path,
            hg_version,
            shlush_donor_id_and_exp_date_df, 
            get_mds_params()['arch_mutations']['mds_irrelevant_shlush_donor_ids'],
        )
        mut_df['ARCH4'] = mut_file_path in arch4_file_paths
        mut_df['mostly_ARCH3'] = mut_file_path in mostly_arch3_file_paths
        mut_dfs.append(mut_df)
    df = pd.concat(mut_dfs, ignore_index=True)

    return df

def get_healthy_atlas_mutation_df(shlush_donor_id_and_exp_date_df):
    mut_file_path = get_mds_params()['arch_mutations']['healthy_atlas_arch_mutation_xlsx_file_path']
    source_file_name = os.path.basename(mut_file_path)
    df = pd.read_excel(mut_file_path)
    df.rename(columns={
        'Individual ID': 'donor_id',
        'Gene': 'gene',
        'Chromosome': 'CHR',
        'Position': 'POS',
        'Reference allele': 'REF',
        'Alternate allele': 'ALT',
        'Average VAF': 'mean_VAF',
    }, inplace=True)
    assert set(shlush_donor_id_and_exp_date_df.columns) == {'donor_id', 'exp_date'}

    irrelevant_shlush_donor_ids = get_mds_params()['arch_mutations']['healthy_atlas_irrelevant_shlush_donor_ids']
    print(f'dropping the following donor_ids because they are irrelevant: {irrelevant_shlush_donor_ids}')
    df = df[~df['donor_id'].isin(irrelevant_shlush_donor_ids)]

    df = generic_utils.merge_preserving_df1_index_and_row_order(df, shlush_donor_id_and_exp_date_df, 
                                                                # how='left',
                                                                )
    df['donor_id'].replace(get_sc_rna_seq_preprocessing_params()['other_donor_id_to_donor_id_i_use'], inplace=True)
    df['source_file_name'] = source_file_name
    return df

def get_agg_donor_small_scale_mut_df(mutation_df):
    print('TODO: make sure to distinguish between zero VAF and nan (not enough DNA data)')
    agg_df = mutation_df.groupby(['donor_id', 'exp_date'])['mean_VAF'].max().reset_index(name='max_mean_VAF')
    donor_date_gene_vaf_df = mutation_df.groupby(['donor_id', 'exp_date', 'gene'])['mean_VAF'].max().reset_index(name='max_mean_VAF')
    for gene in sorted(donor_date_gene_vaf_df['gene'].unique()):
        gene_max_vaf_col = f'{gene}_max_mean_VAF'
        agg_df = generic_utils.merge_preserving_df1_index_and_row_order(agg_df, donor_date_gene_vaf_df.loc[
            donor_date_gene_vaf_df['gene'] == gene, ['donor_id', 'exp_date', 'max_mean_VAF']].rename(
                columns={'max_mean_VAF': gene_max_vaf_col}), on=['donor_id', 'exp_date'], how='left', verbose=False)
        agg_df[gene_max_vaf_col].fillna(0, inplace=True)
        agg_df[f'{gene}_snv'] = agg_df[gene_max_vaf_col] > 0

    donor_gene_mut_count_df = mutation_df[['donor_id', 'exp_date', 'gene', 'CHR', 'POS', 'REF', 'ALT']].drop_duplicates()[
        ['donor_id', 'exp_date', 'gene']].value_counts().reset_index()
    for gene in sorted(donor_gene_mut_count_df['gene'].unique()):
        gene_mut_count_col = f'{gene}_mut_count'
        agg_df = generic_utils.merge_preserving_df1_index_and_row_order(
            agg_df, donor_gene_mut_count_df.loc[
            donor_gene_mut_count_df['gene'] == gene, ['donor_id', 'exp_date', 'count']].rename(
                columns={'count': gene_mut_count_col}), on=['donor_id', 'exp_date'], how='left', verbose=False)
        agg_df[gene_mut_count_col].fillna(0, inplace=True)
    return agg_df

def get_minimal_all_mutation_df():
    return pd.read_csv(get_mds_params()['arch_mutations']['minimal_mutation_df_csv_file_path'])

