import os
import pandas as pd
import numpy as np
from generic import bio_utils
from generic import generic_utils
from generic import gtf_interface_and_utils

HG38_CHROM_DESC_SUFFIX = ':1 REF'
CHR_NAME_WITH_CHR_PREFIX_TO_CHR_NAME_WITHOUT_IT = {
    **{f'chr{x}': str(x) for x in list(range(1,23)) + ['X', 'Y']},
    'chrM': 'MT', 
}

def get_gene_loc_df(hg_genes_gtf_file_path, gene_loc_df_csv_file_path, use_existing_gene_loc_df_csv=False):
    if (not os.path.isfile(gene_loc_df_csv_file_path)) or (not use_existing_gene_loc_df_csv):
        gene_loc_df = gtf_interface_and_utils.get_gene_loc_df(hg_genes_gtf_file_path)
        gene_loc_df.to_csv(gene_loc_df_csv_file_path, index=False)
    
    return pd.read_csv(gene_loc_df_csv_file_path)

def get_relevant_gene_loc_df(c_or_mc_ad, hg_genes_gtf_file_path, gene_loc_df_csv_file_path, remove_chr_prefix_from_chr_name=False):
    gene_loc_df = get_gene_loc_df(
        hg_genes_gtf_file_path=hg_genes_gtf_file_path,
        gene_loc_df_csv_file_path=gene_loc_df_csv_file_path, 
        use_existing_gene_loc_df_csv=True,
    )
    gene_loc_df = generic_utils.merge_preserving_df1_index_and_row_order(
        c_or_mc_ad.var['gene_id'].reset_index().rename(columns={'index': 'gene_name'}), gene_loc_df, how='left')
    # print(gene_loc_df.loc[(~gene_id_missing_mask) & gene_loc_df['seqname'].isna()].head(100))
    missing_gene_mask = gene_loc_df['seqname'].isna()
    missing_gene_names = sorted(gene_loc_df.loc[missing_gene_mask, 'gene_name'].unique())
    num_of_missing_genes = missing_gene_mask.sum()
    print(f'get_relevant_gene_loc_df: returned gene_loc_df is missing {num_of_missing_genes} genes: {missing_gene_names}')
    if remove_chr_prefix_from_chr_name:
        gene_loc_df['seqname'].replace(CHR_NAME_WITH_CHR_PREFIX_TO_CHR_NAME_WITHOUT_IT, inplace=True)
    return gene_loc_df


def get_chrom_pos_df(hg_fasta_file_path, chrom_pos_df_csv_file_path, ordered_chrom_names=None, use_existing_chrom_pos_df_csv=False):
    if (not os.path.isfile(chrom_pos_df_csv_file_path)) or (not use_existing_chrom_pos_df_csv):
        chrom_seq_info_df = bio_utils.get_chrom_seq_info_df_from_fasta_file(hg_fasta_file_path)
        
        flat_dicts = []
        chrom_name_to_pos = {}
        for _, row in chrom_seq_info_df.iterrows():
            chrom_name = row['chrom_name']
            chrom_desc = row['chrom_desc']
            chrom_len = row['chrom_len']
            if chrom_desc.endswith(HG38_CHROM_DESC_SUFFIX):
                # e.g., /dummy/dummy/dummy/raid/human_genome_files/genome.fa (/net/mraid14/export/data/tgdata/users/orenmil/human_genome_files/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa)
                chrom_start_pos, chrom_end_pos = [int(x) for x in chrom_desc.partition(HG38_CHROM_DESC_SUFFIX)[0].split(':')[-2:]]
                # the following is not guaranteed. that's why the description is this way. to not include a prefix or suffix which is millions of Ns.
                # assert chrom_start_pos == 1
                # assert chrom_len == chrom_end_pos
            else:
                
                chrom_start_pos = 1
                chrom_end_pos = chrom_len
                chrom_name = CHR_NAME_WITH_CHR_PREFIX_TO_CHR_NAME_WITHOUT_IT.get(chrom_name, chrom_name)

            # print([chrom_name])
            chrom_name_to_pos[chrom_name] = (chrom_start_pos, chrom_end_pos)
            # print(chrom_start_pos, chrom_end_pos)
            flat_dicts.append({
                'chrom_name': chrom_name,
                'chrom_start_pos': chrom_start_pos, 
                'chrom_end_pos': chrom_end_pos,
            })
        chrom_pos_df = pd.DataFrame(flat_dicts)
        assert chrom_pos_df['chrom_name'].is_unique
        chrom_pos_df['annotated_chrom_len'] = chrom_pos_df['chrom_end_pos'] - chrom_pos_df['chrom_start_pos'] + 1

        if ordered_chrom_names is not None:
            total_prev_bases = 0
            flat_dicts = []
            for chrom_name in ordered_chrom_names:
                chrom_start_pos, chrom_end_pos = chrom_name_to_pos[chrom_name]
                
                unified_pos_diff = total_prev_bases - chrom_start_pos
                total_prev_bases += chrom_end_pos - chrom_start_pos + 1
                flat_dicts.append({
                    'chrom_name': chrom_name,
                    'unified_pos_diff': unified_pos_diff,
                })
            chrom_pos_df = generic_utils.merge_preserving_df1_index_and_row_order(chrom_pos_df, pd.DataFrame(flat_dicts), how='left')
            chrom_pos_df['chrom_unified_start_pos'] = chrom_pos_df['chrom_start_pos'] + chrom_pos_df['unified_pos_diff']
            chrom_pos_df['chrom_unified_end_pos'] = chrom_pos_df['chrom_end_pos'] + chrom_pos_df['unified_pos_diff']
        chrom_pos_df.to_csv(chrom_pos_df_csv_file_path, index=False)
    
    return pd.read_csv(chrom_pos_df_csv_file_path)

def get_bins_with_chroms_as_constraints_and_add_unified_end_pos_bin_i_to_var_df(
        chrom_pos_df, var_df, 
        wanted_bin_size=None, wanted_num_of_genes_in_bin=None, var_df_mask=None,
):
    chrom_pos_df = chrom_pos_df[~chrom_pos_df['unified_pos_diff'].isna()].copy()

    if var_df_mask is None:
        var_df_mask = np.full(len(var_df), True)
    filtered_var_df = var_df[var_df_mask]

    bins = set()
    flat_dicts = []
    for _, row in chrom_pos_df.iterrows():
        chrom_unified_start_pos = row['chrom_unified_start_pos']
        chrom_unified_end_pos = row['chrom_unified_end_pos']
        annotated_chrom_len = row['annotated_chrom_len']
        chrom_name = row['chrom_name']
        # print(f'chrom_name: {chrom_name}')
        
        if wanted_num_of_genes_in_bin is not None:
            assert wanted_bin_size is None
            assert filtered_var_df is not None
            curr_gene_end_positions = sorted(filtered_var_df.loc[filtered_var_df['chrom_name'] == chrom_name, 'unified_end_pos'].unique())
            if not curr_gene_end_positions:
                continue
            
            assert curr_gene_end_positions[0] >= chrom_unified_start_pos
            assert curr_gene_end_positions[-1] <= chrom_unified_end_pos
            
            
            num_of_curr_genes = len(curr_gene_end_positions) # for simplicity, assuming there is exactly one gene for each position. i guess this is mostly true.
            non_edge_bin_indices_in_gene_positions = list(range(0, num_of_curr_genes, wanted_num_of_genes_in_bin))[1:]
            if non_edge_bin_indices_in_gene_positions and (non_edge_bin_indices_in_gene_positions[-1] > (num_of_curr_genes - wanted_num_of_genes_in_bin)):
                non_edge_bin_indices_in_gene_positions = non_edge_bin_indices_in_gene_positions[:-1]

            curr_bins = (
                [chrom_unified_start_pos - 0.5] +
                [curr_gene_end_positions[i] - 0.5 for i in non_edge_bin_indices_in_gene_positions] +
                [chrom_unified_end_pos + 0.5]
            )

            # num_of_genes_in_last_bin = len([x for x in curr_gene_end_positions if x >= curr_bins[-2]])
            # print(f'num_of_genes_in_last_bin: {num_of_genes_in_last_bin}')
            num_of_genes_in_first_bin = len([x for x in curr_gene_end_positions if x < curr_bins[1]])
            # print(f'num_of_genes_in_first_bin: {num_of_genes_in_first_bin}')
            if num_of_genes_in_first_bin == 0:
                print(f'curr_bins: {curr_bins}')
                print(f'chrom_unified_start_pos: {chrom_unified_start_pos}')
                print(f'chrom_unified_end_pos: {chrom_unified_end_pos}')
                print(f'annotated_chrom_len: {annotated_chrom_len}')
                print(f'chrom_name: {chrom_name}')
                raise
            
        if wanted_bin_size is not None:
            assert wanted_num_of_genes_in_bin is None
            num_of_bins = max(1, annotated_chrom_len // wanted_bin_size)
            curr_bins = list(np.linspace(chrom_unified_start_pos - 0.5, chrom_unified_end_pos + 0.5, num_of_bins + 1))


        curr_gene_end_positions_series = pd.Series(curr_gene_end_positions)
        for i, (bin_start, bin_stop) in enumerate(zip(curr_bins[:-1], curr_bins[1:])):
            bin_mask = (curr_gene_end_positions_series >= bin_start) & (curr_gene_end_positions_series < bin_stop)
            flat_dicts.append({
                'chrom_name': chrom_name,
                'bin_start': bin_start,
                'bin_stop': bin_stop,
                'first_gene_unified_end_pos': curr_gene_end_positions_series[bin_mask].min(),
                'last_gene_unified_end_pos': curr_gene_end_positions_series[bin_mask].max(),
                'bin_name': f'{chrom_name}_{i}',
                'num_of_genes_in_bin': bin_mask.sum(),
            })
        bins |= set(curr_bins)

    bins = sorted(bins)
    
    bin_df = pd.DataFrame(flat_dicts)
    
    bin_df['bin_i'] = generic_utils.safe_digitize(bin_df['bin_start'], bins=bins)
    bin_df.sort_values('bin_i', inplace=True)
    bin_df.reset_index(drop=True, inplace=True)

    if 0:
        # for convenience, and not sure i even use this anywhere
        bin_df = generic_utils.merge_preserving_df1_index_and_row_order(bin_df, chrom_pos_df[['chrom_name', 'unified_pos_diff']].rename(
            columns={'unified_pos_diff': 'chrom_unified_pos_diff'}))

    # num_of_bins = len(bins) - 1
    # print(f'num_of_bins: {num_of_bins}')

    var_df.loc[var_df_mask, 'unified_end_pos_bin_i'] = generic_utils.safe_digitize(
        var_df.loc[var_df_mask, 'unified_end_pos'], bins=bins)
    var_df.loc[~var_df_mask, 'unified_end_pos_bin_i'] = np.nan

    return bins, bin_df

def get_var_df(orig_var_df, gene_loc_df, chrom_pos_df):
    # metacells_ad.var['gene_name'] = metacells_ad.var.index
    var_df = generic_utils.merge_preserving_df1_index_and_row_order(
        orig_var_df.rename(columns={'gene': 'gene_name'}), gene_loc_df.rename(columns={'start': 'start_pos_in_chrom', 'seqname': 'chrom_name', 'end': 'end_pos_in_chrom'}), 
        on=['gene_id', 'gene_name']).drop('gene_name', axis=1)
        # left_index=True, right_on='gene_name').drop('gene_name', axis=1)
    var_df = generic_utils.merge_preserving_df1_index_and_row_order(var_df, chrom_pos_df[['chrom_name', 'unified_pos_diff']], how='left')
    
    # assert set(var_df['strand'].unique()) <= {'+', '-'} # not guaranteed if gene_id is missing for some genes...
    forward_strand_mask = var_df['strand'] == '+'
    reverse_strand_mask = var_df['strand'] == '-'
    var_df.loc[forward_strand_mask, 'unified_end_pos'] = var_df.loc[forward_strand_mask, 'end_pos_in_chrom'] + var_df.loc[forward_strand_mask, 'unified_pos_diff']
    var_df.loc[reverse_strand_mask, 'unified_end_pos'] = var_df.loc[reverse_strand_mask, 'start_pos_in_chrom'] + var_df.loc[reverse_strand_mask, 'unified_pos_diff']

    return var_df

# def add_unified_end_pos_bin_i_to_var_df(var_df, bins):
#     relevant_gene_mask = (
#         (~var_df['unified_end_pos'].isna())
#         & (var_df['unified_end_pos'] >= bins[0])
#         & (var_df['unified_end_pos'] <= bins[-1])
#     )
#     var_df.loc[relevant_gene_mask, 'unified_end_pos_bin_i'] = generic_utils.safe_digitize(
#         var_df.loc[relevant_gene_mask, 'unified_end_pos'], bins=bins)

# def add_end_pos_bin_i_to_var_df_such_that_each_bin_contains_one_gene(var_df, gene_mask):
#     var_df.loc[gene_mask, 'unified_end_pos_bin_i'] = np.argsort(var_df.loc[gene_mask, 'unified_end_pos'])

