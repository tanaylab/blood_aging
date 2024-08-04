import os
# import vcf
import importlib
import re
import pandas as pd
import pathlib
from generic import generic_utils


def read_gtf(gtf_file_path):
    num_of_header_rows = generic_utils.get_num_of_header_lines(gtf_file_path, header_line_prefix='#')
    assert num_of_header_rows >= 1

    GTF_COLUMN_NAMES = [
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute',
    ]

    return pd.read_csv(gtf_file_path, names=GTF_COLUMN_NAMES, skiprows=num_of_header_rows, header=None, sep='\t', low_memory=False)

def get_gene_loc_df(gtf_file_path):
    gtf_df = read_gtf(gtf_file_path)

    gene_df = gtf_df[gtf_df['feature'] == 'gene']
    gene_loc_df = gene_df[['seqname', 'start', 'end', 'strand']].copy()
    # gene_loc_df['gene_name'] = gene_df['attribute'].str.split('gene_name "', expand=True)[1].str.partition('"; ', expand=True)[0] # NOTE: unlike gene_id, gene_name is not unique. so dont use it.
    gene_loc_df['gene_id'] = gene_df['attribute'].str.split('gene_id "', expand=True)[1].str.partition('"; ', expand=True)[0]
    assert gene_loc_df['gene_id'].is_unique

    return gene_loc_df.reset_index(drop=True)

def fix_gene_loc_df_to_have_a_single_row_for_each_gene(gene_loc_df):
    genes_with_mulitple_rows_names = list((gene_loc_df[
        # 'gene_name'
        'gene_id'
    ].value_counts() > 1).loc[lambda x: x].index)
    num_of_genes_with_mulitple_rows_names = len(genes_with_mulitple_rows_names)
    print(f'num_of_genes_with_mulitple_rows_names: {num_of_genes_with_mulitple_rows_names}')


    print(f'verifying that for each gene with multiple rows, all rows of the gene are on the same chromosome.')
    for gene_name in genes_with_mulitple_rows_names:
        curr_df = gene_loc_df[gene_loc_df['gene_name'] == gene_name]
        assert len(set(curr_df['seqname'])) == 1
        if 0:
            orig_gene_lens = list(curr_df['end'] - curr_df['start'] + 1)
            max_orig_gene_len = max(orig_gene_lens)
            new_gene_start = curr_df['start'].min()
            new_gene_end = curr_df['end'].max()
            new_gene_len = new_gene_end - new_gene_start + 1
            new_orig_gene_len_ratio = new_gene_len / max_orig_gene_len
            print(f'gene_name: {gene_name}, max_orig_gene_len: {max_orig_gene_len}, new_orig_gene_len_ratio: {new_orig_gene_len_ratio}')


    print('for each gene with multiple rows, filtering to keep only the row with the end in the furthest position (a quite arbitrary choice).')
    fixed_gene_loc_df = gene_loc_df.sort_values('end', ascending=False).drop_duplicates(subset='gene_name', keep='first')
    assert fixed_gene_loc_df['gene_name'].is_unique

    return fixed_gene_loc_df
