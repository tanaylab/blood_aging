import os
import pandas as pd
import pathlib
from generic import generic_utils

LIFTOVER_BIN_PATH = '/dummy/dummy/dummy/raid/software/liftOver'
HG37_TO_HG38_OVER_CHAIN_PATH = '/dummy/dummy/dummy/raid/human_genome_files/hg19ToHg38.over.chain.gz'

OTHER_CHR_TO_LIFTOVER_REPR = {
    # hg19 and GRCh38: https://gatk.broadinstitute.org/hc/en-us/articles/360035890711.
    # maybe could improve this in a smarter way by going over the table in this article?

    # compiled manually simply by searching hg19ToHg38.over.chain. the pattern seems clear, but there seem to not be so many, so i checked them manually.
    'GL000230.1': 'chrUn_gl000230', # https://www.ncbi.nlm.nih.gov/nuccore/GL000230
    'GL000214.1': 'chrUn_gl000214',
    'GL000211.1': 'chrUn_gl000211',
    'GL000232.1': 'chrUn_gl000232',
    'GL000219.1': 'chrUn_gl000219',

    **{str(x): f'chr{x}' for x in range(1,23)},
    **{x: f'chr{x}' for x in range(1,23)},
    'X': 'chrX',
    'Y': 'chrY',
    'MT': 'chrM',

    'GL000192.1': 'chrUn_gl000192', # didn't find it in hg19ToHg38.over.chain. ugh.
}
def get_chr_repr_for_liftover(chr):
    if chr in OTHER_CHR_TO_LIFTOVER_REPR:
        return OTHER_CHR_TO_LIFTOVER_REPR[chr]
    
    print([chr]) # printing in a list allows seeing whether it is a str or an int
    raise RuntimeError('simply add this unknown chr to OTHER_CHR_TO_LIFTOVER_REPR')

LIFTOVER_CHR_REPR_TO_VIREO_REPR = {
    **{f'chr{x}': str(x) for x in list(range(1,23)) + ['X', 'Y']},

    # for the following, I don't really know what they map to in vireo_repr. but verify_all_vireo_chr_reprs_are_known() should let us know.
    'chrM': 'MT', 
    'chrUn_GL000214v1': 'chrUn_GL000214v1',
    'chrUn_GL000219v1': 'chrUn_GL000219v1',
}
def verify_all_vireo_chr_reprs_are_known(reprs):
    all_reprs = set(LIFTOVER_CHR_REPR_TO_VIREO_REPR.values())
    for repr in reprs:
        if repr not in all_reprs:
            print(repr)
            raise RuntimeError('add this vireo repr to LIFTOVER_CHR_REPR_TO_VIREO_REPR')

def get_vireo_chr_repr_for_liftover_repr(chr):
    if chr in LIFTOVER_CHR_REPR_TO_VIREO_REPR:
        return LIFTOVER_CHR_REPR_TO_VIREO_REPR[chr]
    
    print([chr]) # printing in a list allows seeing whether it is a str or an int
    raise RuntimeError('simply add this unknown chr to LIFTOVER_CHR_REPR_TO_VIREO_REPR')


def add_hg38_positions_to_df(df, out_dir_path=None, log_file_path=None, chr_col='Chr', start_col='Start', end_col='End', convert_to_vireo_chr_names=False):
    if ((end_col is not None) and (list(df) != [chr_col, start_col, end_col])) or ((end_col is None) and (list(df) != [chr_col, start_col])):
        err_str = f'df must have exactly the following columns: "{chr_col}", "{start_col}"'
        if end_col is not None:
            err_str += f' and "{end_col}"'
        raise RuntimeError(err_str)
    if out_dir_path is None:
        out_dir_path = 'temp/liftover'
    
    df = df.copy()

    dummy_end_col = end_col is None
    if dummy_end_col:
        assert 'dummy_end' not in df.columns
        end_col = 'dummy_end'
        df[end_col] = df[start_col]

    pathlib.Path(out_dir_path).mkdir(parents=True, exist_ok=True)
    hg37_loci_bed_file_path = os.path.join(out_dir_path, 'hg37_loci.bed')
    hg37_loci_converted_to_hg38_bed_file_path = os.path.join(out_dir_path, 'hg37_loci_converted_to_hg38.bed')
    hg37_loci_that_failed_to_convert_bed_file_path = os.path.join(out_dir_path, 'hg37_loci_that_failed_to_convert.bed')
    df['Chr_for_liftover'] = df[chr_col].apply(get_chr_repr_for_liftover)
    df['uniq_id_to_match_liftover_input_and_output'] = [f'myid{i}' for i in range(len(df))]
    # df[['Chr_for_liftover', start_col, end_col]].to_csv(hg37_loci_bed_file_path, header=False, index=False, sep='\t')
    df[['Chr_for_liftover', start_col, end_col, 'uniq_id_to_match_liftover_input_and_output']].to_csv(hg37_loci_bed_file_path, header=False, index=False, sep='\t')

    # liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed
    generic_utils.run_subprocess(
        [
            LIFTOVER_BIN_PATH,
            hg37_loci_bed_file_path,
            HG37_TO_HG38_OVER_CHAIN_PATH,
            hg37_loci_converted_to_hg38_bed_file_path,
            hg37_loci_that_failed_to_convert_bed_file_path,
        ],
        shell=False,
        log_file_path=log_file_path,
    )

    failed_lines = generic_utils.read_text_file(hg37_loci_that_failed_to_convert_bed_file_path).splitlines()
    failed_lines = [x for x in failed_lines if x != '#Deleted in new']
    failed_lines = [tuple(x.split()) for x in failed_lines]
    failed_df = pd.DataFrame(failed_lines, columns=['Chr_for_liftover', start_col, end_col, 'uniq_id_to_match_liftover_input_and_output'])
    failed_df[start_col] = failed_df[start_col].astype(int)
    failed_df[end_col] = failed_df[end_col].astype(int)
    if not failed_df.empty:
        print(f'failed to convert to hg38:\n{failed_df}')
    else:
        print('all lines were converted successfully')

    filtered_df = df.merge(failed_df, how='outer', indicator=True)
    filtered_df = filtered_df.loc[filtered_df['_merge'] == 'left_only']
    filtered_df.reset_index(drop=True, inplace=True) # This is essential for the pd.concat later to "ignore" row indices.
    filtered_df.drop('_merge', axis=1, inplace=True)
    converted_df = pd.read_csv(
        hg37_loci_converted_to_hg38_bed_file_path, names=['hg38_chr', 'hg38_start', 'hg38_end', 'uniq_id_to_match_liftover_input_and_output'], sep='\t')
    assert len(converted_df) == len(filtered_df)
    assert converted_df['uniq_id_to_match_liftover_input_and_output'].is_unique
    assert filtered_df['uniq_id_to_match_liftover_input_and_output'].is_unique
    orig_filtered_df_len = len(filtered_df)
    filtered_with_converted_df = filtered_df.merge(converted_df)
    assert len(filtered_with_converted_df) == orig_filtered_df_len

    # filtered_with_converted_df = pd.concat([filtered_df, converted_df], axis=1)
    assert not filtered_with_converted_df.isna().any().any()
    # assert not filtered_with_converted_df[[chr_col, start_col, end_col, 'Chr_for_liftover', 'hg38_chr', 'hg38_start', 'hg38_end']].isna().any().any()
    # print(list(filtered_with_converted_df))

    filtered_with_converted_df.drop(['Chr_for_liftover', 'uniq_id_to_match_liftover_input_and_output'], axis=1, inplace=True)
    if dummy_end_col:
        filtered_with_converted_df.drop(columns=[end_col, 'hg38_end'], inplace=True)

    if convert_to_vireo_chr_names:
        filtered_with_converted_df['hg38_chr'].replace(LIFTOVER_CHR_REPR_TO_VIREO_REPR, inplace=True)

    return filtered_with_converted_df

def replace_hg37_with_hg38_coordinates(df, chr_col='Chr', start_col='Start', end_col='End', **kwargs):
    
    orig_cols = [chr_col, start_col]
    if end_col is not None:
        orig_cols.append(end_col)
    with_hg38_df = add_hg38_positions_to_df(df[orig_cols].drop_duplicates(), chr_col=chr_col, start_col=start_col, end_col=end_col, **kwargs)


    assert not ({'hg38_chr', 'hg38_start', 'hg38_end'} & set(df.columns))

    df = generic_utils.merge_preserving_df1_index_and_row_order(df, with_hg38_df, on=orig_cols)

    for orig_col, col38 in [
        (chr_col, 'hg38_chr'),
        (start_col, 'hg38_start'),
        *([] if end_col is None else [(end_col, 'hg38_end')]),
    ]:
        df[f'{orig_col}_hg37'] = df[orig_col]
        df[orig_col] = df[col38]
        df.drop(columns=col38, inplace=True)
    
    return df
    
