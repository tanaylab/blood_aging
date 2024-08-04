from generic import generic_utils
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

def replace_other_donor_ids_with_ones_i_use(
        df, 
        other_donor_id_to_donor_id_i_use, 
        donor_id_column_name='donor_id',
):
    for other_donor_id, donor_id_i_use in other_donor_id_to_donor_id_i_use.items():
        if 0:
            other_donor_id_batch_names = set(df.loc[df[donor_id_column_name] == other_donor_id, 'seq_batch_name'].unique())
            donor_id_batch_names = set(df.loc[df[donor_id_column_name] == donor_id_i_use, 'seq_batch_name'].unique())
            # given current flow of this function, failing this assert wouldn't be a problem, as we would just take the mean of freqs etc for each donor.
            batch_names_with_both_donor_ids = other_donor_id_batch_names & donor_id_batch_names
            assert not batch_names_with_both_donor_ids, str((other_donor_id, donor_id_i_use, batch_names_with_both_donor_ids)) # NOTE: the scenario of sequencing the same donor more than once (but with different ids) in the same batch actually did happen, and isn't a problem, as mentioned above. thus, this whole code block should not run.
        
        df.loc[df[donor_id_column_name] == other_donor_id, donor_id_column_name] = donor_id_i_use

    
def plot_scatter_comparing_two_mip_donors(
        mip_df,
        donor1_id,
        donor2_id,
):
    donor1_df = mip_df[mip_df['sample_name'] == donor1_id]
    donor2_df = mip_df[mip_df['sample_name'] == donor2_id]
    merged_df = generic_utils.merge_and_verify_no_df1_row_was_duplicated_but_allow_losing_df1_rows(
        donor1_df, donor2_df, on=['hg38_chr', 'hg38_start'], suffixes=(f'_{donor1_id}', f'_{donor2_id}'))
    jitter_sd = 0.01
    vaf1_name = f'mean_ref_freq_{donor1_id}'
    vaf2_name = f'mean_ref_freq_{donor2_id}'
    merged_df[vaf1_name] = merged_df[vaf1_name] + np.random.normal(0, jitter_sd, len(merged_df))
    merged_df[vaf2_name] = merged_df[vaf2_name] + np.random.normal(0, jitter_sd, len(merged_df))

    plt.close('all')
    sb.scatterplot(data=merged_df, x=vaf1_name, y=vaf2_name, s=5)
