import os
import sys
import shutil
import importlib
import itertools
import metacells as mc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy
import scipy.spatial.distance
import matplotlib
import seaborn as sb
import matplotlib.pyplot as plt
import pathlib
import logging
import re
import collections
import sklearn
import pickle
import random

# from generic import gene_module_utils
from generic import generic_utils
from mds import mds_analysis_params
# from mds import ipssm_utils
from generic import mc_utils
from mds import arch_mutation_interface_and_utils
from mds import clinical_data_interface_and_utils
from mds import pb_cd34_c_score_threshs
from sc_rna_seq_preprocessing import sc_rna_seq_preprocessing_params

plt.rcParams["patch.force_edgecolor"] = False
plt.rcParams['patch.linewidth'] = 0
plt.rcParams['patch.edgecolor'] = 'none'
# plt.rcParams['scatter.edgecolors'] = 'black' # didnt affect sb.scatterplot

def get_mds_params():
    return mds_analysis_params.MDS_ANALYSIS_PARAMS


def get_sc_rna_seq_preprocessing_params():
    return sc_rna_seq_preprocessing_params.SC_RNA_SEQ_PREPROCESSING_PARAMS


def plot_within_state_sig_marker_heatmaps(mc_ad, out_dir_path, sig_name_to_info, sig_name_to_genes, state_col='state', relevant_state_attr_ordered=['relevant_states']):
    import matplotlib.pyplot as plt
    
    mc_utils.write_expr_and_expr_enrich(mc_ad)
    pathlib.Path(out_dir_path).mkdir(parents=True, exist_ok=True)
    for sig_name, sig_info in sig_name_to_info.items():
        print(f'starting work on {sig_name}')
        for relevant_state_attr in relevant_state_attr_ordered:
            relevant_state_lists = sig_info.get(relevant_state_attr)
            if relevant_state_lists is not None:
                break
        if isinstance(relevant_state_lists[0], str):
            relevant_state_lists = [relevant_state_lists]

        # if sig_name != 'nktdp_samhd1_sig':
        #     continue
        for i, relevant_states in enumerate(relevant_state_lists):
            plt.close('all')
            clustermap_obj = mc_utils.plot_marker_heatmap(
                mc_ad, 
                mc_ad.obs[state_col].isin(relevant_states),
                genes_to_include_anyway=(),
                use_these_genes_instead_of_top_genes=sig_name_to_genes[sig_name],
                show_mc_names=False,
                add_prefixes_to_gene_names=False,
                # cbar_pos=(0.02, 0.92, 0.001, 0.001),
                cbar_pos=None,
            )
            clustermap_obj.fig.suptitle(sig_name)
            clustermap_obj.fig.tight_layout()
            clustermap_obj.savefig(os.path.join(out_dir_path, f'{sig_name}{i+1}.png'))
            plt.close('all')
            # raise

def plot_sig_marker_heatmaps_across_relevant_states(mc_ad, out_dir_path, sig_name_to_genes, sig_names, sig_name_to_relevant_states, state_col='state'):
    import matplotlib.pyplot as plt

    mc_utils.write_expr_and_expr_enrich(mc_ad)
    pathlib.Path(out_dir_path).mkdir(parents=True, exist_ok=True)

    curr_states = list(mc_ad.obs[state_col].unique())
    curr_hspc_states = [x for x in curr_states if x in mds_analysis_params.PB_HSPC_STATE_NAMES]
    curr_myeloid_hspc_states = [x for x in curr_states if x in mds_analysis_params.MYELOID_HSPC_STATE_NAMES_FOUND_IN_NORMAL_PB]

    for sig_name in sig_names:
        # if sig_name != 'higher_in_b_than_pre_b':
        #     continue
        print(f'starting work on {sig_name}')
        genes = sig_name_to_genes[sig_name]
        plt.close('all')
        relevant_states = sig_name_to_relevant_states[sig_name]
        if 'all HSPCs' in relevant_states:
            relevant_states = list((set(relevant_states) - {'all HSPCs'}) | set(curr_hspc_states))
        if 'all myeloid HSPCs' in relevant_states:
            relevant_states = list((set(relevant_states) - {'all myeloid HSPCs'}) | set(curr_myeloid_hspc_states))
        if len(relevant_states) <= 1:
            print(f'skipping {sig_name} because len(relevant_states) <= 1')
            continue
        print(relevant_states)
        clustermap_obj = mc_utils.plot_marker_heatmap(
            mc_ad, 
            mc_ad.obs[state_col].isin(relevant_states),
            equal_contrib_for_each_val_of_col=state_col,
            equal_contrib_for_each_val_of_col_min_mc_count=20, # 240623: i think it was 1 when i generated some of the heatmaps, but this shouldn't change anything, i think
            equal_contrib_for_each_val_of_col_max_mc_count=20,
            genes_to_include_anyway=(),
            use_these_genes_instead_of_top_genes=genes,
            show_mc_names=False,
            add_prefixes_to_gene_names=False,
            # cbar_pos=(0.02, 0.92, 0.001, 0.001),
            cbar_pos=None,
        )
        clustermap_obj.fig.suptitle(sig_name)
        clustermap_obj.fig.tight_layout()
        clustermap_obj.savefig(os.path.join(out_dir_path, f'{sig_name}.png'))
        plt.close('all')
        # break

def plot_sig_state_boxplots(mc_ad, out_dir_path, ordered_states, metadata_cols, sig_name_to_genes, genes_to_ignore=[], state_col='state'):
    import matplotlib.pyplot as plt
    import seaborn as sb

    pathlib.Path(out_dir_path).mkdir(parents=True, exist_ok=True)

    mc_utils.add_cell_scores(
        mc_ad, {k: v for k, v in sig_name_to_genes.items() if k in metadata_cols}, 
        # overwrite_existing_cols=True,
        genes_to_ignore=genes_to_ignore,
    )
    for col in metadata_cols:
        plt.close('all')
        fig, ax = plt.subplots()
        sb.boxplot(
            data=mc_ad.obs,
            y=state_col,
            x=col,
            order=ordered_states,
            palette=mc_utils.get_palette(mc_ad, color_by=state_col),
            showfliers=False,
            # linewidth=0.5,
            ax=ax,
        )
        fig.tight_layout()
        fig.savefig(f'{out_dir_path}/{col}_state_boxplots.png')

def plot_gate_c_state_scatterplots(c_ad, out_dir_path, scatter_name_to_info, palette, name_to_mask):
    scatter_name_to_mask_dict = {}
    for scatter_name, info in scatter_name_to_info.items():
        print(scatter_name)
        x = info['x']
        y = info['y']
        x_thresh = info['x_thresh']
        y_thresh = info['y_thresh']
        background_name = info['background_name']
        get_background_mask_func = info.get('get_background_mask_func', lambda *args: np.full(c_ad.n_obs, True))
        gates_to_target_mask_name = info.get('gates_to_target_mask_name', {})
        skip_plot = info.get('skip_plot', False)

        background_mask = get_background_mask_func(scatter_name_to_mask_dict, c_ad)
        # print(c_ad.obs.loc[background_mask, 'c_state'].value_counts())

        both_low_mask = (c_ad.obs[x] < x_thresh) & (c_ad.obs[y] < y_thresh) & background_mask
        both_high_mask = (c_ad.obs[x] >= x_thresh) & (c_ad.obs[y] >= y_thresh) & background_mask
        low_x_high_y_mask = (c_ad.obs[x] < x_thresh) & (c_ad.obs[y] >= y_thresh) & background_mask
        high_x_low_y_mask = (c_ad.obs[x] >= x_thresh) & (c_ad.obs[y] < y_thresh) & background_mask

        scatter_name_to_mask_dict[scatter_name] = dict(
            background_mask=background_mask,
            both_low=both_low_mask,
            both_high=both_high_mask,
            low_x_high_y=low_x_high_y_mask,
            high_x_low_y=high_x_low_y_mask,
        )
        if skip_plot:
            continue

        gate_to_percentile = dict(
            both_low=(both_low_mask.sum() / background_mask.sum()) * 100,
            both_high=(both_high_mask.sum() / background_mask.sum()) * 100,
            low_x_high_y=(low_x_high_y_mask.sum() / background_mask.sum()) * 100,
            high_x_low_y=(high_x_low_y_mask.sum() / background_mask.sum()) * 100,
        )
        gate_to_percentage_repr = ','.join(f'{k}={v:.1f}%' for k, v in gate_to_percentile.items())

        gates_to_target_percentage = {}
        for gates, target_mask_name in gates_to_target_mask_name.items():
            target_mask = name_to_mask[target_mask_name]
            gates_mask = np.full(c_ad.n_obs, False)
            for gate in gates:
                gates_mask |= scatter_name_to_mask_dict[scatter_name][gate]
            assert not ((~gates_mask) & target_mask).any(), f'{target_mask_name} was found outside {gates}'
            
            gates_to_target_percentage[gates] = target_mask.sum() / gates_mask.sum() * 100
        gates_to_target_percentage_repr = '\n'.join(f'{gates_to_target_mask_name[gates]} is {perc:.1f}% out of {gates}' for gates, perc in gates_to_target_percentage.items())
        if gates_to_target_percentage_repr:
            gates_to_target_percentage_repr = f'\n{gates_to_target_percentage_repr}'

        curr_c_ad_obs = c_ad.obs[background_mask]
        plt.close('all')
        fig, ax = plt.subplots(figsize=(7.5, 6.5))
        
        max_obs_count = int(10e3)
        max_obs_count = np.inf
        sb.scatterplot(
            data=curr_c_ad_obs.sample(n=max_obs_count) if len(curr_c_ad_obs) > max_obs_count else curr_c_ad_obs,
            x=x,
            y=y,
            hue='c_state',
            palette=palette,
            s=5,
            ax=ax,
        )
        ax.axvline(x=x_thresh, color='k', linestyle='-', alpha=0.2)
        ax.axhline(y=y_thresh, color='k', linestyle='-', alpha=0.2)
        ax.get_legend().remove()
        ax_title = f'background: {background_name}\n{gate_to_percentage_repr}{gates_to_target_percentage_repr}'
        ax.set_title(ax_title)
        fig.subplots_adjust(
            top=0.8,
        )
        # fig.tight_layout()
        fig.savefig(f'{out_dir_path}/EDF_9A_{scatter_name}.png', dpi=300)
        
        plt.close('all')


def get_donor_feature_data(c_ad, mc_ad, dataset_name, downsample_cells_target=None, convert_cytopenia_to_ccus_or_icus_according_to_ch=False, skip_ipssm=True):
    cell_freq_epsilon = get_mds_params()['epsilon_for_donor_log_c_state_freq']
    mds_path_dict = mds_analysis_params.get_mds_path_dict(dataset_name=dataset_name)

    # mc_utils.write_downsampled_layer_if_not_exists(c_ad, overwrite_if_existing_downsampled_layer_is_weird=True)
    orig_columns = set(c_ad.obs.columns)
    # c_ad.obs.drop(columns='gmp_l_elane_sig', inplace=True, errors='ignore')
    mc_utils.add_cell_scores(
        c_ad, mds_analysis_params.PB_SIG_NAME_TO_GENES, 
        # overwrite_existing_cols=True,
    )
    print(f'added: {set(c_ad.obs.columns) - orig_columns}')

    clinical_df = clinical_data_interface_and_utils.get_minimal_clinical_data_df()
    c_ad.obs = sc_rna_seq_preprocessing_params.get_df_with_numbered_donor_id(c_ad.obs)
    c_ad.obs = clinical_data_interface_and_utils.get_df_with_diagnosis_and_diagnosis_class_by_exp_name_and_donor_id_and_clinical_df(
        c_ad.obs, clinical_df)

    mc_utils.add_c_state_and_mc_c_state_stats(
        mc_ad, c_ad, 
        cell_state_and_info_list=get_mds_params()['pb_cd34_enriched_cell_state_and_c_info_list'], 
        mask_and_info_list=get_mds_params()['pb_cd34_enriched_mask_and_c_info_list'], 
        cell_type_colors_csv_file_path=get_mds_params()['cell_type_colors_csv_file_path'],
    )

    sc_rna_seq_preprocessing_params.add_exp_and_bleeding_date_according_to_exp_name_and_donor_id(c_ad.obs)

    extra_donor_data_df = c_ad.obs[['donor_id', 'exp_name']].drop_duplicates().astype(str)
    sc_rna_seq_preprocessing_params.add_exp_date_column_to_df(extra_donor_data_df)
    extra_donor_data_df.drop(columns='exp_name', inplace=True)
    extra_donor_data_df.drop_duplicates(inplace=True)
    
    # TODO: would be nicer to add ARCH cols later, after merging illu and ult.
    mutation_df = arch_mutation_interface_and_utils.get_minimal_all_mutation_df()
    agg_donor_small_scale_mut_df = arch_mutation_interface_and_utils.get_agg_donor_small_scale_mut_df(mutation_df)

    extra_donor_data_df = generic_utils.merge_preserving_df1_index_and_row_order(
        extra_donor_data_df, 
        clinical_df,
        on=['donor_id', 'exp_date'],
        how='left',
    )
    extra_donor_data_df = generic_utils.merge_preserving_df1_index_and_row_order(
        extra_donor_data_df, 
        agg_donor_small_scale_mut_df,
        on=['donor_id', 'exp_date'],
        how='left',
    )
    for col in agg_donor_small_scale_mut_df.columns:
        if col.endswith('_snv'):
            extra_donor_data_df[col].fillna(False, inplace=True)
        elif col.endswith('_VAF') or col.endswith('_mut_count'):
            extra_donor_data_df[col].fillna(0, inplace=True)
        elif col not in {'donor_id', 'exp_date'}:
            assert False, col

    extra_donor_data_df.rename(columns={'age': 'donor_age'}, inplace=True)
    extra_donor_data_df.reset_index(drop=True, inplace=True)


    str_col_names = [
        'donor_id', 'exp_name', 'metacell_name', 'projected_type', 'state', 'state_color', 'donor_sex', 'soup_donor_name', 'vireo_donor_name',
        'c_state', 'bleeding_date', 'exp_date',
    ]
    missing_str_col_names = [x for x in str_col_names if x not in c_ad.obs.columns]
    if missing_str_col_names:
        print(f'WARNING: c_ad missing_str_col_names: {missing_str_col_names}')
        str_col_names = [x for x in str_col_names if x not in missing_str_col_names]

    other_col_names = [
        'num_of_non_excluded_umis', 'log_norm_cell_ranger_umi_count', 'donor_age', 'is_ultima',
    ]
    other_col_names += sorted(get_mds_params()['pb_sig_name_to_info'])
        
    missing_other_col_names = [x for x in other_col_names if x not in c_ad.obs.columns]
    if missing_other_col_names:
        print(f'WARNING: c_ad missing_other_col_names: {missing_other_col_names}')
        other_col_names = [x for x in other_col_names if x not in missing_other_col_names]

    c_ad_obs = c_ad.obs[[*str_col_names, *other_col_names]].copy()
    c_ad_obs[str_col_names] = c_ad_obs[str_col_names].astype(str)
    c_ad_obs['log_num_of_non_excluded_umis'] = np.log2(c_ad_obs['num_of_non_excluded_umis'])


    if 'metacell_name' in c_ad_obs.columns:
        outlier_mask = c_ad.obs['metacell_name'] == 'Outliers'
        c_ad_obs.loc[outlier_mask, 'state'] = 'Outliers'
        c_ad_obs.loc[outlier_mask, 'projected_type'] = 'Outliers'
        assert (c_ad_obs['state'] != 'nan').all()
        assert (c_ad_obs['projected_type'] != 'nan').all()


    c_ad_obs = sc_rna_seq_preprocessing_params.get_df_with_numbered_donor_id(
        c_ad_obs, add_bm_according_to_exp_name=True,
        # consider_different_soup_vireo_donors=True,
    )
    c_ad_obs = generic_utils.merge_preserving_df1_index_and_row_order(c_ad_obs, c_ad_obs.groupby('exp_name').size().reset_index(name='exp_final_cell_count'))

    feature_cols = set()
    ds_c_ad_obs = (
        c_ad_obs[generic_utils.sample_mask_per_group(np.full(len(c_ad_obs), True), downsample_cells_target, c_ad_obs['numbered_donor_id'], drop_if_too_few_pos=True)] 
        if downsample_cells_target else c_ad_obs
    )


    if 'state' in ds_c_ad_obs.columns:
        states_missing_from_ordered_state_list = set(ds_c_ad_obs['state'].unique()) - set(mds_analysis_params.ORDERED_CELL_STATE_NAMES)
        assert not states_missing_from_ordered_state_list, f'states_missing_from_ordered_state_list: {states_missing_from_ordered_state_list}'

        state_freq_df = mc_utils.get_cell_state_freq_df(
            c_ad_obs=ds_c_ad_obs,
            group_by_column_names=['numbered_donor_id'],
            cell_state_column_name='state',
            add_count_column=True,
        )
    else:
        state_freq_df = None
    if 'c_state' in ds_c_ad_obs.columns:
        c_state_freq_df = mc_utils.get_cell_state_freq_df(
            c_ad_obs=ds_c_ad_obs,
            group_by_column_names=['numbered_donor_id'],
            cell_state_column_name='c_state',
            add_count_column=True,
        )
    else:
        c_state_freq_df = None
    if 'projected_type' in ds_c_ad_obs.columns:
        projected_type_freq_df = mc_utils.get_cell_state_freq_df(
            c_ad_obs=ds_c_ad_obs,
            group_by_column_names=['numbered_donor_id'],
            cell_state_column_name='projected_type',
        )
    else:
        projected_type_freq_df = None
    
    if 'state' in ds_c_ad_obs.columns:
        state_median_log_num_of_non_excluded_umis_df = ds_c_ad_obs.groupby(['numbered_donor_id', 'state'])['log_num_of_non_excluded_umis'].median().reset_index()
        median_log_norm_cell_ranger_umi_count_df = ds_c_ad_obs.groupby(
            ['numbered_donor_id', 'state'])['log_norm_cell_ranger_umi_count'].median().reset_index()
    if 'c_state' in ds_c_ad_obs.columns:
        c_state_median_log_num_of_non_excluded_umis_df = ds_c_ad_obs.groupby(['numbered_donor_id', 'c_state'])['log_num_of_non_excluded_umis'].median().reset_index()
    
    all_median_log_num_of_non_excluded_umis_df = ds_c_ad_obs.groupby(['numbered_donor_id'])['log_num_of_non_excluded_umis'].median().to_dict()
    donor_info_df = ds_c_ad_obs[['numbered_donor_id', 'bleeding_date', 'donor_id', 'exp_name', 'exp_date']].drop_duplicates()
    numbered_donor_id_to_cell_count = ds_c_ad_obs['numbered_donor_id'].value_counts().to_dict()

    num_of_donors = ds_c_ad_obs['numbered_donor_id'].nunique()
    flat_dicts = []
    for i, numbered_donor_id in enumerate(sorted(ds_c_ad_obs['numbered_donor_id'].unique())):
        # print(f'{numbered_donor_id} ({i}/{num_of_donors})')
        if 'projected_type' in ds_c_ad_obs.columns:
            projected_composition_dict = {
                f'projected___{projected_type}': freq 
                for projected_type, freq in projected_type_freq_df.loc[projected_type_freq_df['numbered_donor_id'] == numbered_donor_id, ['projected_type', 'freq']].to_records(index=False)
            }
        else:
            projected_composition_dict = {}
                
        if 'state' in ds_c_ad_obs.columns:
            state_median_log_num_of_non_excluded_umis_dict = {
                f'{state}_median_log_num_of_non_excluded_umis': v 
                for state, v in state_median_log_num_of_non_excluded_umis_df.loc[state_median_log_num_of_non_excluded_umis_df['numbered_donor_id'] == numbered_donor_id, ['state', 'log_num_of_non_excluded_umis']].to_records(index=False)
            }
            state_median_log_norm_cell_ranger_umi_count_dict = {
                f'{state}_median_log_norm_cell_ranger_umi_count': v 
                for state, v in median_log_norm_cell_ranger_umi_count_df.loc[median_log_norm_cell_ranger_umi_count_df['numbered_donor_id'] == numbered_donor_id, ['state', 'log_norm_cell_ranger_umi_count']].to_records(index=False)
            }
            hspc_composition_dict = {
                state: freq 
                for state, freq in state_freq_df.loc[state_freq_df['numbered_donor_id'] == numbered_donor_id, ['state', 'freq']].to_records(index=False)
            }
            hspc_composition_count_dict = {
                f'{state}_count': count 
                for state, count in state_freq_df.loc[state_freq_df['numbered_donor_id'] == numbered_donor_id, ['state', 'count']].to_records(index=False)
            }
        else:
            state_median_log_num_of_non_excluded_umis_dict = {}
            state_median_log_norm_cell_ranger_umi_count_dict = {}
            hspc_composition_dict = {}
            hspc_composition_count_dict = {}
        
        if 'c_state' in ds_c_ad_obs.columns:
            c_state_median_log_num_of_non_excluded_umis_dict = {
                f'c_{state}_median_log_num_of_non_excluded_umis': v 
                for state, v in c_state_median_log_num_of_non_excluded_umis_df.loc[c_state_median_log_num_of_non_excluded_umis_df['numbered_donor_id'] == numbered_donor_id, ['c_state', 'log_num_of_non_excluded_umis']].to_records(index=False)
            }
            hspc_c_composition_dict = {
                f'c_{state}': freq 
                for state, freq in c_state_freq_df.loc[c_state_freq_df['numbered_donor_id'] == numbered_donor_id, ['c_state', 'freq']].to_records(index=False)
            }
            hspc_c_composition_count_dict = {
                f'c_{state}_count': count 
                for state, count in c_state_freq_df.loc[c_state_freq_df['numbered_donor_id'] == numbered_donor_id, ['c_state', 'count']].to_records(index=False)
            }
        else:
            c_state_median_log_num_of_non_excluded_umis_dict = {}
            hspc_c_composition_dict = {}
            hspc_c_composition_count_dict = {}
        
        donor_row = donor_info_df.loc[donor_info_df['numbered_donor_id'] == numbered_donor_id, ['bleeding_date', 'donor_id', 'exp_name', 'exp_date']]
        assert len(donor_row) == 1
        donor_info = donor_row.iloc[0].to_dict()
        cell_count = numbered_donor_id_to_cell_count[numbered_donor_id]
        flat_dicts.append({
            'numbered_donor_id': numbered_donor_id,
            'cell_count': cell_count,
            'median_log_num_of_non_excluded_umis': all_median_log_num_of_non_excluded_umis_df[numbered_donor_id],
            **donor_info,
            **projected_composition_dict,
            **hspc_composition_dict,
            **hspc_c_composition_dict,
            **hspc_composition_count_dict,
            **hspc_c_composition_count_dict,
            **state_median_log_num_of_non_excluded_umis_dict,
            **state_median_log_norm_cell_ranger_umi_count_dict,
            **c_state_median_log_num_of_non_excluded_umis_dict,
        })

    
    df = pd.DataFrame(flat_dicts)
    
    for column_name in {
        *(state_freq_df['state'] if (state_freq_df is not None) else []),
        *([f'{x}_count' for x in state_freq_df['state']] if (state_freq_df is not None) else []),
        *([f'c_{x}' for x in c_state_freq_df['c_state']] if (c_state_freq_df is not None) else []),
        *([f'c_{x}_count' for x in c_state_freq_df['c_state']] if (c_state_freq_df is not None) else []),
        *([f'projected___{x}' for x in projected_type_freq_df['projected_type']] if (projected_type_freq_df is not None) else []),
    } & set(df.columns):
        df[column_name].fillna(0, inplace=True)

    for prefix in [
        '',
        'c_',
    ]:
        total_hspc_col_name = f'{prefix}HSPC'
        total_non_hspc_col_name = f'{prefix}non_HSPC'

        curr_ordered_hspc_state_col_names = [f'{prefix}{x}' for x in mds_analysis_params.ORDERED_PB_HSPC_STATE_NAMES]
        ordered_observed_hspc_state_col_names = [x for x in curr_ordered_hspc_state_col_names if x in df.columns]
        assert ordered_observed_hspc_state_col_names
        df[total_hspc_col_name] = df[ordered_observed_hspc_state_col_names].sum(axis=1)
        assert (df[total_hspc_col_name] <= 1).all()
        curr_ordered_non_hspc_state_col_names = [f'{prefix}{x}' for x in mds_analysis_params.ORDERED_PB_NON_HSPC_STATE_NAMES]
        ordered_observed_non_hspc_state_col_names = [x for x in curr_ordered_non_hspc_state_col_names if x in df.columns]
        df[total_non_hspc_col_name] = df[ordered_observed_non_hspc_state_col_names].sum(axis=1)
        assert (df[total_non_hspc_col_name] <= 1).all()
        # assert ((df[total_hspc_col_name] <= 1) | np.isclose(df[total_hspc_col_name], 1)).all()

        any_hspc_mask = df[total_hspc_col_name] > 0
        donors_without_any_hspcs_numbered_donor_ids = sorted(df.loc[~any_hspc_mask, 'numbered_donor_id'].unique())

        # i guess this should always hold. nope. but why should i care?
        # assert not donors_without_any_hspcs_numbered_donor_ids, f'donors_without_any_hspcs_numbered_donor_ids: {donors_without_any_hspcs_numbered_donor_ids}'
        for state_col in ordered_observed_hspc_state_col_names:
            df.loc[any_hspc_mask, state_col] /= df.loc[any_hspc_mask, total_hspc_col_name]
            # df.loc[any_hspc_mask, state_col] = np.log2(df.loc[any_hspc_mask, state_col])
            log_state_col = f'log_{state_col}'
            df[log_state_col] = np.log2(df[state_col] + cell_freq_epsilon)
            if prefix == 'c_':
                feature_cols.add(log_state_col)
        
        any_non_hspc_mask = df[total_non_hspc_col_name] > 0
        donors_with_only_hspcs_numbered_donor_ids = sorted(df.loc[~any_non_hspc_mask, 'numbered_donor_id'].unique())
        
        for state_col in ordered_observed_non_hspc_state_col_names:
            df.loc[any_non_hspc_mask, state_col] /= df.loc[any_non_hspc_mask, total_non_hspc_col_name]
            assert (df.loc[~any_non_hspc_mask, state_col] == 0).all()
    
    
    assert set(df.columns) & set(extra_donor_data_df.columns) == {'donor_id', 'exp_date'}
    
    
    def filter_df_by_enough_cells_per_numbered_donor(curr_c_ad_obs, min_count=10):
        enough_cells_numbered_donor_ids = (curr_c_ad_obs['numbered_donor_id'].value_counts() >= min_count).loc[lambda x: x].index
        return curr_c_ad_obs[curr_c_ad_obs['numbered_donor_id'].isin(enough_cells_numbered_donor_ids)]

    if 'c_state' in c_ad_obs.columns:
        pb_sig_default_min_donor_cell_count = get_mds_params()['pb_sig_default_min_donor_cell_count']
        pb_sig_c_states_repr_to_min_donor_cell_count = get_mds_params()['pb_sig_c_states_repr_to_min_donor_cell_count']
        pb_weak_sig_min_donor_cell_count = get_mds_params()['pb_weak_sig_min_donor_cell_count']
        pb_sig_quantiles = get_mds_params()['pb_sig_quantiles']
        cols_before_adding_quantiles = set(df.columns)
        c_state_repr_to_c_ad_obs_for_strong_sigs = {}
        c_state_repr_to_c_ad_obs_for_weak_sigs = {}
        for sig_name, sig_info in get_mds_params()['pb_sig_name_to_info'].items():
            print(f'starting work on {sig_name}')
            relevant_state_lists = sig_info['relevant_states']
            if isinstance(relevant_state_lists[0], str):
                relevant_state_lists = [relevant_state_lists]

            c_state_sets_to_calc_statistics = sig_info.get('c_state_sets_to_calc_statistics', [set(x) for x in relevant_state_lists])
            for c_states in c_state_sets_to_calc_statistics:
                assert isinstance(c_states, set)
                assert not (c_states & {'HSC', 'GMP-E', 'MPP', 'MEBEMP-E', 'CLP-M', 'CLP-L'})
                c_states_repr = '_'.join(sorted(c_states))
                curr_min_donor_cell_count = pb_sig_c_states_repr_to_min_donor_cell_count.get(c_states_repr, pb_sig_default_min_donor_cell_count)
                print(f'{c_states_repr}: curr_min_donor_cell_count={curr_min_donor_cell_count}')

                if c_states_repr not in c_state_repr_to_c_ad_obs_for_strong_sigs:
                    c_state_repr_to_c_ad_obs_for_strong_sigs[c_states_repr] = filter_df_by_enough_cells_per_numbered_donor(
                        c_ad_obs[c_ad_obs['c_state'].isin(list(c_states))], min_count=curr_min_donor_cell_count,
                    )
                curr_c_ad_obs = c_state_repr_to_c_ad_obs_for_strong_sigs[c_states_repr]
                # curr_df = curr_c_ad_obs.groupby('numbered_donor_id')[sig_name].median().reset_index(name=f'{s_phase_sig}_{c_states_repr}_')
                curr_df = curr_c_ad_obs.groupby('numbered_donor_id')[sig_name].quantile(pb_sig_quantiles).unstack().reset_index().rename(
                    columns={x: f'{sig_name}_{c_states_repr}_{x}' for x in pb_sig_quantiles})
                df = generic_utils.merge_preserving_df1_index_and_row_order(df, curr_df, how='left')
                
                if c_states_repr not in c_state_repr_to_c_ad_obs_for_weak_sigs:
                    c_state_repr_to_c_ad_obs_for_weak_sigs[c_states_repr] = filter_df_by_enough_cells_per_numbered_donor(
                        c_ad_obs[c_ad_obs['c_state'].isin(list(c_states))], min_count=pb_weak_sig_min_donor_cell_count,
                    )
                curr_c_ad_obs = c_state_repr_to_c_ad_obs_for_weak_sigs[c_states_repr]
                curr_df = curr_c_ad_obs.groupby('numbered_donor_id')[sig_name].mean().reset_index(name=f'{sig_name}_{c_states_repr}_mean')
                df = generic_utils.merge_preserving_df1_index_and_row_order(df, curr_df, how='left')
        feature_cols |= set(df.columns) - cols_before_adding_quantiles
    

    ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
        df, extra_donor_data_df, on=['donor_id', 'exp_date'], how='left')
    
    col_names = [
        'donor_id', 'exp_name', 
        'donor_sex', 
        'donor_age', 
        'is_ultima', 'bleeding_date', 'exp_date',
        'numbered_donor_id', # 240205: previously discarded this one. no idea why.
        'numbered_donor_id_without_tech_rep_suffix', 
        'exp_final_cell_count',
    ]
    missing_col_names = [x for x in col_names if x not in ds_c_ad_obs.columns]
    if missing_col_names:
        print(f'WARNING: ext_donor_feature_df missing_col_names: {missing_col_names}')
    col_names = [x for x in col_names if x not in missing_col_names]
    
    

    ds_c_ad_obs['donor_age'] = ds_c_ad_obs['donor_age'].astype(float) # to fix an error in the next line (caused by attempting to merge on donor_age which is float64 in one df and object in the other)
    # def is_allowed_conflict(col, vec1, vec2):
    #     if col == 'donor_age':
    #         return (np.abs(vec1 - vec2) <= 1) | (vec2 == -3) # as i used -3 as a dummy value for some time, which was a bad call.
    #     return vec1 == vec2
    ext_donor_feature_df = generic_utils.replace_vals_by_join(
        ext_donor_feature_df, ds_c_ad_obs.drop_duplicates(subset=['donor_id', 'exp_name'])[col_names], 
        cols_to_join_on=['donor_id', 'exp_name', 'numbered_donor_id', 'exp_date', 'bleeding_date'],
        # is_allowed_conflict=is_allowed_conflict,
    )

    # TODO: would be nicer to add CNA cols later, after merging illu and ult.
    donor_agg_cna_df_csv_file_path = get_mds_params()['donor_agg_cna_df_csv_file_path']
    if os.path.isfile(donor_agg_cna_df_csv_file_path):
        donor_agg_cna_df = pd.read_csv(donor_agg_cna_df_csv_file_path)
        all_exp_names_mask = donor_agg_cna_df['exp_name'] == 'all'

        donor_agg_cna_df_all_exp_names = donor_agg_cna_df.loc[
            all_exp_names_mask, donor_agg_cna_df.columns[donor_agg_cna_df.columns != 'exp_name']]
        donor_agg_cna_df_all_exp_names = donor_agg_cna_df_all_exp_names.merge(ext_donor_feature_df[['donor_id', 'exp_name']])
        donor_agg_cna_df = pd.concat([donor_agg_cna_df[~all_exp_names_mask], donor_agg_cna_df_all_exp_names], ignore_index=True)

        # print(f'donor_agg_cna_df:\n{donor_agg_cna_df}')

        ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
            ext_donor_feature_df, 
            donor_agg_cna_df, 
            how='left',
            on=['donor_id', 'exp_name'],
        )

        for col in donor_agg_cna_df.columns:
            if col.startswith('has_del'):
                ext_donor_feature_df[col].fillna(False, inplace=True)
            elif col == 'ipssr_cytogenetic_risk_group':
                ext_donor_feature_df[col].fillna('Good', inplace=True)
            elif col == 'cna_count':
                ext_donor_feature_df[col].fillna(0, inplace=True)
            elif col not in {'donor_id', 'exp_name'}:
                assert False, col

    else:
        donor_agg_cna_df = None
        ext_donor_feature_df['cna_count'] = np.nan
        ext_donor_feature_df['ipssr_cytogenetic_risk_group'] = np.nan


    

    ext_donor_feature_df = (
        clinical_data_interface_and_utils.get_df_with_diagnosis_and_diagnosis_class_by_exp_name_and_donor_id_and_clinical_df(
            ext_donor_feature_df, clinical_df))

    ext_donor_feature_df['diagnosis_class'].fillna('nan', inplace=True)

    valid_age_mask = ext_donor_feature_df['donor_age'] >= 0
    ext_donor_feature_df.loc[valid_age_mask, 'log_donor_age'] = np.log2(ext_donor_feature_df.loc[valid_age_mask, 'donor_age'])

    sex_known_mask = ext_donor_feature_df['donor_sex'].isin(['male', 'female'])
    ext_donor_feature_df.loc[sex_known_mask, 'is_female_as_int'] = (ext_donor_feature_df.loc[sex_known_mask, 'donor_sex'] == 'female').astype(int)

    if 'max_mean_VAF' in ext_donor_feature_df.columns:
        ext_donor_feature_df['CH'] = (
            (ext_donor_feature_df['max_mean_VAF'] >= 0.02) # 0.02 according to WHO 2022 (https://www.nature.com/articles/s41375-022-01613-1, 2022)
            | (ext_donor_feature_df['cna_count'] > 0)
        )

    if convert_cytopenia_to_ccus_or_icus_according_to_ch:
        ext_donor_feature_df.loc[
            (ext_donor_feature_df['diagnosis_class'] == 'cytopenia') & ext_donor_feature_df['CH'], 'diagnosis_class'] = 'CCUS'
        ext_donor_feature_df.loc[
            (ext_donor_feature_df['diagnosis_class'] == 'cytopenia') & (~ext_donor_feature_df['CH']), 'diagnosis_class'] = 'ICUS'

    ext_donor_feature_df['sex_and_diagnosis_class'] = (
        ext_donor_feature_df['donor_sex'].astype(str) + ',' + ext_donor_feature_df['diagnosis_class'].astype(str)
    )

    ext_donor_feature_df['MF'] = ext_donor_feature_df['diagnosis'].isin(clinical_data_interface_and_utils.MF_DIAGNOSES)
    ext_donor_feature_df['MF?'] = ext_donor_feature_df['diagnosis'].isin(clinical_data_interface_and_utils.MF_SUSPECT_DIAGNOSES)
    ext_donor_feature_df['ET'] = ext_donor_feature_df['diagnosis'].isin(clinical_data_interface_and_utils.ET_DIAGNOSES)
    ext_donor_feature_df['CMML'] = ext_donor_feature_df['diagnosis'].isin(clinical_data_interface_and_utils.CMML_DIAGNOSES)
    ext_donor_feature_df['PV'] = ext_donor_feature_df['diagnosis'] == 'PV'

    ext_donor_feature_df.loc[ext_donor_feature_df['MF'], 'mpn_diagnosis_class'] = 'MF'
    ext_donor_feature_df.loc[ext_donor_feature_df['ET'], 'mpn_diagnosis_class'] = 'ET'
    ext_donor_feature_df.loc[ext_donor_feature_df['PV'], 'mpn_diagnosis_class'] = 'PV'
    ext_donor_feature_df.loc[ext_donor_feature_df['MF?'], 'mpn_diagnosis_class'] = 'MF?'
    ext_donor_feature_df.loc[ext_donor_feature_df['diagnosis_class'] == 'normal', 'mpn_diagnosis_class'] = 'normal'

    mpn_missing_mpn_diagnosis_class_mask = ext_donor_feature_df['mpn_diagnosis_class'].isna() & (ext_donor_feature_df['diagnosis_class'] == 'MPN')
    mpn_diagnoses_missing_mpn_diagnosis_class = ext_donor_feature_df.loc[mpn_missing_mpn_diagnosis_class_mask, 'diagnosis'].drop_duplicates()
    assert mpn_diagnoses_missing_mpn_diagnosis_class.isin(clinical_data_interface_and_utils.NON_MPN_CYTOSIS).all(), str(sorted(mpn_diagnoses_missing_mpn_diagnosis_class))

    bm_mask = ext_donor_feature_df['exp_name'].str.contains('_bm_')
    assert (bm_mask == (ext_donor_feature_df['numbered_donor_id'].str.contains('_bm'))).all()
    ext_donor_feature_df['is_bm'] = bm_mask

            
    preprocessing_out_exp_df_csv_file_path = mds_path_dict['preprocessing_out_exp_df_csv_file_path']
    if os.path.isfile(preprocessing_out_exp_df_csv_file_path):
        ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
            ext_donor_feature_df, pd.read_csv()[['exp_name', 'num_of_cells']].rename(columns={'num_of_cells': 'exp_cell_count'}), how='left')
        ext_donor_feature_df['fraction_of_exp'] = ext_donor_feature_df['cell_count'] / ext_donor_feature_df['exp_cell_count']
    if 'HSPC' in ext_donor_feature_df.columns:
        ext_donor_feature_df['HSPC_count'] = ext_donor_feature_df['cell_count'] * ext_donor_feature_df['HSPC']
        ext_donor_feature_df['log_HSPC_count'] = np.log2(ext_donor_feature_df['HSPC_count'] + 1e-4)
        ext_donor_feature_df['non_HSPC_count'] = ext_donor_feature_df['cell_count'] * (1 - ext_donor_feature_df['HSPC'])
        ext_donor_feature_df['log_non_HSPC_count'] = np.log2(ext_donor_feature_df['non_HSPC_count'] + 1e-4)

    if 'c_HSC_MPP_median_log_num_of_non_excluded_umis' in ext_donor_feature_df.columns:
        # c_HSC_MPP_median_log_num_of_non_excluded_umis because on average, HSC_MPP is the most common c_state.
        print('\n\nAdding c_state_HSC_MPP_median_umi_count_log_ratio columns\n\n')
        possible_c_state_to_median_col = {x: f'c_{x}_median_log_num_of_non_excluded_umis' for x in mds_analysis_params.PB_HSPC_STATE_NAMES}
        curr_c_state_to_median_col = {k: v for k, v in possible_c_state_to_median_col.items() if v in ext_donor_feature_df.columns}
        for c_state, median_umi_count_col in curr_c_state_to_median_col.items():
            if c_state != 'HSC_MPP':
                ext_donor_feature_df[f'c_{c_state}_HSC_MPP_median_umi_count_log_ratio'] = (
                    ext_donor_feature_df[median_umi_count_col] - 
                    ext_donor_feature_df[f'c_HSC_MPP_median_log_num_of_non_excluded_umis']
                )


    if 'HSPC' in ext_donor_feature_df.columns:
        ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
            ext_donor_feature_df, 
            ext_donor_feature_df.groupby('exp_name')['HSPC_count'].sum().reset_index(
                name='exp_curr_total_HSPC_count'),
            verbose=False,
        )
        ext_donor_feature_df['HSPC_fraction_of_curr_exp_HSPC'] = (
            ext_donor_feature_df['HSPC_count'] / ext_donor_feature_df['exp_curr_total_HSPC_count'])

        ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
            ext_donor_feature_df, 
            ext_donor_feature_df.groupby('exp_name')['HSPC'].median().reset_index(name='exp_median_HSPC'),
            verbose=False,
        )
        ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
            ext_donor_feature_df, 
            ext_donor_feature_df.groupby('exp_name')['HSPC'].min().reset_index(name='exp_min_HSPC'),
            verbose=False,
        )
    ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
        ext_donor_feature_df, 
        ext_donor_feature_df.groupby('exp_name')['donor_id'].nunique().reset_index(name='curr_num_of_donors_in_exp'),
        verbose=False,
    )

    exp_donor_hpsc_count_estimation_df = pd.read_csv(mds_path_dict['presumably_hspc_count_df_csv_file_path'])
    ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
        ext_donor_feature_df,
        exp_donor_hpsc_count_estimation_df.rename(columns={'presumably_hspc_count': 'hspc_count_estimated_after_cell_exclusion'}),
    )

    ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
        ext_donor_feature_df,
        c_ad.obs.loc[c_ad.obs['not_dc_nkt_monocyte_endothel_b_lamp3'], ['donor_id', 'exp_name']].value_counts().reset_index(name='final_c_hspc_count'),
    )
        
    assert ext_donor_feature_df['numbered_donor_id'].is_unique
    ext_donor_feature_df['biased_composition_due_to_discarded_low_umi_count_barcodes'] = ext_donor_feature_df['exp_name'].isin(
        sc_rna_seq_preprocessing_params.get_biased_composition_due_to_discarded_low_umi_count_barcodes_exp_names())    


    ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
        ext_donor_feature_df,
        mc_utils.get_c_mask_freq_per_obs_column(
            c_ad, 
            c_ad.obs['CLP-E'],
            ['donor_id', 'exp_name'],
            # background_c_mask=c_ad.obs['HSC_CLP_traj'],
            background_c_mask=c_ad.obs['not_dc_nkt_monocyte_endothel_b_lamp3'], # 240110: maybe better than HSC_CLP_traj. maybe makes more sense because this background is less arbitrary?
        )[['donor_id', 'exp_name', 'freq']].rename(columns={'freq': 'old_c_clp_e'}),
        on=['donor_id', 'exp_name'],
    )
    ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(
        ext_donor_feature_df,
        mc_utils.get_c_mask_freq_per_obs_column(
            c_ad, 
            (
                c_ad.obs['c_state'].isin(['CLP', 'HSC_MPP']) &
                (c_ad.obs['higher_in_clp_m_than_all_myeloid_hspcs'] > pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_CLP_E_LOW_THRESH) &
                (c_ad.obs['higher_in_clp_m_than_all_myeloid_hspcs'] < pb_cd34_c_score_threshs.HIGHER_IN_CLP_M_THAN_ALL_MYELOID_HSPCS_CLP_E_HIGH_THRESH)),
            ['donor_id', 'exp_name'],
            background_c_mask=c_ad.obs['not_dc_nkt_monocyte_endothel_b_lamp3'],
        )[['donor_id', 'exp_name', 'freq']].rename(columns={'freq': 'c_clp_e'}),
        on=['donor_id', 'exp_name'],
    )
    ext_donor_feature_df['log_c_clp_e'] = np.log2(ext_donor_feature_df['c_clp_e'] + cell_freq_epsilon)

    for state1, state2 in mds_analysis_params.NEIGHBOR_STATE_PAIRS:
        log_ratio_col = f'log_ratio_c_{state2}_{state1}'
        ext_donor_feature_df[log_ratio_col] = ext_donor_feature_df[f'log_c_{state2}'] - ext_donor_feature_df[f'log_c_{state1}']
    ext_donor_feature_df['complex_karyo'] = ext_donor_feature_df['cna_count'] >= 3
    binary_snv_cols = [x for x in ext_donor_feature_df.columns if x.endswith('_snv')]
    ext_donor_feature_df['snv_gene_count'] = ext_donor_feature_df[binary_snv_cols].sum(axis=1)

    # TODO: would be nicer to add ipssm cols later, after merging illu and ult.
    if not skip_ipssm:
        raise NotImplementedError('did not update this code, so it currently fails')
        ipssm_df = ipssm_utils.get_ipssm_final_cols_and_scores_df(ext_donor_feature_df, mc_model_paths['ipssm_df_csv'], mc_model_paths['ipssm_res_df_csv'])
        ext_donor_feature_df = generic_utils.merge_preserving_df1_index_and_row_order(ext_donor_feature_df, ipssm_df)


    ext_donor_feature_df.drop(columns='numbered_donor_id', inplace=True) # numbered_donor_id is dangerous because it might not be the same if you take different tech reps... 

    return dict(
        ext_donor_feature_df=ext_donor_feature_df,
        feature_cols=feature_cols,
    )
      

def write_donor_diff_expr(
        c_ad, mc_ad, root_out_dir_path, 
        diagnosis_classes_to_min_donor_count, # this is including bio reps
        releavant_c_states=None,
        per_donor_rather_than_per_donor_sample=False, 
        min_cell_count=20, 
):
    out_dir_path = os.path.join(root_out_dir_path, 'per_donor' if per_donor_rather_than_per_donor_sample else 'per_donor_sample')

    mc_utils.write_expr_and_expr_enrich(mc_ad)
    mc_utils.add_c_state_and_mc_c_state_stats(
        mc_ad, c_ad, 
        cell_state_and_info_list=get_mds_params()['pb_cd34_enriched_cell_state_and_c_info_list'], 
        mask_and_info_list=get_mds_params()['pb_cd34_enriched_mask_and_c_info_list'], 
        cell_type_colors_csv_file_path=get_mds_params()['cell_type_colors_csv_file_path'],
        only_add_c_state=True,
    )

    minimal_c_ad_obs = c_ad.obs[['donor_id', 'exp_name', 'bleeding_date', 'diagnosis', 'is_ultima', 'num_of_non_excluded_umis', 'is_bm']].copy()
    minimal_c_ad_obs['donor_id_and_exp_name'] = minimal_c_ad_obs['donor_id'].astype(str) + '__' + minimal_c_ad_obs['exp_name'].astype(str)

    if per_donor_rather_than_per_donor_sample:
        mask_of_cells_allowing_pooling_all_cells_per_donor = sc_rna_seq_preprocessing_params.get_mask_of_cells_allowing_pooling_all_pb_cells_per_donor(
            minimal_c_ad_obs)
        curr_donor_ids = sorted(minimal_c_ad_obs.loc[mask_of_cells_allowing_pooling_all_cells_per_donor, 'donor_id'].unique())
        donor_or_donor_exp_to_c_mask = {
            x: ((minimal_c_ad_obs['donor_id'] == x) & mask_of_cells_allowing_pooling_all_cells_per_donor)
            for x in curr_donor_ids
        }
    else:
        donor_or_donor_exp_to_c_mask = {x: minimal_c_ad_obs['donor_id_and_exp_name'] == x for x in minimal_c_ad_obs['donor_id_and_exp_name'].unique()}

    clinical_df = clinical_data_interface_and_utils.get_minimal_clinical_data_df()
    c_ad.obs = sc_rna_seq_preprocessing_params.get_df_with_numbered_donor_id(c_ad.obs)
    c_ad.obs = clinical_data_interface_and_utils.get_df_with_diagnosis_and_diagnosis_class_by_exp_name_and_donor_id_and_clinical_df(
        c_ad.obs, clinical_df)
    misc_numbered_donor_info = sc_rna_seq_preprocessing_params.get_misc_numbered_donor_info(c_ad)
    numbered_donor_df = misc_numbered_donor_info['numbered_donor_df']

    plt.close('all')
    
    observed_c_states = set(c_ad.obs['c_state'].unique())
    if releavant_c_states is None:
        releavant_c_states = observed_c_states
    else:
        releavant_c_states = set(releavant_c_states) & observed_c_states
    releavant_c_states = sorted(releavant_c_states)

    c_state_to_c_state_repr = {x: generic_utils.replace_special_chars_with_underscores(x) for x in releavant_c_states}
    assert len(set(c_state_to_c_state_repr.values())) == len(c_state_to_c_state_repr)
    states_repr_to_state_set_info = {
        c_state_repr: {
            'c_state': c_state,
            'states': None,
            'c_mask': c_ad.obs['c_state'] == c_state,
        }
        for c_state, c_state_repr in c_state_to_c_state_repr.items()
    }
    
    states_reprs = []
    for states_repr, state_set_info in states_repr_to_state_set_info.items():
        # if states_repr != 'all_states': 
        #     continue
        # print(f'states_repr: {states_repr}')

        # continue
        states = state_set_info['states']
        c_state = state_set_info['c_state']
        # raise
        if 'c_mask' in state_set_info:
            general_c_mask = state_set_info['c_mask']
        elif states == 'all_states':
            general_c_mask = np.full(c_ad.n_obs, True)
        elif states == 'all_HSPCs':
            general_c_mask = c_ad.obs['state'].isin(mds_analysis_params.PB_HSPC_STATE_NAMES)
        else:
            general_c_mask = c_ad.obs['state'].isin(states)
        

        donor_cell_count_df = c_ad.obs.loc[general_c_mask, ['exp_name', 'donor_id']].value_counts().reset_index(name='donor_cell_count')
        donor_cell_count_df = generic_utils.merge_preserving_df1_index_and_row_order(donor_cell_count_df, numbered_donor_df[['exp_name', 'donor_id', 'diagnosis_class', 'numbered_donor_id_without_tech_rep_suffix']], on=['exp_name', 'donor_id'])
        filtered_donor_cell_count_df = donor_cell_count_df[donor_cell_count_df['donor_cell_count'] >= min_cell_count]
        
        enough_donor_samples = True
        for diagnosis_classes, min_donor_count in diagnosis_classes_to_min_donor_count.items():
            donor_count = filtered_donor_cell_count_df.loc[
                filtered_donor_cell_count_df['diagnosis_class'].isin(list(diagnosis_classes)), 'numbered_donor_id_without_tech_rep_suffix'].nunique()
            # print(f'{diagnosis_classes}: {donor_count}')
            if donor_count < min_donor_count:
                enough_donor_samples = False
                break
        if not enough_donor_samples:
            continue
        
        states_reprs.append(states_repr)
        print(f'states_repr: {states_repr}')
        
        # continue

        curr_out_dir_path = os.path.join(out_dir_path, f'{states_repr}')
        cache_dir_path = os.path.join(curr_out_dir_path, f'{states_repr}_cache')
        pathlib.Path(curr_out_dir_path).mkdir(parents=True, exist_ok=True)
        pathlib.Path(cache_dir_path).mkdir(parents=True, exist_ok=True)

        curr_log_ratio_df_csv_file_path = os.path.join(curr_out_dir_path, f'log_ratio_df.csv')
        curr_c_mc_log_ratio_df_csv_file_path = os.path.join(curr_out_dir_path, f'c_mc_log_ratio_df.csv')
        curr_cell_count_df_csv_file_path = os.path.join(curr_out_dir_path, f'cell_count_df.csv')
        curr_proj_expr_df_csv_file_path = os.path.join(curr_out_dir_path, f'proj_expr_df.csv')
        curr_c_expr_df_csv_file_path = os.path.join(curr_out_dir_path, f'c_expr_df.csv')

        c_proj_log_ratio_df, cell_count_df = mc_utils.get_c_proj_log_ratio_df(
            c_ad, mc_ad, name_to_c_mask=donor_or_donor_exp_to_c_mask, general_c_mask=general_c_mask, out_dir_path=curr_out_dir_path,
            # use_existing_df_csv_file=True,
            min_num_of_cells=min_cell_count,
            cache_dir_path=cache_dir_path,
            log_ratio_df_csv_file_path=curr_log_ratio_df_csv_file_path,
            cell_count_df_csv_file_path=curr_cell_count_df_csv_file_path,
            c_mc_log_ratio_df_csv_file_path=curr_c_mc_log_ratio_df_csv_file_path,
            proj_expr_df_csv_file_path=curr_proj_expr_df_csv_file_path,
            c_expr_df_csv_file_path=curr_c_expr_df_csv_file_path,
            # epsilon_for_c_log=1/7, #  230824: actually i think 1/7 is worse even now (at least for N227 and N229), both when comparing c_expr to mc_expr and c_expr to proj_expr (mc_expr and proj_expr are not affected by this change)... 230813: remove this when we move to an atlas that was calculated later (than today - i don't know when exactly the default value of metacell_umis_regularization (which had a different name previously) changed from 1/7 to 1/16 (collect_metacells uses it, and accordingly the script that converts old MC h5ad files to new ones, as it calls collect_metacells).
        )
        # c_expr_df = pd.read_csv(curr_c_expr_df_csv_file_path, index_col=0)
    generic_utils.write_text_file(os.path.join(out_dir_path, 'states_reprs.txt'), '\n'.join(states_reprs))

def write_donor_diff_expr_mw_df(
        donor_diff_expr_dir_path,
        numbered_donor_df,
        max_bh_fdr_corrected_mw_pval=0.1,
        obs_expr_quantiles=[0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 0.99],
        min_99_percentile_expr=-16,
):
    
    states_reprs = generic_utils.read_text_file(os.path.join(donor_diff_expr_dir_path, 'states_reprs.txt')).splitlines()
    res_dfs = []
    for c_states_repr in states_reprs:
        # if c_states_repr != 'BEMP':
        #     continue
        c_states_dir_path = os.path.join(donor_diff_expr_dir_path, c_states_repr)
        log_ratio_df = pd.read_csv(os.path.join(c_states_dir_path, 'log_ratio_df.csv'))
        c_expr_df = pd.read_csv(os.path.join(c_states_dir_path, 'c_expr_df.csv'))
        cell_count_df = pd.read_csv(os.path.join(c_states_dir_path, 'cell_count_df.csv'))
        cell_count_df = cell_count_df[cell_count_df['name'].isin(log_ratio_df.columns)]
        cell_count_df[['donor_id', 'exp_name']] = cell_count_df['name'].str.split('__', expand=True)
        cell_count_df = cell_count_df.sort_values('cell_count', ascending=False).drop_duplicates(subset='donor_id', keep='first')
        assert cell_count_df['donor_id'].is_unique
        cell_count_df = generic_utils.merge_preserving_df1_index_and_row_order(cell_count_df, numbered_donor_df[['donor_id', 'exp_name', 'diagnosis_class']])
        donor_exps = list(cell_count_df['name'])
        log_ratio_df = log_ratio_df[['gene'] + donor_exps]
        c_expr_df = c_expr_df[['gene'] + donor_exps]

        res_df = c_expr_df[donor_exps].quantile(obs_expr_quantiles, axis=1).T.rename(columns={q: f'q{q}' for q in obs_expr_quantiles})
        res_df.insert(0, 'gene', c_expr_df['gene'])
        res_df = res_df[res_df['q0.99'] >= min_99_percentile_expr]
        log_ratio_df = log_ratio_df[log_ratio_df['gene'].isin(res_df['gene'])]
        

        flat_dicts = []

        mds_donor_exps = sorted(cell_count_df.loc[cell_count_df['diagnosis_class'].isin(['MDS', 'MDS/MPN']), 'name'])
        non_mds_donor_exps = sorted(cell_count_df.loc[cell_count_df['diagnosis_class'].isin(['normal', 'cytopenia']), 'name'])

        log_ratio_df.set_index('gene', inplace=True)
        assert set(donor_exps) == set(mds_donor_exps) | set(non_mds_donor_exps)
        mds_df = log_ratio_df[mds_donor_exps]
        non_mds_df = log_ratio_df[non_mds_donor_exps]
        gene_count = len(log_ratio_df)
        for i, gene in enumerate(log_ratio_df.index):
            if i % 1000 == 500:
                print(f'{gene} ({i+1}/{gene_count})')        
                # break
            mw_res = generic_utils.perform_mw_test(mds_df.loc[gene], non_mds_df.loc[gene], alternative='two-sided')
            flat_dicts.append({
                'gene': gene,
                'mw_pval': mw_res['pvalue'],
                'mw_u_minus_mean_u': mw_res['u_minus_mean_u'],
            })
        mw_df = pd.DataFrame(flat_dicts)
        res_df = generic_utils.merge_preserving_df1_index_and_row_order(res_df, mw_df)
        res_df['bh_fdr_corrected_mw_pval'] = scipy.stats.false_discovery_control(res_df['mw_pval'])
        res_df = res_df[res_df['bh_fdr_corrected_mw_pval'] <= max_bh_fdr_corrected_mw_pval]
        res_df['c_states_repr'] = c_states_repr
        res_dfs.append(res_df)

    all_res_df = pd.concat(res_dfs, ignore_index=True)
    all_res_df.sort_values('bh_fdr_corrected_mw_pval', inplace=True)
    all_res_df.to_csv(os.path.join(donor_diff_expr_dir_path, 'mw_df.csv'), index=False)

def write_donor_diff_exp_bio_reps(
        out_dir_path1,
        out_dir_path2,
        all_bio_rep_corr_df_csv_file_path,
        bio_rep_frac_for_downsample=0.9,
        downsample_reps=40,
):
    state_reprs
        
    bio_rep_corr_dfs = []

    for states_repr in states_repr:
        out_dir_path1 = f'/dummy/dummy/dummy/raid/mds/240222_pb_exps_potentially_containing_mds_or_cytopenia_illu/intermediate_output/mc_models/final_all_except_other_exp_condition/numbered_donor_cell_and_proj_mean_expr_dfs/per_donor_{states_repr}'
        out_dir_path2 = f'/dummy/dummy/dummy/raid/mds/240222_pb_exps_potentially_containing_mds_or_cytopenia_ult/intermediate_output/mc_models/final_all_except_other_exp_condition/numbered_donor_cell_and_proj_mean_expr_dfs/per_donor_{states_repr}'

        curr_log_ratio_df_csv_file_path1 = os.path.join(out_dir_path1, 'log_ratio_df.csv')
        curr_log_ratio_df_csv_file_path2 = os.path.join(out_dir_path2, 'log_ratio_df.csv')

        print(f'{curr_log_ratio_df_csv_file_path1}: {generic_utils.get_file_last_modification_time_str_repr(curr_log_ratio_df_csv_file_path1)}')
        print(f'{curr_log_ratio_df_csv_file_path2}: {generic_utils.get_file_last_modification_time_str_repr(curr_log_ratio_df_csv_file_path2)}')

        curr_log_ratio_df1 = pd.read_csv(curr_log_ratio_df_csv_file_path1, index_col=0)
        curr_log_ratio_df2 = pd.read_csv(curr_log_ratio_df_csv_file_path2, index_col=0)

        common_genes = sorted(set(curr_log_ratio_df1.index) & set(curr_log_ratio_df2.index))
        print(f'len(common_genes): {len(common_genes)}')
        filtered_df1 = curr_log_ratio_df1.loc[common_genes]
        filtered_df2 = curr_log_ratio_df2.loc[common_genes]
        combined_df = pd.concat([filtered_df1, filtered_df2], axis=1)

        curr_bio_rep_df = sc_rna_seq_preprocessing_params.get_earliest_and_latest_bio_rep_df(numbered_donor_df[numbered_donor_df['numbered_donor_id'].isin(combined_df.columns)])

        curr_bio_rep_count = len(curr_bio_rep_df)
        print(f'curr_bio_rep_count: {curr_bio_rep_count}')

        curr_bio_rep_frac_count = int(bio_rep_frac_for_downsample * curr_bio_rep_count)

        bio_rep_corr_df = pd.DataFrame({'gene': combined_df.index})
        random.seed(0)
        random_desc_and_pair_indices = [
            (f'pair_random_{bio_rep_frac_for_downsample}_{i+1}', random.sample(range(len(curr_bio_rep_df)), curr_bio_rep_frac_count))
            for i in range(downsample_reps)
        ]

        for desc, pair_indices in [
            ('all', range(len(curr_bio_rep_df))),
            *random_desc_and_pair_indices,
        ]:
            rand_bio_rep_df = curr_bio_rep_df.iloc[pair_indices]
            
            corr_col = f'corr_{desc}'
            assert corr_col not in bio_rep_corr_df

            early_df = combined_df[rand_bio_rep_df['numbered_donor_id_early']]
            late_df = combined_df[rand_bio_rep_df['numbered_donor_id_late']]
            corrs = mc.ut.pairs_corrcoef_rows(
                mc.ut.to_layout(mc.ut.to_numpy_matrix(early_df), layout='row_major'), 
                mc.ut.to_layout(mc.ut.to_numpy_matrix(late_df), layout='row_major'), 
                reproducible=True,
            )
            bio_rep_corr_df[corr_col] = corrs

        bio_rep_corr_df['max_rand_frac_corr'] = bio_rep_corr_df[[f'corr_{x[0]}' for x in random_desc_and_pair_indices]].max(axis=1)    
        bio_rep_corr_df = bio_rep_corr_df[[x for x in bio_rep_corr_df.columns if 'pair_random' not in x]]
        bio_rep_corr_df['state'] = c_state
        if 0:
            plt.close('all')
            # sb.histplot(bio_rep_corr_df['corr_all'])
            sb.histplot(bio_rep_corr_df['max_rand_frac_corr'])
        if 0:
            plt.close('all')
            sb.scatterplot(bio_rep_corr_df, x='corr_all', y='max_rand_frac_corr')
            print(((bio_rep_corr_df['max_rand_frac_corr'] - bio_rep_corr_df['corr_all']) > 0).value_counts())
        bio_rep_corr_dfs.append(bio_rep_corr_df)

        # raise RuntimeError('yay')

    all_bio_rep_corr_df = pd.concat(bio_rep_corr_dfs, ignore_index=True)
    all_bio_rep_corr_df.sort_values('max_rand_frac_corr', ascending=False, inplace=True)

    if 0:
        all_bio_rep_corr_df.to_csv(all_bio_rep_corr_df_csv_file_path, index=False)


def get_all_misc_numbered_donor_info():
    dataset_and_mc_model_names = [
        ('illu_mds', 'final_normal_pb_atlas'),
        ('illu_mds', 'final_mds_cyto_normal_excluding_atlas'),
        ('ult_mds', 'final_mds_cyto_normal_excluding_atlas'),
    ]
    mc_model_paths_dicts = [mds_analysis_params.get_mc_model_paths(x[1], x[0]) for x in dataset_and_mc_model_names]
    all_misc_numbered_donor_info = sc_rna_seq_preprocessing_params.combine_misc_numbered_donor_infos([
        generic_utils.read_object_from_pickle_file(x['misc_numbered_donor_info_pickle']) for x in mc_model_paths_dicts
    ])
    return all_misc_numbered_donor_info

def get_all_numbered_donor_df():
    return get_all_misc_numbered_donor_info()['numbered_donor_df']

def identify_inconsistent_arch_and_get_non_arch4_donor_bleeding_date_to_exclude():
    all_numbered_donor_df = get_all_numbered_donor_df()
    assert (all_numbered_donor_df['exp_date'] == all_numbered_donor_df['bleeding_date']).all()
    
    mutation_df = arch_mutation_interface_and_utils.get_minimal_all_mutation_df()
    presumably_non_amplicon_mut_df = mutation_df[~mutation_df['source_file_name'].isin(['240621_CALR_without_MPN.xlsx', '240621_SRSF2.xlsx'])]
    # mostly_arch3_mut_df = presumably_non_amplicon_mut_df.loc[presumably_non_amplicon_mut_df['mostly_ARCH3'], ['donor_id', 'exp_date', 'gene']].drop_duplicates()
    # arch4_mut_df = presumably_non_amplicon_mut_df.loc[presumably_non_amplicon_mut_df['ARCH4'], ['donor_id', 'exp_date', 'gene']].drop_duplicates()
    all_arch4_df = pd.read_csv(get_sc_rna_seq_preprocessing_params()['donor_table_paths']['arch4_donor_id_and_exp_date_df_csv'])
    any_mut_donor_ids = sorted(set(presumably_non_amplicon_mut_df['donor_id']))
    multi_sample_donor_df = all_numbered_donor_df[['donor_id', 'bleeding_date']].drop_duplicates()['donor_id'].value_counts().reset_index(name='sample_count')
    multi_sample_donor_df = multi_sample_donor_df[multi_sample_donor_df['sample_count'] > 1]
    multi_sample_donor_ids = sorted(set(multi_sample_donor_df['donor_id']))
    any_mut_multi_sample_donor_ids = sorted(set(multi_sample_donor_ids) & set(any_mut_donor_ids))

    inconsistent_donor_id_to_exp_date_to_mutated_genes = {}
    non_arch4_inconsistent_donor_id_and_bleeding_date_flat_dicts = []
    print('\ninconsistent ARCH:')
    for donor_id in any_mut_multi_sample_donor_ids:
        # if donor_id != 'N48':
        #     continue
        curr_mut_df = presumably_non_amplicon_mut_df[presumably_non_amplicon_mut_df['donor_id'] == donor_id]
        exp_date_to_mutated_genes = {}
        genes = None
        inconsistent_muts = False
        all_dates = set(all_numbered_donor_df.loc[all_numbered_donor_df['donor_id'] == donor_id, 'bleeding_date']) # ugh. no bleeding date for mutations...
        mut_dates = set(curr_mut_df['exp_date'])
        if not mut_dates <= all_dates:
            # print(f'{donor_id}: mut_dates={mut_dates}, all_dates={all_dates}')
            all_dates |= mut_dates
        for exp_date in sorted(all_dates):
            curr_date_df = curr_mut_df[curr_mut_df['exp_date'] == exp_date]
            exp_date_to_mutated_genes[exp_date] = sorted(set(curr_date_df['gene']))
            if genes is None:
                genes = exp_date_to_mutated_genes[exp_date]
            else:
                if genes != exp_date_to_mutated_genes[exp_date]:
                    inconsistent_muts = True
        if inconsistent_muts:
            arch4_dates = set(all_arch4_df.loc[all_arch4_df['donor_id'] == donor_id, 'exp_date'])
            non_arch4_dates = set(exp_date_to_mutated_genes) - arch4_dates
            non_arch4_inconsistent_donor_id_and_bleeding_date_flat_dicts.extend(
                {'donor_id': donor_id, 'bleeding_date': x} for x in non_arch4_dates)
            exp_date_to_mutated_genes_for_print = {f'(4){k}' if k in arch4_dates else k: v for k,v in exp_date_to_mutated_genes.items()}
            print(f'{donor_id}: {exp_date_to_mutated_genes_for_print}')
            inconsistent_donor_id_to_exp_date_to_mutated_genes[donor_id] = exp_date_to_mutated_genes
    return pd.DataFrame(non_arch4_inconsistent_donor_id_and_bleeding_date_flat_dicts)

def add_non_ARCH4_sample_in_donor_with_any_inconsistent_ARCH_col(numbered_donor_df):
    non_arch4_donor_bleeding_date_to_exclude_df = identify_inconsistent_arch_and_get_non_arch4_donor_bleeding_date_to_exclude()
    numbered_donor_df['non_ARCH4_sample_in_donor_with_any_inconsistent_ARCH'] = generic_utils.get_mask_of_df1_rows_that_inner_join_will_keep(
        numbered_donor_df, non_arch4_donor_bleeding_date_to_exclude_df)
    discarded_due_to_inconsistent_arch_donor_ids = set(numbered_donor_df['donor_id']) - set(numbered_donor_df.loc[~numbered_donor_df['non_ARCH4_sample_in_donor_with_any_inconsistent_ARCH'], 'donor_id'])
    assert not discarded_due_to_inconsistent_arch_donor_ids, f'discarded_due_to_inconsistent_arch_donor_ids: {discarded_due_to_inconsistent_arch_donor_ids}'

def assert_multi_sample_donors_have_explicit_exp_names_in_cna_infos():
    all_numbered_donor_df = get_all_numbered_donor_df()    
    multi_sample_donor_df = all_numbered_donor_df[['donor_id', 'bleeding_date']].drop_duplicates()['donor_id'].value_counts().reset_index(name='sample_count')
    multi_sample_donor_df = multi_sample_donor_df[multi_sample_donor_df['sample_count'] > 1]
    multi_sample_donor_ids = sorted(set(multi_sample_donor_df['donor_id']))

    donor_agg_cna_df_csv_file_path = get_mds_params()['donor_agg_cna_df_csv_file_path']
    donor_agg_cna_df = pd.read_csv(donor_agg_cna_df_csv_file_path)
    multi_sample_donor_agg_cna_df = donor_agg_cna_df[donor_agg_cna_df['donor_id'].isin(multi_sample_donor_ids)]
    print(f'len(multi_sample_donor_agg_cna_df): {len(multi_sample_donor_agg_cna_df)}')
    # print(f'multi_sample_donor_agg_cna_df: {multi_sample_donor_agg_cna_df}')
    all_mask = multi_sample_donor_agg_cna_df['exp_name'] == 'all'
    if all_mask.any():
        multi_sample_donor_agg_cna_df_with_all_exp_name = multi_sample_donor_agg_cna_df[all_mask]
        print(f'multi_sample_donor_agg_cna_df_with_all_exp_name:\n{multi_sample_donor_agg_cna_df_with_all_exp_name}')
        raise RuntimeError('shouldnt have exp_name=all')

def get_fig4_and_xgboost_basic_allowed_donor_mask(ext_donor_feature_df):
    return ~(
        ext_donor_feature_df['complex_karyo'] | 
        ext_donor_feature_df['non_ARCH4_sample_in_donor_with_any_inconsistent_ARCH'] |
        ((ext_donor_feature_df['final_c_hspc_count'] < get_mds_params()['min_hspc_count_for_composition_analysis'])) |
        ext_donor_feature_df['biased_composition_due_to_discarded_low_umi_count_barcodes']
    )

def get_df_with_usage_in_analysis_masks(all_ext_donor_feature_df):
    usage_in_analysis_masks_df_csv_file_path = get_mds_params()['usage_in_analysis_masks_df_csv_file_path']
    if os.path.isfile(usage_in_analysis_masks_df_csv_file_path):
        return generic_utils.merge_preserving_df1_index_and_row_order(
            all_ext_donor_feature_df, pd.read_csv(usage_in_analysis_masks_df_csv_file_path))
    
    all_ext_donor_feature_df = all_ext_donor_feature_df.copy()

    mask_names = [
        'in_orig_ref_model',
        'in_fig4_ref_model',
        'in_mds_cyto_groups_analysis',
        'used_to_fit_hspc_compos_knn',
        'in_train_set',
        'in_test_set',
        'in_blast_vs_clp_e_scatter',
        'in_within_state_sig_bio_rep_analysis',
    ]
    # print(f'mask_names:\n{mask_names}')

    
    n1_donor_df_csv_file_path = get_sc_rna_seq_preprocessing_params()['donor_table_paths']['n1_donor_id_df_csv']
    if os.path.isfile(n1_donor_df_csv_file_path):
        n1_donor_df = pd.read_csv(n1_donor_df_csv_file_path)
    else:
        n1_c_ad = ad.read(get_mds_params()['nimrod_148_atlas']['c_ad_file_path'])
        n1_donor_df = n1_c_ad.obs[['exp_name', 'indiv_id']].rename(columns={'indiv_id': 'nimrod_compatible_donor_id'}).drop_duplicates().reset_index(drop=True)
        n1_donor_df.to_csv(n1_donor_df_csv_file_path, index=False)
    n1_donor_df['donor_id'] = n1_donor_df['nimrod_compatible_donor_id'].astype(str).replace(get_sc_rna_seq_preprocessing_params()['other_donor_id_to_donor_id_i_use'])
    sc_rna_seq_preprocessing_params.add_exp_and_bleeding_date_according_to_exp_name_and_donor_id(n1_donor_df)
    n1_donor_df.loc[n1_donor_df['exp_name'] == 'demux_27_02_23_1', 'bleeding_date'] = '27.02.23' 
    n1_donor_df.loc[n1_donor_df['exp_name'] == 'demux_28_02_23_1', 'bleeding_date'] = '28.02.23' 

    assert {'nimrod_compatible_donor_id', 'bleeding_date'} <= set(all_ext_donor_feature_df.columns)
    all_ext_donor_feature_df['in_orig_ref_model'] = generic_utils.get_mask_of_df1_rows_that_inner_join_will_keep(
        all_ext_donor_feature_df, 
        n1_donor_df[['nimrod_compatible_donor_id', 'bleeding_date']].drop_duplicates(),
    )

    orig_ref_healthy_donor_ids_i_skipped_due_to_suspected_contamination = {'N350','N351','N352','N353'}
    orig_ref_donor_ids_i_skipped_due_to_other_disease = {'N135', 'N297'}
    orig_ref_donor_ids_i_skipped_due_to_abnormally_high_counts = {'N132', 'N177', 'N75', 'N83', 'N96'}
    orig_ref_donor_ids_i_skipped = (
        orig_ref_healthy_donor_ids_i_skipped_due_to_suspected_contamination
        | orig_ref_donor_ids_i_skipped_due_to_other_disease
        | orig_ref_donor_ids_i_skipped_due_to_abnormally_high_counts
    )
    assert all_ext_donor_feature_df.loc[all_ext_donor_feature_df['in_orig_ref_model'], 'nimrod_compatible_donor_id'].nunique() == (
        148 - len(orig_ref_donor_ids_i_skipped))
    # to make sure this is as expected, despite swapped labels
    assert set(all_ext_donor_feature_df.loc[
        all_ext_donor_feature_df['in_orig_ref_model'] & 
        all_ext_donor_feature_df['bleeding_date'].isin(['27.02.23', '28.02.23']), 
        'donor_id'
    ]) == {'N354', 'N355', 'N356', 'N357', 'N260', 'N307', 'N362', 'N363'} 
    
    all_ext_donor_feature_df['in_fig4_ref_model'] = all_ext_donor_feature_df['dataset_name'] == 'ref'
    
    curr_feature_df = all_ext_donor_feature_df.copy()
    curr_feature_df = curr_feature_df[get_fig4_and_xgboost_basic_allowed_donor_mask(curr_feature_df)]

    curr_feature_df = curr_feature_df[curr_feature_df['numbered_donor_id'].isin(
        sc_rna_seq_preprocessing_params.get_best_or_random_numbered_donor_id_per_donor_or_donor_sample(
            curr_feature_df, single_nd_id_per_donor=True,
            sort_by_cols=[
                # NOTE: better also to always prefer ARCH4, as unless we reanilize ARCH3, it is possible that it was overwritten by ARCH4 result (and if the ARCH4 result was no-mutations, then i would not know that ARCH3 was overwritten).
                'known_treatment', 
                'cell_count', # NOTE: ugh. that's a (not too bad) bug here. cell_count is the count in the mc_model. so for donors with very high count, this is lower due to downsampling to 3*median in MC model, though that's what i actually use in the end, so maybe even better this way. also, it doesn't make sense to order by cell_count - should order by HSPC count. but nevermind, the effect on the results is negligible.
                # 'final_c_hspc_count', # this would have been better than using cell_count, see above.
            ], 
            sort_by_ascending=[True, False],
        ))]
    
    all_ext_donor_feature_df['in_mds_cyto_groups_analysis'] = all_ext_donor_feature_df['numbered_donor_id'].isin(curr_feature_df['numbered_donor_id'])
    fig4b_donor_count = all_ext_donor_feature_df['in_mds_cyto_groups_analysis'].sum()
    assert fig4b_donor_count == 191, f"fig4b_donor_count: {fig4b_donor_count}" # 240711: a sanity check that nothing (major) changed recently. things could change without the 191 changing...

    all_ext_donor_feature_df['used_to_fit_hspc_compos_knn'] = (
        all_ext_donor_feature_df['in_mds_cyto_groups_analysis'] & (all_ext_donor_feature_df['dataset_name'] == 'ref'))
    assert all_ext_donor_feature_df['used_to_fit_hspc_compos_knn'].sum() == 70 # 240720: a sanity check that nothing (major) changed recently.

    all_ext_donor_feature_df['in_train_set'] = (
        all_ext_donor_feature_df['in_mds_cyto_groups_analysis'] & (all_ext_donor_feature_df['dataset_name'] != 'ref') 
        & all_ext_donor_feature_df['is_ultima']
    )
    all_ext_donor_feature_df['in_test_set'] = (
        all_ext_donor_feature_df['in_mds_cyto_groups_analysis'] & (all_ext_donor_feature_df['dataset_name'] != 'ref') 
        & ~all_ext_donor_feature_df['is_ultima']
        & (all_ext_donor_feature_df['diagnosis_class'] != 'normal')
    )
    
    assert not (
        set(all_ext_donor_feature_df.loc[all_ext_donor_feature_df['in_train_set'], 'donor_id']) & 
        set(all_ext_donor_feature_df.loc[all_ext_donor_feature_df['in_test_set'], 'donor_id'])
    )
    assert all_ext_donor_feature_df.loc[all_ext_donor_feature_df['in_train_set'], 'donor_id'].is_unique
    assert all_ext_donor_feature_df.loc[all_ext_donor_feature_df['in_test_set'], 'donor_id'].is_unique

    # all_ext_donor_feature_df.to_csv('temp/all_ext_donor_feature_df.csv', index=False)
    print('\ntrain:\n' + str(all_ext_donor_feature_df.loc[all_ext_donor_feature_df['in_train_set'], 'diagnosis_class'].value_counts()))
    print('\ntest:\n' + str(all_ext_donor_feature_df.loc[all_ext_donor_feature_df['in_test_set'], 'diagnosis_class'].value_counts()))

    multiple_diagnosis_class_donor_ids = list((all_ext_donor_feature_df.drop_duplicates(subset='numbered_donor_id_without_tech_rep_suffix')[[
        'donor_id', 'diagnosis_class']].drop_duplicates()['donor_id'].value_counts() > 1).loc[lambda x: x].index)
    multiple_diagnosis_class_donor_df = all_ext_donor_feature_df.loc[all_ext_donor_feature_df['donor_id'].isin(multiple_diagnosis_class_donor_ids), ['donor_id', 'exp_name', 'diagnosis', 'diagnosis_class']].drop_duplicates().sort_values('donor_id')
    if not multiple_diagnosis_class_donor_df.empty:
        print(f'WARNING: multiple_diagnosis_class_donor_df:\n{multiple_diagnosis_class_donor_df}')

    curr_df = all_ext_donor_feature_df.copy()
    curr_df['bm_facs_blast_10x_sample_approx_months_dist'] = curr_df['bm_facs_blast_10x_sample_days_dist'] // 30
    curr_df = curr_df[
        # as the median of c_clp_e_freq across normal is 0.034, a total of 500 cells seems reasonable, as we expect to sample 17 c_clp_e cells.
        (curr_df['final_c_hspc_count'] >= get_mds_params()['min_hspc_count_for_composition_analysis'])
        & ~curr_df['biased_composition_due_to_discarded_low_umi_count_barcodes']
        & (curr_df['bm_facs_blast_10x_sample_days_dist'] <= get_mds_params()['max_bm_facs_blast_10x_sample_days_dist'])
    ].copy()

    all_ext_donor_feature_df['in_blast_vs_clp_e_scatter'] = all_ext_donor_feature_df['numbered_donor_id'].isin(
        sc_rna_seq_preprocessing_params.get_best_or_random_numbered_donor_id_per_donor_or_donor_sample(
            curr_df, single_nd_id_per_donor=True,
            sort_by_cols=['bm_facs_blast_10x_sample_approx_months_dist', 'final_c_hspc_count'],
            sort_by_ascending=[True, False],
            # sort_by_cols=['final_c_hspc_count'], sort_by_ascending=[False], # pretty much the same, but it seems important to prioritize time distance between scRNA and BM FACS.
        )
    )

    curr_df = all_ext_donor_feature_df.copy()
    curr_df = curr_df[curr_df['numbered_donor_id'].isin(
        sc_rna_seq_preprocessing_params.get_best_or_random_numbered_donor_id_per_donor_or_donor_sample(
            curr_df[
                ~curr_df['complex_karyo']
                # & ~curr_df['known_treatment'] # losing too many samples...
            ]
        )
    )]
    bio_rep_df = sc_rna_seq_preprocessing_params.get_earliest_and_latest_bio_rep_df(curr_df)
    chosen_bio_rep_samples_df = pd.concat([
        pd.DataFrame({'donor_id': bio_rep_df['donor_id'], 'exp_name': bio_rep_df[f'exp_name_{x}']}) 
        for x in ('early', 'late')
    ], ignore_index=True)
    all_ext_donor_feature_df['in_within_state_sig_bio_rep_analysis'] = (
        generic_utils.get_mask_of_df1_rows_that_inner_join_will_keep(all_ext_donor_feature_df, chosen_bio_rep_samples_df))
    
    for mask in mask_names:
        if mask not in {'in_within_state_sig_bio_rep_analysis', 'in_orig_ref_model', 'in_fig4_ref_model'}:
            # print(mask)
            assert all_ext_donor_feature_df.loc[all_ext_donor_feature_df[mask], 'donor_id'].is_unique, f'all_ext_donor_feature_df.loc[{mask}, "donor_id"] is not unique'
    
    all_ext_donor_feature_df[['nimrod_compatible_donor_id', 'scrna_exp_id', *mask_names]].to_csv(usage_in_analysis_masks_df_csv_file_path, index=False)
    return all_ext_donor_feature_df

def add_unreliable_arch_cols_and_remove_vaf_vals(all_ext_donor_feature_df):
    all_ext_donor_feature_df['bleeding_date_as_date'] = pd.to_datetime(
        all_ext_donor_feature_df['bleeding_date'], format='%d.%m.%y', errors='raise')

    all_ext_donor_feature_df['ARCH4'] = (
        all_ext_donor_feature_df['bleeding_date_as_date'] >= clinical_data_interface_and_utils.FIRST_ARCH4_DATE_AS_DATE)
    all_ext_donor_feature_df['CH_detected_in_ARCH4'] = (
        all_ext_donor_feature_df['max_mean_VAF'] > 0) & all_ext_donor_feature_df['ARCH4']
    has_arch4_ch_df = all_ext_donor_feature_df.groupby('donor_id')['CH_detected_in_ARCH4'].any().reset_index(name='donor_has_ARCH4_sample_with_CH')
    all_ext_donor_feature_df['donor_has_ARCH4_sample_with_CH'] = generic_utils.merge_preserving_df1_index_and_row_order(all_ext_donor_feature_df, has_arch4_ch_df)['donor_has_ARCH4_sample_with_CH']
    all_ext_donor_feature_df['ARCH3_and_donor_has_ARCH4_sample_with_CH'] = (
        ~all_ext_donor_feature_df['ARCH4'] & all_ext_donor_feature_df['donor_has_ARCH4_sample_with_CH'])
    
    all_ext_donor_feature_df.loc[all_ext_donor_feature_df['ARCH3_and_donor_has_ARCH4_sample_with_CH'], 'max_mean_VAF'] = np.nan
    mutation_df = arch_mutation_interface_and_utils.get_minimal_all_mutation_df()
    mut_genes = sorted(mutation_df['gene'].drop_duplicates())
    all_ext_donor_feature_df.loc[
        all_ext_donor_feature_df['ARCH3_and_donor_has_ARCH4_sample_with_CH'], [f'{x}_max_mean_VAF' for x in mut_genes]] = np.nan

def get_unreliable_arch_blood_sample_ids(all_ext_donor_feature_df):
    return sorted(all_ext_donor_feature_df.loc[all_ext_donor_feature_df[
        'ARCH3_and_donor_has_ARCH4_sample_with_CH'], 'blood_sample_id'].drop_duplicates())

def add_mean_near_neighbor_dist(all_ext_donor_feature_df):
    curr_feature_df = all_ext_donor_feature_df.copy()
    
    min_hspc_count = get_mds_params()['min_hspc_count_for_composition_analysis']
    epsilon_for_cell_freq = get_mds_params()['epsilon_for_donor_log_c_state_freq']
    nearest_neighbor_count = get_mds_params()['nearest_neighbor_count_for_composition_analysis']

    ordered_pb_hspc_state_names = [f'c_{x}' for x in mds_analysis_params.ORDERED_PB_HSPC_STATE_NAMES]
    observed_hspc_state_names = sorted(set(curr_feature_df.columns) & set(ordered_pb_hspc_state_names))
    log_hspc_freq_cols = [f'log_{x}' for x in observed_hspc_state_names]
    curr_feature_df[log_hspc_freq_cols] = np.log2(curr_feature_df[observed_hspc_state_names] + epsilon_for_cell_freq)

    relevant_mask_to_calc_mean_nn_dist = (
        (curr_feature_df['final_c_hspc_count'] >= min_hspc_count) & 
        ~curr_feature_df['biased_composition_due_to_discarded_low_umi_count_barcodes']
    )
    assert not (curr_feature_df['in_mds_cyto_groups_analysis'] & ~relevant_mask_to_calc_mean_nn_dist).any()
    curr_feature_df = curr_feature_df[relevant_mask_to_calc_mean_nn_dist]

    unfiltered_curr_feature_df = curr_feature_df.copy()
    curr_feature_df = curr_feature_df[curr_feature_df['in_mds_cyto_groups_analysis']]

    knn = sklearn.neighbors.NearestNeighbors(n_neighbors=nearest_neighbor_count)
    feature_df_for_knn = curr_feature_df.loc[curr_feature_df['used_to_fit_hspc_compos_knn'], log_hspc_freq_cols]
    knn.fit(feature_df_for_knn)

    dists_from_nearest_neighbors = knn.kneighbors(unfiltered_curr_feature_df[log_hspc_freq_cols])[0]
    unfiltered_curr_feature_df['mean_near_neighbor_dist'] = dists_from_nearest_neighbors.mean(axis=1)
    
    # set to nan because there is a bias for the atlas samples - as they have one neighbor with a distance of 0 by definition.
    unfiltered_curr_feature_df.loc[unfiltered_curr_feature_df['used_to_fit_hspc_compos_knn'], 'mean_near_neighbor_dist'] = np.nan

    all_ext_donor_feature_df['mean_near_neighbor_dist'] = generic_utils.merge_preserving_df1_index_and_row_order(
        all_ext_donor_feature_df, unfiltered_curr_feature_df[['mean_near_neighbor_dist', 'numbered_donor_id']],
        how='left',
    )['mean_near_neighbor_dist'].to_numpy()



def write_and_get_cell_filtering_count_df(dataset_name, use_existing=False):
    mds_path_dict = mds_analysis_params.get_mds_path_dict(dataset_name=dataset_name)
    cell_filtering_count_df_csv_file_path = mds_path_dict['cell_filtering_count_df_csv_file_path']
    if use_existing and os.path.isfile(cell_filtering_count_df_csv_file_path):
        count_df = pd.read_csv(cell_filtering_count_df_csv_file_path)
    else:
        c_ad_file_path = mds_path_dict['assigned_and_no_doublet_enirhced_mc_c_ad_file_path']
        if not os.path.isfile(c_ad_file_path):
            raise RuntimeError(f'{c_ad_file_path} was not found.')
        print(f'the last modification time of the input files are:')
        print(f'{c_ad_file_path}: {generic_utils.get_file_last_modification_time_str_repr(c_ad_file_path)}')
        print('please make sure all are as expected.')

        preprocessing_out_exp_df_csv_file_path = mds_path_dict['preprocessing_out_exp_df_csv_file_path']
        curr_exp_df = pd.read_csv(preprocessing_out_exp_df_csv_file_path)
        c_ad = ad.read_h5ad(c_ad_file_path)
        single_donor_exp_names = sc_rna_seq_preprocessing_params.get_single_donor_exp_names()

        flat_dicts = []
        for _, row in curr_exp_df[~curr_exp_df['donor_id_to_num_of_cells'].isna()].iterrows():
            exp_name = row['exp_name']
            if exp_name in single_donor_exp_names:
                # print(f'NOTE: {exp_name} is in single_donor_exp_names')
                assert ',' not in row['ids_in_exp']
                flat_dicts.append({
                    'exp_name': exp_name,
                    'donor_id': row['ids_in_exp'],
                    'orig_count': row['num_of_cells'],
                })
            else:
                donor_id_to_num_of_cells = row['donor_id_to_num_of_cells']

                for donor_id_to_num_of_cells_repr in donor_id_to_num_of_cells.split(';'):
                    donor_id, _, num_of_cells = donor_id_to_num_of_cells_repr.partition(':')
                    num_of_cells = int(num_of_cells)
                    flat_dicts.append({
                        'exp_name': exp_name,
                        'donor_id': donor_id,
                        'orig_count': num_of_cells,
                    })

        count_df = pd.DataFrame(flat_dicts)
        count_df = count_df[count_df['exp_name'].isin(c_ad.obs['exp_name'].unique())]

        for df_path in [
            os.path.join(mds_path_dict['preprocessing_out_dir_file_path'], 'mouse_human_doublet_discarded_df.csv'),
            os.path.join(mds_path_dict['preprocessing_out_dir_file_path'], 'whitelist_discarded_df.csv'),
            os.path.join(mds_path_dict['preprocessing_out_dir_file_path'], 'under_manual_exp_min_umi_count_df.csv'),
            mds_path_dict['too_few_non_excluded_umis_count_df_csv_file_path'],
            mds_path_dict['too_high_mito_and_malat1_count_df_csv_file_path'],
            mds_path_dict['too_low_ribo_prot_count_df_csv_file_path'],
            mds_path_dict['hba_hbb_contaminated_count_df_csv_file_path'],
            mds_path_dict['neut_or_neut_contaminated_count_df_csv_file_path'],
            mds_path_dict['platelet_or_platelet_contaminated_count_df_csv_file_path'],
            mds_path_dict['too_low_norm_nuclear_expr_count_df_csv_file_path'],
            mds_path_dict['mebemp_monocyte_doublet_count_df_csv_file_path'],
            mds_path_dict['enriched_genotype_doublet_mc_count_df_csv_file_path'],
            mds_path_dict['high_log_norm_umi_count_count_df_csv_file_path'],
            mds_path_dict['dc_nkt_monocyte_endothel_b_count_df_csv_file_path'],
            mds_path_dict['presumably_hspc_count_df_csv_file_path'],
        ]:
            curr_count_df = pd.read_csv(df_path)
            new_cols = set(curr_count_df.columns) - {'exp_name', 'donor_id'}
            assert len(new_cols) == 1
            new_col = next(iter(new_cols))
            if curr_count_df['donor_id'].isna().all():
                count_df[new_col] = 0
            else:
                count_df = generic_utils.merge_preserving_df1_index_and_row_order(count_df, curr_count_df, how='left', on=['exp_name', 'donor_id'])
                count_df[new_col].fillna(0, inplace=True)
            

        count_cols = [x for x in count_df.columns if x not in {'orig_count', 'exp_name', 'donor_id'}]
        count_df = count_df[count_df[[x for x in count_cols if x not in {'not_whitelist_barcode_count', 'under_manual_exp_min_umi_count_count'}]].sum(axis=1) > 0]
        counts_dont_match_mask = count_df['orig_count'] != count_df[count_cols].sum(axis=1)
        assert not counts_dont_match_mask.any(), f'count_df[counts_dont_match_mask]:\n{count_df[counts_dont_match_mask]}'
        if 0:
            count_df[counts_dont_match_mask]['orig_count'] - count_df[counts_dont_match_mask][count_cols].sum(axis=1)

        c_ad.obs[['exp_name', 'donor_id']] = c_ad.obs[['exp_name', 'donor_id']].astype(str)
        count_df = generic_utils.merge_preserving_df1_index_and_row_order(
            count_df,
            c_ad.obs[c_ad.obs['not_dc_nkt_monocyte_endothel_b_lamp3']].groupby(['exp_name', 'donor_id'])['log_num_of_non_excluded_umis'].median().reset_index(name='presumably_hspc_median_log_num_of_non_excluded_umis')
        )
        count_df['manual_biased_composition_type'] = count_df['exp_name'].map(sc_rna_seq_preprocessing_params.get_exp_name_to_manual_biased_composition_type()).to_numpy()
        assert not count_df['manual_biased_composition_type'].isna().any(), f"exps missing from type_to_exp_name_to_forced_num_of_filtered_barcodes_and_min_umi_count: {sorted(count_df.loc[count_df['manual_biased_composition_type'].isna(), 'exp_name'].unique())}"
        count_df.sort_values('presumably_hspc_median_log_num_of_non_excluded_umis', inplace=True)
        count_df.to_csv(cell_filtering_count_df_csv_file_path, index=False)
    
    count_cols = list(x for x in count_df.columns if x not in {
        'orig_count', 'exp_name', 'donor_id', 'presumably_hspc_median_log_num_of_non_excluded_umis', 'manual_biased_composition_type'})
    assert (count_df['orig_count'] == count_df[count_cols].sum(axis=1)).all(), 'already did this sanity check, so this one should never fail'
    return count_df, count_cols

def generate_scrna_exp_id_df(all_ext_donor_feature_df):
    # print(all_ext_donor_feature_df[['blood_sample_id', 'exp_name']].sort_values(['blood_sample_id', 'exp_name']))
    scrna_exp_id_df = all_ext_donor_feature_df[['blood_sample_id', 'exp_name']].sort_values(['blood_sample_id', 'exp_name'])['exp_name'].drop_duplicates().reset_index(drop=True).to_frame()
    scrna_exp_id_df['scrna_exp_id'] = scrna_exp_id_df.index + 1
    return scrna_exp_id_df

def generate_blood_sample_id_df(all_ext_donor_feature_df):
    all_ext_donor_feature_df['bleeding_date_as_date'] = pd.to_datetime(
        all_ext_donor_feature_df['bleeding_date'], format='%d.%m.%y', errors='raise')
    all_ext_donor_feature_df.sort_values(['bleeding_date_as_date', 'donor_id'], inplace=True)
    blood_sample_id_df = all_ext_donor_feature_df[['donor_id', 'bleeding_date']].drop_duplicates().reset_index(drop=True)
    blood_sample_id_df['blood_sample_id'] = blood_sample_id_df.index + 1
    return blood_sample_id_df
