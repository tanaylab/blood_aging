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
from generic import hg38_utils
from mds import arch_mutation_interface_and_utils
from mds import clinical_data_interface_and_utils
from sc_rna_seq_preprocessing import sc_rna_seq_preprocessing_params

plt.rcParams["patch.force_edgecolor"] = False
plt.rcParams['patch.linewidth'] = 0
plt.rcParams['patch.edgecolor'] = 'none'
# plt.rcParams['scatter.edgecolors'] = 'black' # didnt affect sb.scatterplot

def get_mds_params():
    return mds_analysis_params.MDS_ANALYSIS_PARAMS


def get_sc_rna_seq_preprocessing_params():
    return sc_rna_seq_preprocessing_params.SC_RNA_SEQ_PREPROCESSING_PARAMS

def get_karyotype_estimation_params():
    return get_mds_params()['karyotype_estimation']


MC_STRICT_THRESHOLD_NAMES = [
    'min_median_cluster_median_projected_fold_threshold',
    'max_median_cluster_median_projected_fold_threshold',
]
MC_PARTIAL_THRESHOLD_NAMES = [
    'min_median_cluster_median_projected_fold_partial_threshold',
    'max_median_cluster_median_projected_fold_partial_threshold',
]
ALL_MC_THRESHOLD_NAMES = MC_STRICT_THRESHOLD_NAMES + MC_PARTIAL_THRESHOLD_NAMES
CELL_THRESHOLD_NAMES = [
    'max_mean_cell_cluster_mean_projected_fold_threshold',
    'min_mean_cell_cluster_mean_projected_fold_threshold',
]

ALL_CELL_AND_MC_THRESHOLD_NAMES = CELL_THRESHOLD_NAMES + ALL_MC_THRESHOLD_NAMES



def get_info_for_karyo_estimation(unified_mc_model_name, unified_mc_ad_file_path, unified_c_ad_file_path, dataset_name):
    mds_path_dict = mds_analysis_params.get_mds_path_dict(dataset_name=dataset_name)

    chrom_pos_df = hg38_utils.get_chrom_pos_df(
        get_mds_params()['hg_fasta_file_path'], mds_path_dict['chrom_pos_df_csv'], 
        ordered_chrom_names=get_karyotype_estimation_params()['ordered_chrom_names_for_karyotype_estimation'],
        use_existing_chrom_pos_df_csv=True,
    )
        
    if unified_c_ad_file_path is not None:
        unified_c_ad = ad.read_h5ad(unified_c_ad_file_path)
    else:
        unified_c_ad = None
    unified_mc_ad = ad.read_h5ad(unified_mc_ad_file_path)

    # if ('gene_id' not in list(unified_mc_ad.var)) and (unified_c_ad_file_path is not None):
    if 'gene_id' not in list(unified_mc_ad.var):
        # this is needed because gene_id is unique while gene_name might not be.
        mc_utils.add_gene_id_to_mc_ad(unified_mc_ad, c_ad=unified_c_ad)
    unified_mc_ad_var_names = unified_mc_ad.var_names.copy()
    gene_loc_df = hg38_utils.get_relevant_gene_loc_df(
        unified_mc_ad, 
        hg_genes_gtf_file_path=get_mds_params()['hg_genes_gtf_file_path'], 
        gene_loc_df_csv_file_path=get_mds_params()['gene_loc_df_csv_file_path'],
        remove_chr_prefix_from_chr_name=True,
    )
    assert unified_mc_ad_var_names.isin(gene_loc_df['gene_name']).all()

    if unified_c_ad_file_path is not None:
        donor_id_to_exp_sharing_donor_ids = mc_utils.get_donor_id_to_exp_sharing_donor_ids(unified_c_ad)
        
        unified_mc_model_minimal_donor_exp_info = mds_analysis_params.get_mc_model_minimal_donor_exp_info(unified_mc_model_name, dataset_name)
        donor_ids = unified_mc_model_minimal_donor_exp_info['donor_ids']
        total_num_of_donors = len(donor_ids)
    else:
        donor_id_to_exp_sharing_donor_ids = None


    if 'expr' not in unified_mc_ad.layers:
        mc_utils.write_expr_and_expr_enrich(unified_mc_ad)


    atlas_to_use = get_karyotype_estimation_params()['atlas_to_use']
    print(f'atlas_to_use: {atlas_to_use}')
        
    atlas_params = get_mds_params()[atlas_to_use]
    atlas_mc_model_paths = mds_analysis_params.get_mc_model_paths(
        os.path.basename(atlas_params['mc_model_dir_path']), atlas_params['dataset_name'])
    atlas_misc_numbered_donor_info = generic_utils.read_object_from_pickle_file(atlas_mc_model_paths['misc_numbered_donor_info_pickle'])
    atlas_numbered_donor_df = atlas_misc_numbered_donor_info['numbered_donor_df']
    atlas_donor_ids = set(atlas_numbered_donor_df['donor_id'].unique())

    karyotype_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [(0, 'blue'), (0.45, 'white'), (0.55, 'white'), (1, 'red')])
    epsilon_for_log_projected_fold = get_karyotype_estimation_params()['epsilon_for_log_projected_fold']

    if 'proj_expr' not in unified_mc_ad.layers:
        mc_utils.write_proj_expr(unified_mc_ad)
    log_norm_expr_ratio = unified_mc_ad.layers['expr'] - unified_mc_ad.layers['proj_expr']
    
    max_abs_log_norm_expression_ratio_20_80_percentiles = np.maximum(
        np.quantile(log_norm_expr_ratio, 0.8, axis=0),
        -np.quantile(log_norm_expr_ratio, 0.2, axis=0),
    )
    # donor_id_to_all_cna_infos = mds_analysis_params.get_donor_id_to_all_cna_infos()

    var_df = hg38_utils.get_var_df(unified_mc_ad.var, gene_loc_df, chrom_pos_df)


    plt.close('all')

    too_often_high_abs_projected_fold = max_abs_log_norm_expression_ratio_20_80_percentiles > get_mds_params()[
        'karyotype_estimation']['max_max_abs_log_norm_expression_ratio_20_80_percentiles']

    low_log_norm_expression_thresh = get_karyotype_estimation_params()['low_log_norm_expression_thresh']

    low_log_norm_expression = unified_mc_ad.layers['expr'] < low_log_norm_expression_thresh
    low_log_norm_projected_expression = unified_mc_ad.layers['proj_expr'] < low_log_norm_expression_thresh

    # TODO: why logical_and and not logical_or? if it is low in the unified model (query), then it is bad because we won't be able to distinguish between different copy numbers. if it is only low in the projection, then it might mask deletions as it would always give a signal for a higher than expected expression...
    low_log_norm_expression_and_projected_expression = np.logical_and(low_log_norm_expression, low_log_norm_projected_expression)

    # TODO: maybe only consider presumably healthy metacells here??
    fraction_of_mcs_with_low_log_norm_expression_and_projected_expression_per_gene = low_log_norm_expression_and_projected_expression.mean(axis=0)

    if 0:
        plt.close('all')
        fig, ax = plt.subplots()
        sb.histplot(fraction_of_mcs_with_low_log_norm_expression_and_projected_expression_per_gene, bins=100, ax=ax)
        sb.histplot((unified_mc_ad.layers['expr'] < low_log_norm_expression_thresh).mean(axis=0), bins=100, ax=ax)
        raise

    max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression = get_karyotype_estimation_params()[
        'max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression']

    too_often_pretty_silent_in_both_query_and_atlas_gene_mask = (
        fraction_of_mcs_with_low_log_norm_expression_and_projected_expression_per_gene > 
        max_fraction_of_mcs_with_low_log_norm_expression_and_projected_expression
    )

    always_pretty_silent_in_query_or_atlas_mask = low_log_norm_expression.all(axis=0) | low_log_norm_projected_expression.all(axis=0)

    genes_to_ignore_names = mds_analysis_params.get_karyotype_estimation_genes_to_ignore_names(
        get_mds_params()['noisy_gene_names']['sex_diff_expressed'], 
        ignore_problematic_ultima_illumina_genes=get_karyotype_estimation_params()['ignore_problematic_ultima_illumina_genes'],
    )

    ignore_gene_mask = (
        unified_mc_ad.var_names.isin(genes_to_ignore_names) 
        | too_often_pretty_silent_in_both_query_and_atlas_gene_mask 
        | var_df['unified_end_pos'].isna()
        # | unified_mc_ad.var['ignored_gene'] # TODO: do we want that???
    )

    if total_num_of_donors >= get_karyotype_estimation_params()['min_total_num_of_donors_for_extra_ignored_genes']:
        num_of_ignored_genes_before = ignore_gene_mask.sum()
        ignore_gene_mask |= (
            always_pretty_silent_in_query_or_atlas_mask
            | too_often_high_abs_projected_fold 
        )
        num_of_ignored_genes_after = ignore_gene_mask.sum()
        curr_num_of_ignored_genes = num_of_ignored_genes_after - num_of_ignored_genes_before
        print(f'total_num_of_donors is high enough ({total_num_of_donors}), so we assume that no CNA appears in all donors, allowing us to add ({curr_num_of_ignored_genes}) extra ignored genes.')

    unified_mc_model_paths = mds_analysis_params.get_mc_model_paths(unified_mc_model_name, dataset_name)
    generic_utils.write_text_file(
        unified_mc_model_paths['genes_excluded_from_cna_due_to_consistent_obs_exp_diff'], 
        '\n'.join(unified_mc_ad.var_names[too_often_high_abs_projected_fold]),
    )
    generic_utils.write_text_file(
        unified_mc_model_paths['genes_excluded_from_cna_due_to_low_expr'], 
        '\n'.join(unified_mc_ad.var_names[too_often_pretty_silent_in_both_query_and_atlas_gene_mask | always_pretty_silent_in_query_or_atlas_mask]),
    )
    generic_utils.write_text_file(
        unified_mc_model_paths['genes_excluded_from_cna_due_to_missing_position_or_on_y'], 
        '\n'.join(unified_mc_ad.var_names[var_df['unified_end_pos'].isna()]),
    )

    gene_mask = np.array(~ignore_gene_mask)
    assert (unified_mc_ad.var_names == var_df.index).all()
    gene_names = list(unified_mc_ad.var_names[gene_mask])

    var_df_mask = gene_mask
    wanted_num_of_genes_in_bin = get_karyotype_estimation_params()['wanted_num_of_genes_in_bin']
    bins, bin_df = hg38_utils.get_bins_with_chroms_as_constraints_and_add_unified_end_pos_bin_i_to_var_df(
        chrom_pos_df,
        var_df,
        wanted_num_of_genes_in_bin=wanted_num_of_genes_in_bin, 
        var_df_mask=var_df_mask,
    )
    num_of_bins = len(bins) - 1
    print(f'num_of_bins: {num_of_bins}')

    if 'ultima_illumina_technical_repeats' in get_mds_params():
        gene_read_gc_df = pd.read_csv(get_mds_params()['ultima_illumina_technical_repeats']['gene_read_gc_df_csv_file_path'])
        gene_read_gc_df.rename(columns={'gene_name': 'gene'}, inplace=True)
    else:
        gene_read_gc_df = None

    num_of_gc_bins = get_karyotype_estimation_params().get('num_of_gc_bins', None)

    return dict(
        chrom_pos_df=chrom_pos_df,
        var_df=var_df,
        gene_mask=gene_mask,
        unified_c_ad=unified_c_ad,
        donor_id_to_exp_sharing_donor_ids=donor_id_to_exp_sharing_donor_ids,
        num_of_gc_bins=num_of_gc_bins,
        gene_read_gc_df=gene_read_gc_df,
        bins=bins,
        bin_df=bin_df,
        gene_names=gene_names,
        atlas_donor_ids=atlas_donor_ids,
        karyotype_cmap=karyotype_cmap,
        epsilon_for_log_projected_fold=epsilon_for_log_projected_fold,
        atlas_to_use=atlas_to_use,
        unified_mc_ad_var_names=unified_mc_ad_var_names,
        wanted_num_of_genes_in_bin=wanted_num_of_genes_in_bin,
        var_df_mask=var_df_mask,
    )



        
def get_bin_projected_fold_ordered_by_loc(
        projected_fold_mat,
        var_df,
        bins,
        gene_mask,
        wanted_num_of_genes_in_bin,
        bin_df,
        mean_or_median='median',
        truncate_x_tick_labels=False,
):
    assert not np.isnan(projected_fold_mat).any()

    num_of_bins = len(bins) - 1
    mean_projected_fold_columns = []
    for bin_i in range(num_of_bins):
            
        bin_gene_mask = var_df['unified_end_pos_bin_i'] == bin_i
        curr_gene_mask = bin_gene_mask & gene_mask
        num_of_genes = curr_gene_mask.sum()
        if num_of_genes:
            curr_gene_mask_after_exclusion = curr_gene_mask[gene_mask]
            
            if curr_gene_mask_after_exclusion.sum() == 1:
                curr_column = mc.ut.to_numpy_vector(projected_fold_mat[:, curr_gene_mask_after_exclusion])
            else:
                if mean_or_median == 'mean':
                    # curr_column = np.nanmean(projected_fold_mat[:, curr_gene_mask_after_exclusion], axis=1)
                    # asserted earlier that there are no nans
                    curr_column = mc.ut.to_numpy_vector(np.mean(projected_fold_mat[:, curr_gene_mask_after_exclusion], axis=1))
                else:
                    assert mean_or_median == 'median'
                    # curr_column = mc.ut.to_numpy_vector(np.nanmedian(projected_fold_mat[:, curr_gene_mask_after_exclusion], axis=1))
                    # asserted earlier that there are no nans
                    curr_column = mc.ut.to_numpy_vector(np.median(projected_fold_mat[:, curr_gene_mask_after_exclusion], axis=1))
            
        else:
            print(f'num_of_genes == 0, bin_i: {bin_i}')
            curr_column = np.full(projected_fold_mat.shape[0], 0)

        mean_projected_fold_columns.append(curr_column)
    mean_projected_fold_mat = np.vstack(mean_projected_fold_columns).transpose()

    nan_mask = np.isnan(mean_projected_fold_mat)
    assert not nan_mask.any()
    if 0:
        num_of_nan_vals_in_mean_projected_fold_mat = nan_mask.sum()
        if num_of_nan_vals_in_mean_projected_fold_mat > 0:
            print(f'WARNING! replacing {num_of_nan_vals_in_mean_projected_fold_mat} nan values in mean_projected_fold_mat (which is displayed in the heatmap) with zeros. (mean_projected_fold_mat.shape: {mean_projected_fold_mat.shape}). this means that {num_of_nan_vals_in_mean_projected_fold_mat} of the WHITE rectangles in the heatmap are missing values, NOT mean_projected_fold=0. beware.')
            mean_projected_fold_mat[nan_mask] = 0
    
    if wanted_num_of_genes_in_bin == 1:
        bin_gene_reprs = [','.join(var_df[(var_df['unified_end_pos_bin_i'] == x) & gene_mask].index) for x in range(num_of_bins)]
        if truncate_x_tick_labels:
            bin_gene_reprs = [x if (len(x) < 20) else (x[:15] + '[...]') for x in bin_gene_reprs]
        xticklabels = bin_gene_reprs
        # xticklabels = generic_utils.dilute_str_list(bin_gene_reprs, 1)
    else:
        xticklabels = [x[:-len('_0')] if x.endswith('_0') else '' for x in list(bin_df['bin_name'])]

    # This must hold for the manipulation of xticklabels later to be valid. 
    assert len(xticklabels) == mean_projected_fold_mat.shape[1]

    return mean_projected_fold_mat, xticklabels

def get_donor_clone_and_cna_info(donor_id):
    if donor_id in get_karyotype_estimation_params()['donor_id_to_clone_and_cna_info']:
        return get_karyotype_estimation_params()['donor_id_to_clone_and_cna_info'][donor_id]
    return None

def plot_multiple_gene_bin_clustermap(
        donor_id,
        mc_ad,
        projected_fold_mat,
        bins,
        bin_df,
        var_df,
        gene_mask,
        wanted_num_of_genes_in_bin,
        karyotype_cmap,
        mc_model_name,
        exp_name,
):
    bin_projected_fold_mat, xticklabels = get_bin_projected_fold_ordered_by_loc(
        projected_fold_mat=projected_fold_mat, var_df=var_df, bins=bins, gene_mask=gene_mask,
        wanted_num_of_genes_in_bin=wanted_num_of_genes_in_bin, bin_df=bin_df,
        # mean_or_median='mean',
    )

    # mc_correlations = mc.ut.corrcoef(
    #     matrix=mc.ut.to_layout(bin_projected_fold_mat, layout='row_major')
    #     per='row',
    #     reproducible=True,
    # )

    clustering_params = get_karyotype_estimation_params()['clustering_cells_or_metacells']
    mc_linkage_mat = scipy.cluster.hierarchy.linkage(
        bin_projected_fold_mat,
        method=clustering_params['linkage_method'], 
        metric=clustering_params['linkage_metric'],
    )

    donor_clone_and_cna_info = get_donor_clone_and_cna_info(donor_id)
    if (donor_clone_and_cna_info is not None) and ('num_of_mc_clusters_for_identifying_clones' in donor_clone_and_cna_info):
        num_of_mc_clusters_for_identifying_clones = donor_clone_and_cna_info['num_of_mc_clusters_for_identifying_clones']
        assert 'mc_clusters_for_identifying_clones_threshold_dist_quantile' not in donor_clone_and_cna_info
        fcluster_t = num_of_mc_clusters_for_identifying_clones
        if isinstance(fcluster_t, dict):
            fcluster_t = fcluster_t[(mc_model_name, exp_name)]
        fcluster_criterion = 'maxclust'
    # else:
    #     num_of_mc_clusters_for_identifying_clones = mc_ad.n_obs // clustering_params['default_num_of_mcs_per_cluster_for_identifying_clones']
    #     num_of_mc_clusters_for_identifying_clones = max(num_of_mc_clusters_for_identifying_clones, clustering_params['min_num_of_clusters_for_identifying_clones'])
    else:
        if (donor_clone_and_cna_info is not None) and ('mc_clusters_for_identifying_clones_threshold_dist_quantile' in donor_clone_and_cna_info):
            threshold_dist_quantile = donor_clone_and_cna_info['mc_clusters_for_identifying_clones_threshold_dist_quantile']
            if isinstance(threshold_dist_quantile, dict):
                threshold_dist_quantile = threshold_dist_quantile[mc_model_name]
        else:
            threshold_dist_quantile = clustering_params['default_mc_clusters_for_identifying_clones_threshold_dist_quantile']
        
        dist_thresh = np.quantile(mc_linkage_mat[:,2], threshold_dist_quantile)
        fcluster_t = dist_thresh
        fcluster_criterion = 'distance'
    # raise
    cluster_i_of_mcs = scipy.cluster.hierarchy.fcluster(
        Z=mc_linkage_mat,
        t=fcluster_t,
        criterion=fcluster_criterion,
        # t=num_of_mc_clusters_for_identifying_clones,
        # criterion='maxclust',
    ) - 1

    cluster_indices = sorted(set(cluster_i_of_mcs))
    cluster_i_to_color = {
        k: v
        # for k,v in zip(list(range(num_of_mc_clusters_for_identifying_clones)), generic_utils.get_n_colors(num_of_mc_clusters_for_identifying_clones))
        for k,v in zip(cluster_indices, generic_utils.get_n_colors(len(cluster_indices)))
    }

    if 'state_color' in mc_ad.obs.columns:
        state_color_of_mcs = mc_ad.obs['state_color']
    else:
        print('NOTE: no state_color column in mc_ad.obs. using black instead.')
        state_color_of_mcs = np.full(mc_ad.n_obs, 'black')
    
    if 'projected_type_color' in mc_ad.obs.columns:
        projected_type_color_of_mcs = mc_ad.obs['projected_type_color']
    else:
        print('NOTE: no projected_type_color column in mc_ad.obs. using black instead.')
        projected_type_color_of_mcs = np.full(mc_ad.n_obs, 'black')
        
    row_colors = [
        state_color_of_mcs,
        projected_type_color_of_mcs,
        generic_utils.get_num_colors(mc_ad.obs['projected_correlation'], vmin=0.75, vmax=1),
        [cluster_i_to_color[x] for x in cluster_i_of_mcs],
    ]

    row_cluster = True

    clustermap_obj = sb.clustermap(
        bin_projected_fold_mat,
        # cmap='RdBu_r',
        cmap=karyotype_cmap,
        # row_cluster=False,
        row_cluster=row_cluster,
        row_linkage=mc_linkage_mat,
        # row_cluster=True,
        col_cluster=False,
        # row_colors=row_colors,
        row_colors=row_colors,
        # cbar_kws=None,
        # yticklabels=yticklabels,
        yticklabels=False,
        xticklabels=xticklabels,
        vmin=-1.5,
        vmax=1.5,
        figsize=(13,10),
        # cbar_pos=None, # if we do this, then clustermap_obj.ax_cbar gives another ax.
        # **curr_clustermap_kwargs,
    )

    heatmap_ax = clustermap_obj.ax_heatmap
    # heatmap_ax.set_title(f'row order: {order_cells_or_metacells_by_desc}')

    if wanted_num_of_genes_in_bin > 1:
        xtick_diffs = np.diff(heatmap_ax.get_xticks())
        bin_size_on_x = xtick_diffs.mean()
        assert np.isclose(xtick_diffs, bin_size_on_x).all()

        chrom_name_and_heatmap_ax_x_of_first_chrom_ticks = [
            (x.get_text(), x.get_position()[0]) for x in heatmap_ax.get_xticklabels() if x.get_text() != '']

        half_bin_size_on_x = bin_size_on_x / 2

        left_xlim, right_xlim = heatmap_ax.get_xlim()
        chrom_name_and_left_and_right_heatmap_ax_xs = [
            (chrom_name, bin_x - half_bin_size_on_x, next_bin_x - half_bin_size_on_x) for (chrom_name, bin_x), (_, next_bin_x)
            in zip(chrom_name_and_heatmap_ax_x_of_first_chrom_ticks, chrom_name_and_heatmap_ax_x_of_first_chrom_ticks[1:] + [('dummy', right_xlim + half_bin_size_on_x)])
        ]

        # heatmap_ax.axvline(left_xlim, color='black')
        new_xticks_and_labels = []
        for chrom_name, bin_left_x, bin_right_x in chrom_name_and_left_and_right_heatmap_ax_xs:
            if not np.isclose(bin_right_x, right_xlim):
                heatmap_ax.axvline(bin_right_x, color='black', linewidth=0.7)
            bin_middle_x = np.mean([bin_left_x, bin_right_x])
            new_xticks_and_labels.append((bin_middle_x, chrom_name))

        assert pd.Series([x[0] for x in new_xticks_and_labels]).is_monotonic_increasing

        heatmap_ax.set_xticks([x[0] for x in new_xticks_and_labels])
        heatmap_ax.set_xticklabels(
            [x[1] for x in new_xticks_and_labels],
            # rotation=45, ha='right', 
            # rotation_mode='anchor', 
            # fontsize='small',
        )
    generic_utils.make_all_spines_visible(heatmap_ax)

    return cluster_i_of_mcs, cluster_i_to_color, clustermap_obj

def get_valid_moving_median_mask(single_gene_bin_df, num_of_bins_truncated_due_to_bin_size_on_each_edge):
    valid_moving_median_mask = np.full(len(single_gene_bin_df), True)
    for chrom_name in sorted(single_gene_bin_df['chrom_name'].unique()):
        chrom_mask = single_gene_bin_df['chrom_name'] == chrom_name
        min_bin = single_gene_bin_df.loc[chrom_mask, 'bin_i'].min()
        max_bin = single_gene_bin_df.loc[chrom_mask, 'bin_i'].max()
        curr_invalid_bins_mask = (
            (single_gene_bin_df['bin_i'] >= min_bin)
            & (single_gene_bin_df['bin_i'] < (min_bin + num_of_bins_truncated_due_to_bin_size_on_each_edge))
        ) | (
            (single_gene_bin_df['bin_i'] > (max_bin - num_of_bins_truncated_due_to_bin_size_on_each_edge))
            & (single_gene_bin_df['bin_i'] <= max_bin)
        )
        expected_num_of_invalid_bins = num_of_bins_truncated_due_to_bin_size_on_each_edge * 2
        curr_num_of_gene_bins = chrom_mask.sum()
        all_chrom_is_invalid = curr_num_of_gene_bins <= expected_num_of_invalid_bins
        if all_chrom_is_invalid:
            print(f'ATTENTION: the number of gene bins in {chrom_name} ({curr_num_of_gene_bins}) is smaller than the moving window size')
        assert (curr_invalid_bins_mask.sum() == expected_num_of_invalid_bins) | all_chrom_is_invalid
        valid_moving_median_mask[curr_invalid_bins_mask] = False
    return valid_moving_median_mask

def plot_moving_median_and_get_cluster_i_to_clone_name(
        donor_id,
        exp_sharing_donor_ids,
        cluster_i_of_mcs,
        single_gene_bin_df,
        moving_window_size,
        valid_moving_median_mask,
        cluster_i_to_color,
        single_gene_bin_projected_fold_mat,
        single_gene_xticklabels,
        var_df,
        cells_or_metacells,
        mean_or_median,
        clusters_to_scatter_plot_indices=None,
        final_clusters_to_plot_indices=None,
        print_non_filtered_cluster_i_to_cna_median_cluster_median_projected_fold=False,
        print_non_filtered_cluster_i_to_cna_names=False,
        show_individual_genes_in_cluster_moving_median_plot=False,
        quantiles_to_plot_hlines=(0.01, 0.99),
        skip_cnas_with_skip_cna_attr=True,
        specific_genes_to_show_names=None,
        show_all_small_clusters_combined=True,
):
    cells_or_mcs = 'mcs' if cells_or_metacells == 'metacells' else 'cells'
    
    if not np.issubdtype(cluster_i_of_mcs.dtype, np.integer):
        orig_cluster_i_of_mcs = cluster_i_of_mcs.copy()
        orig_cluster_is = sorted(set(orig_cluster_i_of_mcs))
        # cluster_i_to_orig_cluster_i = {i: x for i, x in enumerate(orig_cluster_is)}
        orig_cluster_i_to_cluster_i = {x: i for i, x in enumerate(orig_cluster_is)}
        cluster_i_of_mcs = np.array([orig_cluster_i_to_cluster_i[x] for x in orig_cluster_i_of_mcs])
        orig_cluster_i_to_color = cluster_i_to_color.copy()
        print(f'orig_cluster_i_to_cluster_i: {orig_cluster_i_to_cluster_i}')
        cluster_i_to_color = {orig_cluster_i_to_cluster_i[k]: v for k, v in orig_cluster_i_to_color.items() if k in orig_cluster_i_to_cluster_i}
    assert np.issubdtype(cluster_i_of_mcs.dtype, np.integer), f'cluster_i_of_mcs must be integers. cluster_i_of_mcs: {cluster_i_of_mcs}'

    if specific_genes_to_show_names:
        specific_genes_unified_end_positions = var_df.loc[list(specific_genes_to_show_names), 'unified_end_pos']

    chrom_name_to_cna_infos = {}
    cna_infos = get_karyotype_estimation_params()['seemingly_batchy_cna_infos'].copy()
    for exp_sharing_donor_id in exp_sharing_donor_ids:
        curr_donor_clone_and_cna_info = get_donor_clone_and_cna_info(exp_sharing_donor_id)
        if (curr_donor_clone_and_cna_info is not None):
            if ('cna_infos' in curr_donor_clone_and_cna_info):
                curr_cna_infos = curr_donor_clone_and_cna_info['cna_infos']
                if skip_cnas_with_skip_cna_attr:
                    curr_cna_infos = [x for x in curr_cna_infos if 'skip_cna' not in x]
                # NOTE: taking also cna_infos with a skip_cna attr.
                cna_infos.extend(curr_cna_infos)

    chrom_name_to_cna_infos = mds_analysis_params.get_chrom_name_to_cna_infos(cna_infos)

    donor_clone_and_cna_info = get_donor_clone_and_cna_info(donor_id)
    if (donor_clone_and_cna_info is not None) and (f'min_num_of_{cells_or_mcs}_in_cluster_to_show_in_cluster_moving_median_plot' in donor_clone_and_cna_info):
        min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot = donor_clone_and_cna_info[f'min_num_of_{cells_or_mcs}_in_cluster_to_show_in_cluster_moving_median_plot']
    else:
        min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot = get_karyotype_estimation_params()[
            f'default_min_num_of_{cells_or_mcs}_in_cluster_to_show_in_cluster_moving_median_plot']

    cluster_counts = pd.Series(cluster_i_of_mcs).value_counts()
    small_cluster_indices = sorted(cluster_counts[cluster_counts < min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot].index)
    clusters_to_plot_indices = sorted(cluster_counts[cluster_counts >= min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot].index)
    # print(f'clusters_to_plot_indices: {clusters_to_plot_indices}')
    fixed_cluster_i_of_mcs = cluster_i_of_mcs.copy()
    cluster_indices = sorted(set(fixed_cluster_i_of_mcs))
    fixed_cluster_indices = cluster_indices.copy()
    combined_small_clusters_cluster_i = None
    # if small_cluster_indices and (donor_clone_and_cna_info is not None) and ('consider_small_clusters_combined' in donor_clone_and_cna_info):
    if show_all_small_clusters_combined and small_cluster_indices:
        print(f'cluster_indices: {cluster_indices}')
        combined_small_clusters_cluster_i = max(cluster_indices) + 1
        fixed_cluster_i_of_mcs[np.isin(fixed_cluster_i_of_mcs, small_cluster_indices)] = combined_small_clusters_cluster_i
        fixed_cluster_indices = sorted(set(fixed_cluster_i_of_mcs))
        cluster_i_to_color[combined_small_clusters_cluster_i] = 'black'
        clusters_to_plot_indices.append(combined_small_clusters_cluster_i)

    if clusters_to_scatter_plot_indices is None:
        clusters_to_scatter_plot_indices = clusters_to_plot_indices

    assert single_gene_bin_projected_fold_mat.shape[1] == len(single_gene_bin_df)
    dfs = []
    mean_or_median_func = np.mean if mean_or_median == 'mean' else np.median
    cluster_mean_or_median_column_name = f'cluster_{mean_or_median}_projected_fold'
    moving_mean_or_median_column_name = f'moving_median_{cluster_mean_or_median_column_name}'
    for cluster_i in fixed_cluster_indices:
        # if cluster_i != 0:
        #     continue
        # print(f'cluster_i: {cluster_i}')
        cluster_mask = fixed_cluster_i_of_mcs == cluster_i
        num_of_mcs_in_cluster = cluster_mask.sum()
        # if (num_of_mcs_in_cluster < min_num_of_mcs_in_cluster_to_show_in_cluster_moving_median_plot) & (cluster_i != combined_small_clusters_cluster_i):
        #     continue
        # print(cluster_i)
        assert cluster_mask.any()
        df = pd.DataFrame({
            cluster_mean_or_median_column_name: mean_or_median_func(single_gene_bin_projected_fold_mat[cluster_mask, :], axis=0),
            'genes_repr': single_gene_xticklabels,
            'mean_unified_end_pos': single_gene_bin_df['mean_unified_end_pos'],
            'chrom_name': single_gene_bin_df['chrom_name'],
            'bin_i': single_gene_bin_df['bin_i'],
        })
        df['cluster_i'] = cluster_i
        df['cluster_i_color'] = cluster_i_to_color[cluster_i]
        df['num_of_mcs_in_cluster'] = num_of_mcs_in_cluster

        moving_window_mat = np.lib.stride_tricks.sliding_window_view(df[cluster_mean_or_median_column_name], moving_window_size)
        dummy_median_vals_on_each_edge = [np.nan] * (moving_window_size // 2)
        df[moving_mean_or_median_column_name] = (
            dummy_median_vals_on_each_edge + 
            list(
                # np.median(
                # np.mean(
                mean_or_median_func(
                    moving_window_mat, axis=1)) + 
            dummy_median_vals_on_each_edge)
        df.loc[~valid_moving_median_mask, moving_mean_or_median_column_name] = np.nan

        dfs.append(df)

    all_moving_median_cluster_df = pd.concat(dfs, ignore_index=True)

    plot_mask = all_moving_median_cluster_df['cluster_i'].isin(clusters_to_plot_indices)
    if final_clusters_to_plot_indices:
        plot_mask = plot_mask & all_moving_median_cluster_df['cluster_i'].isin(final_clusters_to_plot_indices)
    scatter_plot_mask = all_moving_median_cluster_df['cluster_i'].isin(clusters_to_scatter_plot_indices)

    max_min_y_lim = -1.1
    min_max_y_lim = 0.68
    if show_individual_genes_in_cluster_moving_median_plot:
        y_lims = (min(all_moving_median_cluster_df[cluster_mean_or_median_column_name].min(), max_min_y_lim), max(all_moving_median_cluster_df[cluster_mean_or_median_column_name].max(), min_max_y_lim))
    else:
        y_lims = (min(all_moving_median_cluster_df[moving_mean_or_median_column_name].min(), max_min_y_lim), max(
            all_moving_median_cluster_df[moving_mean_or_median_column_name].max(), min_max_y_lim))

    chrom_names_of_fig_rows = [
        [str(x) for x in range(1, 6)], 
        [str(x) for x in range(6, 13)], 
        [str(x) for x in range(13, 23)] + ['X'], 
    ]
    nrows = len(chrom_names_of_fig_rows)

    fig, axes = plt.subplots(figsize=(15, 7), nrows=nrows, sharey=True, sharex=False)
    cluster_i_to_cna_names = collections.defaultdict(set)
    for curr_chrom_names, ax in zip(chrom_names_of_fig_rows, axes):
        curr_single_gene_bin_df = single_gene_bin_df[single_gene_bin_df['chrom_name'].isin(curr_chrom_names)]

        if show_individual_genes_in_cluster_moving_median_plot:
            sb.scatterplot(
                data=all_moving_median_cluster_df[scatter_plot_mask],
                x='mean_unified_end_pos', 
                y=cluster_mean_or_median_column_name, 
                hue='cluster_i',
                palette=cluster_i_to_color,
                # label=f'cluster {cluster_i}',
                legend=False,
                
                # for general audience?
                # s=6,
                # alpha=0.3,

                # for identifying CNA edges
                s=8,
                alpha=1,

                ax=ax,
            )
        
        assert curr_single_gene_bin_df['bin_i'].is_monotonic_increasing
        prev_chrom_max_unified_end_pos = None
        for chrom_name in curr_single_gene_bin_df['chrom_name'].drop_duplicates():
            chrom_mask = all_moving_median_cluster_df['chrom_name'] == chrom_name

            chrom_min_unified_end_pos = curr_single_gene_bin_df.loc[chrom_mask, 'mean_unified_end_pos'].min()
            chrom_max_unified_end_pos = curr_single_gene_bin_df.loc[chrom_mask, 'mean_unified_end_pos'].max()

            if specific_genes_to_show_names:
                curr_specific_genes_mask = (
                    (specific_genes_unified_end_positions >= chrom_min_unified_end_pos) &
                    (specific_genes_unified_end_positions <= chrom_max_unified_end_pos)
                )
                for specific_gene_unified_end_pos in specific_genes_unified_end_positions[curr_specific_genes_mask]:
                    ax.axvline(specific_gene_unified_end_pos, color='red', linestyle=':', alpha=0.7)

            ax.text(
                np.mean((chrom_min_unified_end_pos, chrom_max_unified_end_pos)), -0.02, chrom_name, ha='center', va='top',
                transform=matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes),
            )
        
            if prev_chrom_max_unified_end_pos is not None:
                ax.axvline(np.mean([prev_chrom_max_unified_end_pos, chrom_min_unified_end_pos]), color='black', linewidth=0.5)
            # sb.scatterplot(
            sb.lineplot(
                data=all_moving_median_cluster_df[chrom_mask & plot_mask],
                x='mean_unified_end_pos', 
                # x='bin_i', # seems worse to me.
                y=moving_mean_or_median_column_name, 
                hue='cluster_i',
                palette=cluster_i_to_color,
                # label=f'cluster {cluster_i}',
                legend=False,
                # s=10,
                ax=ax,
            )

            if chrom_name in chrom_name_to_cna_infos:
                for cna_info in chrom_name_to_cna_infos[chrom_name]:
                    cna_name = cna_info['name']
                    print(f'cna_name: {cna_name}')
                    first_gene_unified_end_pos = var_df.loc[cna_info['first_gene'], 'unified_end_pos'] if 'first_gene' in cna_info else chrom_min_unified_end_pos
                    last_gene_unified_end_pos = var_df.loc[cna_info['last_gene'], 'unified_end_pos'] if 'last_gene' in cna_info else chrom_max_unified_end_pos
                    ax.axvspan(
                        first_gene_unified_end_pos, 
                        last_gene_unified_end_pos, 
                        color='grey', 
                        alpha=0.1,
                    )

                    # NOTE: probably really not important, but note that taking the median of each row and then the median of the n_rows medians might give a different result than taking the median of each column and then the median of the n_cols medians. for example:
                    # mat = np.array([
                    #     [0, 1, 1, 1, 0],
                    #     [0, 3, 3, 3, 0],
                    #     [2, 1, 2, 3, 2],
                    # ])
                    # print(np.median(np.median(mat, axis=1))) # 2
                    # print(np.median(np.median(mat, axis=0))) # 1
                    
                    temp_grouped_df = all_moving_median_cluster_df.loc[
                        (all_moving_median_cluster_df['mean_unified_end_pos'] >= first_gene_unified_end_pos) & 
                        (all_moving_median_cluster_df['mean_unified_end_pos'] <= last_gene_unified_end_pos)
                    ].groupby('cluster_i')[cluster_mean_or_median_column_name]
                    if mean_or_median == 'mean':
                        cluster_i_to_cna_median_cluster_median_projected_fold = temp_grouped_df.mean().to_dict()
                    else:
                        cluster_i_to_cna_median_cluster_median_projected_fold = temp_grouped_df.median().to_dict()
                    if print_non_filtered_cluster_i_to_cna_median_cluster_median_projected_fold:
                        print(f'cluster_i_to_cna_median_cluster_median_projected_fold: {cluster_i_to_cna_median_cluster_median_projected_fold}')

                    filtered_cluster_i_to_cna_median_cluster_median_projected_fold = {
                        k:v for k,v in cluster_i_to_cna_median_cluster_median_projected_fold.items() if k in clusters_to_plot_indices}
                    print(f'filtered_cluster_i_to_cna_median_cluster_median_projected_fold: {filtered_cluster_i_to_cna_median_cluster_median_projected_fold}')
                    for cluster_i, cna_median_cluster_median_projected_fold in cluster_i_to_cna_median_cluster_median_projected_fold.items():
                        if (
                            (
                                ('min_median_cluster_median_projected_fold_threshold' in cna_info) and
                                (cna_median_cluster_median_projected_fold >= cna_info['min_median_cluster_median_projected_fold_threshold'])
                            ) | 
                            (
                                ('max_median_cluster_median_projected_fold_threshold' in cna_info) and
                                (cna_median_cluster_median_projected_fold <= cna_info['max_median_cluster_median_projected_fold_threshold'])
                            )
                        ):
                            cluster_i_to_cna_names[cluster_i].add(cna_info['name'])


            prev_chrom_max_unified_end_pos = chrom_max_unified_end_pos
        
        ax.axhline(0, color='grey', linestyle='-', alpha=0.7)
        ax.axhline(0.58, color='grey', linestyle='--', alpha=0.3)
        ax.axhline(-1, color='grey', linestyle='--', alpha=0.3)
        for quantile in quantiles_to_plot_hlines:
            val = all_moving_median_cluster_df[moving_mean_or_median_column_name].quantile(quantile)
            ax.axhline(val, color='grey', linestyle='-', alpha=0.3)

        
        x_lims = (curr_single_gene_bin_df['mean_unified_end_pos'].min(), curr_single_gene_bin_df['mean_unified_end_pos'].max())
        ax.set_xlim(x_lims)
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            labelbottom=False,
        )
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        # ax.grid(axis='y', alpha=0.3)

        
    axes[-1].set_ylim(y_lims)
    axes[-1].set_xlabel('chromosomal coordinate', labelpad=18)
    axes[1].set_ylabel('moving_median_cluster_median_normalized_projected_fold')

    cluster_i_to_cna_names = dict(cluster_i_to_cna_names) # I don't want a defaultdict moving around.
    if print_non_filtered_cluster_i_to_cna_names:
        print(f'cluster_i_to_cna_names: {cluster_i_to_cna_names}')
    filtered_cluster_i_to_cna_names = {k: v for k,v in cluster_i_to_cna_names.items() if k in clusters_to_plot_indices}
    print(f'filtered_cluster_i_to_cna_names: {filtered_cluster_i_to_cna_names}')

    cluster_i_to_clone_name = {}
    if (donor_clone_and_cna_info is not None) and ('clone_name_to_cna_names' in donor_clone_and_cna_info) and ('cna_infos' in donor_clone_and_cna_info):
        all_donor_cna_names = {x['name'] for x in donor_clone_and_cna_info['cna_infos'] if 'skip_cna' not in x}
        clone_name_to_cna_names = donor_clone_and_cna_info['clone_name_to_cna_names']
        for cluster_i, cna_names in cluster_i_to_cna_names.items():
            for clone_name, clone_cna_names in clone_name_to_cna_names.items():
                if (cna_names & all_donor_cna_names) == set(clone_cna_names):
                    cluster_i_to_clone_name[cluster_i] = clone_name
        print(f'cluster_i_to_clone_name: {cluster_i_to_clone_name}')
    
    return fig, axes, all_moving_median_cluster_df, cluster_i_to_clone_name, cluster_i_to_cna_names




def get_cell_attrs_df_for_write_cna_hists_and_dfs(
        donor_id,
        c_ad,
        mc_ad,
        cluster_i_of_mcs,
        cluster_i_to_color,
        cna_median_normalized_projected_fold_df,
):
    cna_infos = mds_analysis_params.get_donor_id_to_cna_infos()[donor_id]
    cell_attrs_df = generic_utils.merge_preserving_df1_index_and_row_order(
        c_ad.obs[['num_of_non_excluded_umis', 'metacell_i']], 
        pd.DataFrame({
            'mc_cluster_i': cluster_i_of_mcs, 
            'mc_cluster_color': [cluster_i_to_color[x] for x in cluster_i_of_mcs], 
            'metacell_i': mc_ad.obs['metacell_i'].to_numpy(),
        }),
    )

    # cna_median_normalized_projected_fold_df_copy = cna_median_normalized_projected_fold_df.copy()
    # for cna_name in .rename(columns={
    #     x: f'{x}_detected' for x in cna_median_normalized_projected_fold_df.columns})

    extra_cell_attrs_df = mc_ad.obs[['metacell_i']].copy()
    for cna_info in cna_infos:
        cna_name = cna_info['name']
        if set(cna_info) & set(MC_STRICT_THRESHOLD_NAMES):
            if 'min_median_cluster_median_projected_fold_threshold' in cna_info:
                new_column = cna_median_normalized_projected_fold_df[cna_name] >= cna_info['min_median_cluster_median_projected_fold_threshold']
            else:
                new_column = cna_median_normalized_projected_fold_df[cna_name] <= cna_info['max_median_cluster_median_projected_fold_threshold']
            extra_cell_attrs_df[f'{cna_name}_detected_in_mc'] = new_column.to_numpy()
        if set(cna_info) & set(MC_PARTIAL_THRESHOLD_NAMES):
            if 'min_median_cluster_median_projected_fold_partial_threshold' in cna_info:
                new_column = cna_median_normalized_projected_fold_df[cna_name] >= cna_info['min_median_cluster_median_projected_fold_partial_threshold']
            else:
                new_column = cna_median_normalized_projected_fold_df[cna_name] <= cna_info['max_median_cluster_median_projected_fold_partial_threshold']
            extra_cell_attrs_df[f'{cna_name}_partially_detected_in_mc'] = new_column.to_numpy()
    cell_attrs_df = generic_utils.merge_preserving_df1_index_and_row_order(cell_attrs_df, extra_cell_attrs_df)

    return cell_attrs_df

def write_cna_hists_and_dfs(
        donor_id,
        chrom_pos_df,
        var_df,
        single_gene_bin_df,
        out_dir_path,
        single_gene_bin_projected_fold_mat,
        mean_or_median,
        cells_or_metacells,
        cell_attrs_df=None,
        c_ad=None,
        mc_ad=None,
        epsilon_for_log_norm_expr=1e-5,
):
    
    donor_clone_and_cna_info = get_donor_clone_and_cna_info(donor_id)
    cna_infos = [x for x in donor_clone_and_cna_info['cna_infos'] if 'skip_cna' not in x]
    cell_or_mc = 'cell' if cells_or_metacells == 'cells' else 'mc'

    if cells_or_metacells == 'cells':
        mc_cluster_out_dir_path = os.path.join(out_dir_path, 'mc_cluster_hists')
        pathlib.Path(mc_cluster_out_dir_path).mkdir(parents=True, exist_ok=True)
        high_umi_count_out_dir_path = os.path.join(out_dir_path, 'high_umi_count_hists')
        pathlib.Path(high_umi_count_out_dir_path).mkdir(parents=True, exist_ok=True)
        uncorrected_hists_out_dir_path = os.path.join(out_dir_path, 'uncorrected_c_ad_X_vs_projected_hists')
        pathlib.Path(uncorrected_hists_out_dir_path).mkdir(parents=True, exist_ok=True)

    hist_masks_and_titles_and_hues_and_extra_descs = [
        (np.full(single_gene_bin_projected_fold_mat.shape[0], True), '', False, None, out_dir_path, False),
    ]
    if cell_attrs_df is not None:
        assert cells_or_metacells == 'cells'
        for mc_cluster_i in sorted(cell_attrs_df['mc_cluster_i'].unique()):
            cell_mask = cell_attrs_df['mc_cluster_i'] == mc_cluster_i
            mc_cluster_color = cell_attrs_df.loc[cell_mask, 'mc_cluster_color']
            assert mc_cluster_color.nunique() == 1
            mc_cluster_color = mc_cluster_color.iloc[0]
            num_of_mcs_in_cluster = cell_attrs_df.loc[cell_mask, 'metacell_i'].nunique()
            # print(f'num_of_mcs_in_cluster: {num_of_mcs_in_cluster}')
            if (donor_clone_and_cna_info is not None) and ('min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist' in donor_clone_and_cna_info):
                min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist = donor_clone_and_cna_info['min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist']
            else:
                min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist = get_karyotype_estimation_params()['default_min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist']
            
            if num_of_mcs_in_cluster >= min_num_of_mcs_in_cluster_to_show_its_cna_cell_hist:
                hist_masks_and_titles_and_hues_and_extra_descs.append(
                    (cell_mask, f'mc_cluster_{mc_cluster_i}_{mc_cluster_color}', True, None, mc_cluster_out_dir_path, False))

        for min_quantile in [0.5, 0.8, 0.9]:
            min_num_of_non_excluded_umis = cell_attrs_df['num_of_non_excluded_umis'].quantile(min_quantile)
            cell_mask = cell_attrs_df['num_of_non_excluded_umis'] >= min_num_of_non_excluded_umis
            hist_masks_and_titles_and_hues_and_extra_descs.append(
                (cell_mask, f'total_umis_min_quantile_{min_quantile}', False, f'>= {min_num_of_non_excluded_umis:.1f}', high_umi_count_out_dir_path, False))
        
        for attr_column_name in cell_attrs_df.select_dtypes(include=[bool]).columns:
            curr_out_dir_path = os.path.join(out_dir_path, f'{attr_column_name}_hists')
            pathlib.Path(curr_out_dir_path).mkdir(parents=True, exist_ok=True)
            cell_mask = cell_attrs_df[attr_column_name]
            hist_masks_and_titles_and_hues_and_extra_descs.append((cell_mask, attr_column_name, True, None, curr_out_dir_path, True))


    cna_name_to_cna_median_normalized_projected_fold = {}

    hist_x_column_name = f'cna_{mean_or_median}_normalized_projected_fold'
    for cna_info in cna_infos:
        cna_name = cna_info['name']
        extended_cna_info = get_extended_cna_info(cna_info, chrom_pos_df, var_df)

        first_gene_unified_end_pos = extended_cna_info['first_gene_unified_end_pos']
        last_gene_unified_end_pos = extended_cna_info['last_gene_unified_end_pos']
        
        if c_ad is not None:
            assert mc_ad is not None
            assert cell_or_mc == 'cell'
            cna_gene_names = var_df[(var_df['unified_end_pos'] >= first_gene_unified_end_pos) & (var_df['unified_end_pos'] <= last_gene_unified_end_pos)].index
            cna_gene_names = list(set(cna_gene_names) & set(c_ad.var_names) & set(mc_ad.var_names))
            c_ad_cna_gene_mask = c_ad.var_names.isin(cna_gene_names)
            cna_c_df = pd.DataFrame({
                'cna_region_expr': np.log2((mc.ut.to_numpy_vector(c_ad.X[:, c_ad_cna_gene_mask].sum(axis=1)) / c_ad.obs['num_of_non_excluded_umis']) + epsilon_for_log_norm_expr),
                'metacell_name': c_ad.obs['metacell_name'],
            })
            # TODO: discard outliers. in the current model, i had no outliers, i think...
            # cna_c_df = cna_c_df[cna_c_df]
            
            mc_ad_cna_gene_mask = mc_ad.var_names.isin(cna_gene_names)
            cna_c_df = generic_utils.merge_preserving_df1_index_and_row_order(
                cna_c_df, 
                pd.DataFrame({
                    'projected_cna_region_expr': np.log2(mc.ut.to_numpy_vector(mc_ad.layers['projected_fraction'][:, mc_ad_cna_gene_mask].sum(axis=1)) + epsilon_for_log_norm_expr),
                    'metacell_name': mc_ad.obs['metacell_name'],
                }),
            )
            cna_c_df['cna_region_expr_log_ratio'] = cna_c_df['cna_region_expr'] - cna_c_df['projected_cna_region_expr']

            fig, ax = plt.subplots()
            sb.histplot(cna_c_df['cna_region_expr_log_ratio'], ax=ax)
            for x in (0, 0.5, -0.5):
                ax.axvline(x, color='red', linestyle='--', alpha=0.3)
            
            ax_title = f'{cna_name}_uncorrected_X'
            ax.set_title(ax_title)
            fig.savefig(os.path.join(uncorrected_hists_out_dir_path, f'{ax_title}_cna_region_expr_log_ratio_hist.png'))
            # raise



        cna_gene_mask = (
            (single_gene_bin_df['mean_unified_end_pos'] >= first_gene_unified_end_pos) &
            (single_gene_bin_df['mean_unified_end_pos'] <= last_gene_unified_end_pos)
        )
        # print(cna_name, cna_gene_mask.sum())
        assert single_gene_bin_df['mean_unified_end_pos'].is_monotonic_increasing
        mean_or_median_func = np.median if (mean_or_median == 'median') else np.mean
        cna_median_normalized_projected_fold_of_mcs = mean_or_median_func(single_gene_bin_projected_fold_mat[:, cna_gene_mask], axis=1)
        cna_name_to_cna_median_normalized_projected_fold[cna_name] = cna_median_normalized_projected_fold_of_mcs


        for curr_mask, curr_title, use_hue, extra_desc, curr_out_dir_path, require_cna_name_contained_in_title in hist_masks_and_titles_and_hues_and_extra_descs:
            if require_cna_name_contained_in_title:
                if cna_name not in curr_title:
                    continue
            curr_df = pd.DataFrame({
                hist_x_column_name: cna_median_normalized_projected_fold_of_mcs,
                'mask': curr_mask,
            })
            
            plt.close('all')
            fig, ax = plt.subplots()
            histplot_kwargs = dict(
                x=hist_x_column_name,
                ax=ax,
                bins=cna_info.get('mc_median_normalized_projected_fold_hist_bins', 'auto'),
            )
            
            if use_hue:
                histplot_kwargs = {
                    **histplot_kwargs,
                    **dict(
                        data=curr_df,
                        hue='mask',
                        stat='proportion',
                        common_norm=False,
                    ),
                }
            else:
                histplot_kwargs = {
                    **histplot_kwargs,
                    **dict(data=curr_df[curr_mask]),
                }

            sb.histplot(**histplot_kwargs)
            for thresh_column_name in ALL_MC_THRESHOLD_NAMES if (cells_or_metacells == 'metacells') else CELL_THRESHOLD_NAMES:
                if thresh_column_name in cna_info:
                    thresh = cna_info[thresh_column_name]
                    ax.axvline(thresh, color='red', linestyle='--', alpha=0.3)
            
            ax_title = curr_title
            if ax_title:
                ax_title += ' '
            ax_title += f'{cna_name}\n'
            if extra_desc is not None:
                ax_title += f' ({extra_desc})'
            if (cells_or_metacells == 'cells') and (~curr_mask).any():
                num_of_true_cells = curr_mask.sum()
                ax_title += f' (#true_cells={num_of_true_cells})'

            ax.set_title(ax_title)
            if require_cna_name_contained_in_title:
                fig_file_name = f'{curr_title}_{mean_or_median}_normalized_projected_fold_hist.png'
            else:
                fig_file_name = f'{cna_name}_{curr_title}_{mean_or_median}_normalized_projected_fold_hist.png'

            fig.savefig(os.path.join(curr_out_dir_path, fig_file_name))

    cna_median_normalized_projected_fold_df = pd.DataFrame(cna_name_to_cna_median_normalized_projected_fold)
    return cna_median_normalized_projected_fold_df

def get_clone_name_to_color(donor_clone_and_cna_info):
    clone_names = sorted(donor_clone_and_cna_info['clone_name_to_cna_names_and_detection_levels'])
    num_of_clones = len(clone_names)
    clone_name_to_color = {
        k: v
        for k,v in zip(clone_names, generic_utils.get_n_colors(num_of_clones))
    }
    assert 'grey' not in clone_name_to_color.values()
    clone_name_to_color['nan'] = 'grey'
    
    return clone_name_to_color

def calculate_cna_clusters_and_assign_clones(
        cna_median_normalized_projected_fold_df,
        curr_cells_or_metacells_ad,
        donor_id,
        karyotype_cmap,
        donor_exp_out_dir_path,
        cells_or_metacells,
        mean_or_median,
        cnas_to_use_for_clustering_names=None,
):
    donor_clone_and_cna_info = get_donor_clone_and_cna_info(donor_id)
    cna_name_to_cna_info = mds_analysis_params.get_cna_name_to_cna_info(donor_id)
    cna_median_normalized_projected_fold_df = cna_median_normalized_projected_fold_df.copy() # because we are editing this df, and the caller is not expecting that.
    
    if cnas_to_use_for_clustering_names is None:
        cnas_to_use_for_clustering_names = [
            x for x in list(cna_median_normalized_projected_fold_df) 
            if set(cna_name_to_cna_info[x]) & set(CELL_THRESHOLD_NAMES if cells_or_metacells == 'cells' else ALL_MC_THRESHOLD_NAMES)]
    
    if not cnas_to_use_for_clustering_names:
        return None

    filtered_cna_median_normalized_projected_fold_df = cna_median_normalized_projected_fold_df[cnas_to_use_for_clustering_names]


    cell_or_mc = 'mc' if cells_or_metacells == 'metacells' else 'cell'

    single_cna = filtered_cna_median_normalized_projected_fold_df.shape[1] == 1

    if single_cna:
        if 'single_cna_median_projected_fold_bins' in donor_clone_and_cna_info:
            bins = donor_clone_and_cna_info.get('single_cna_median_projected_fold_bins')
        else:
            relevant_threshs = (
                (
                    # NOTE: putting MC_STRICT_THRESHOLD_NAMES before MC_PARTIAL_THRESHOLD_NAMES here means that later we would prefer to use the strict threshold, if it exists.
                    MC_STRICT_THRESHOLD_NAMES + MC_PARTIAL_THRESHOLD_NAMES
                ) if (cells_or_metacells == 'metacells') else CELL_THRESHOLD_NAMES
            )
            single_cna_info = [x for x in donor_clone_and_cna_info['cna_infos'] if x['name'] in cnas_to_use_for_clustering_names]
            assert len(single_cna_info) == 1
            single_cna_info = single_cna_info[0]
            bins = None
            for thresh in relevant_threshs:
                if thresh in single_cna_info:
                    bins = (-np.inf, single_cna_info[thresh], np.inf)
                    # raise RuntimeError(thresh + '\n\n' + str(bins))
                    break
            assert bins is not None
        cna_cluster_i_of_mcs = generic_utils.safe_digitize(
            filtered_cna_median_normalized_projected_fold_df.iloc[:, 0], 
            bins=bins,
        )
        sort_indices = np.argsort(filtered_cna_median_normalized_projected_fold_df.iloc[:, 0])
        row_cluster = False
    else:
        clustering_params = get_karyotype_estimation_params()['clustering_cells_or_metacells']
        mc_cna_linkage_mat = scipy.cluster.hierarchy.linkage(
            filtered_cna_median_normalized_projected_fold_df.to_numpy(),
            method=clustering_params['linkage_method'], 
            metric=donor_clone_and_cna_info.get('cna_cluster_linkage_metric', clustering_params['linkage_metric']),
        )

        # if (donor_clone_and_cna_info is not None) and ('num_of_mc_clusters_for_identifying_clones' in donor_clone_and_cna_info):
        #     num_of_mc_clusters_for_identifying_clones = donor_clone_and_cna_info['num_of_mc_clusters_for_identifying_clones']
        # else:
        #     num_of_mc_clusters_for_identifying_clones = mc_ad.n_obs // clustering_params['default_num_of_mcs_per_cluster_for_identifying_clones']
        #     num_of_mc_clusters_for_identifying_clones = max(num_of_mc_clusters_for_identifying_clones, clustering_params['min_num_of_clusters_for_identifying_clones'])

        if (donor_clone_and_cna_info is not None) and (f'{cell_or_mc}_cna_clusters_threshold_dist_quantile' in donor_clone_and_cna_info):
            threshold_dist_quantile = donor_clone_and_cna_info[f'{cell_or_mc}_cna_clusters_threshold_dist_quantile']
        else:
            threshold_dist_quantile = clustering_params['default_mc_cna_clusters_threshold_dist_quantile']
        dist_thresh = np.quantile(mc_cna_linkage_mat[:,2], threshold_dist_quantile)
        cna_cluster_i_of_mcs = scipy.cluster.hierarchy.fcluster(
            Z=mc_cna_linkage_mat,
            t=dist_thresh,
            criterion='distance',
            # t=donor_clone_and_cna_info.get(f'num_of_{cell_or_mc}_cna_clusters', clustering_params['default_num_of_mc_cna_clusters']),
            # criterion='maxclust',
        ) - 1

        sort_indices = list(range(curr_cells_or_metacells_ad.n_obs))
        row_cluster = True

        # row_colors.append()

    cna_cluster_indices = sorted(set(cna_cluster_i_of_mcs))
    num_of_mc_cna_clusters = len(cna_cluster_indices)
    cna_cluster_i_to_color = {
        k: v
        for k,v in zip(cna_cluster_indices, generic_utils.get_n_colors(num_of_mc_cna_clusters))
    }


    cna_infos = mds_analysis_params.get_donor_id_to_cna_infos()[donor_id]

    cna_median_normalized_projected_fold_df_for_clustermap = cna_median_normalized_projected_fold_df.copy()
    cna_median_normalized_projected_fold_df['cna_cluster_i'] = cna_cluster_i_of_mcs
    df_grouped_by_cna_cluster = cna_median_normalized_projected_fold_df.groupby('cna_cluster_i')
    cna_cluster_df = df_grouped_by_cna_cluster.median() if (mean_or_median == 'median') else df_grouped_by_cna_cluster.mean()
    cna_detected_column_names = []
    for cna_info in cna_infos:
        cna_name = cna_info['name']
        cna_detected_column_name = f'{cna_name}_detected'

        # NOTE: assignment of 'partially' before 'yes' means that 'yes' will overwrite 'partially', as it should.
        if {'max_median_cluster_median_projected_fold_threshold', 'max_median_cluster_median_projected_fold_partial_threshold'} & set(cna_info):
            if 'max_median_cluster_median_projected_fold_partial_threshold' in cna_info:
                cna_cluster_df.loc[cna_cluster_df[cna_name] <= cna_info['max_median_cluster_median_projected_fold_partial_threshold'], cna_detected_column_name] = 'partially'
            if 'max_median_cluster_median_projected_fold_threshold' in cna_info:
                cna_cluster_df.loc[cna_cluster_df[cna_name] <= cna_info['max_median_cluster_median_projected_fold_threshold'], cna_detected_column_name] = 'yes'
        elif {'min_median_cluster_median_projected_fold_threshold', 'min_median_cluster_median_projected_fold_partial_threshold'} & set(cna_info):
            if 'min_median_cluster_median_projected_fold_partial_threshold' in cna_info:
                cna_cluster_df.loc[cna_cluster_df[cna_name] >= cna_info['min_median_cluster_median_projected_fold_partial_threshold'], cna_detected_column_name] = 'partially'
            if 'min_median_cluster_median_projected_fold_threshold' in cna_info:
                cna_cluster_df.loc[cna_cluster_df[cna_name] >= cna_info['min_median_cluster_median_projected_fold_threshold'], cna_detected_column_name] = 'yes'
        
        column_names = set(cna_cluster_df)
        if cna_detected_column_name in column_names:
            cna_cluster_df[cna_detected_column_name].fillna('no', inplace=True)
            cna_detected_column_names.append(cna_detected_column_name)


    curr_cells_or_metacells_ad.obs['cna_cluster_i'] = cna_cluster_i_of_mcs
    mc_obs_grouped_by_cna_cluster = curr_cells_or_metacells_ad.obs.groupby('cna_cluster_i')
    cluster_size_df = mc_obs_grouped_by_cna_cluster.size().reset_index(name=f'num_of_{cells_or_metacells}')
    if cells_or_metacells == 'metacells':
        assert 'num_of_cells' not in list(cluster_size_df)
        cluster_size_df = generic_utils.merge_preserving_df1_index_and_row_order(cluster_size_df, mc_obs_grouped_by_cna_cluster['grouped'].sum().reset_index(name='cell_count'))
        cluster_size_df['fraction_of_cells_with_metacells'] = cluster_size_df['cell_count'] / cluster_size_df['cell_count'].sum()

    cna_cluster_df = generic_utils.merge_preserving_df1_index_and_row_order(cluster_size_df, cna_cluster_df.reset_index())

    row_colors = [
        curr_cells_or_metacells_ad.obs['state_color'],
    ]
    if cells_or_metacells == 'metacells':
        row_colors.append(generic_utils.get_num_colors(curr_cells_or_metacells_ad.obs['projected_correlation'], vmin=0.75, vmax=1))
    # if cells_or_metacells == 'cells':
    #     row_colors.append(generic_utils.get_num_colors(curr_cells_or_metacells_ad.obs['projected_correlation'], vmin=-4e-4, vmax=4e-4))
    row_colors.append([cna_cluster_i_to_color[x] for x in cna_cluster_i_of_mcs])
    if cells_or_metacells == 'cells':
        row_colors.append(generic_utils.get_num_colors(np.log2(curr_cells_or_metacells_ad.obs['num_of_non_excluded_umis']), vmin=9, vmax=12))

    if 'clone_name_to_cna_names_and_detection_levels' in donor_clone_and_cna_info and donor_clone_and_cna_info['clone_name_to_cna_names_and_detection_levels']:
        donor_clone_and_cna_info = get_donor_clone_and_cna_info(donor_id)
        cna_cluster_df.insert(1, 'clone_name', 'nan')
        for clone_name, cna_names_and_detection_levels in donor_clone_and_cna_info['clone_name_to_cna_names_and_detection_levels'].items():
            clone_mask = np.full(len(cna_cluster_df), True)
            for cna_name, detection_level in cna_names_and_detection_levels:
                cna_detected_column_name = f'{cna_name}_detected'
                assert cna_detected_column_name in list(cna_cluster_df), 'a cna in clone_name_to_cna_names_and_detection_levels is missing from cna_infos'
                clone_mask &= cna_cluster_df[cna_detected_column_name] == detection_level
            assert (cna_cluster_df.loc[clone_mask, 'clone_name'] == 'nan').all()
            cna_cluster_df.loc[clone_mask, 'clone_name'] = clone_name

        clone_name_of_mcs = generic_utils.merge_preserving_df1_index_and_row_order(
            cna_median_normalized_projected_fold_df[['cna_cluster_i']], cna_cluster_df[['cna_cluster_i', 'clone_name']])['clone_name'].to_numpy()
        
        clone_name_to_color = get_clone_name_to_color(donor_clone_and_cna_info)
        row_colors.append([clone_name_to_color[x] for x in clone_name_of_mcs])
    else:
        clone_name_of_mcs = clone_name_to_color = None

    abs_vmin_vmax = donor_clone_and_cna_info.get('cna_clustermap_abs_vmin_vmax', 1.5) 

    clustermap_kwargs = dict(
        # data=cna_median_normalized_projected_fold_df_for_clustermap,
        data=cna_median_normalized_projected_fold_df_for_clustermap.iloc[sort_indices],
        # data=cna_median_normalized_projected_fold_df_for_clustermap[list(cna_median_normalized_projected_fold_df_for_clustermap)[:1]],
        # cmap='RdBu_r',
        cmap=karyotype_cmap,
        # row_cluster=False,
        row_cluster=row_cluster,
        col_cluster=False,
        # row_colors=row_colors,
        row_colors=[
            np.array(x)[sort_indices] for x in row_colors
        ],
        # cbar_kws=None,
        # yticklabels=yticklabels,
        yticklabels=False,
        xticklabels=list(cna_median_normalized_projected_fold_df_for_clustermap),
        vmin=-abs_vmin_vmax,
        vmax=abs_vmin_vmax,
        figsize=(donor_clone_and_cna_info.get('cna_clustermap_width', 5), donor_clone_and_cna_info.get('cna_clustermap_height', 10)),
        # cbar_pos=None, # if we do this, then clustermap_obj.ax_cbar gives another ax.
        # **curr_clustermap_kwargs,
        colors_ratio=donor_clone_and_cna_info.get('cna_clustermap_colors_ratio', 0.15),
    )
    if not single_cna:
        clustermap_kwargs = {
            **clustermap_kwargs, 
            **dict(
                row_linkage=mc_cna_linkage_mat,
            ),
        }

    plt.close('all')

    clustermap_obj = sb.clustermap(**clustermap_kwargs)
    # clustermap_obj = sb.clustermap(cna_median_normalized_projected_fold_df.iloc[sort_indices], row_cluster=False, col_cluster=False)
    clustermap_obj.savefig(os.path.join(donor_exp_out_dir_path, 'cna_clustermap.png'))

    return cna_cluster_df, cna_cluster_i_of_mcs, cna_cluster_i_to_color, clone_name_of_mcs, clone_name_to_color


def plot_mc_vs_mean_cell_normalized_projected_fold_scatter_per_cna(
        cna_median_normalized_projected_fold_df,
        cells_cna_mean_normalized_projected_fold_df,
        c_ad,
        mc_ad,
        out_dir_path,
):
    cells_cna_df = cells_cna_mean_normalized_projected_fold_df.copy()
    cna_mc_df = cna_median_normalized_projected_fold_df.copy()
    cna_names = list(cna_median_normalized_projected_fold_df)

    cells_cna_df['metacell_i'] = c_ad.obs['metacell_i'].astype(int).to_numpy() # without to_numpy(), it doesnt work. it tries to use the index, i guess?
    cna_mc_df['metacell_i'] = mc_ad.obs['metacell_i'].astype(int).to_numpy()
    cna_mc_df['state'] = mc_ad.obs['state'].to_numpy()

    cna_mc_df = generic_utils.merge_preserving_df1_index_and_row_order(
        cna_mc_df, 
        cells_cna_df.groupby('metacell_i').mean().reset_index().rename(columns={x: f'cell_mean_{x}' for x in list(cells_cna_df) if x != 'metacell_i'}),
    )
    
    plt.close('all')
    for cna_name in cna_names:
        fig, ax = plt.subplots()
        sb.scatterplot(
            data=cna_mc_df,
            x=f'cell_mean_{cna_name}',
            y=cna_name,
            hue='state',
            palette=mc_utils.get_palette(mc_ad, color_by='state'),
            ax=ax,
            legend=False,
        )
        generic_utils.plot_y_equals_x_line_on_ax(ax)
        fig.savefig(os.path.join(out_dir_path, f'{cna_name}_mc_vs_cell_mean_normalized_projected_fold_scatter.png'))

def write_cell_projected_fold_df_csv(
        cna_median_normalized_projected_fold_df,
        cells_cna_mean_normalized_projected_fold_df,
        c_ad,
        mc_ad,
        out_dir_path,
):
    cell_projected_fold_df = cells_cna_mean_normalized_projected_fold_df.copy()
    cell_projected_fold_df['metacell_i'] = c_ad.obs['metacell_i'].astype(int).to_numpy()
    cell_projected_fold_df.index = c_ad.obs_names.to_numpy()
    mc_projected_fold_df = cna_median_normalized_projected_fold_df.rename(
        columns={x: f'mc_{x}' for x in list(cna_median_normalized_projected_fold_df)})
    mc_projected_fold_df['metacell_i'] = mc_ad.obs['metacell_i'].astype(int).to_numpy()
    cell_projected_fold_df = generic_utils.merge_preserving_df1_index_and_row_order(cell_projected_fold_df, mc_projected_fold_df)
    cell_projected_fold_df.drop(columns=['metacell_i'], inplace=True)
    # cell_projected_fold_df['total_num_of_umis_in_cna_region'] = mc.ut.to_numpy_vector(c_ad.X[:, var_df_cna_only_mask].sum(axis=1))
    
    cell_projected_fold_df_csv_file_path = os.path.join(out_dir_path, 'cell_projected_fold_df.csv')
    cell_projected_fold_df.to_csv(cell_projected_fold_df_csv_file_path)

def plot_cell_cna_state_boxplots(
        cells_cna_mean_normalized_projected_fold_df,
        c_ad,
        out_dir_path,
):
    df = cells_cna_mean_normalized_projected_fold_df.copy()
    df['state'] = c_ad.obs['state'].to_numpy()
    curr_cell_states = set(df['state'])
    ordered_types = [x for x in mds_analysis_params.ORDERED_CELL_STATE_NAMES if x in curr_cell_states]

    print(df['state'].value_counts())

    plt.close('all')
    for cna_name in list(cells_cna_mean_normalized_projected_fold_df):
        fig, ax = plt.subplots(figsize=(10, 5))
        # sb.swarmplot(
        # sb.stripplot(
        # sb.violinplot(
        sb.boxplot(
            data=df, 
            # data=df[df[cna_name] <= 0], 
            x='state', 
            y=cna_name, 
            hue='state', 
            palette=mc_utils.get_palette(c_ad, color_by='state'),
            order=ordered_types,
            dodge=False,
            # legend=False,
            fliersize=3,
            flierprops=dict(alpha=0.5),
            ax=ax,
            notch=True,
            # medianprops=dict(color={"color": "coral"}
        )
        ax.set_xlabel(None)
        ax.set_ylabel('mean_normalized_projected_fold')
        ax.axhline(0, linestyle='--', color='black', alpha=0.5)
        ax.get_legend().remove()
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45, ha='right', 
            rotation_mode='anchor', 
            # fontsize='small',
        )
        ax.set_title(cna_name)
        fig.tight_layout()

        fig_file_name = f'{cna_name}_state_boxplots.png'
        fig.savefig(os.path.join(out_dir_path, fig_file_name))
        
def get_extended_cna_info(cna_info, chrom_pos_df, var_df):
    extended_cna_info = cna_info.copy()
    chrom_row = chrom_pos_df.loc[chrom_pos_df['chrom_name'] == cna_info['chrom']].iloc[0]
    extended_cna_info['first_gene_unified_end_pos'] = var_df.loc[
        cna_info['first_gene'], 'unified_end_pos'] if 'first_gene' in cna_info else chrom_row['chrom_unified_start_pos']
    extended_cna_info['last_gene_unified_end_pos'] = var_df.loc[
        cna_info['last_gene'], 'unified_end_pos'] if 'last_gene' in cna_info else chrom_row['chrom_unified_end_pos']
    extended_cna_info['chrom_size'] = chrom_row['annotated_chrom_len']
    return extended_cna_info

chrom_name_to_centromere_approx_pos = {
    # just by https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38
    '11': int(53e6),
    '5': int(48.5e6),
    '12': int(35.5e6),
    '20': int(28e6),
    '7': int(60e6),
    '17': int(25e6), # https://mds-risk-model.com/ requires this
    '3': int(92e6),
}

def get_genes_in_region(var_df, chrom_name, p_arm_only=False, q_arm_only=False):
    gene_mask = var_df['chrom_name'] == chrom_name
    if q_arm_only:
        assert not p_arm_only
        gene_mask &= var_df['end_pos_in_chrom'] >= chrom_name_to_centromere_approx_pos[chrom_name]
    if p_arm_only:
        assert not q_arm_only
        gene_mask &= var_df['start_pos_in_chrom'] <= chrom_name_to_centromere_approx_pos[chrom_name]
    return var_df[gene_mask].index
        
def get_cna_attrs_df(
        donor_id,
        chrom_pos_df,
        var_df,
        gene_mask,
):
    donor_clone_and_cna_info = get_donor_clone_and_cna_info(donor_id)
    if (donor_clone_and_cna_info is None) or ('cna_infos' not in donor_clone_and_cna_info):
        return None
    cna_infos = [x for x in donor_clone_and_cna_info['cna_infos'] if 'skip_cna' not in x]

    flat_dicts = []
    if 'donor_general_manual_comment' in donor_clone_and_cna_info:
        flat_dicts.append({'donor_general_manual_comment': donor_clone_and_cna_info['donor_general_manual_comment']})
    
    any_cna_info_contains_exp_name = False
    for cna_info in cna_infos:
        extended_cna_info = get_extended_cna_info(cna_info, chrom_pos_df, var_df)
        cna_size = extended_cna_info['last_gene_unified_end_pos'] - extended_cna_info['first_gene_unified_end_pos'] + 1
        whole_chrom_cna = cna_size == extended_cna_info['chrom_size']
        cna_gene_mask = (var_df['chrom_name'] == cna_info['chrom']) if whole_chrom_cna else (
            (var_df['unified_end_pos'] >= extended_cna_info['first_gene_unified_end_pos']) & 
            (var_df['unified_end_pos'] <= extended_cna_info['last_gene_unified_end_pos'])
        )
        cna_genes = list(var_df[cna_gene_mask].index)
        first_gene = cna_genes[0]
        last_gene = cna_genes[-1]
        filtered_cna_genes = list(var_df[cna_gene_mask & gene_mask].index)
        
        cna_genes_repr = '' if whole_chrom_cna else ', '.join(cna_genes)
        num_of_genes = len(cna_genes)
        
        filtered_cna_genes_repr = ', '.join(filtered_cna_genes)
        num_of_filtered_genes = len(filtered_cna_genes)

        manual_comment = cna_info.get('manual_comment', '')
        manual_comment_for_blood_aging_paper = cna_info.get('manual_comment_for_blood_aging_paper', '')

        flat_dict = {
            'name': cna_info['name'],
            'chrom': cna_info['chrom'],
            'whole_chrom': whole_chrom_cna,
            'first_gene': first_gene,
            'last_gene': last_gene,
            'genes': cna_genes_repr,
            'filtered_genes': filtered_cna_genes_repr,
            'num_of_genes': num_of_genes,
            'num_of_filtered_genes': num_of_filtered_genes,
            'manual_comment': manual_comment,
            'manual_comment_for_blood_aging_paper': manual_comment_for_blood_aging_paper,
        }

        if cna_info['chrom'] == '5':
            flat_dict['contains_CSNK1A1'] = 'CSNK1A1' in cna_genes

        if 'exp_names' in extended_cna_info:
            any_cna_info_contains_exp_name = True
        else:
            assert not any_cna_info_contains_exp_name, 'if a donor has multiple samples, exp_names should exist for every CNA'
        exp_names = extended_cna_info.get('exp_names', ['all'])

        for exp_name in exp_names:
            curr_flat_dict = {'exp_name': exp_name, **flat_dict}
            if exp_name != 'all':
                curr_flat_dict['name'] += f'_{exp_name}'
            flat_dicts.append(curr_flat_dict)

    if not flat_dicts:
        return None
    cna_attrs_df = pd.DataFrame(flat_dicts)

    if 'donor_general_manual_comment' not in list(cna_attrs_df.columns):
        cna_attrs_df['donor_general_manual_comment'] = np.nan

    del_mask = cna_attrs_df['name'].str.contains('_del_')
    dup_mask = cna_attrs_df['name'].str.contains('_dup_')
    assert not (del_mask & dup_mask).any()
    neither_dul_or_dup_mask = ~(del_mask | dup_mask) & cna_attrs_df['donor_general_manual_comment'].isna()
    assert not neither_dul_or_dup_mask.any(), f'cna_attrs_df[neither_dul_or_dup_mask]:\n{cna_attrs_df[neither_dul_or_dup_mask]}'
    cna_attrs_df['is_del'] = del_mask
    cna_attrs_df['is_dup'] = dup_mask
    
    # print((cna_attrs_df['donor_general_manual_comment'] == cna_attrs_df['name'].isna()).value_counts())
    assert (cna_attrs_df['donor_general_manual_comment'].notna() == cna_attrs_df['name'].isna()).all()

    assert (cna_attrs_df.loc[cna_attrs_df['donor_general_manual_comment'].isna(), 'first_gene']).isin(var_df.index).all()
    assert (cna_attrs_df.loc[cna_attrs_df['donor_general_manual_comment'].isna(), 'last_gene']).isin(var_df.index).all()
    
    for chrom_name in chrom_name_to_centromere_approx_pos:
        # print('verify in genome browser:', chrom_name, get_genes_in_region(var_df, chrom_name, p_arm_only=True)) # done verifying
        for del_or_dup in ['del', 'dup']:
            col = f'is_{del_or_dup}_{chrom_name}p'
            cna_attrs_df[col] = cna_attrs_df['donor_general_manual_comment'].isna() & (cna_attrs_df['whole_chrom'] == False) & cna_attrs_df[f'is_{del_or_dup}'] & cna_attrs_df['first_gene'].isin(get_genes_in_region(var_df, chrom_name, p_arm_only=True))
            assert not ((cna_attrs_df[col] == True) & (cna_attrs_df['chrom'] != chrom_name)).any()
    
            col = f'is_{del_or_dup}_{chrom_name}q'
            cna_attrs_df[col] = cna_attrs_df['donor_general_manual_comment'].isna() & (cna_attrs_df['whole_chrom'] == False) & cna_attrs_df[f'is_{del_or_dup}'] & cna_attrs_df['last_gene'].isin(get_genes_in_region(var_df, chrom_name, q_arm_only=True))
            assert not ((cna_attrs_df[col] == True) & (cna_attrs_df['chrom'] != chrom_name)).any()
    for chrom_name in get_karyotype_estimation_params()['ordered_chrom_names_for_karyotype_estimation']:
        for del_or_dup in ['del', 'dup']:
            col = f'is_{del_or_dup}_{chrom_name}'
            cna_attrs_df[col] = (
                cna_attrs_df['donor_general_manual_comment'].isna() & (cna_attrs_df['whole_chrom'] == True) & 
                cna_attrs_df[f'is_{del_or_dup}'] & (cna_attrs_df['chrom'] == chrom_name)
            )

    cna_attrs_df = generic_utils.merge_preserving_df1_index_and_row_order(
        cna_attrs_df,
        var_df[['end_pos_in_chrom']].reset_index().rename(columns={'index': 'first_gene', 'end_pos_in_chrom': 'first_gene_end_pos_in_chrom'}),
        how='left',
    )
    assert (cna_attrs_df['first_gene'].isna() == cna_attrs_df['first_gene_end_pos_in_chrom'].isna()).all()
    
    cna_attrs_df = generic_utils.merge_preserving_df1_index_and_row_order(
        cna_attrs_df,
        var_df[['end_pos_in_chrom']].reset_index().rename(columns={'index': 'last_gene', 'end_pos_in_chrom': 'last_gene_end_pos_in_chrom'}),
        how='left',
    )
    assert (cna_attrs_df['last_gene'].isna() == cna_attrs_df['last_gene_end_pos_in_chrom'].isna()).all()


    return cna_attrs_df

def get_ipssr_cytogenetic_risk_group(donor_cna_attrs_df):
    # https://www.mds-foundation.org/ipss-r-calculator/
    # https://ashpublications.org/blood/article/120/12/2454/30571/Revised-International-Prognostic-Scoring-System
    donor_cna_attrs_df = donor_cna_attrs_df[donor_cna_attrs_df['donor_general_manual_comment'].isna()].copy()
    cna_count = len(donor_cna_attrs_df)
    cna_row = donor_cna_attrs_df.iloc[0]
    if (cna_count == 0):
        return 'Good'
    elif (cna_count == 1):
        if cna_row['is_del_11q']:
            return 'Very Good'
        if cna_row[['is_del_5q', 'is_del_12p', 'is_del_20q']].any():
            return 'Good'
        if cna_row[['is_del_7q', 'is_dup_8', 'is_dup_19']].any():
            return 'Intermediate'
        if cna_row[['is_del_7', 'is_del_3q']].any():
            return 'Poor'
        return 'Intermediate' 
    elif (cna_count == 2):
        if cna_row['is_del_5q'].any(axis=0):
            return 'Good'
        if cna_row[['is_del_7q', 'is_del_7']].any(axis=None):
            return 'Poor'
        return 'Intermediate' 
    elif (cna_count == 3):
        return 'Poor'
    elif (cna_count >= 4):
        return 'Very Poor'
    else:
        assert False


def get_ipssr_cytogenetic_risk_group_df(all_cna_attrs_df):
    all_cna_attrs_df = all_cna_attrs_df[all_cna_attrs_df['donor_general_manual_comment'].isna()].copy()
    flat_dicts = []
    for (donor_id, exp_name), curr_cna_df in all_cna_attrs_df.groupby(['donor_id', 'exp_name']):
        flat_dicts.append({
            'donor_id': donor_id,
            'exp_name': exp_name,
            'ipssr_cytogenetic_risk_group': get_ipssr_cytogenetic_risk_group(curr_cna_df),
        })
        
    return pd.DataFrame(flat_dicts)

def estimate_karyo(info_for_karyo, dataset_name, correct_by_gc, stop_after_mc_moving_median=True, forced_donor_ids=[], forced_mc_model_names=[], allow_donor_also_in_atlas=False, stop_after_writing_cna_hists_and_dfs=False, dpi=100, out_dir_path=None):
    chrom_pos_df = info_for_karyo['chrom_pos_df']
    var_df = info_for_karyo['var_df']
    gene_mask = info_for_karyo['gene_mask']
    unified_c_ad = info_for_karyo['unified_c_ad']
    donor_id_to_exp_sharing_donor_ids = info_for_karyo['donor_id_to_exp_sharing_donor_ids']
    num_of_gc_bins = info_for_karyo['num_of_gc_bins']
    gene_read_gc_df = info_for_karyo['gene_read_gc_df']
    bins = info_for_karyo['bins']
    bin_df = info_for_karyo['bin_df']
    # gene_names = info_for_karyo['gene_names']
    atlas_donor_ids = info_for_karyo['atlas_donor_ids']
    karyotype_cmap = info_for_karyo['karyotype_cmap']
    epsilon_for_log_projected_fold = info_for_karyo['epsilon_for_log_projected_fold']
    atlas_to_use = info_for_karyo['atlas_to_use']
    unified_mc_ad_var_names = info_for_karyo['unified_mc_ad_var_names']
    wanted_num_of_genes_in_bin = info_for_karyo['wanted_num_of_genes_in_bin']
    var_df_mask = info_for_karyo['var_df_mask']
    
    correct_by_gc_repr = '_corrected_by_gc' if correct_by_gc else ''

    mds_path_dict = mds_analysis_params.get_mds_path_dict(dataset_name)
    if out_dir_path is None:
        out_dir_path = mds_path_dict['karyotype_estimation_out_dir_path']

    examined_donor_exp_out_dir_path_flat_dicts = []
    # mc_model_names = [x for x in mds_analysis_params.get_all_existing_mc_model_names() if x.startswith('final_only_')]

    all_existing_mc_model_names = mds_analysis_params.get_all_existing_mc_model_names(mds_path_dict)
    mc_model_names = [x[1] for x in mds_analysis_params.get_single_donor_single_batch_mc_model_types_and_names(c_ad=unified_c_ad)]
    if forced_mc_model_names:
        mc_model_names = sorted(set(forced_mc_model_names) & set(mc_model_names))



    for mc_model_name in mc_model_names:
        if (all_existing_mc_model_names is not None) and (mc_model_name not in all_existing_mc_model_names):
            print(f'skipping {mc_model_name}, as it does not exist')
            continue

        mc_model_paths = mds_analysis_params.get_mc_model_paths(mc_model_name, dataset_name=dataset_name)
        if (not os.path.isfile(mc_model_paths['donor_df_before_computing_metacells_csv'])) or (not os.path.isfile(mc_model_paths['cells_with_metacell_attrs_ad'])):
            continue
        
        mc_model_minimal_donor_exp_info = mds_analysis_params.get_mc_model_minimal_donor_exp_info(mc_model_name, dataset_name=dataset_name)
        curr_donor_ids = mc_model_minimal_donor_exp_info['donor_ids']

        if len(curr_donor_ids) > 1:
            continue
        
        donor_id_ = mc_model_minimal_donor_exp_info['donor_id']
        if forced_donor_ids and (donor_id_ not in forced_donor_ids):
            continue
        # print(f'mc_model_minimal_donor_exp_info: {mc_model_minimal_donor_exp_info}')
        
        print(mc_model_name)

        exp_names_repr = mc_model_minimal_donor_exp_info.get('exp_names_repr')
        exp_names = mc_model_minimal_donor_exp_info.get('exp_names')
        exp_name = exp_names[0] if (len(exp_names) == 1) else None

        # print('exp_names_repr')
        # print(exp_names_repr)
        # print(exp_names)
        # continue

        if (donor_id_ in atlas_donor_ids) and (not allow_donor_also_in_atlas):
            raise RuntimeError('must not use projection of a donor to an atlas with the data of that donor. if the metacells with deletions/duplications are projected to themselves, we would miss the deletion/duplication')
        
        donor_id = donor_id_
        print(f'donor_id: {donor_id}')
        
        # if not mc_model_name.endswith(donor_id):
        #     continue
        print(f'mc_model_name: {mc_model_name}')

        # continue

        donor_exp_out_dir_path = mds_analysis_params.get_karyotype_donor_exp_out_dir_path(donor_id, exp_names_repr, corrected_by_gc=correct_by_gc, mds_path_dict=mds_path_dict)
        assert donor_id not in examined_donor_exp_out_dir_path_flat_dicts
            
        examined_donor_exp_out_dir_path_flat_dicts.append({
            'donor_id': donor_id,
            'exp_names_repr': exp_names_repr,
            'out_dir_path': donor_exp_out_dir_path,
        })


        
        # continue

        file_path_attr_prefix = ''
        if atlas_to_use:
            atlas_params = get_mds_params()[atlas_to_use]
            if 'file_path_attr_prefix' in atlas_params:
                file_path_attr_prefix = atlas_params['file_path_attr_prefix']
                
        
        curr_mc_ad_file_path = mc_model_paths[f'{file_path_attr_prefix}metacells_with_projection_ad']
        print(f'curr_mc_ad_file_path: {generic_utils.get_file_last_modification_time_str_repr(curr_mc_ad_file_path)}')
        curr_mc_ad = ad.read_h5ad(curr_mc_ad_file_path)
        
        if 'gene_names_file_path' in get_karyotype_estimation_params():
            curr_gene_names = generic_utils.read_text_file(get_karyotype_estimation_params()['gene_names_file_path']).split('\n')
            curr_gene_mask = curr_mc_ad.var_names.isin(curr_gene_names)
            curr_gene_names = curr_mc_ad.var_names[curr_gene_mask]
        else:
            assert (curr_mc_ad.var_names == unified_mc_ad_var_names).all()
            assert (curr_mc_ad.var_names == var_df.index).all()
            curr_gene_mask = gene_mask

        if 0:
            # TODO: remove this after calculating MCs again - here only because i changed my assigned cell types
            # curr_mc_ad.obs.drop(['state', 'state_color'], axis=1, inplace=True)
            mc_utils.add_state(
                curr_mc_ad, 
                cell_state_and_info_list=get_mds_params()['pb_cd34_enriched_cell_state_and_info_list'],
                cell_type_colors_csv_file_path=get_mds_params()['cell_type_colors_csv_file_path'],
                silently_overwrite_state_and_state_color_columns=True,
                default_state_column_name='projected_type',
                # verbose=False,
            )

        if correct_by_gc:
            plt.close('all')
            gc_correction_out_dir_path = os.path.join(donor_exp_out_dir_path, 'gc_correction')
            pathlib.Path(gc_correction_out_dir_path).mkdir(parents=True, exist_ok=True)

            correction_common_kwargs = dict(
                num_of_gc_bins=num_of_gc_bins,
                gene_read_gc_df=gene_read_gc_df,
                epsilon_to_add_for_log=epsilon_for_log_projected_fold,
                out_fig_dir_path=gc_correction_out_dir_path,
                use_median_for_gc_bin_correction=get_karyotype_estimation_params()['use_median_for_gc_bin_correction'],
            )
            # IIUC, for get_gene_gc_correction_df it is only important that the weight we give each metacell is according to the number of cells it has, because we assume the sequencing technology doesn't differentiate between UMIs of different cells (it just sees a big soup of DNA molecules). essentially, we want to give each UMI the same weight, and if we used the fraction for each metacell, it would give each metacell the same weight, thus giving more weight to UMIs belonging to small metacells...
            projected_umi_count_mat = mc_utils.get_umi_count_mat(curr_mc_ad, mc_x=curr_mc_ad.layers['projected_fraction'], gene_mask=curr_gene_mask)
            gc_correction_df = mc_utils.get_gene_gc_correction_df(
                mc_utils.get_umi_count_mat(curr_mc_ad, gene_mask=curr_gene_mask),
                projected_umi_count_mat,
                curr_mc_ad.var_names[curr_gene_mask],
                before_or_after_correction='before',
                **correction_common_kwargs,
            )

            corrected_x = curr_mc_ad.X[:, curr_gene_mask] / gc_correction_df['gc_bin_correction_factor'].to_numpy()[np.newaxis, :]
            
            # NOTE: if i understand correctly, this doesnt really matter because anyway later we always normalize each metacell to the total number of UMIs of the curr_gene_mask genes.
            # corrected_x *= (curr_mc_ad.X[:, curr_gene_mask].sum(axis=1) / corrected_x.sum(axis=1))[:, np.newaxis]
            
            gc_correction_df_after_correction = mc_utils.get_gene_gc_correction_df(
                mc_utils.get_umi_count_mat(curr_mc_ad, mc_x=corrected_x),
                projected_umi_count_mat, 
                curr_mc_ad.var_names[curr_gene_mask], 
                before_or_after_correction='after',
                **correction_common_kwargs,
            )
        else:
            corrected_x = curr_mc_ad.X[:, curr_gene_mask]

        plt.close('all')

        norm_observed_expr = corrected_x / corrected_x.sum(axis=1)
        if scipy.sparse.issparse(norm_observed_expr):
            norm_observed_expr = mc.ut.to_numpy_matrix(norm_observed_expr)
        # print(epsilon_for_log_projected_fold)
        log_norm_observed_expr = np.log2(norm_observed_expr + epsilon_for_log_projected_fold)
        projected_frac_before_renorm = curr_mc_ad.layers['projected_fraction'][:, curr_gene_mask]
        norm_projected_expr = projected_frac_before_renorm / projected_frac_before_renorm.sum(axis=1)
        if scipy.sparse.issparse(norm_projected_expr):
            norm_projected_expr = mc.ut.to_numpy_matrix(norm_projected_expr)
        log_norm_projected_expr = np.log2(norm_projected_expr + epsilon_for_log_projected_fold)
        projected_fold_mat = log_norm_observed_expr - log_norm_projected_expr
            
        print(f'projected_fold_mat.min(): {projected_fold_mat.min()}')
        print(f'projected_fold_mat.max(): {projected_fold_mat.max()}')
        
        if get_karyotype_estimation_params()['normalize_each_mc_projected_fold_row_by_median']:
            # bin_projected_fold_mat -= np.median(bin_projected_fold_mat, axis=1)[:, np.newaxis]
            projected_fold_mat -= np.median(projected_fold_mat, axis=1)[:, np.newaxis]
            # projected_fold_mat -= np.median(projected_fold_mat, axis=1)

        cluster_i_of_mcs, cluster_i_to_color, clustermap_obj = plot_multiple_gene_bin_clustermap(
            donor_id=donor_id,
            mc_ad=curr_mc_ad,
            projected_fold_mat=projected_fold_mat,
            bins=bins,
            bin_df=bin_df,
            var_df=var_df,
            gene_mask=curr_gene_mask,
            wanted_num_of_genes_in_bin=wanted_num_of_genes_in_bin,
            karyotype_cmap=karyotype_cmap,
            mc_model_name=mc_model_name,
            exp_name=exp_name,
        )
        # raise

        fig_file_name = None
        if exp_names:
            for exp_name in exp_names:
                fig_file_name = f'at_least_in_part_{exp_name}__{donor_id}_wanted_{wanted_num_of_genes_in_bin}_genes_in_bin_mc_karyotype_heatmap_{file_path_attr_prefix}{correct_by_gc_repr}.png'
                clustermap_obj.savefig(os.path.join(out_dir_path, fig_file_name))
        else:
            print('WARNING: exp_names is None or empty. make sure this makes sense.') 
            assert re.match(r'final_only_\d+', mc_model_name) 
            fig_file_name = f'{donor_id}_wanted_{wanted_num_of_genes_in_bin}_genes_in_bin_mc_karyotype_heatmap_{file_path_attr_prefix}{correct_by_gc_repr}.png'
        clustermap_obj.savefig(os.path.join(donor_exp_out_dir_path, fig_file_name), dpi=dpi)
        
        plt.close('all')
        
        # raise

        moving_window_size = get_karyotype_estimation_params()['moving_window_size']
        assert moving_window_size % 2 == 1
        num_of_bins_truncated_due_to_bin_size_on_each_edge = moving_window_size // 2

        curr_var_df = var_df.copy() # because we are going to edit it now.
        single_gene_bins, single_gene_bin_df = hg38_utils.get_bins_with_chroms_as_constraints_and_add_unified_end_pos_bin_i_to_var_df(
            chrom_pos_df,
            curr_var_df,
            # wanted_bin_size=wanted_bin_size,
            wanted_num_of_genes_in_bin=1, 
            var_df_mask=var_df_mask,
        )
        single_gene_bin_df['mean_unified_end_pos'] = single_gene_bin_df[['first_gene_unified_end_pos', 'last_gene_unified_end_pos']].mean(axis=1)
        single_gene_bin_projected_fold_mat, single_gene_xticklabels = get_bin_projected_fold_ordered_by_loc(
            projected_fold_mat=projected_fold_mat, var_df=curr_var_df, bins=single_gene_bins, gene_mask=curr_gene_mask,
            wanted_num_of_genes_in_bin=1, bin_df=single_gene_bin_df,
            truncate_x_tick_labels=True,
            # mean_or_median='mean',
        )

        valid_moving_median_mask = get_valid_moving_median_mask(single_gene_bin_df, num_of_bins_truncated_due_to_bin_size_on_each_edge)
        
        exp_sharing_donor_ids = donor_id_to_exp_sharing_donor_ids[donor_id] if donor_id_to_exp_sharing_donor_ids else [donor_id]
        print(f'exp_sharing_donor_ids: {exp_sharing_donor_ids}')
        print(pd.Series(cluster_i_of_mcs).value_counts().iloc[:5])

        fig, axes, first_mc_all_moving_median_cluster_df, cluster_i_to_clone_name, cluster_i_to_cna_names = plot_moving_median_and_get_cluster_i_to_clone_name(
            donor_id=donor_id,
            exp_sharing_donor_ids=exp_sharing_donor_ids,
            cluster_i_of_mcs=cluster_i_of_mcs,
            single_gene_bin_df=single_gene_bin_df,
            moving_window_size=moving_window_size,
            valid_moving_median_mask=valid_moving_median_mask,
            cluster_i_to_color=cluster_i_to_color,
            single_gene_bin_projected_fold_mat=single_gene_bin_projected_fold_mat,
            single_gene_xticklabels=single_gene_xticklabels,
            var_df=curr_var_df,
            cells_or_metacells='metacells',
            mean_or_median='median',
            # clusters_to_scatter_plot_indices=[13],
            # print_non_filtered_cluster_i_to_cna_median_cluster_median_projected_fold=True,
            # print_non_filtered_cluster_i_to_cna_names=True,
            # show_individual_genes_in_cluster_moving_median_plot=True,
        )
            
        fig.savefig(os.path.join(donor_exp_out_dir_path, f'moving_median_cluster_median_normalized_projected_fold.png'))
        if stop_after_mc_moving_median:
            continue

        plt.close('all')

        mc_hists_out_dir_path = os.path.join(donor_exp_out_dir_path, 'mc_hists')
        pathlib.Path(mc_hists_out_dir_path).mkdir(parents=True, exist_ok=True)

        cna_infos = mds_analysis_params.get_donor_id_to_cna_infos().get(donor_id, None)
        if cna_infos:
            cna_median_normalized_projected_fold_df = write_cna_hists_and_dfs(
                donor_id=donor_id,
                chrom_pos_df=chrom_pos_df,
                var_df=curr_var_df,
                single_gene_bin_df=single_gene_bin_df,
                single_gene_bin_projected_fold_mat=single_gene_bin_projected_fold_mat,
                out_dir_path=mc_hists_out_dir_path,
                mean_or_median='median',
                cells_or_metacells='metacells',
            )
            cna_median_normalized_projected_fold_df.to_csv(os.path.join(mc_hists_out_dir_path, 'cna_median_normalized_projected_fold_df.csv'), index=False)
            if stop_after_writing_cna_hists_and_dfs:
                continue

            calculate_cna_clusters_and_assign_clones_res = calculate_cna_clusters_and_assign_clones(
                cna_median_normalized_projected_fold_df=cna_median_normalized_projected_fold_df,
                curr_cells_or_metacells_ad=curr_mc_ad,
                donor_id=donor_id,
                karyotype_cmap=karyotype_cmap,
                donor_exp_out_dir_path=donor_exp_out_dir_path,
                cells_or_metacells='metacells',
                mean_or_median='median',
            )
            if calculate_cna_clusters_and_assign_clones_res is not None:
                cna_cluster_df, cna_cluster_i_of_mcs, cna_cluster_i_to_color, clone_name_of_mcs, clone_name_to_color = calculate_cna_clusters_and_assign_clones_res
            else:
                cna_cluster_df = cna_cluster_i_of_mcs = cna_cluster_i_to_color = clone_name_of_mcs = clone_name_to_color = None

            curr_c_ad_file_path = mc_model_paths[f'{file_path_attr_prefix}cells_with_metacell_attrs_ad']
            print(f'curr_c_ad_file_path: {generic_utils.get_file_last_modification_time_str_repr(curr_c_ad_file_path)}')
            curr_c_ad = ad.read_h5ad(curr_c_ad_file_path)

            orig_num_of_cells = curr_c_ad.n_obs
            curr_c_ad = mc.ut.slice(curr_c_ad, obs=~curr_c_ad.obs['metacell_i'].isna())
            log_num_of_non_excluded_umis = np.log2(mc.ut.to_numpy_vector(curr_c_ad.X.sum(axis=1)))

            num_of_outlier_cells = orig_num_of_cells - curr_c_ad.n_obs
            outlier_cell_percentage = num_of_outlier_cells / orig_num_of_cells * 100
            print(f'{num_of_outlier_cells} ({outlier_cell_percentage:.1f}%) outlier cells were discarded')
            curr_c_ad.obs['metacell_i'] = curr_c_ad.obs['metacell_i'].astype(int)

            # log_norm_projected_expr[curr_c_ad.obs['metacell_i'].astype(int), :].shape

            cells_corrected_x = curr_c_ad.X[:, curr_gene_mask]
            if correct_by_gc:
                cells_corrected_x /= gc_correction_df['gc_bin_correction_factor'].to_numpy()[np.newaxis, :]
            
            cells_corrected_x = mc.ut.to_numpy_matrix(cells_corrected_x) # to make things the same as with corrected_x, which is an np.ndarray, rather than cells_corrected_x, which was an np.matrix (until now). (e.g., np.newaxis will be required similarly to corrected_x)

            cells_norm_observed_expr = cells_corrected_x / cells_corrected_x.sum(axis=1)[:, np.newaxis]
            cells_log_norm_observed_expr = np.log2(cells_norm_observed_expr + epsilon_for_log_projected_fold)
            cells_projected_fold_mat = cells_log_norm_observed_expr - log_norm_projected_expr[curr_c_ad.obs['metacell_i'].astype(int), :]
            if get_karyotype_estimation_params()['normalize_each_cell_projected_fold_row_by_mean']:
                # plt.close('all')
                # sb.histplot(np.median(cells_projected_fold_mat, axis=1))
                # sb.histplot(np.mean(cells_projected_fold_mat, axis=1))
                # NOTE: as i later use the mean, i think it makes sense to also normalize by the mean.
                cells_projected_fold_mat -= np.mean(cells_projected_fold_mat, axis=1)[:, np.newaxis]
                # cells_projected_fold_mat -= np.median(cells_projected_fold_mat, axis=1)[:, np.newaxis]

            cells_single_gene_bin_projected_fold_mat, _ = get_bin_projected_fold_ordered_by_loc(
                projected_fold_mat=cells_projected_fold_mat, var_df=curr_var_df, bins=single_gene_bins, gene_mask=curr_gene_mask,
                wanted_num_of_genes_in_bin=1, bin_df=single_gene_bin_df,
                truncate_x_tick_labels=True,
                # mean_or_median='mean',
            )

            cells_donor_exp_out_dir_path = os.path.join(donor_exp_out_dir_path, 'cells')
            pathlib.Path(cells_donor_exp_out_dir_path).mkdir(parents=True, exist_ok=True)

            cell_attrs_df = get_cell_attrs_df_for_write_cna_hists_and_dfs(
                donor_id=donor_id,
                c_ad=curr_c_ad,
                mc_ad=curr_mc_ad,
                cluster_i_of_mcs=cluster_i_of_mcs,
                cluster_i_to_color=cluster_i_to_color,
                cna_median_normalized_projected_fold_df=cna_median_normalized_projected_fold_df,
            )

            cells_cna_mean_normalized_projected_fold_df = write_cna_hists_and_dfs(
                donor_id=donor_id,
                chrom_pos_df=chrom_pos_df,
                var_df=curr_var_df,
                single_gene_bin_df=single_gene_bin_df,
                single_gene_bin_projected_fold_mat=cells_single_gene_bin_projected_fold_mat,
                out_dir_path=cells_donor_exp_out_dir_path,
                mean_or_median='mean',
                cells_or_metacells='cells',
                cell_attrs_df=cell_attrs_df,
                c_ad=curr_c_ad,
                mc_ad=curr_mc_ad,
            )

            write_cell_projected_fold_df_csv(
                cna_median_normalized_projected_fold_df=cna_median_normalized_projected_fold_df,
                cells_cna_mean_normalized_projected_fold_df=cells_cna_mean_normalized_projected_fold_df,
                c_ad=curr_c_ad,
                mc_ad=curr_mc_ad,
                out_dir_path=donor_exp_out_dir_path,
            )
            
            mc_vs_cell_mean_out_dir_path = os.path.join(cells_donor_exp_out_dir_path, 'mc_vs_cell_mean_scatters')
            pathlib.Path(mc_vs_cell_mean_out_dir_path).mkdir(parents=True, exist_ok=True)
            plot_mc_vs_mean_cell_normalized_projected_fold_scatter_per_cna(
                cna_median_normalized_projected_fold_df=cna_median_normalized_projected_fold_df,
                cells_cna_mean_normalized_projected_fold_df=cells_cna_mean_normalized_projected_fold_df,
                c_ad=curr_c_ad,
                mc_ad=curr_mc_ad,
                out_dir_path=mc_vs_cell_mean_out_dir_path,
            )

            state_boxplots_out_dir_path = os.path.join(cells_donor_exp_out_dir_path, 'state_boxplots')
            pathlib.Path(state_boxplots_out_dir_path).mkdir(parents=True, exist_ok=True)
            plot_cell_cna_state_boxplots(
                cells_cna_mean_normalized_projected_fold_df=cells_cna_mean_normalized_projected_fold_df,
                c_ad=curr_c_ad,
                out_dir_path=state_boxplots_out_dir_path,
            )

            cnas_with_cell_threshold_names = [
                x['name'] for x in cna_infos 
                if {'max_mean_cell_cluster_mean_projected_fold_threshold', 'min_mean_cell_cluster_mean_projected_fold_threshold'} & set(x)
            ]
            if cnas_with_cell_threshold_names:
                cells_cna_cluster_df, cna_cluster_i_of_cells, cells_cna_cluster_i_to_color, clone_name_of_cells, _ = calculate_cna_clusters_and_assign_clones(
                    cna_median_normalized_projected_fold_df=cells_cna_mean_normalized_projected_fold_df,
                    curr_cells_or_metacells_ad=curr_c_ad,
                    donor_id=donor_id,
                    karyotype_cmap=karyotype_cmap,
                    donor_exp_out_dir_path=cells_donor_exp_out_dir_path,
                    cells_or_metacells='cells',
                    mean_or_median='mean',
                    cnas_to_use_for_clustering_names=cnas_with_cell_threshold_names,
                )

            # raise
            for curr_cluster_i_of_mcs, curr_cluster_i_to_color, curr_single_gene_bin_projected_fold_mat, cells_or_metacells, curr_mean_or_median, fig_title in (
                # (cluster_i_of_mcs, cluster_i_to_color, single_gene_bin_projected_fold_mat, 'metacells', 'median', 'orig_clusters'),
                (cna_cluster_i_of_mcs, cna_cluster_i_to_color, single_gene_bin_projected_fold_mat, 'metacells', 'median', 'cna_clusters'),
                # (clone_name_of_mcs, clone_name_to_color, single_gene_bin_projected_fold_mat, 'metacells', 'median', 'final_clones'),
                *(
                    (
                        (cna_cluster_i_of_cells, cells_cna_cluster_i_to_color, cells_single_gene_bin_projected_fold_mat, 'cells', 'mean', 'cell_cna_clusters'),
                        # (clone_name_of_cells, clone_name_to_color, cells_single_gene_bin_projected_fold_mat, 'cells', 'mean', 'final_cell_clones'),
                    ) if cnas_with_cell_threshold_names else ()
                )
            ):
                if (curr_cluster_i_of_mcs is None) or (curr_cluster_i_to_color is None):
                    continue
                # print(pd.Series(curr_cluster_i_of_mcs).value_counts().iloc[:15])
                fig, axes, curr_all_moving_median_cluster_df, cluster_i_to_clone_name, cluster_i_to_cna_names = plot_moving_median_and_get_cluster_i_to_clone_name(
                    donor_id=donor_id,
                    exp_sharing_donor_ids=exp_sharing_donor_ids,
                    cluster_i_of_mcs=curr_cluster_i_of_mcs,
                    single_gene_bin_df=single_gene_bin_df,
                    moving_window_size=moving_window_size,
                    valid_moving_median_mask=valid_moving_median_mask,
                    cluster_i_to_color=curr_cluster_i_to_color,
                    single_gene_bin_projected_fold_mat=curr_single_gene_bin_projected_fold_mat,
                    single_gene_xticklabels=single_gene_xticklabels,
                    var_df=curr_var_df,
                    cells_or_metacells=cells_or_metacells,
                    mean_or_median=curr_mean_or_median,
                    # clusters_to_scatter_plot_indices=[0,1],
                    # print_non_filtered_cluster_i_to_cna_median_cluster_median_projected_fold=True,
                    # print_non_filtered_cluster_i_to_cna_names=True,
                    # show_individual_genes_in_cluster_moving_median_plot=True,
                    # quantiles_to_plot_hlines=(0.01, 0.09),
                )
                fig.suptitle(fig_title)
                fig.savefig(os.path.join(donor_exp_out_dir_path, f'{fig_title}_moving_median_cluster_median_normalized_projected_fold.png'))
                plt.close('all')
        # break

def get_all_cna_attrs_df_and_donor_agg_cna_df(info_for_karyo):
    dfs = []

    donors_with_cnas_ids = set(mds_analysis_params.get_donor_id_to_cna_infos())

    for donor_id in sorted(donors_with_cnas_ids):
        print(donor_id)
        curr_df = get_cna_attrs_df(
            donor_id,
            chrom_pos_df=info_for_karyo['chrom_pos_df'],
            var_df=info_for_karyo['var_df'],
            gene_mask=info_for_karyo['gene_mask'],
        )
        if curr_df is not None:
            curr_df.insert(0, 'donor_id', donor_id)
            dfs.append(curr_df)
    all_cna_attrs_df = pd.concat(dfs, ignore_index=True)
    assert all_cna_attrs_df['name'].is_unique
    if 'contains_CSNK1A1' in list(all_cna_attrs_df):
        all_cna_attrs_df['contains_CSNK1A1'].fillna(False, inplace=True)

    donor_agg_cna_df = get_ipssr_cytogenetic_risk_group_df(all_cna_attrs_df)

    agg_col_to_is_mut_cols = {
        'has_del_17_or_del_17p': ['is_del_17', 'is_del_17p'],
        'has_del_7_or_del_7q': ['is_del_7', 'is_del_7q'],
        'has_del_5q': ['is_del_5q'],
    }

    for col, mut_cols in agg_col_to_is_mut_cols.items():
        curr_donor_exp_df = all_cna_attrs_df.loc[all_cna_attrs_df[mut_cols].any(axis=1), ['donor_id', 'exp_name']]
        donor_agg_cna_df[col] = generic_utils.get_mask_of_df1_rows_that_inner_join_will_keep(donor_agg_cna_df, curr_donor_exp_df)

    donor_agg_cna_df = generic_utils.merge_preserving_df1_index_and_row_order(
        donor_agg_cna_df,
        all_cna_attrs_df[all_cna_attrs_df['donor_general_manual_comment'].isna()][['donor_id', 'exp_name']].value_counts().reset_index(name='cna_count'),
    )
    return all_cna_attrs_df, donor_agg_cna_df
        



