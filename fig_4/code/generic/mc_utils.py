import os
import logging
import metacells as mc
import traceback
import numpy as np
import re
import math
import contextlib


# import scanpy as sc  # type: ignore
import sys
import random
import pathlib
import anndata as ad
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib
import scipy.cluster.hierarchy
import scipy.spatial.distance
import scipy.stats
import scipy.sparse
import sklearn.cluster
import collections.abc

from generic import generic_utils

# print(mc.__file__)
# print(mc.__version__)

CELLS_OBS_COLUMN_NAMES_ADDED_BY_METACELL_AND_DOWNSTREAM_ANALYSIS = (
    'grouped', 'total_umis', '__zeros_downsample_umis',
    'metacells_rare_gene_module', 'rare_metacell', 'metacell_i',
    'metacell_name', 'most_similar',
)

MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG = 1e-5
C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG = 2**(-13) # TODO: should it be lower? not sure at all...


def plot_cell_umi_total_count_dist_per_cell_attr(
        data_ad,
        cell_attr_name,
        violinplot_kwargs=dict(),
        y_upper_lim=None,
):
    umis_of_cells = mc.ut.get_o_numpy(data_ad, name='__x__', sum=True)
    num_of_cells = data_ad.n_obs
    umis_of_cells_with_cell_attr_df = pd.DataFrame({cell_attr_name: data_ad.obs[cell_attr_name], 'cell_umi_total_count': umis_of_cells})
    assert len(umis_of_cells_with_cell_attr_df) == num_of_cells

    fig1, ax1 = plt.subplots(figsize=(14, 5))
    sb.violinplot(
        data=umis_of_cells_with_cell_attr_df, 
        x=cell_attr_name, 
        y='cell_umi_total_count', 
        ax=ax1,
        width=1,
        **violinplot_kwargs,
    )
    if y_upper_lim:
        ax1.set_ylim(0, y_upper_lim)
    else:
        ax1.set_ylim(0, np.quantile(umis_of_cells, 0.99))
    ax1.set_xticklabels(
        ax1.get_xticklabels(),
        rotation=45, ha='right', 
        rotation_mode='anchor', 
        # fontsize='small',
    )
    return (fig1, ax1)



def get_total_umi_count_of_genes(data_ad):
    if 'total_umi_count' in data_ad.var:
        print("ATTENTION: using existing data_ad.var['total_umi_count']")
    else:
        print("computing and storing data_ad.var['total_umi_count']")
        total_umi_count_of_genes = mc.ut.get_v_numpy(data_ad, '__x__', sum=True)
        mc.ut.set_v_data(data_ad, 'total_umi_count', total_umi_count_of_genes)
    
    return data_ad.var['total_umi_count']

def get_total_umi_count_of_genes_for_downsampled_cells(data_ad):
    if 'downsampled_cells_total_umi_count' in data_ad.var:
        print("ATTENTION: using existing data_ad.var['downsampled_cells_total_umi_count']")
    else:
        print("computing and storing data_ad.var['total_umi_count']")
        total_umi_count_of_genes = mc.ut.get_v_numpy(data_ad, 'downsampled', sum=True)
        mc.ut.set_v_data(data_ad, 'downsampled_cells_total_umi_count', total_umi_count_of_genes)
    
    return data_ad.var['downsampled_cells_total_umi_count']

def get_umi_counts_normalized_by_total_cell_umis(data_ad, return_row_major=True):
    if 'normalized_X_row_major' in data_ad.layers:
        print('ATTENTION: using existing normalized_X_row_major')
    else:
        print('computing and storing normalized_X_row_major and normalized_X_col_major')
        mc.ut.set_vo_data(data_ad, 'normalized_X_row_major', mc.ut.fraction_by(data_ad.X, by='row'))
        mc.ut.set_vo_data(data_ad, 'normalized_X_col_major', mc.ut.to_layout(data_ad.layers['normalized_X_row_major'], layout='column_major'))
    
    if return_row_major:
        return data_ad.layers['normalized_X_row_major']
    return data_ad.layers['normalized_X_col_major']

def old_plot_gene_umi_normalized_count_dist_per_cell_attr(
        data_ad,
        gene_name,
        cell_attr_name,
        violinplot_kwargs=dict(),
        y_upper_lim=None,
):
    num_of_cells = data_ad.n_obs
    umi_counts_normalized_by_total_cell_umis = get_umi_counts_normalized_by_total_cell_umis(data_ad, return_row_major=False)
    gene_umi_normalized_counts = mc.ut.to_numpy_vector(umi_counts_normalized_by_total_cell_umis[:,data_ad.var_names == gene_name])

    y_axis_label = f'cell_{gene_name}_umi_normalized_count'
    umis_of_cells_with_cell_attr_df = pd.DataFrame({cell_attr_name: data_ad.obs[cell_attr_name], y_axis_label: gene_umi_normalized_counts})
    assert len(umis_of_cells_with_cell_attr_df) == num_of_cells

    fig1, ax1 = plt.subplots(figsize=(14, 5))
    # sb.boxplot(
    sb.violinplot(
        data=umis_of_cells_with_cell_attr_df, 
        x=cell_attr_name, 
        y=y_axis_label, 
        ax=ax1,
        width=1,
        **violinplot_kwargs,
    )
    if y_upper_lim:
        ax1.set_ylim(0, y_upper_lim)
    else:
        ax1.set_ylim(0, np.quantile(gene_umi_normalized_counts, 0.99))

    ax1.set_xticklabels(
        ax1.get_xticklabels(),
        rotation=45, ha='right', 
        rotation_mode='anchor', 
        # fontsize='small',
    )
    return (fig1, ax1)

def get_mw_for_cell_attr_vals_vs_other_cell_attr_vals(
        data_ad,
        value_of_cells,
        cell_attr_name,
        cell_attr_vals1,
        cell_attr_vals2,
):
    cell_attr_vals1 = set(cell_attr_vals1)
    cell_attr_vals2 = set(cell_attr_vals2)
    assert not (set(cell_attr_vals1) & set(cell_attr_vals2))

    cells1_mask = data_ad.obs[cell_attr_name].isin(cell_attr_vals1)
    cells2_mask = data_ad.obs[cell_attr_name].isin(cell_attr_vals2)

    cells1_values = value_of_cells[cells1_mask]
    cells2_values = value_of_cells[cells2_mask]

    u, pval = scipy.stats.mannwhitneyu(
        cells1_values,
        cells2_values,
        # use_continuity=use_continuity,
        alternative='two-sided',
    )
    return u, pval

# def plot_(
#         c_ad,
#         traj_pos_obs_col_name,
#         traj_pos_min_truncate_val,
#         traj_pos_max_truncate_val,
#         bin_edge_to_tick_label_func=lambda x: generic_utils.remove_redundant_trailing_zeros(f'{x:.1f}'),
# ):
    



    


def plot_heatmap_hist_per_cell_attr(
        attr_val_of_cells,
        val_of_cells,
        list_of_attr_val_to_color=None,
        ordered_attr_vals=None,
        row_cluster=False,
        normalize_bin_counts=False, # TODO: change default to True? sounds more sensible now (230706)
        log_bin_counts=False,
        min_count=None,
        bins=None,
        cell_mask=None,
        dilute_x_labels=True,
        dilute_y_labels=True,
        left_y_ticks=True,
        add_num_of_cells_row_color=False,
        obj_desc=None,
        list_of_heatmap_vline_x_and_kwargs=list(),
        cumulative=False,
        cmap_name='rocket_r', # hot_r
        x_axis_in_log=False,
        epsilon_for_log=None,
        cbar_kwargs={},
        num_of_cells_cbar_kwargs={},
        clustermap_kwargs={},
        vmin=None,
        vmax=None,
        # log_x_axis_ticks=None,
        x_axis_target_num_of_strs=20,
        num_of_spaces_to_append_to_y_ticklabels=0,
        min_truncate_val=None,
        max_truncate_val=None,
        bin_edge_xticks=True,
        bin_edge_to_tick_label_func=lambda x: generic_utils.remove_redundant_trailing_zeros(f'{x:.1f}'),
        also_return_attr_val_to_bin_counts=False,
        attr_val_to_yticklabel=None,
):
    attr_val_of_cells = np.array(attr_val_of_cells) # to avoid annoying category types and stuff like that
    val_of_cells = np.array(val_of_cells) # to avoid annoying category types and stuff like that

    if min_truncate_val is not None:
        val_of_cells = np.maximum(val_of_cells, min_truncate_val)
    if max_truncate_val is not None:
        val_of_cells = np.minimum(val_of_cells, max_truncate_val)

    if ordered_attr_vals is not None:
        assert not row_cluster
        assert set(attr_val_of_cells) <= set(ordered_attr_vals)
        assert len(ordered_attr_vals) == len(set(ordered_attr_vals))

    if cell_mask is not None:
        attr_val_of_cells = attr_val_of_cells[cell_mask]
        val_of_cells = val_of_cells[cell_mask]

    df = pd.DataFrame({
        'attr': attr_val_of_cells,
        'val': val_of_cells,
    })
    if x_axis_in_log:
        assert (df['val'] >= 0).all()
        if not (df['val'] > 0).all():
            assert epsilon_for_log is not None
            df['val'] += epsilon_for_log
        df['val'] = np.log2(df['val'])

    bins = generic_utils.get_bins(bins, df['val'])
    
    bin_means = (bins[:-1] + bins[1:]) / 2
    
    # print(f'bins: {bins}')
    # print(f'bin_means: {bin_means}')

    num_of_decimal_digits_for_bin_mean_round = max(0, math.ceil(-np.log10(bins[1] - bins[0])))
    
    if ordered_attr_vals is None:
        if row_cluster:
            ordered_attr_vals = list(df['attr'].unique()) # arbitrary order
        else:
            ordered_attr_vals = list(df.groupby('attr')['val'].mean().sort_values().index)

    attr_val_to_bin_counts = {}
    attr_vals_skipped_due_to_low_count = []
    for attr_val, curr_df in df.groupby('attr'):
        # print(f'attr_val: {attr_val}')
        num_of_vals = len(curr_df)
        if min_count is not None:
            if num_of_vals < min_count:
                attr_vals_skipped_due_to_low_count.append(attr_val)
                continue
        bin_counts = np.histogram(curr_df['val'], bins=bins)[0]

        if cumulative:
            bin_counts = bin_counts.cumsum()
            bin_counts += (curr_df['val'] < bins[0]).sum()
            assert (bin_counts[-1] + (curr_df['val'] > bins[-1]).sum()) == num_of_vals
        
        if normalize_bin_counts:
            bin_counts = bin_counts / num_of_vals
        # print(bin_counts[-1])
        num_of_dim_of_bin_counts = len(bin_counts.shape)
        if num_of_dim_of_bin_counts != 1:
            # return bin_counts
            raise RuntimeError(f'num_of_dim_of_bin_counts: {num_of_dim_of_bin_counts}')
        attr_val_to_bin_counts[attr_val] = bin_counts


    if attr_vals_skipped_due_to_low_count:
        print(f'ATTENTION: attr_vals_skipped_due_to_low_count: {attr_vals_skipped_due_to_low_count}')
    
    ordered_attr_vals = [x for x in ordered_attr_vals if x in attr_val_to_bin_counts]


    mat = np.vstack([
        attr_val_to_bin_counts[attr_val]
        for attr_val in ordered_attr_vals
    ])
    if log_bin_counts:
        if mat.min() > 0:
            epsilon_for_log = 0
        else:
            epsilon_for_log = mat[mat > 0].min() / 2
        mat = np.log2(mat + epsilon_for_log)
    # print(next(iter(attr_val_to_bin_counts.values())))
    # raise
    # print(mat)

    if list_of_attr_val_to_color is None:
        row_colors = None
    else:
        row_colors = [
            [attr_val_to_color[x] for x in ordered_attr_vals]
            for attr_val_to_color in list_of_attr_val_to_color
        ]
    if add_num_of_cells_row_color:
        attr_val_to_count = df['attr'].value_counts().to_dict()
        
        cell_count_cmap_name = 'binary'
        cell_count_colors, cell_count_color_norm = generic_utils.get_num_colors(
            [attr_val_to_count[x] for x in ordered_attr_vals], colormap=cell_count_cmap_name, also_return_norm=True)
        if row_colors is None:
            row_colors = [cell_count_colors]
        else:
            row_colors.append(cell_count_colors)

    yticklabels = [attr_val_to_yticklabel[x] for x in ordered_attr_vals] if attr_val_to_yticklabel else ordered_attr_vals
    if dilute_y_labels and (len(yticklabels) > 50):
        dilution_factor = math.ceil(len(yticklabels) / 50)
        yticklabels = generic_utils.dilute_str_list(yticklabels, dilution_factor)

    if num_of_spaces_to_append_to_y_ticklabels > 0:
        yticklabels = [x + ' ' * num_of_spaces_to_append_to_y_ticklabels for x in yticklabels]

    if bin_edge_xticks:
        xticklabels = True
    else:
        xticklabels = [generic_utils.remove_redundant_trailing_zeros(str(np.round(x, num_of_decimal_digits_for_bin_mean_round))) for x in bin_means]
        # print(xticklabels)
        # raise
        if dilute_x_labels and (len(xticklabels) > x_axis_target_num_of_strs):
            # dilution_factor = math.ceil(len(xticklabels) / 20)
            # xticklabels = generic_utils.dilute_str_list(xticklabels, dilution_factor)
            xticklabels = generic_utils.dilute_str_list(xticklabels, target_num_of_strs=x_axis_target_num_of_strs)
        # else:
        #     xticklabels = True

    cbar_label = 'fraction' if normalize_bin_counts else 'count'
    if obj_desc:
        cbar_label = f'fraction of {obj_desc}s' if normalize_bin_counts else f'#{obj_desc}s'
    if log_bin_counts:
        cbar_label = f'log2({cbar_label})'

    # print(f'mat.shape: {mat.shape}')
    # raise
    # if x_axis_in_log:
    #     assert not dilute_x_labels # i think?
    #     # assert not bin_edge_xticks
    #     mat = np.hstack((np.full((mat.shape[0], 1), np.nan), mat)) # what the hell is this doing?
    #     # xticklabels = [''] + xticklabels

    curr_clustermap_kwargs = dict(
        dendrogram_ratio=(0.07, 0.15),
    )
    for k,v in clustermap_kwargs.items():
        curr_clustermap_kwargs[k] = v
    
    if vmin is None:
        vmin = np.nanmin(mat)
    if vmax is None:
        vmax = np.nanmax(mat)
    clustermap_obj = sb.clustermap(
        mat,
        cmap=cmap_name,
        row_cluster=row_cluster,
        col_cluster=False,
        row_colors=row_colors,
        # cbar_kws=None,
        yticklabels=yticklabels,
        xticklabels=xticklabels,
        vmin=vmin,
        vmax=vmax,
        # cbar_pos=None, # if we do this, then clustermap_obj.ax_cbar gives another ax.
        **curr_clustermap_kwargs,
    )
    heatmap_ax = clustermap_obj.ax_heatmap
    if left_y_ticks:
        heatmap_ax.yaxis.set_ticks_position("left")
    
    generic_utils.make_all_spines_and_x_and_y_axes_invisible(clustermap_obj.ax_cbar)
    clustermap_obj.ax_cbar.clear()

    curr_cbar_kwargs = dict(
        # location='top',
        # anchor=(0.5, 0.9),
        orientation='horizontal',
        # ax=clustermap_obj.cax,
        # ax=clustermap_obj.fig.axes[1],
        aspect=15,
        fraction=0.5, # no idea why this works, but it does. maybe because the default location is 'bottom', and then this determines the fraction of the ax which is considered as the bottom.
        shrink=0.7,
    )
    for k,v in cbar_kwargs.items():
        curr_cbar_kwargs[k] = v
    clustermap_obj.fig.colorbar(
        matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap_name), 
        ax=clustermap_obj.ax_col_dendrogram,
        label=cbar_label,
        **curr_cbar_kwargs,
    )

    if add_num_of_cells_row_color:
        curr_num_of_cells_cbar_kwargs = dict(
            fraction=0.9, # no idea why this works, but it does.
            shrink=0.8,
            aspect=15,
            # location='left',
            # anchor=(0.5, 0.9),
            # orientation='horizontal',
            # pad=0.7,
            # panchor=(0.5, 0.4),
        )
        for k,v in num_of_cells_cbar_kwargs.items():
            curr_num_of_cells_cbar_kwargs[k] = v
        clustermap_obj.fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=cell_count_color_norm, cmap=cell_count_cmap_name), 
            ax=clustermap_obj.ax_cbar,
            label=f'#{obj_desc}s',
            **curr_num_of_cells_cbar_kwargs,
        )

    real_x_and_heatmap_ax_x_of_ticks_with_labels = [
        (float(x.get_text()), x.get_position()[0]) for x in heatmap_ax.get_xticklabels() if x.get_text() != '']
    
    for vline_x, heatmap_vline_kwargs in list_of_heatmap_vline_x_and_kwargs:
        # assert np.isclose(vline_x, bins).any(), f'vline_x must be equal to one of the bin edges'
        vline_heatmap_ax_x = np.interp(
            vline_x, [x[0] for x in real_x_and_heatmap_ax_x_of_ticks_with_labels], [x[1] for x in real_x_and_heatmap_ax_x_of_ticks_with_labels])
        heatmap_ax.axvline(vline_heatmap_ax_x, **heatmap_vline_kwargs)

    if bin_edge_xticks:
        generic_utils.set_heatmap_bin_edge_xticks(
            heatmap_ax, bins,
            bin_edge_to_tick_label_func=bin_edge_to_tick_label_func,
            dilute_x_labels=dilute_x_labels,
            x_axis_target_num_of_strs=x_axis_target_num_of_strs,
            min_truncate_val=min_truncate_val,
            max_truncate_val=max_truncate_val,
        )

    # if x_axis_in_log:
    #     # assert not bin_edge_xticks
    #     heatmap_ax.set_xscale('log', base=2)
    #     heatmap_ax.set_xlim(1, heatmap_ax.get_xlim()[1])

    #     if 0:
    #         real_x_of_min_heatmap_ax_x_tick = real_x_and_heatmap_ax_x_of_ticks_with_labels[0][0]
    #         real_x_of_max_heatmap_ax_x_tick = real_x_and_heatmap_ax_x_of_ticks_with_labels[-1][0]

    #         new_ticks = [x for x in log_x_axis_ticks if real_x_of_min_heatmap_ax_x_tick <= x <= real_x_of_max_heatmap_ax_x_tick]
    #         if bins_are_consecutive_integers: # because in this case, other ticks are meaningless.
    #             assert all(np.isclose(x, bin_means).any() for x in new_ticks)
    #         real_x_and_heatmap_ax_x_of_new_ticks = list(zip(
    #             new_ticks, 
    #             np.interp(
    #                 new_ticks, 
    #                 [x[0] for x in real_x_and_heatmap_ax_x_of_ticks_with_labels], 
    #                 [x[1] for x in real_x_and_heatmap_ax_x_of_ticks_with_labels],
    #             ),
    #         ))

    #         heatmap_ax.set_xticks([x[1] for x in real_x_and_heatmap_ax_x_of_new_ticks])
    #         heatmap_ax.get_xaxis().set_major_formatter(matplotlib.ticker.FixedFormatter(
    #             [generic_utils.remove_redundant_trailing_zeros(str(x[0])) for x in real_x_and_heatmap_ax_x_of_new_ticks]))

    # if row_colors is not None:
    #     clustermap_obj.fig.tight_layout() # this creates distance between the row colors and the matrix...
    if also_return_attr_val_to_bin_counts:
        return (clustermap_obj, attr_val_to_bin_counts)
    return clustermap_obj


def plot_dist_per_cell_attr(
        attr_value_of_cells,
        value_of_cells,
        value_description,
        cell_attr_name,
        cell_attr_value=None,
        also_plot_for_all=False,
        y_log_scale=False,
        x_upper_lim=None,
        x_lower_lim=None,
        x_upper_lim_qauntile=0.99,
        # ignore_zeros=False, # NOTE: I think this is not ok, as zeros should not be ignored.
        cell_attr_val_prefix=None,
        ordered_cell_attr_vals_to_show=None,
        bins=None,
        histplot_kwargs=dict(
            # stat='density',
            # color='tab:blue',

            # element='step', fill=False,
        ),
        extra_histplot_kwargs=dict(),
        fig_size=(5,7),
        y_labelpad=None,
        return_also_cell_attr_vals=False,
        order_by_mean_val=True,
):
    assert len(attr_value_of_cells) == len(value_of_cells)
    if bins is not None:
        histplot_kwargs['bins'] = bins

    cell_attr_vals_df = pd.DataFrame({cell_attr_name: attr_value_of_cells, value_description: value_of_cells})

    if cell_attr_val_prefix:
        assert ordered_cell_attr_vals_to_show is None
        cell_attr_vals_df = cell_attr_vals_df.loc[cell_attr_vals_df[cell_attr_name].str.startswith(cell_attr_val_prefix)]
        assert not cell_attr_vals_df.empty
    if ordered_cell_attr_vals_to_show:
        assert cell_attr_val_prefix is None
        cell_attr_vals_df = cell_attr_vals_df.loc[cell_attr_vals_df[cell_attr_name].isin(ordered_cell_attr_vals_to_show)]
        assert not cell_attr_vals_df.empty


    if x_upper_lim is None:
        x_upper_lim = np.quantile(cell_attr_vals_df[value_description], x_upper_lim_qauntile)
    if x_lower_lim is None:
        x_lower_lim = cell_attr_vals_df[value_description].min()

    if cell_attr_value:
        fig1, ax1 = plt.subplots(figsize=(7, 5))
        ax_or_axes_to_return = ax1
        cell_attr_vals_df['cell_is_of_specified_attr_val'] = cell_attr_vals_df[cell_attr_name] == cell_attr_value

        for cell_is_of_specified_attr_val, curr_group_df in cell_attr_vals_df.groupby('cell_is_of_specified_attr_val'):
            sb.histplot(
                curr_group_df[value_description], 
                ax=ax1,
                label=(cell_attr_value if cell_is_of_specified_attr_val else 'rest') + f' ({len(curr_group_df[value_description])})',
                **histplot_kwargs,
                **extra_histplot_kwargs,
            )
        
        if y_log_scale:
            ax1.set_yscale('log')
            ax1.set_ylabel('log(Density)')
        else:
            ax1.set_ylabel('Density')
        ax1.set_yticks([])
        ax1.legend()
        
        ax1.set_xlim(x_lower_lim, x_upper_lim)
        
    else:
        num_of_cell_attr_values = cell_attr_vals_df[cell_attr_name].nunique()
        fig1, axes = plt.subplots(nrows=num_of_cell_attr_values + also_plot_for_all, sharex=True,figsize=fig_size)
        if not isinstance(axes, collections.abc.Iterable):
            axes = [axes]
        ax_or_axes_to_return = axes

        if also_plot_for_all:
            sb.histplot(
                cell_attr_vals_df[value_description], 
                ax=axes[0],
                label=f'all ({len(curr_group_df[value_description])})',
                **histplot_kwargs,
                **extra_histplot_kwargs,
            )
            if y_log_scale:
                ax1.set_yscale('log')
                ax1.set_ylabel('log(Density) all')
            else:
                ax1.set_ylabel('Density all')
            axes = axes[1:]

        if ordered_cell_attr_vals_to_show is None:
            ordered_cell_attr_vals_to_show = cell_attr_vals_df[cell_attr_name].unique()
            if order_by_mean_val:
                print(f'sorting {cell_attr_name} by mean {value_description}')
                sorted_means = cell_attr_vals_df.groupby(cell_attr_name)[value_description].mean().sort_values()
                # print(f'sorted_means:\n{sorted_means}')
                ordered_cell_attr_vals_to_show = list(sorted_means.index)
        for ax1, curr_cell_attr_val in zip(axes, ordered_cell_attr_vals_to_show):
            curr_group_df = cell_attr_vals_df.loc[cell_attr_vals_df[cell_attr_name] == curr_cell_attr_val]
            # print(f'curr_cell_attr_val: {curr_cell_attr_val}')
            sb.histplot(
                curr_group_df[value_description], 
                ax=ax1,
                label=curr_cell_attr_val + f' ({len(curr_group_df[value_description])})',
                **histplot_kwargs,
                **extra_histplot_kwargs,
            )
            # ax1.set_yticks([])
            if y_log_scale:
                ax1.set_yscale('log')
            ax1.set_ylabel(
                curr_cell_attr_val, rotation=0, 
                labelpad=y_labelpad,
            )
            ax1.yaxis.set_label_position("right")

            ax1.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                labelbottom=False
            ) 
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
        ax1.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=True,      # ticks along the bottom edge are off
            labelbottom=True
        )


            
        
        
        if x_upper_lim:
            ax1.set_xlim(x_lower_lim, x_upper_lim)
        else:
            ax1.set_xlim(x_lower_lim, np.quantile(value_of_cells, x_upper_lim_qauntile))
        

    
    if return_also_cell_attr_vals:
        return (fig1, ax_or_axes_to_return, ordered_cell_attr_vals_to_show)
    else:
        return (fig1, ax_or_axes_to_return)

def get_filtered_gene_names(data_ad, gene_names):
    # gene_names_missing_from_ad = set(gene_names) - set(data_ad.var_names)
    gene_names_missing_from_ad = [x for x in gene_names if x not in data_ad.var_names]

    if gene_names_missing_from_ad:
        print(f'WARNING: the following genes are missing from the ad: {gene_names_missing_from_ad}. continuing without them.')
        return [x for x in gene_names if x in data_ad.var_names]
    return gene_names


def get_gene_names(data_ad, gene_or_gene_names, skip_missing_genes=False):
    if isinstance(gene_or_gene_names, str):
        gene_names = [gene_or_gene_names]
    else:
        gene_names = gene_or_gene_names
    
    if skip_missing_genes:
        gene_names = get_filtered_gene_names(data_ad, gene_names)
    else:
        missing_gene_names = set(gene_names) - set(data_ad.var_names)
        assert not missing_gene_names, f'the following received gene names are missing from data_ad.var_names: {missing_gene_names}'
        if 0:
            # 231220: don't know i wrote this
            missing_genes_mask = data_ad.var_names.isin(list(missing_gene_names))
            assert not missing_genes_mask.any(), f'{np.array(gene_names)[missing_genes_mask]} not in data_ad.var_names'
    
    generic_utils.assert_iter_is_unique(gene_names)
    # assert len(gene_names) == len(set(gene_names)), str(sorted((pd.Series(gene_names).value_counts() > 1).loc[lambda x: x].index))
    return gene_names

def get_gene_indices(data_ad, gene_names, skip_missing_genes=False):
    gene_names = get_gene_names(data_ad, gene_names, skip_missing_genes=skip_missing_genes)
    var_names = list(data_ad.var_names)
    return [var_names.index(x) for x in gene_names]

def get_gene_mask_allowing_none_gene_names(data_ad, gene_names):
    if gene_names is not None:
        return get_gene_mask(data_ad, gene_names)
    return np.full(data_ad.n_vars, True)

def get_gene_mask(data_ad, gene_names, skip_missing_genes=False):
    gene_names = get_gene_names(data_ad, gene_names, skip_missing_genes=skip_missing_genes)
    return data_ad.var_names.isin(gene_names)

def get_gene_index(data_ad, gene_name):
    return list(data_ad.var_names).index(gene_name)

def plot_gene_umi_count_dist_per_cell_attr_for_downsampled_cells(
        data_ad,
        cell_attr_name,
        gene_names=None,
        gene_name_patterns=None,
        x_upper_lim=None,
        x_upper_lim_qauntile=0.99,
        also_plot_for_all=False,
        y_log_scale=False,
        min_gene_total_umi_count_for_downsampled_cells=1e5,
        # ignore_zeros=False, # NOTE: I think this is not ok, as zeros should not be ignored.
        extra_histplot_kwargs=dict(),
        cell_attr_val_prefix=None,
        cell_attr_vals=None,
        max_bin=None,
        **histplot_kwargs,
):
    all_genes_mask = mc.tl.find_named_genes(data_ad, names=gene_names, patterns=gene_name_patterns)
    all_gene_indices = np.where(all_genes_mask)[0]
    total_umi_count_of_genes = get_total_umi_count_of_genes_for_downsampled_cells(data_ad)
    
    figs_and_axes = []
    for gene_index in all_gene_indices:
        gene_name = data_ad.var_names[gene_index]
        gene_total_umi_count = total_umi_count_of_genes[gene_index]
        if min_gene_total_umi_count_for_downsampled_cells and gene_total_umi_count < min_gene_total_umi_count_for_downsampled_cells:
            print(f'not plotting for {gene_name} because its total_umi_count is {gene_total_umi_count} (min={min_gene_total_umi_count_for_downsampled_cells})')
            continue
        
        count_col_name = f'downsampled_cell_{gene_name}_count'
        count_of_cells = mc.ut.to_numpy_vector(data_ad.layers['downsampled'][:,gene_index])
        
        max_count = int(np.ceil(count_of_cells.max()))
        if max_bin is None:
            bins = list(range(max_count + 1))
        else:
            bins = list(range(min(max_bin, max_count + 1)))
            if bins[-1] < max_count:
                bins.append(max_count)
        fig1, ax1 = plot_dist_per_cell_attr(
            attr_value_of_cells=data_ad.obs[cell_attr_name],
            value_of_cells=count_of_cells,
            value_description=count_col_name,
            cell_attr_name=cell_attr_name,
            also_plot_for_all=also_plot_for_all,
            y_log_scale=y_log_scale,
            x_upper_lim=x_upper_lim,
            x_upper_lim_qauntile=x_upper_lim_qauntile,
            extra_histplot_kwargs=extra_histplot_kwargs,
            cell_attr_val_prefix=cell_attr_val_prefix,
            ordered_cell_attr_vals_to_show=cell_attr_vals,
            bins=bins,
            **histplot_kwargs,
        )
        figs_and_axes.append((fig1, ax1))


    return figs_and_axes

def sample_cells_and_clean_excluded_genes(
        input_file_path_raw_data,
        cell_fraction_to_sample,
        num_of_cells_to_sample,
        min_num_of_umis_per_gene,
        min_total_umis_per_cell,
        max_total_umis_per_cell,
        excluded_gene_names,
        excluded_gene_name_patterns,
        max_excluded_gene_umi_fraction_per_cell,
        intermediate_out_dir_path,
        output_file_path_clean_data,
        use_existing_clean_gene_data_file_if_exists=False,
        fix_var_names_and_obs_names=False,
):
    raw_data_file_name = os.path.basename(input_file_path_raw_data)
    
    pathlib.Path(intermediate_out_dir_path).mkdir(parents=True, exist_ok=True)
    
    log_file_path = os.path.join(intermediate_out_dir_path, 'mc_log.txt')
    log_file = open(log_file_path, 'w')
    try:
        # raise RuntimeError('err messageaaaa')
        mc.ut.setup_logger(level=logging.DEBUG, long_level_names=True, time=True, to=log_file)
    # except RuntimeError as err:
    except AssertionError as err:
        assert 'mcs/utilities/logging.py, line 255 in setup_logger' in str(traceback.extract_tb(sys.exc_info()[2]))


    fixed_raw_data_file_path = os.path.join(intermediate_out_dir_path, f'{raw_data_file_name}.fixed.h5ad')
    sampled_data_file_path = os.path.join(intermediate_out_dir_path, f'{raw_data_file_name}.sample.h5ad')
    data_after_cleaning_genes_file_path = os.path.join(intermediate_out_dir_path, f'{raw_data_file_name}.sample_after_cleaning_genes.h5ad')

    if use_existing_clean_gene_data_file_if_exists and os.path.isfile(data_after_cleaning_genes_file_path):
        sampled_data = ad.read_h5ad(data_after_cleaning_genes_file_path)

    else:
        raw_data = ad.read_h5ad(input_file_path_raw_data)
        mc.ut.set_name(raw_data, 'raw_UMIs')

        assert raw_data.obs_names[0] != 0
        assert raw_data.var_names[0] != 0

        orig_num_of_cells = raw_data.shape[0]
        print(f'orig_num_of_cells: {orig_num_of_cells}')
        # one must be None and one must not be None
        assert (cell_fraction_to_sample is None) ^ (num_of_cells_to_sample is None)
        if cell_fraction_to_sample is not None:
            assert 0 < cell_fraction_to_sample <= 1
            sample_size = orig_num_of_cells * cell_fraction_to_sample
        else:
            assert num_of_cells_to_sample is not None
            assert 1 <= num_of_cells_to_sample <= orig_num_of_cells
            sample_size = num_of_cells_to_sample
        print(f'sample_size: {sample_size}')
        
        

        np.random.seed(0)
        # https://anndata.readthedocs.io/en/latest/
        sampled_cell_indices = np.random.choice(raw_data.n_obs, sample_size, replace=False)
        sampled_data = ad.AnnData(
            raw_data[sampled_cell_indices,:].X,
            obs=raw_data[sampled_cell_indices,:].obs,
            var=raw_data.var,
        )

        sampled_data.write(sampled_data_file_path)

        print(f'excluded_gene_names: {excluded_gene_names}')
        print(f'excluded_gene_name_patterns: {excluded_gene_name_patterns}')
        
        # NOTE: this function assumes that <AnnDataArg>.var_names has proper gene names.
        mc.pl.analyze_clean_genes(
            sampled_data,
            excluded_gene_names=excluded_gene_names,
            excluded_gene_patterns=excluded_gene_name_patterns,
            properly_sampled_min_gene_total=min_num_of_umis_per_gene,
            random_seed=123456,
        )
        # mc.tl.find_properly_sampled_genes(sampled_data, min_gene_total=min_num_of_umis_per_gene)
        mc.pl.pick_clean_genes(sampled_data)

        sampled_data.write(data_after_cleaning_genes_file_path)

    final_excluded_gene_names = sampled_data.var_names[sampled_data.var['excluded_gene']]
    print(f'final_excluded_gene_names: {final_excluded_gene_names}')


    total_umis_of_cells = mc.ut.get_o_numpy(sampled_data, name='__x__', sum=True)

    fig1, ax1 = plt.subplots(figsize=(8, 5))

    plot = sb.histplot(total_umis_of_cells, ax=ax1)
    plot.set(xlabel='UMIs', ylabel='Density', yticks=[])
    if min_total_umis_per_cell is not None:
        ax1.axvline(x=min_total_umis_per_cell, color='red', linestyle=':')
    if max_total_umis_per_cell is not None:
        ax1.axvline(x=max_total_umis_per_cell, color='red', linestyle=':')
    # ax.set_xscale('log')
    ax1.set_xlim(0, 11000)

    if min_total_umis_per_cell is None or max_total_umis_per_cell is None:
        raise RuntimeError('specify min_total_umis_per_cell and max_total_umis_per_cell according to the distribution...')


    fig2, ax2 = plt.subplots(figsize=(8, 5))

    excluded_genes_data = mc.tl.filter_data(sampled_data, var_masks=['~clean_gene'])[0]
    excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
    excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells
    plot = sb.histplot(excluded_fraction_of_umis_of_cells, ax=ax2)
    plot.set(xlabel='excluded gene UMIs fraction', ylabel='Density', yticks=[])
    if max_excluded_gene_umi_fraction_per_cell is not None:
        ax2.axvline(x=max_excluded_gene_umi_fraction_per_cell, color='red', linestyle=':')
    # ax.set_xscale('log')

    if max_excluded_gene_umi_fraction_per_cell is None:
        raise RuntimeError('specify max_excluded_gene_umi_fraction_per_cell according to the distribution...')

    mc.pl.analyze_clean_cells(
        sampled_data,
        properly_sampled_min_cell_total=min_total_umis_per_cell,
        properly_sampled_max_cell_total=max_total_umis_per_cell,
        properly_sampled_max_excluded_genes_fraction=max_excluded_gene_umi_fraction_per_cell,
    )

    mc.pl.pick_clean_cells(sampled_data)

    clean_data = mc.pl.extract_clean_data(sampled_data)

    clean_data.write(output_file_path_clean_data)

    
def get_most_prevalent_genes_rho_rho_and_linkage(
        data_ad,
        approx_num_of_most_prevalent_genes,
        downsample_min_samples=750, # the default on mcs v0.9.0
        downsample_min_cell_quantile=0.05, # the default on mcs v0.9.0
        downsample_max_cell_quantile=0.5, # the default on mcs v0.9.0
):
    mc.tl.downsample_cells(
        data_ad,
        '__x__',
        downsample_min_samples=downsample_min_samples,
        downsample_min_cell_quantile=downsample_min_cell_quantile,
        downsample_max_cell_quantile=downsample_max_cell_quantile,
        random_seed=1234,
    )
    downsample_total_umis_of_genes = mc.ut.get_v_numpy(data_ad, name='downsampled', sum=True)
    downsample_total_umis = downsample_total_umis_of_genes.sum()
    downsample_normalized_total_umis_of_genes = downsample_total_umis_of_genes / downsample_total_umis

    sorted_downsample_normalized_total_umis_of_genes = sorted(downsample_normalized_total_umis_of_genes, reverse=True)

    top_20_genes_mask = downsample_normalized_total_umis_of_genes >= sorted_downsample_normalized_total_umis_of_genes[19]
    top_20_gene_names = list(data_ad.var_names[top_20_genes_mask])
    print(f'top_20_gene_names: {top_20_gene_names}')

    most_prevalent_gene_threshold = sorted_downsample_normalized_total_umis_of_genes[approx_num_of_most_prevalent_genes - 1]
    print(f'most_prevalent_gene_threshold: {most_prevalent_gene_threshold}')
    most_prevalent_genes_mask = downsample_normalized_total_umis_of_genes >= most_prevalent_gene_threshold

    downsampled_x_col_major = mc.ut.to_layout(data_ad.layers['downsampled'], layout='column_major')[:,most_prevalent_genes_mask]
    downsample_gene_rho = mc.ut.corrcoef(downsampled_x_col_major, per='column', reproducible=True)
    downsample_gene_rho_rho = mc.ut.corrcoef(downsample_gene_rho, per=None, reproducible=True)
    downsample_gene_rho_rho_linkage = scipy.cluster.hierarchy.ward(downsample_gene_rho_rho)
    return (most_prevalent_genes_mask, downsample_gene_rho_rho, downsample_gene_rho_rho_linkage) 

def get_umis_of_cells_for_given_gene_mask(data_ad, genes_mask):
    return mc.ut.to_numpy_vector(mc.ut.to_layout(mc.ut.get_vo_proper(data_ad, layout='column_major')[:,genes_mask], layout='row_major').sum(axis=1))

def get_umis_of_cells_for_given_gene_names(data_ad, gene_names):
    genes_mask = mc.tl.find_named_genes(data_ad, names=gene_names)
    return get_umis_of_cells_for_given_gene_mask(data_ad, genes_mask)

def add_gene_id_to_mc_ad(mc_ad, c_ad=None, feature_df=None):
    # assert 'gene_id' not in list(mc_ad.var) # 230702 wasn't true in 230702 metacell version (not sure about previous versions)
    if c_ad is not None:
        assert mc_ad.n_vars == c_ad.n_vars
        if 'gene_id' not in c_ad.var.columns:
            print('WARNING: c_ad.var does not have gene_id column. add_gene_id_to_mc_ad DOES NOTHING.')
            return
        mc_ad.var = generic_utils.merge_preserving_df1_index_and_row_order(
            mc_ad.var, c_ad.var['gene_id'], left_index=True, right_index=True)
    else:
        feature_df = feature_df.copy()
        duplicated_gene_names = sorted(feature_df['gene'][feature_df['gene'].duplicated()].to_numpy())
        assert len(duplicated_gene_names) == len(set(duplicated_gene_names))
        num_of_gene_names_with_multiple_gene_ids = len(duplicated_gene_names)
        print(f'discarding the following {num_of_gene_names_with_multiple_gene_ids} gene names as they have multiple gene ids: {duplicated_gene_names}')
        feature_df = feature_df[~feature_df['gene'].isin(duplicated_gene_names)]
        mc_ad.var = generic_utils.merge_preserving_df1_index_and_row_order(
            mc_ad.var, feature_df[['gene_id', 'gene']], left_index=True, right_on='gene', how='left')

def get_corr_df(
        mc_ad, gene_names=None, obs_column_name=None, vec_to_find_genes_corr_with=None,
        gene_mask=None, mc_mask=None, corrs_with_obs_columns=False, other_gene_min_expr_range=1,
        verbose=True, ret_df=True,
):
    if mc_mask is None:
        mc_mask = np.full(mc_ad.n_obs, True)
    if gene_mask is None:
        gene_mask = np.full(mc_ad.n_vars, True)

    if other_gene_min_expr_range is not None:
        expr_range_df = get_gene_expr_range_df(mc_ad, mc_mask)
        gene_mask &= mc_ad.var_names.isin(expr_range_df.loc[expr_range_df['exp_range'] >= other_gene_min_expr_range, 'gene'])
    
    if vec_to_find_genes_corr_with is None:
        if gene_names is not None:
            assert obs_column_name is None
            vec_to_find_genes_corr_with = get_genes_expr(mc_ad, gene_names, c_or_mc_mask=mc_mask)
        else:
            assert obs_column_name is not None
            assert gene_names is None
            vec_to_find_genes_corr_with = mc.ut.to_numpy_vector(mc_ad.obs.loc[mc_mask, obs_column_name])
    vec_to_find_genes_corr_with = vec_to_find_genes_corr_with.astype(np.float32)

    if verbose:
        print(f'pd.Series(vec_to_find_genes_corr_with).describe():\n{pd.Series(vec_to_find_genes_corr_with).describe()}')
    # print(f'vec_to_find_genes_corr_with.dtype: {vec_to_find_genes_corr_with.dtype}')

    if corrs_with_obs_columns:
        names = [x for x in list(mc_ad.obs) if ((mc_ad.obs[x].dtype.kind in 'biuf') and (mc_ad.obs[x].nunique() > 1))]
        mat = mc_ad.obs.loc[mc_mask, names].to_numpy().T.astype(np.float32)
    else:
        expr = mc_ad.layers['expr'][np.ix_(mc_mask, gene_mask)]
        names = list(mc_ad.var_names[gene_mask])
        mat = expr.T
    corrs = mc.ut.to_numpy_vector(mc.ut.cross_corrcoef_rows(
        mc.ut.to_layout(mc.ut.to_numpy_matrix(vec_to_find_genes_corr_with[np.newaxis,:]), layout='row_major'),
        mc.ut.to_layout(mat, layout='row_major'),
        reproducible=True,
    ))
    corrs = pd.Series(corrs, index=names).sort_values(ascending=False)
    if ret_df:
        return pd.DataFrame({'corr': corrs, 'gene': corrs.index})
    return corrs

def get_gene_gene_corr_mat_df(
        c_or_mc_ad, genes=None, c_or_mc_mask=None, 
        also_plot_clustermap_and_return_clustermap_obj_and_ordered_gene_names=False,
        max_obs_count=int(100e3),
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        **clustermap_kwargs,
):
    ad_is_mc_ad = is_mc_ad(c_or_mc_ad)
    
    if c_or_mc_mask is None:
        c_or_mc_mask = np.full(c_or_mc_ad.n_obs, True)
    
    obs_count = c_or_mc_mask.sum()
    assert obs_count <= max_obs_count, f'obs_count={obs_count} > max_obs_count={max_obs_count}. increase max_obs_count if you want to proceed.'
    print(f'obs_count: {obs_count}')

    gene_mask = np.full(c_or_mc_ad.n_vars, True) if (genes is None) else c_or_mc_ad.var_names.isin(genes)

    
    if ad_is_mc_ad:
        write_expr_and_expr_enrich(c_or_mc_ad)
        expr = c_or_mc_ad.layers['expr'][np.ix_(c_or_mc_mask, gene_mask)]
    else:
        write_downsampled_layer_if_not_exists(c_or_mc_ad)
        expr = np.log2(
            (mc.ut.to_numpy_matrix(c_or_mc_ad.layers['downsampled'][np.ix_(c_or_mc_mask, gene_mask)]) / c_or_mc_ad.uns['downsample_samples'])
            + c_epsilon_to_add_to_fraction_before_log
        )
    
    corr_mat = mc.ut.corrcoef(mc.ut.to_layout(expr, layout='column_major'), per='column', reproducible=True)
    gene_names = c_or_mc_ad.var_names[gene_mask]
    corr_mat_df = pd.DataFrame(corr_mat, index=gene_names, columns=gene_names)

    if also_plot_clustermap_and_return_clustermap_obj_and_ordered_gene_names:
        if 'selected_gene' in c_or_mc_ad.var.columns:
            selected_gene_names = sorted(c_or_mc_ad.var_names[c_or_mc_ad.var['selected_gene']])
            tick_labels = [x if x in selected_gene_names else f'**{x}' for x in gene_names]
        else:
            tick_labels = gene_names

        if 'vmin' not in clustermap_kwargs:
            clustermap_kwargs['vmin'] = -1
        if 'vmax' not in clustermap_kwargs:
            clustermap_kwargs['vmax'] = 1
        if 'cmap' not in clustermap_kwargs:
            clustermap_kwargs['cmap'] = 'bwr'
        if 'metric' not in clustermap_kwargs:
            clustermap_kwargs['metric'] = 'correlation'
        clustermap_obj = sb.clustermap(
            corr_mat_df, 
            xticklabels=tick_labels, 
            yticklabels=tick_labels,
            **clustermap_kwargs,
        )
        ordered_gene_indices = clustermap_obj.dendrogram_row.reordered_ind
        ordered_gene_names = list(corr_mat_df.index[ordered_gene_indices])
        return corr_mat_df, clustermap_obj, ordered_gene_names

    return corr_mat_df

def get_gene_modules(
        linkage_mat,
        clustering_threshold=None,
        clustering_criterion='distance',
):
    if clustering_criterion == 'inconsistent':
        # clustering_inconsistency_threshold = 1.05
        # clustering_inconsistency_threshold = 0.9
        clustering_inconsistency_threshold = clustering_threshold
        fig1, ax1 = plt.subplots(figsize=(8, 5))
        gene_inconsistency_mat = scipy.cluster.hierarchy.inconsistent(linkage_mat)
        plot = sb.histplot(gene_inconsistency_mat[:,3], ax=ax1)
        plot.set(
            xlabel='inconsistency coefficient', 
            ylabel='Density', 
            yticks=[],
        )
        if clustering_inconsistency_threshold is None:
            raise RuntimeError('specify clustering_inconsistency_threshold according to the distribution...')
        else:
            ax1.axvline(x=clustering_inconsistency_threshold, color='red', linestyle=':')

        module_of_genes = scipy.cluster.hierarchy.fcluster(
            Z=linkage_mat,
            t=clustering_inconsistency_threshold,
            criterion='inconsistent',
            R=gene_inconsistency_mat,
        )
    else:
        assert clustering_criterion == 'distance'
        fig1, ax1 = plt.subplots(figsize=(8, 5))
        linkage_distances = linkage_mat[:,2]
        filtered_linkage_distances_for_plot_only = linkage_distances[linkage_distances < 100]
        # clustering_distance_threshold = 9
        clustering_distance_threshold = clustering_threshold
        plot = sb.histplot(filtered_linkage_distances_for_plot_only, ax=ax1)
        plot.set(
            xlabel='cluster distance', 
            ylabel='Density', 
            yticks=[],
        )
        # ax1.set_xlim(0, 25)
        if clustering_distance_threshold is None:
            raise RuntimeError('specify clustering_distance_threshold according to the distribution...')
        else:
            ax1.axvline(x=clustering_distance_threshold, color='red', linestyle=':')

        module_of_genes = scipy.cluster.hierarchy.fcluster(
            Z=linkage_mat,
            t=clustering_distance_threshold,
            criterion='distance',
        )

    num_of_modules = len(set(module_of_genes)) - 1
    print(f'num_of_modules: {num_of_modules}')
    cluster_sizes = np.unique(module_of_genes, return_counts=True)[1]
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    plot = sb.histplot(cluster_sizes, ax=ax2, bins=list(range(cluster_sizes.max() + 1)))
    plot.set(
        xlabel='cluster size', 
        ylabel='Density', 
        yticks=[],
    )
    return module_of_genes

def get_and_plot_gene_modules(
        gene_names,
        module_of_genes,
        genes_of_interest,
        similarity_of_genes,
        genes_similarity_linkage,
        max_num_of_gene_modules_to_plot=20,
        also_plot_all_modules_together=True,
):
    ordered_gene_indices = scipy.cluster.hierarchy.leaves_list(genes_similarity_linkage)
    del genes_similarity_linkage # can't use it after i reorder
    
    # reorder genes according to position in dendogram
    assert not isinstance(gene_names, set)
    gene_names = np.array(gene_names)[ordered_gene_indices]
    
    module_of_genes = module_of_genes[ordered_gene_indices]
    similarity_of_genes = np.array(similarity_of_genes)[ordered_gene_indices,:][:,ordered_gene_indices]


    # for some reason, np.isin doesn't like sets.
    module_indices = sorted(set(module_of_genes[np.isin(gene_names, list(genes_of_interest))]) - {-1})
    # print(f'module_indices (of modules containing any gene of interest): {module_indices}')

    gene_modules = {}
    num_of_modules_plotted = 0
    all_modules_gene_mask = np.full(len(module_of_genes), False)
    for module_index in module_indices:
        module_genes_mask = module_of_genes == module_index
        all_modules_gene_mask |= module_genes_mask
        if isinstance(similarity_of_genes, pd.DataFrame):
            similarity_of_module = similarity_of_genes.loc[module_genes_mask, module_genes_mask]
            similarity_of_module = similarity_of_module.to_numpy()
        else:
            similarity_of_module = similarity_of_genes[np.ix_(module_genes_mask, module_genes_mask)]
        curr_module_gene_names = gene_names[module_genes_mask]
        curr_module_genes_of_interest_names = [x for x in curr_module_gene_names if x in genes_of_interest]
        mod_name = f'mod_{curr_module_genes_of_interest_names[0]}'
        curr_module_genes_of_interest_names = set(curr_module_genes_of_interest_names)
        assert mod_name not in gene_modules
        gene_modules[mod_name] = list(curr_module_gene_names)
        # print(f'list(curr_module_gene_names): {list(curr_module_gene_names)}')

        assert num_of_modules_plotted <= max_num_of_gene_modules_to_plot
        if num_of_modules_plotted < max_num_of_gene_modules_to_plot:
            gene_names_with_stars = [
                '(*) ' + name if name in curr_module_genes_of_interest_names else name
                for name in curr_module_gene_names
            ]
            # print(similarity_of_module)
            similarity_of_module_for_heatmap = pd.DataFrame(data=similarity_of_module, columns=gene_names_with_stars, index=gene_names_with_stars)
            fig, ax = plt.subplots(figsize=(8, 5))
            sb.heatmap(similarity_of_module_for_heatmap, vmin=0, vmax=1, xticklabels=True, yticklabels=True, ax=ax, cmap="YlGnBu")
            ax.set_title(f'Gene Module {module_index}')
            num_of_modules_plotted += 1
        else:
            print(f'not plotting because already plotted max_num_of_gene_modules_to_plot ({max_num_of_gene_modules_to_plot})')
            # return
    
    if also_plot_all_modules_together:
        all_gene_names = gene_names[all_modules_gene_mask]
        # if isinstance(similarity_of_genes, pd.DataFrame):
        #     similarity_of_module = similarity_of_genes.loc[module_genes_mask, module_genes_mask]
        # else:
        similarity_of_all = similarity_of_genes[np.ix_(all_modules_gene_mask, all_modules_gene_mask)]
        similarity_of_all_for_heatmap = pd.DataFrame(similarity_of_all, columns=all_gene_names, index=all_gene_names)
        fig, ax = plt.subplots(figsize=(8, 5))
        sb.heatmap(similarity_of_all_for_heatmap, vmin=0, vmax=1, xticklabels=True, yticklabels=True, ax=ax, cmap="YlGnBu")
        ax.set_title('all modules')


    gene_names_in_dendogram_order = gene_names
    return gene_modules, gene_names_in_dendogram_order, ax



def filter_and_plot_gene_modules(
        module_of_genes,
        data_ad=None,
        gene_module_indices=None,
        max_num_of_gene_modules_to_plot=20,
        genes_of_interest=set(),
        show_gene_modules_only_for_these_genes=None,
        correlation_percentile=65,
        min_correlation_percentile_val=0.5,
        min_module_size=5,
        max_module_size=100,
        return_after_done_plotting=False,
        similarity_of_genes=None,
        gene_names=None,
):
    if similarity_of_genes is None:
        similarity_of_genes = mc.ut.get_vv_frame(data_ad, 'var_similarity')
    if gene_names is None:
        gene_names = data_ad.var_names

    if isinstance(similarity_of_genes, pd.DataFrame):
        similarity_of_genes = similarity_of_genes.to_numpy()

    if gene_module_indices is None:
        gene_module_indices = list(range(max(module_of_genes) + 1))
    
    if show_gene_modules_only_for_these_genes:
        gene_module_indices = set(gene_module_indices) & set(
            module_of_genes[gene_names.isin(show_gene_modules_only_for_these_genes)])
        gene_module_indices = sorted(gene_module_indices)

    print(similarity_of_genes.shape)
    filtered_gene_modules = []
    num_of_modules_plotted = 0
    for gene_module in gene_module_indices:
        # print(f'\ngene_module: {gene_module}')
        module_genes_mask = module_of_genes == gene_module
        similarity_of_module = similarity_of_genes[np.ix_(module_genes_mask, module_genes_mask)]
        curr_module_gene_names = gene_names[module_genes_mask]
        num_of_genes_in_module = len(curr_module_gene_names)
        if (num_of_genes_in_module < min_module_size) or (num_of_genes_in_module > max_module_size):
            if num_of_genes_in_module > max_module_size:
                print(f'num_of_genes_in_module: {num_of_genes_in_module}')
            continue
        # print(type(module_genes_mask), module_genes_mask.dtype)
        # print(similarity_of_module.shape)
        # raise

        percentile_val = np.percentile(similarity_of_module, correlation_percentile, axis=None)
        print(f'Gene Module {gene_module} ({num_of_genes_in_module} genes): {correlation_percentile}th percentile: {percentile_val}')
        if percentile_val < min_correlation_percentile_val:
            continue

        curr_module_genes_of_interest_names = {x for x in curr_module_gene_names if x in genes_of_interest}
        # print(f'curr_module_genes_of_interest_names: {curr_module_genes_of_interest_names}')
        print(f'genes in module {gene_module}: {curr_module_gene_names}')
        if curr_module_genes_of_interest_names:
            print(f'suspect genes in module {gene_module}: {list(curr_module_genes_of_interest_names)}')
        
        # for suspect_gene_name in sorted(curr_module_genes_of_interest_names)[:3]:
        #     print(f'{suspect_gene_name} similarity describe: {pd.Series(list(similarity_of_module[suspect_gene_name])).describe()}')

        gene_names_with_stars = [
            '(*) ' + name if name in curr_module_genes_of_interest_names else name
            for name in curr_module_gene_names
        ]
        similarity_of_module_for_heatmap = pd.DataFrame(similarity_of_module, columns=gene_names_with_stars, index=gene_names_with_stars)
        assert (similarity_of_module_for_heatmap.columns == similarity_of_module_for_heatmap.index).all()

        assert num_of_modules_plotted <= max_num_of_gene_modules_to_plot
        if num_of_modules_plotted < max_num_of_gene_modules_to_plot:
            fig, ax = plt.subplots(figsize=(8, 5))
            sb.heatmap(similarity_of_module_for_heatmap, vmin=0, vmax=1, xticklabels=True, yticklabels=True, ax=ax, cmap="YlGnBu")
            ax.set_title(f'Gene Module {gene_module}')
            num_of_modules_plotted += 1
        elif return_after_done_plotting:
            return

        filtered_gene_modules.append(curr_module_gene_names)
    return filtered_gene_modules


def choose_lateral_genes(
        data_ad,
        suspect_gene_names,
        suspect_gene_name_patterns,
        use_existing_related_genes_similarity_and_related_genes_module_if_exist=False,
):
    raise RuntimeError('old, i think. i guess i shouldnt use it')
    suspect_genes_mask = mc.tl.find_named_genes(data_ad, names=suspect_gene_names, patterns=suspect_gene_name_patterns)
    suspect_gene_names = sorted(data_ad.var_names[suspect_genes_mask])
    print(f'suspect_gene_names: {suspect_gene_names}')

    if (
        (not use_existing_related_genes_similarity_and_related_genes_module_if_exist) or 
        ('varp' not in dir(data_ad)) or ('related_genes_similarity' not in data_ad.varp) or
        ('related_genes_module' not in data_ad.var)
    ):
        mc.pl.relate_genes(
            data_ad, 
            random_seed=123456,
            # min_genes_of_modules=20,
            # genes_similarity_method='pearson', # default is 'repeated_pearson', IIUC, and it makes sense, IIUC.
        )

    module_of_genes = data_ad.var['related_genes_module']
    suspect_gene_modules = {x for x in module_of_genes[suspect_genes_mask] if x != -1}
    print(f'len(suspect_gene_modules): {len(suspect_gene_modules)}')

    filter_and_plot_gene_modules(
        data_ad=data_ad,
        module_of_genes=module_of_genes,
        gene_module_indices=suspect_gene_modules,
        genes_of_interest=suspect_gene_names,
        correlation_percentile=80,
        min_correlation_percentile_val=0.2,
    )
    

def write_gene_relative_variance(mc_ad, epsilon_to_add_to_mc_gene_mean_expression):
    relative_var = mc_ad.X.var(axis=0) / (mc_ad.X.mean(axis=0) + epsilon_to_add_to_mc_gene_mean_expression)
    mc.ut.set_v_data(mc_ad, name='relative_var', data=relative_var)
    
def get_norm_expression(c_ad=None, ad_x=None, gene_mask=None, mask_genes_after_normalizing=False):
    # NOTE: makes no sense for mc_ad, as now X is norm expr.
    if ad_x is None:
        ad_x = c_ad.X
    else:
        assert c_ad is None
    
    if mask_genes_after_normalizing:
        sums = ad_x.sum(axis=1)
    
    if gene_mask is not None:
        ad_x = ad_x[:, gene_mask]
    
    if not mask_genes_after_normalizing:
        sums = ad_x.sum(axis=1)
        
    assert sums.max() >= 2 # i.e., makes sense to be here only if X is UMI count.
    if len(sums.shape) == 1: # 230309: for some reason, today for mc_ad it was true, and for c_ad it wasn't.
        sums = sums[:, np.newaxis]
    return ad_x / sums

def get_log_norm_expression(c_or_mc_ad=None, ad_x=None, gene_mask=None, epsilon_for_log=1e-5, mask_genes_after_normalizing=False):
    norm_expr = get_norm_expression(c_or_mc_ad, ad_x, gene_mask, mask_genes_after_normalizing=mask_genes_after_normalizing)
    if scipy.sparse.issparse(norm_expr):
        print(f'(get_log_norm_expression) converting to dense matrix (shape={norm_expr.shape})')
        norm_expr = mc.ut.to_numpy_matrix(norm_expr)
    return np.log2(norm_expr + epsilon_for_log)

def get_genes_umi_frac_by_gene_mask(c_or_mc_ad, gene_mask, mc_or_cell_mask=None, layer_name=None):
    if (layer_name is None) or (layer_name == 'X'):
        mat = c_or_mc_ad.X
    else:
        assert layer_name in {'downsampled'}, f'{layer_name} doesnt seem to make sense in this context (because why take the fraction?)'
        mat = c_or_mc_ad.layers[layer_name]
    
    if mc_or_cell_mask is None:
        mc_or_cell_mask = np.full(c_or_mc_ad.n_obs, True)
    mat = mat[mc_or_cell_mask, :]
    return mc.ut.to_numpy_vector(mat[:, gene_mask].sum(axis=1) / mat.sum(axis=1))

def get_genes_umi_frac(c_or_mc_ad, gene_names, mc_or_cell_mask=None, skip_missing_genes=False, layer_name=None):
    gene_mask = get_gene_mask(c_or_mc_ad, gene_names, skip_missing_genes=skip_missing_genes)
    # print(f'get_genes_umi_frac: {layer_name}')
    return get_genes_umi_frac_by_gene_mask(c_or_mc_ad, gene_mask, mc_or_cell_mask=mc_or_cell_mask, layer_name=layer_name)

def get_genes_umi_count_by_gene_mask(c_or_mc_ad, gene_mask, mc_or_cell_mask=None, layer_name=None):
    if (layer_name is None) or (layer_name == 'X'):
        mat = c_or_mc_ad.X
        assert not is_mc_ad(c_or_mc_ad)
    else:
        assert layer_name in {'downsampled'}, 'only X or downsampled seem sensible in this context'
        mat = c_or_mc_ad.layers[layer_name]
    
    if mc_or_cell_mask is None:
        mc_or_cell_mask = np.full(c_or_mc_ad.n_obs, True)
    mat = mat[mc_or_cell_mask, :]
    return mc.ut.to_numpy_vector(mat[:, gene_mask].sum(axis=1))

def get_genes_umi_count(c_or_mc_ad, gene_names, mc_or_cell_mask=None, skip_missing_genes=False, layer_name=None):
    gene_mask = get_gene_mask(c_or_mc_ad, gene_names, skip_missing_genes=skip_missing_genes)
    return get_genes_umi_count_by_gene_mask(c_or_mc_ad, gene_mask, mc_or_cell_mask=mc_or_cell_mask, layer_name=layer_name)

def write_downsampled_layer_if_not_exists(c_ad, random_seed=1234, overwrite_if_existing_downsampled_layer_is_weird=False):
    assert not is_mc_ad(c_ad)
    if overwrite_if_existing_downsampled_layer_is_weird and ('downsampled' in c_ad.layers):
        weird_existing = pd.Series(mc.ut.to_numpy_vector(c_ad.layers['downsampled'].sum(axis=1))).value_counts(normalize=True).iloc[0] < 0.9
        if weird_existing:
            print('write_downsampled_layer_if_not_exists: existing downsampled layer is weird. overwriting.')
            del c_ad.layers['downsampled']

    if 'downsampled' not in c_ad.layers:
        print('write_downsampled_layer_if_not_exists: calling mc.tl.downsample_cells()')
        mc.tl.downsample_cells(c_ad, random_seed=random_seed)

def write_umi_count_obs_col_if_not_exists(c_ad):
    if 'umi_count' not in c_ad.obs.columns:
        print('write_umi_count_obs_col_if_not_exists: writing umi_count obs column')
        c_ad.obs['umi_count'] = mc.ut.to_numpy_vector(c_ad.X.sum(axis=1))

def get_downsampled_c_mask(c_ad):
    write_umi_count_obs_col_if_not_exists(c_ad)
    return (c_ad.obs['umi_count'] >= c_ad.uns['downsample_samples']).to_numpy()

def get_genes_expr(
        c_or_mc_ad, gene_names, c_or_mc_mask=None, 
        mc_epsilon_to_add_to_fraction_before_log=MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        skip_missing_genes=False, 
        layer_name=None, ret_df=False, mc_obs_columns_to_add=['state'], return_umi_frac=False, return_umi_count=False,
):
    ad_is_mc = is_mc_ad(c_or_mc_ad)
    if ad_is_mc:
        epsilon_for_log = mc_epsilon_to_add_to_fraction_before_log
    else:
        epsilon_for_log = c_epsilon_to_add_to_fraction_before_log
    
    
    # if layer_name is None:
    #     if ad_is_mc:
    #         layer_name = 'expr' if 'expr' in c_or_mc_ad.layers else 'X'
    #     else:
    #         # layer_name = 'downsampled'
    #         layer_name = 'X'
    
    if layer_name == 'downsampled':
        write_downsampled_layer_if_not_exists(c_or_mc_ad)

    if c_or_mc_mask is None:
        c_or_mc_mask = np.full(c_or_mc_ad.n_obs, True)
    gene_mask = get_gene_mask(c_or_mc_ad, gene_names, skip_missing_genes=skip_missing_genes)
    # print('gene_mask.sum()')
    # print(gene_mask.sum())
    if (layer_name is None) and (gene_mask.sum() == 1) and ('expr' in c_or_mc_ad.layers):
        if return_umi_frac:
            raise NotImplementedError('sorry')
        print('the expr layer already exists. using it as get_genes_expr received a single gene')
        # raise
        expr_vec = mc.ut.to_numpy_vector(c_or_mc_ad.layers['expr'][np.ix_(c_or_mc_mask, gene_mask)])
    # raise
    elif return_umi_count:
        assert (layer_name in ('downsampled', 'X')) or (layer_name is None)
        assert not is_mc_ad(c_or_mc_ad)
        expr_vec = get_genes_umi_count(
            c_or_mc_ad, gene_names, mc_or_cell_mask=c_or_mc_mask, skip_missing_genes=skip_missing_genes, layer_name=layer_name,
        )
    else:
        umi_frac_vec = get_genes_umi_frac(
            c_or_mc_ad, gene_names, mc_or_cell_mask=c_or_mc_mask, skip_missing_genes=skip_missing_genes, layer_name=layer_name,
        )
        if return_umi_frac:
            expr_vec = umi_frac_vec
        else:
            expr_vec = np.log2(umi_frac_vec + epsilon_for_log)
    
    if ret_df:
        return pd.DataFrame({
            'name': c_or_mc_ad.obs_names[c_or_mc_mask], 
            'expr': expr_vec,
            **{
                x: c_or_mc_ad.obs.loc[c_or_mc_mask, x]
                for x in mc_obs_columns_to_add if x in c_or_mc_ad.obs.columns
            },
        }).sort_values('expr', ascending=False)
    return expr_vec


def get_genes_pooled_total_norm_expression(c_or_mc_ad, cell_or_mc_mask, gene_names, downsample_layer_name='downsampled'):
    if downsample_layer_name is None:
        pooled_x = c_or_mc_ad.X[cell_or_mc_mask, :].sum(axis=0)
    else:
        pooled_x = c_or_mc_ad.layers[downsample_layer_name][cell_or_mc_mask, :].sum(axis=0)
    pooled_x = mc.ut.to_numpy_vector(pooled_x)
    
    gene_mask = get_gene_mask(c_or_mc_ad, gene_names)
    return pooled_x[gene_mask].sum() / pooled_x.sum()

def get_genes_pooled_log_total_norm_expression(
        c_or_mc_ad, cell_or_mc_mask, gene_names, downsample_layer_name='downsampled', epsilon_for_log=1e-5,
):
    return np.log2(get_genes_pooled_total_norm_expression(
        c_or_mc_ad, cell_or_mc_mask, gene_names, downsample_layer_name=downsample_layer_name) + epsilon_for_log)

def get_gene_or_genes_log_norm_or_downsampled_expression_and_desc(
        c_or_mc_ad, gene_name_or_names, 
        mc_epsilon_to_add_to_fraction_before_log=MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG, 
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        jitter_sd=0, layer_name='downsampled',
        skip_missing_genes=False, mc_or_cell_mask=None, return_umi_frac=False,
):
    # if 'expr' in c_or_mc_ad.layers:
    #     mat_type = 'expr'
    #     mat = c_or_mc_ad.layers['expr']
    # elif layer_name is None:
    #     mat = c_or_mc_ad.X
    #     mat_type = 'X'
    # else:
    #     mat = c_or_mc_ad.layers[layer_name]
    #     # c_or_mc_ad is actually c_ad. let's verify that
    #     assert 'barcode' in list(c_or_mc_ad.obs)
    #     mat_type = 'downsample'

    
    # if mat_type in ('expr', 'X'):
    vals = get_genes_expr(
        c_or_mc_ad, gene_name_or_names, 
        mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log, 
        c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log, 
        skip_missing_genes=skip_missing_genes, 
        layer_name=layer_name, c_or_mc_mask=mc_or_cell_mask, return_umi_frac=return_umi_frac,
    )
    
    # else:
    #     assert mat_type == 'downsample'
    #     vals = mc.ut.to_numpy_vector(mat[:,get_gene_indices(c_or_mc_ad, gene_name_or_names, skip_missing_genes=skip_missing_genes)].sum(axis=1))
    
    if mc_or_cell_mask is None:
        mc_or_cell_mask = np.full(c_or_mc_ad.n_obs, True)
    
    if jitter_sd:
        jitter_vals = np.random.normal(0, jitter_sd, mc_or_cell_mask.sum())
        vals += jitter_vals
    return (vals, get_gene_desc(c_or_mc_ad, gene_name_or_names, skip_missing_genes=skip_missing_genes))


def get_diff_exp_df(
        c_or_mc_ad=None, c_or_mc_mask1=None, c_or_mc_mask2=None, layer_name=None, 
        max_num_of_obs=int(2e3), epsilon_for_log=1e-5, epsilon_for_c_log=1/16, use_geomean=True,
        expr1_df=None, expr2_df=None, allow_non_disjoint_masks=False, sample_mask_if_too_many_cells=True,
        add_min_expr_columns=False,
):
    arg2_name_and_list = None
    if isinstance(c_or_mc_mask2, list):
        arg2_name_and_list = ('c_or_mc_mask2', c_or_mc_mask2)
    if isinstance(expr2_df, list):
        arg2_name_and_list = ('expr2_df', expr2_df)
    if arg2_name_and_list:
        arg2_name, arg2_list = arg2_name_and_list
        df = None
        for i, curr_arg2 in enumerate(arg2_list):
            curr_df = get_diff_exp_df(
                c_or_mc_ad=c_or_mc_ad, c_or_mc_mask1=c_or_mc_mask1, layer_name=layer_name, 
                max_num_of_obs=max_num_of_obs, epsilon_for_log=epsilon_for_log, epsilon_for_c_log=epsilon_for_c_log, use_geomean=use_geomean,
                expr1_df=expr1_df, allow_non_disjoint_masks=allow_non_disjoint_masks, sample_mask_if_too_many_cells=sample_mask_if_too_many_cells,
                add_min_expr_columns=add_min_expr_columns, **{arg2_name: curr_arg2},
            )
            curr_df.rename(columns={x: f'{x}_{i+1}' for x in set(curr_df.columns) - {'expr1', 'gene'}}, inplace=True)
            if df is None:
                df = curr_df
            else:
                df = generic_utils.merge_preserving_df1_index_and_row_order(df, curr_df)
        log_ratio_cols = [x for x in df.columns if x.startswith('log_ratio')]
        expr2_cols = [x for x in df.columns if x.startswith('expr2_')]
        df['max_log_ratio'] = df[log_ratio_cols].max(axis=1)
        df['min_log_ratio'] = df[log_ratio_cols].min(axis=1)
        df['max_expr2'] = df[expr2_cols].max(axis=1)
        df['min_expr2'] = df[expr2_cols].min(axis=1)
        df['same_sign_for_all_log_ratios'] = np.sign(df['max_log_ratio']) == np.sign(df['min_log_ratio'])
        df.loc[df['same_sign_for_all_log_ratios'], 'same_sign_min_abs_log_ratio'] = df.loc[
            df['same_sign_for_all_log_ratios'], ['max_log_ratio', 'min_log_ratio']].abs().min(axis=1)
        
        df.insert(1, 'same_sign_min_abs_log_ratio', df.pop('same_sign_min_abs_log_ratio'))
        df.sort_values('same_sign_min_abs_log_ratio', ascending=False, inplace=True)
        return df



    if (c_or_mc_mask1 is not None) and (c_or_mc_mask2 is not None):
        if (c_or_mc_mask1 & c_or_mc_mask2).any():
            if allow_non_disjoint_masks:
                print('warning: masks are not disjoint')
            else:
                raise RuntimeError('masks must be disjoint')


    get_mean_log_total_norm_expression_common_kwargs = dict(
        c_or_mc_ad=c_or_mc_ad, layer_name=layer_name, 
        max_num_of_obs=max_num_of_obs,
        epsilon_for_log=epsilon_for_log,
        epsilon_for_c_log=epsilon_for_c_log,
        use_geomean=use_geomean,
    )
    if expr1_df is None:
        assert c_or_mc_mask1 is not None
        
        if sample_mask_if_too_many_cells and (c_or_mc_mask1.sum() > max_num_of_obs):
            print(f'downsampling mask1 to {max_num_of_obs}')
            c_or_mc_mask1 = generic_utils.sample_mask(c_or_mc_mask1, max_num_of_obs, allow_not_sampling_due_to_too_few_pos=True)
        
        cell_count1 = c_or_mc_mask1.sum()
        print(f'cell_count1: {cell_count1}')
        expr1 = get_mean_expr(c_or_mc_mask=c_or_mc_mask1, **get_mean_log_total_norm_expression_common_kwargs)
        expr1_df = pd.DataFrame({
            'gene': c_or_mc_ad.var_names,
            'expr1': expr1,
        })
    else:
        assert set(expr1_df.columns) >= {'gene', 'expr'}, str(expr1_df.columns)
        expr1_df = expr1_df[['gene', 'expr']].rename(columns={'expr': 'expr1'})

    if (c_or_mc_mask2 is None) and (expr2_df is None):
        print('cell_or_mc_mask2 is None. assigning cell_or_mc_mask2=~cell_or_mc_mask1')
        c_or_mc_mask2 = ~c_or_mc_mask1
    
    if expr2_df is None:
        assert c_or_mc_mask2 is not None

        if sample_mask_if_too_many_cells and (c_or_mc_mask2.sum() > max_num_of_obs):
            print(f'downsampling mask2 to {max_num_of_obs}')
            c_or_mc_mask2 = generic_utils.sample_mask(c_or_mc_mask2, max_num_of_obs, allow_not_sampling_due_to_too_few_pos=True)
        
        cell_count2 = c_or_mc_mask2.sum()
        print(f'cell_count2: {cell_count2}')
        expr2 = get_mean_expr(c_or_mc_mask=c_or_mc_mask2, **get_mean_log_total_norm_expression_common_kwargs)
        expr2_df = pd.DataFrame({
            'gene': c_or_mc_ad.var_names,
            'expr2': expr2,
        })
    else:
        assert set(expr2_df.columns) >= {'gene', 'expr'}, str(expr2_df.columns)
        expr2_df = expr2_df[['gene', 'expr']].rename(columns={'expr': 'expr2'})
    

    df = generic_utils.inner_join_preserving_df1_row_order(expr1_df, expr2_df)

    df['log_ratio'] = (df['expr1'] - df['expr2'])
    df['abs_log_ratio'] = df['log_ratio'].abs()
    df.sort_values('abs_log_ratio', ascending=False, inplace=True)

    if add_min_expr_columns:
        df['min_expr'] = df[['expr1', 'expr2']].min(axis=1)
        df['min_expr_is_zero'] = np.isclose(df['min_expr'], np.log2(1e-5))
    return df

def get_3_way_diff_exp_df(
        c_or_mc_ad=None, c_or_mc_mask=None, c_or_mc_mask1=None, c_or_mc_mask2=None, 
        expr_df=None, expr1_df=None, expr2_df=None, show_scatter=True,
        min_log_ratio_for_plot_gene_names=1,
        genes_to_plot_names_anyway=None,
        **get_diff_exp_df_kwargs,
):
    # if (c_or_mc_mask is not None):
    #     c_or_mc_mask = mc.ut.to_numpy_vector(c_or_mc_mask)
    # if (c_or_mc_mask1 is not None):
    #     c_or_mc_mask1 = mc.ut.to_numpy_vector(c_or_mc_mask1)
    # if (c_or_mc_mask2 is not None):
    #     c_or_mc_mask2 = mc.ut.to_numpy_vector(c_or_mc_mask2)

    if (c_or_mc_mask is not None) and (c_or_mc_mask1 is not None):
        df1 = get_diff_exp_df(c_or_mc_ad, c_or_mc_mask, c_or_mc_mask1, **get_diff_exp_df_kwargs)
    elif (c_or_mc_mask is not None):
        df1 = get_diff_exp_df(c_or_mc_ad, c_or_mc_mask, expr2_df=expr1_df, **get_diff_exp_df_kwargs)
    elif (c_or_mc_mask1 is not None):
        df1 = get_diff_exp_df(c_or_mc_ad, expr1_df=expr_df, c_or_mc_mask2=c_or_mc_mask1, **get_diff_exp_df_kwargs)
    else:
        df1 = get_diff_exp_df(c_or_mc_ad, expr1_df=expr_df, expr2_df=expr1_df, **get_diff_exp_df_kwargs)
    
    if (c_or_mc_mask is not None) and (c_or_mc_mask2 is not None):
        df2 = get_diff_exp_df(c_or_mc_ad, c_or_mc_mask, c_or_mc_mask2, **get_diff_exp_df_kwargs)
    elif (c_or_mc_mask is not None):
        df2 = get_diff_exp_df(c_or_mc_ad, c_or_mc_mask, expr2_df=expr2_df, **get_diff_exp_df_kwargs)
    elif (c_or_mc_mask2 is not None):
        df2 = get_diff_exp_df(c_or_mc_ad, expr1_df=expr_df, c_or_mc_mask2=c_or_mc_mask2, **get_diff_exp_df_kwargs)
    else:
        df2 = get_diff_exp_df(c_or_mc_ad, expr1_df=expr_df, expr2_df=expr2_df, **get_diff_exp_df_kwargs)

    df1.rename(columns={'expr1': 'expr', 'expr2': 'expr1'}, inplace=True)
    df2.rename(columns={'expr1': 'expr'}, inplace=True)

    df = generic_utils.merge_preserving_df1_index_and_row_order(df1, df2, on=['gene', 'expr'], suffixes=('1', '2'))

    if show_scatter:
        fig, ax = plt.subplots()
        sb.scatterplot(
            data=df,
            x='log_ratio1',
            y='log_ratio2',
            s=8,
            ax=ax,
        )
        ax.axhline(0, color='grey', alpha=0.5)
        ax.axvline(0, color='grey', alpha=0.5)
        
        genes_to_plot_names_mask = (
            ((df['log_ratio1'] > min_log_ratio_for_plot_gene_names) & (df['log_ratio2'] > min_log_ratio_for_plot_gene_names))
            | ((df['log_ratio1'] < -min_log_ratio_for_plot_gene_names) & (df['log_ratio2'] < -min_log_ratio_for_plot_gene_names))
            | df['gene'].isin(['XIST', 'RPS4Y1'])
        )
        if genes_to_plot_names_anyway is not None:
            genes_to_plot_names_mask |= df['gene'].isin(genes_to_plot_names_anyway)

        for i, row in df[genes_to_plot_names_mask].iterrows():
            ax.text(row['log_ratio1'], row['log_ratio2'], row['gene'])


    return df

def get_pairwise_diff_exp_df(
        c_or_mc_ad=None, c_or_mc_mask_pairs=None, expr_df_pairs=None,
        show_scatter=True,
        min_abs_log_ratio_for_plot_gene_names=1,
        genes_to_plot_names_anyway=None,
        **get_diff_exp_df_kwargs,
):
    dfs = []
    if c_or_mc_mask_pairs:
        assert expr_df_pairs is None
        for c_or_mc_mask1, c_or_mc_mask2 in c_or_mc_mask_pairs:
            dfs.append(get_diff_exp_df(c_or_mc_ad, c_or_mc_mask1, c_or_mc_mask2, **get_diff_exp_df_kwargs))
    elif expr_df_pairs:
        assert c_or_mc_mask_pairs is None
        for expr1_df, expr2_df in expr_df_pairs:
            dfs.append(get_diff_exp_df(expr1_df=expr1_df, expr2_df=expr2_df, **get_diff_exp_df_kwargs))
    else:
        assert False, 'either c_or_mc_mask_pairs or expr_df_pairs must be provided'
    
    all_genes = set()
    for df in dfs:
        all_genes |= set(df['gene'])

    df = pd.DataFrame({'gene': sorted(all_genes)})
    for i, curr_df in enumerate(dfs):
        curr_df.rename(columns={
            'expr1': f'expr1{i + 1}', 'expr2': f'expr2{i + 1}', 
            'log_ratio': f'log_ratio{i + 1}', 'abs_log_ratio': f'abs_log_ratio{i + 1}',
        }, inplace=True)
        df = generic_utils.merge_preserving_df1_index_and_row_order(df, curr_df, how='left')
    if show_scatter:
        fig, ax = plt.subplots()
        sb.scatterplot(
            data=df,
            x=f'log_ratio1',
            y=f'log_ratio2',
            s=8,
            ax=ax,
        )
        generic_utils.plot_y_equals_x_line_on_ax(ax)
        # ax.axhline(0, color='grey', alpha=0.5)
        # ax.axvline(0, color='grey', alpha=0.5)
        
        genes_to_plot_names_mask = (
            (
                (df[f'abs_log_ratio1'] > min_abs_log_ratio_for_plot_gene_names)
                & (df[f'abs_log_ratio2'] > min_abs_log_ratio_for_plot_gene_names)
            )
            | df['gene'].isin(['XIST', 'RPS4Y1'])
        )
        if genes_to_plot_names_anyway is not None:
            genes_to_plot_names_mask |= df['gene'].isin(genes_to_plot_names_anyway)

        for i, row in df[genes_to_plot_names_mask].iterrows():
            ax.text(row[f'log_ratio1'], row[f'log_ratio2'], row['gene'])
    
    log_ratio_cols = [x for x in df.columns if x.startswith('log_ratio')]
    abs_log_ratio_cols = [x for x in df.columns if x.startswith('abs_log_ratio')]
    df['max_log_ratio'] = df[log_ratio_cols].max(axis=1)
    df['min_log_ratio'] = df[log_ratio_cols].min(axis=1)
    df['max_abs_log_ratio'] = df[abs_log_ratio_cols].max(axis=1)
    df['min_abs_log_ratio'] = df[abs_log_ratio_cols].min(axis=1)
    df['same_sign_for_all_log_ratios'] = np.sign(df['max_log_ratio']) == np.sign(df['min_log_ratio'])
    df.loc[df['same_sign_for_all_log_ratios'], 'same_sign_min_abs_log_ratio'] = df.loc[
        df['same_sign_for_all_log_ratios'], 'min_abs_log_ratio']
    
    df.insert(1, 'same_sign_min_abs_log_ratio', df.pop('same_sign_min_abs_log_ratio'))
    df.sort_values('same_sign_min_abs_log_ratio', ascending=False, inplace=True)
    return df

def get_leave_one_out_diff_exp_df(
        c_or_mc_ad, genes, c_or_mc_mask=None, max_num_of_obs=int(5e3), sort_by_quarter_log_ratios=False,
        sample_max_obs_count=True, genes_to_also_plot=None, bin_count_for_plot=3, layer_name='downsampled',
        plot_rest_min_vs_rest_non_min=False,
):
    print('NOTE: you probably want to use generic_utils.get_equal_contrib_mask() to generate c_or_mc_mask')
    
    if c_or_mc_mask is None:
        c_or_mc_mask = np.full(c_or_mc_ad.n_obs, True)

    if genes_to_also_plot:
        assert set(genes_to_also_plot) <= set(genes)

    if layer_name == 'downsampled':
        write_downsampled_layer_if_not_exists(c_or_mc_ad)
        c_or_mc_mask &= get_downsampled_c_mask(c_or_mc_ad)
    
    num_of_obs = c_or_mc_mask.sum()
    assert num_of_obs > 0, 'num_of_obs=0'
    if num_of_obs > max_num_of_obs:
        if sample_max_obs_count:
            print(f'NOTE: sampling c_or_mc_mask to max_num_of_obs ({max_num_of_obs})')
            c_or_mc_mask = generic_utils.sample_mask(c_or_mc_mask, max_num_of_obs)
        else:            
            raise RuntimeError(f'num_of_obs > max_num_of_obs ({num_of_obs} > {max_num_of_obs}). increase max_num_of_obs to proceed anyway.')
    print(f'num_of_obs: {num_of_obs}')

    flat_dicts = []
    for gene in genes:
        curr_genes = [x for x in genes if x != gene]
        rest_expr = get_genes_expr(c_or_mc_ad, curr_genes, c_or_mc_mask=c_or_mc_mask, layer_name=layer_name)
        rest_expr_q1, rest_expr_median, rest_expr_q3 = np.quantile(rest_expr, [0.25, 0.5, 0.75])
        curr_high_half_mask = generic_utils.get_mask_filtered_by_submask(c_or_mc_mask, rest_expr > rest_expr_median)
        curr_low_half_mask = c_or_mc_mask & ~curr_high_half_mask
        curr_high_quarter_mask = generic_utils.get_mask_filtered_by_submask(c_or_mc_mask, rest_expr > rest_expr_q3)
        curr_low_quarter_mask = generic_utils.get_mask_filtered_by_submask(c_or_mc_mask, rest_expr <= rest_expr_q1)
        
        assert curr_high_half_mask.any(), str(pd.Series(rest_expr).value_counts(normalize=True))
        assert curr_low_half_mask.any(), str(pd.Series(rest_expr).value_counts(normalize=True))
        assert curr_high_quarter_mask.any(), str(pd.Series(rest_expr).value_counts(normalize=True))
        assert curr_low_quarter_mask.any(), str(pd.Series(rest_expr).value_counts(normalize=True))
        
        rest_high_half_expr = get_mean_expr(c_or_mc_ad, curr_high_half_mask, genes=gene, verbose=False, scale_final_geomean_to_sum_to_1=False)
        assert len(rest_high_half_expr) == 1
        rest_high_half_expr = rest_high_half_expr[0]
        
        rest_low_half_expr = get_mean_expr(c_or_mc_ad, curr_low_half_mask, genes=gene, verbose=False, scale_final_geomean_to_sum_to_1=False)
        assert len(rest_low_half_expr) == 1
        rest_low_half_expr = rest_low_half_expr[0]
        
        rest_high_quarter_expr = get_mean_expr(c_or_mc_ad, curr_high_quarter_mask, genes=gene, verbose=False, scale_final_geomean_to_sum_to_1=False)
        assert len(rest_high_quarter_expr) == 1
        rest_high_quarter_expr = rest_high_quarter_expr[0]
        
        # print(f'curr_low_quarter_mask.sum(): {curr_low_quarter_mask.sum()}')
        rest_low_quarter_expr = get_mean_expr(c_or_mc_ad, curr_low_quarter_mask, genes=gene, verbose=False, scale_final_geomean_to_sum_to_1=False)
        assert len(rest_low_quarter_expr) == 1
        rest_low_quarter_expr = rest_low_quarter_expr[0]
        
        flat_dicts.append({
            'gene': gene,
            'rest_high_half_expr': rest_high_half_expr,
            'rest_low_half_expr': rest_low_half_expr,
            'rest_high_quarter_expr': rest_high_quarter_expr,
            'rest_low_quarter_expr': rest_low_quarter_expr,
            'rest_expr_q1': rest_expr_q1,
            'rest_expr_median': rest_expr_median,
            'rest_expr_q3': rest_expr_q3,
        })

        if genes_to_also_plot and (gene in genes_to_also_plot):
            print(f'plotting {gene}')
            fig, ax = plt.subplots()
            if plot_rest_min_vs_rest_non_min:
                rest_expr_min = rest_expr.min()
                # if rest_expr_min != rest_expr_median:
                #     print('rest_expr_min != rest_expr_median')
                rest_expr_min2 = rest_expr[rest_expr > rest_expr_min].min()
                curr_hues = generic_utils.get_bin_median_and_repr_vecs(rest_expr, bins=[-np.inf, rest_expr_min2, np.inf])
            else:
                curr_hues = generic_utils.get_bin_median_and_repr_vecs(rest_expr, num_of_bins=bin_count_for_plot)
            plot_gene_hist(
                c_or_mc_ad,
                gene,
                c_or_mc_mask=c_or_mc_mask,
                hue_named_series=pd.Series(curr_hues[0], name=f'bin_median_genes_except_{gene}'),
                ax=ax,
            )
            if plot_rest_min_vs_rest_non_min:
                ax.set_title(ax.get_title() + ' (min vs non-min)')



    df = pd.DataFrame(flat_dicts)
    df.insert(1, 'half_log_ratio', df['rest_high_half_expr'] - df['rest_low_half_expr'])
    df.insert(2, 'quarter_log_ratio', df['rest_high_quarter_expr'] - df['rest_low_quarter_expr'])
    df.sort_values('quarter_log_ratio' if sort_by_quarter_log_ratios else 'half_log_ratio', ascending=False, inplace=True)
    return df
    

def get_mean_expr(
        c_or_mc_ad, c_or_mc_mask, layer_name=None, max_num_of_obs=int(2e3),
        epsilon_for_log=1e-5, 
        ret_df=False,
        sort_df_by_expr_rest_max_expr_log_ratio=None,
        background_mask_for_rest_max_expr=None,
        min_max_expr_to_keep_in_df=-np.inf,
        genes=None,
        use_geomean=True, epsilon_for_c_log=1/16,
        scale_final_geomean_to_sum_to_1=True, # 240325: only now changed to True (as mc.pl.collect_metacells() does that...)
        also_return_iqr=False,
        mc_ad=None,
        verbose=True,
        sample_max_obs_count=True,
):
    num_of_obs = c_or_mc_mask.sum()
    assert num_of_obs > 0, 'num_of_obs=0'
    if num_of_obs > max_num_of_obs:
        if sample_max_obs_count:
            print(f'NOTE: sampling c_or_mc_mask to max_num_of_obs ({max_num_of_obs})')
            c_or_mc_mask = generic_utils.sample_mask(c_or_mc_mask, max_num_of_obs)
        else:
            raise RuntimeError(f'num_of_obs > max_num_of_obs ({num_of_obs} > {max_num_of_obs}). increase max_num_of_obs to proceed anyway.')
    elif verbose:
        print(f'num_of_obs: {num_of_obs}')
    ad_is_mc = is_mc_ad(c_or_mc_ad)

    if also_return_iqr:
        if layer_name != 'expr':
            raise NotImplementedError('sorry')

    if layer_name is None:
        if 'expr' in c_or_mc_ad.layers:
            layer_name = 'expr'
        # elif 'downsampled' in c_or_mc_ad.layers:
        #     layer_name = 'downsampled'
        else:
            layer_name = 'X'
    assert not (sort_df_by_expr_rest_max_expr_log_ratio and (layer_name != 'expr'))
    gene_mask = get_gene_mask(c_or_mc_ad, genes) if genes else np.full(c_or_mc_ad.n_vars, True)
    if layer_name == 'expr':
        expr = mc.ut.to_numpy_matrix(c_or_mc_ad.layers[layer_name][np.ix_(c_or_mc_mask, gene_mask)])
        mean_expr = mc.ut.to_numpy_vector(expr.mean(axis=0))
        if also_return_iqr:
            expr_iqr = scipy.stats.iqr(expr, axis=0)
            
    else:
        if layer_name == 'X':
            x = c_or_mc_ad.X
        else:
            x = c_or_mc_ad.layers[layer_name]

        curr_x = mc.ut.to_numpy_matrix(x[np.ix_(c_or_mc_mask, gene_mask)])
        # print(f'ad_is_mc: {ad_is_mc}')
        if ad_is_mc:
            x_norm = curr_x
        else:
            umi_count_of_cells = mc.ut.to_numpy_vector(x[c_or_mc_mask, :].sum(axis=1))
            x_norm = curr_x / umi_count_of_cells[:, np.newaxis]
        if use_geomean:
            if ad_is_mc:
                expr = np.log2(x_norm + epsilon_for_log)
                mean_expr = mc.ut.to_numpy_vector(expr.mean(axis=0))
            else:
                mean_expr = np.log2(get_norm_umi_count_geomean(
                    curr_x, epsilon_for_log=epsilon_for_c_log, umi_count_of_cells=umi_count_of_cells, 
                    scale_final_geomean_to_sum_to_1=scale_final_geomean_to_sum_to_1,
                ) + epsilon_for_log)
        else:
            raise RuntimeError('why not use geomean? sounds like you want to pool cells, not take their mean if you dont use geomean, i think...')

    
    if ret_df:
        df = pd.DataFrame({
            'gene': c_or_mc_ad.var_names[gene_mask],
            'expr': mean_expr,
        })
        if also_return_iqr:
            df['expr_iqr'] = expr_iqr
        if ((layer_name == 'expr') or ((not ad_is_mc) and (mc_ad is not None))) and (
            sort_df_by_expr_rest_max_expr_log_ratio or (sort_df_by_expr_rest_max_expr_log_ratio is None)
        ):
            if ad_is_mc:
                rest_mask = ~c_or_mc_mask
                mc_ad = c_or_mc_ad
            else:
                rest_mask = np.full(mc_ad.n_obs, True)
            if background_mask_for_rest_max_expr is not None:
                rest_mask &= background_mask_for_rest_max_expr
            df['rest_max_expr'] = mc.ut.to_numpy_vector(mc_ad.layers['expr'][np.ix_(rest_mask, gene_mask)].max(axis=0))
            df = df[df['rest_max_expr'] >= min_max_expr_to_keep_in_df]
            df['expr_rest_max_expr_log_ratio'] = df['expr'] - df['rest_max_expr']
            df.sort_values('expr_rest_max_expr_log_ratio', ascending=False, inplace=True)
        else:
            df.sort_values('expr', ascending=False, inplace=True)
        return df
    return mean_expr

def write_state_expr_and_state_expr_enrich(mc_ad, state_col_name='state', overwrite_state_to_median_expr=False, overwrite_state_expr_enrich=False):
    states = sorted(mc_ad.obs[state_col_name].unique())
    assert 'expr' in mc_ad.layers
    if ('state_expr_enrich' not in mc_ad.uns) or overwrite_state_to_median_expr:
        mc_ad.uns['state_to_median_expr'] = {
            state: np.median(mc_ad.layers['expr'][mc_ad.obs[state_col_name] == state, :], axis=0)
            for state in states
        }
    if ('state_expr_enrich' not in mc_ad.layers) or overwrite_state_expr_enrich:
        mc_ad.layers['state_expr_enrich'] = np.full(mc_ad.X.shape, np.nan)
        for state in states:
            # print(state)
            state_mask = mc_ad.obs[state_col_name] == state
            mc_ad.layers['state_expr_enrich'][state_mask, :] = mc_ad.layers['expr'][state_mask, :] - mc_ad.uns['state_to_median_expr'][state][np.newaxis, :]

def get_mc_state_expr_enrich_df(mc_ad, mc_name):
    df = pd.DataFrame({
        'gene': mc_ad.var_names,
        'state_expr_enrich': mc.ut.to_numpy_vector(mc_ad.layers['state_expr_enrich'][mc_ad.obs_names == mc_name, :]),
    })
    df['abs_state_expr_enrich'] = df['state_expr_enrich'].abs()
    df.sort_values('abs_state_expr_enrich', ascending=False, inplace=True)
    return df

def get_state_outlier_mc_df(mc_ad, mc_mask=None, min_abs_state_expr_enrich=1, min_state_mc_count=4):
    write_state_expr_and_state_expr_enrich(
        mc_ad,
        overwrite_state_to_median_expr=True,
        overwrite_state_expr_enrich=True,
    )

    if mc_mask is None:
        mc_mask = np.full(mc_ad.n_obs, True)

    states = list((mc_ad.obs.loc[mc_mask, 'state'].value_counts() >= min_state_mc_count).loc[lambda x: x].index)
    mc_mask = mc_mask & mc_ad.obs['state'].isin(states)

    df = pd.DataFrame({
        'metacell_name': mc_ad.obs_names[mc_mask],
        'state': mc_ad.obs.loc[mc_mask, 'state'],
        'high_abs_state_expr_enrich_count': (np.abs(mc_ad.layers['state_expr_enrich'][mc_mask, :]) > min_abs_state_expr_enrich).sum(axis=1),
    })
    
    df = generic_utils.merge_preserving_df1_index_and_row_order(
        df, 
        df.groupby('state')['high_abs_state_expr_enrich_count'].median().reset_index(
            name='state_median_high_abs_state_expr_enrich_count')
    )
    df['norm_high_abs_state_expr_enrich_count'] = (df['high_abs_state_expr_enrich_count'] + 1) / (df['state_median_high_abs_state_expr_enrich_count'] + 1)

    df.sort_values('norm_high_abs_state_expr_enrich_count', ascending=False, inplace=True)
    return df

def get_mc_sig_var_mean_df(c_ad, mc_ad, sig_name_to_genes, max_sig_inner_fold_count=None, epsilon_for_log=1e-5, overwrite_existing_c_ad_sig_cols=True): 
    if max_sig_inner_fold_count is None:
        print('choose max_sig_inner_fold_count using the returned df')
        return mc_ad.var[mc_ad.var['significant_inner_folds_count'] > 5].sort_values('significant_inner_folds_count', ascending=False)

    genes_to_ignore = list(mc_ad.var_names[mc_ad.var['significant_inner_folds_count'] > max_sig_inner_fold_count])

    if 'downsampled' not in c_ad.layers:
        mc.tl.downsample_cells(c_ad, random_seed=1234)
    sig_name_to_col_name = add_cell_scores(
        c_ad, sig_name_to_genes, 
        overwrite_existing_cols=overwrite_existing_c_ad_sig_cols, 
        genes_to_ignore=genes_to_ignore, 
        layer_name='downsampled', add_as_umi_count=True,
    )

    dfs = []
    c_ad_obs = c_ad.obs[c_ad.obs['metacell_name'] != 'Outliers'].copy()
    c_ad_obs[['metacell_name', 'state']] = c_ad_obs[['metacell_name', 'state']].astype(str)
    for sig_name, col_name in sig_name_to_col_name.items():
        print(sig_name)
        df = c_ad_obs.groupby(['metacell_name', 'state'])[col_name].agg([np.nanmean, np.nanvar]).reset_index().rename(
            columns={'nanmean': 'ds_umi_count_mean', 'nanvar': 'ds_umi_count_var'})
        df['var_mean_log_ratio'] = np.log2(df['ds_umi_count_var'] + 1) - np.log2(df['ds_umi_count_mean'] + 1)
        df['sig_name'] = sig_name
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    df['ds_umi_count_mean_expr'] = np.log2((df['ds_umi_count_mean'] / c_ad.uns['downsample_samples']) + epsilon_for_log)

    df = generic_utils.merge_preserving_df1_index_and_row_order(df, df.groupby(['state', 'sig_name'])['var_mean_log_ratio'].median().reset_index(
        name='state_median_var_mean_log_ratio'))
    df['norm_var_mean_log_ratio'] = df['var_mean_log_ratio'] - df['state_median_var_mean_log_ratio']
    df.sort_values('norm_var_mean_log_ratio', ascending=False, inplace=True)
    df.insert(2, 'norm_var_mean_log_ratio', df.pop('norm_var_mean_log_ratio'))
    return df

def plot_mc_sig_var_mean_scatter(mc_ad, mc_sig_var_mean_df, sig_name, states=None, legend=False, **scatterplot_kwargs):
    fig, ax = plt.subplots()
    df = mc_sig_var_mean_df[mc_sig_var_mean_df['sig_name'] == sig_name]
    if states is not None:
        df = df[df['state'].isin(states)]

    sb.scatterplot(
        df,
        x='ds_umi_count_mean_expr',
        y='var_mean_log_ratio',
        hue='state',
        palette=get_palette(mc_ad, color_by='state'),
        legend=legend,
        linewidth=0.3,
        edgecolor='black',
        ax=ax,
        **scatterplot_kwargs
    )
    ax.set_title(sig_name)
    return fig, ax

def old_get_state_outlier_mc_df(mc_ad, mc_mask=None, min_log_ratio=1, min_high_low_ratio_gene_count=1):
    raise RuntimeError('you probably want get_state_outlier_mc_df')
    # NOTE: two similar outliers will hide each other. ugh.
    if mc_mask is None:
        mc_mask = np.full(mc_ad.n_obs, True)
    num_of_mcs = mc_mask.sum()
    assert num_of_mcs, 'empty mc_mask'
    state_to_mc_mask = {
        state: mc_ad.obs['state'] == state
        for state in mc_ad.obs.loc[mc_mask, 'state'].unique()
    }
    for state, curr_mc_mask in state_to_mc_mask.items():
        if curr_mc_mask.sum() == 1:
            print(f'skipping {state} because it has a single metacell')
            mc_mask = mc_mask & (mc_ad.obs['state'] != state)
    rows = []
    for i, mc_name in enumerate(mc_ad.obs_names[mc_mask]):
        if i % 20 == 19:
            print(f'get_state_outlier_mc_df {i}/{num_of_mcs}')
        curr_mc_mask = mc_ad.obs_names == mc_name
        curr_mc_state = mc_ad.obs.loc[mc_name, 'state']
        curr_df = get_mean_expr(
            mc_ad, curr_mc_mask, background_mask_for_rest_max_expr=state_to_mc_mask[curr_mc_state], ret_df=True,
            verbose=False,
        )
        curr_df = curr_df[curr_df['expr_rest_max_expr_log_ratio'] >= min_log_ratio]
        high_low_ratio_gene_count = len(curr_df)
        if high_low_ratio_gene_count >= min_high_low_ratio_gene_count:
            curr_row = curr_df.iloc[0].copy()
            curr_row['high_low_ratio_gene_count'] = high_low_ratio_gene_count
            curr_row['metacell_name'] = mc_name
            curr_row['state'] = curr_mc_state
            rows.append(curr_row)
    if not rows:
        print('didnt find any outlier metacell')
        return
    df = pd.DataFrame(rows)
    # df.sort_values('expr_rest_max_expr_log_ratio', ascending=False, inplace=True)
    df.sort_values('high_low_ratio_gene_count', ascending=False, inplace=True)
    return df

def get_pooled_total_norm_expression(c_or_mc_ad, cell_or_mc_mask, downsample_layer_name='downsampled'):
    if downsample_layer_name is None:
        pooled_x = c_or_mc_ad.X[cell_or_mc_mask, :].sum(axis=0)
    else:
        pooled_x = c_or_mc_ad.layers[downsample_layer_name][cell_or_mc_mask, :].sum(axis=0)
    pooled_x = mc.ut.to_numpy_vector(pooled_x)
    
    return pooled_x / pooled_x.sum()

def get_pooled_log_total_norm_expression(
        c_or_mc_ad, cell_or_mc_mask, downsample_layer_name='downsampled', epsilon_for_log=1e-5, return_expr_df=False,
):
    expr = np.log2(get_pooled_total_norm_expression(
        c_or_mc_ad, cell_or_mc_mask, downsample_layer_name=downsample_layer_name) + epsilon_for_log)
    if return_expr_df:
        return pd.DataFrame({
            'gene': c_or_mc_ad.var_names,
            'expr': expr,
        })
    return expr

def old_naive_get_mean_expr(c_ad, c_mask, gene_mask=None):
    if gene_mask is None:
        gene_mask = np.full(c_ad.n_vars, True)
    
    return mc.ut.to_numpy_vector(c_or_mc_ad.layers[layer_name][cell_or_mc_mask, :][:, gene_mask].mean(axis=0))

def get_expr_corr_with_mcs_df(
        mc_ad, expr_vec=None, expr_df=None, gene_mask=None, 
        # expr_vec_already_filtered=False, 
        mc_obs_columns_to_add=['state', 'projected_type'],
        num_of_top_corr_mcs_to_print_diff_exp_df=2,
        diff_exp_min_abs_log_ratio=1,
        show_only_up_or_down_regulated_genes=None,
):
    if gene_mask is None:
        # gene_mask = mc_ad.var['significant_gene'].to_numpy()
        gene_mask = mc_ad.var['selected_gene'].to_numpy()


    if expr_vec is None:
        gene_names = mc_ad.var_names[gene_mask]
        genes_missing_from_expr_df = set(gene_names) - set(expr_df['gene'])
        if genes_missing_from_expr_df:
            print(f'NOTE: genes_missing_form_expr_df ({len(genes_missing_from_expr_df)}): {genes_missing_from_expr_df}')
            gene_mask &= mc_ad.var_names.isin(expr_df['gene'])
        expr_vec = mc.ut.to_numpy_vector(generic_utils.merge_preserving_df1_index_and_row_order(
            pd.DataFrame({'gene': mc_ad.var_names[gene_mask]}),
            expr_df,
        )['expr'].to_numpy())
    else:
        if expr_vec.ndim != 1:
            expr_vec = mc.ut.to_numpy_vector(expr_vec)[gene_mask]
    expr = expr_vec[np.newaxis, :]

    orig_expr_df = expr_df
    expr_df = pd.DataFrame({'expr': expr_vec, 'gene': mc_ad.var_names[gene_mask]})
    generic_utils.merge_preserving_df1_index_and_row_order(expr_df, orig_expr_df[['gene', 'expr']]) # a sanity check
    
    assert not np.isnan(expr).any()
    expr_of_mcs = mc_ad.layers['expr'][:, gene_mask]
    corr_df = pd.DataFrame({
        'metacell_name': mc_ad.obs_names,
        'corr': mc.ut.to_numpy_vector(mc.ut.cross_corrcoef_rows(
            mc.ut.to_layout(expr, layout='row_major'),
            mc.ut.to_layout(expr_of_mcs, layout='row_major'),
            reproducible=True,
        )),
    })
    corr_df.sort_values('corr', ascending=False, inplace=True)
    
    mc_obs_columns_to_add = [x for x in mc_obs_columns_to_add if x in mc_ad.obs.columns]
    if 'metacell_name' in mc_ad.obs.columns:
        curr_mc_ad_obs = mc_ad.obs[mc_obs_columns_to_add + ['metacell_name']]
    else:
        curr_mc_ad_obs = mc_ad.obs[mc_obs_columns_to_add].copy()
        curr_mc_ad_obs['metacell_name'] = mc_ad.obs_names

    corr_df = generic_utils.merge_preserving_df1_index_and_row_order(
        corr_df,
        curr_mc_ad_obs,
    )

    for top_mc_i in range(num_of_top_corr_mcs_to_print_diff_exp_df):
        mc_name = corr_df.iloc[top_mc_i]['metacell_name']
        diff_exp_df = get_diff_exp_df(mc_ad, expr1_df=expr_df, c_or_mc_mask2=mc_ad.obs['metacell_name'] == mc_name)
        # print(diff_exp_df[diff_exp_df['gene'] == 'CAMP'])
        diff_exp_df = diff_exp_df[diff_exp_df['abs_log_ratio'] >= diff_exp_min_abs_log_ratio]
        if show_only_up_or_down_regulated_genes == 'up':
            diff_exp_df = diff_exp_df[diff_exp_df['log_ratio'] > 0]
        elif show_only_up_or_down_regulated_genes == 'down':
            diff_exp_df = diff_exp_df[diff_exp_df['log_ratio'] < 0]
        else:
            assert show_only_up_or_down_regulated_genes is None
        print(f'top mc {top_mc_i} ({mc_name}) diff exp:')
        print(diff_exp_df)

    return corr_df

def get_pooled_cells_expr_corr_with_mcs_df(mc_ad, c_ad, c_mask, layer_name='X', gene_names=None):
    pooled_expr_df = get_mean_expr(
        c_ad, c_mask, layer_name=layer_name,
        ret_df=True,
    )
    pooled_expr = generic_utils.merge_preserving_df1_index_and_row_order(mc_ad.var_names.to_frame(name='gene'), pooled_expr_df, how='left')['expr'].to_numpy()
    if gene_names is None:
        gene_mask = None
    else:
        gene_mask = get_gene_mask(mc_ad, gene_names)
    return get_expr_corr_with_mcs_df(mc_ad, pooled_expr, gene_mask=gene_mask)

def get_mean_mc_expr_corr_with_mcs_df(mc_ad, mc_mask, gene_mask=None, name_to_add_to_mc_ad_obs=None, layer_name='X'):
    mean_mc_expr = get_mean_expr(mc_ad, mc_mask, layer_name=layer_name)
    df = get_expr_corr_with_mcs_df(mc_ad, mean_mc_expr, gene_mask=gene_mask)
    if name_to_add_to_mc_ad_obs is not None:
        if name_to_add_to_mc_ad_obs in mc_ad.obs.columns:
            mc_ad.obs.drop(columns=name_to_add_to_mc_ad_obs, inplace=True)
        mc_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(mc_ad.obs, df.rename(columns={'corr': name_to_add_to_mc_ad_obs}))
        plot_manifold_umap(mc_ad)
        plot_manifold_umap(mc_ad, color_by=name_to_add_to_mc_ad_obs)

    return df



def get_total_norm_expression_of_cells(c_ad, gene_names):
    raise RuntimeError('move to use get_genes_umi_frac')
    gene_indices = get_gene_indices(c_ad, gene_names)
    return mc.ut.to_numpy_vector(c_ad.X[:, gene_indices].sum(axis=1) / c_ad.X.sum(axis=1))


def get_norm_projected_expression(mc_ad):
    return mc_ad.layers['projected'] / mc_ad.layers['projected'].sum(axis=1)[:, np.newaxis]

def get_log_norm_expression_zm_score(mc_ad, metacell_mask, metacell_of_interest_name_or_i=None):
    medians = np.median(mc_ad.layers['expr'][metacell_mask, :], axis=0)
    stds = mc_ad.layers['expr'][metacell_mask, :].std(axis=0)
    with generic_utils.allow_numpy_division_by_zero_context_manager():
        zm_scores = (mc_ad.layers['expr'][metacell_mask, :] - medians) / stds
    if metacell_of_interest_name_or_i:
        mc_of_interest_i = get_metacell_i(metacell_of_interest_name_or_i)
        with generic_utils.allow_numpy_division_by_zero_context_manager():
            mc_of_interest_zm_scores = (mc_ad.layers['expr'][mc_of_interest_i, :] - medians) / stds
        return zm_scores, mc_of_interest_zm_scores
    return zm_scores

def get_log_norm_expression_cell_type_zm_score(mc_ad, cell_state_column_name):
    cell_type_zm_scores = np.full(mc_ad.X.shape, np.nan)
    for cell_type in mc_ad.obs[cell_state_column_name].unique():
        metacell_mask = mc_ad.obs['projected_type'] == cell_type
        cell_type_zm_scores[metacell_mask, :] = get_log_norm_expression_zm_score(mc_ad, metacell_mask=metacell_mask)
    return cell_type_zm_scores

def write_expr_and_expr_enrich(mc_ad, epsilon_to_add_for_log=1e-5):
    if {'expr', 'expr_enrich'} <= set(mc_ad.layers):
        print('skipping writing expr and expr_enrich as they already exist')
        return
    # if not allow_earlier_mc_version: # this makes no sense because this isn't the version used to run mc.pl.collect_metacells.
    #     assert mc_ad.uns['metacells_algorithm'] in {'metacells.0.9.0-dev.1', 'metacells.0.9.0', 'metacells.0.9.1-dev', 'metacells.0.9.1', 'metacells.0.9.4dev', 'metacells.0.9.4'}
    expr = np.log2(mc.ut.to_numpy_matrix(mc_ad.X) + epsilon_to_add_for_log)
    expr_enrich = expr - np.median(expr, axis=0)
    mc.ut.set_vo_data(mc_ad, name='expr', data=expr)
    mc.ut.set_vo_data(mc_ad, name='expr_enrich', data=expr_enrich)

def write_proj_expr(mc_ad, epsilon_to_add_for_log=1e-5, silent=False):
    if 'proj_expr' in mc_ad.layers:
        if not silent:
            print('skipping writing proj_expr as it already exists')
        return
    assert mc_ad.uns['metacells_algorithm'] in {'metacells.0.9.0-dev.1', 'metacells.0.9.0', 'metacells.0.9.1-dev', 'metacells.0.9.1', 'metacells.0.9.4dev', 'metacells.0.9.4'}
    proj_expr = np.log2(mc.ut.to_numpy_matrix(mc_ad.layers['projected_fraction']) + epsilon_to_add_for_log)
    mc.ut.set_vo_data(mc_ad, name='proj_expr', data=proj_expr)
    

def get_log_norm_expression_and_log_enrichment_and_zm_score_and_norm_expression_std(
        mc_ad, epsilon_to_add_to_mc_gene_norm_expression, 
        metacell_mask_to_calculate_median_and_std=None,
):
    norm_expression = get_norm_expression(mc_ad)
    
    if metacell_mask_to_calculate_median_and_std is None:
        metacell_mask_to_calculate_median_and_std = np.full(mc_ad.n_obs, True)
    median_norm_expression = np.median(norm_expression[metacell_mask_to_calculate_median_and_std, :], axis=0)
    norm_expression_std = norm_expression[metacell_mask_to_calculate_median_and_std, :].std(axis=0)

    if 'projected' in mc_ad.layers:
        norm_projected_expression = get_norm_projected_expression(mc_ad)
        log_norm_projected_expression = np.log2(norm_projected_expression + epsilon_to_add_to_mc_gene_norm_expression) 
    else:
        log_norm_projected_expression = None

    expr = np.log2(norm_expression + epsilon_to_add_to_mc_gene_norm_expression) 
    log_enrichment = expr - np.median(expr, axis=0)

    with generic_utils.allow_numpy_division_by_zero_context_manager():
        zm_score = (norm_expression - median_norm_expression) / norm_expression_std

    return expr, log_norm_projected_expression, log_enrichment, zm_score, norm_expression_std

def write_log_norm_expression_and_log_enrichment_and_zm_score_and_norm_expression_std(mc_ad, epsilon_to_add_to_mc_gene_norm_expression):
    expr, log_norm_projected_expression, log_enrichment, zm_score, norm_expression_std = (
        get_log_norm_expression_and_log_enrichment_and_zm_score_and_norm_expression_std(mc_ad, epsilon_to_add_to_mc_gene_norm_expression))

    mc.ut.set_vo_data(
        mc_ad,
        name='expr',
        data=expr,
    )
    if log_norm_projected_expression is not None:
        mc.ut.set_vo_data(
            mc_ad,
            name='log_norm_projected_expression',
            data=log_norm_projected_expression,
        )
    mc.ut.set_vo_data(
        mc_ad,
        name='log_enrichment',
        data=log_enrichment,
    )
    mc.ut.set_vo_data(
        mc_ad,
        name='zm_score',
        data=zm_score,
    )
    mc.ut.set_v_data(
        mc_ad,
        name='norm_expression_std',
        data=norm_expression_std,
    )
    
    for cell_state_column_name in ('cell_type', 'projected_type'):
        if cell_state_column_name in list(mc_ad.obs):
            projected_type_zm_score = get_log_norm_expression_cell_type_zm_score(mc_ad, cell_state_column_name=cell_state_column_name)
            mc.ut.set_vo_data(
                mc_ad,
                name=f'{cell_state_column_name}_zm_score',
                data=projected_type_zm_score,
            )

    

def get_gene_module_log_norm_expression_and_log_enrichment(mc_ad, gene_module_names, epsilon_to_add_to_mc_gene_norm_expression):
    gene_module_mask = mc_ad.var_names.isin(gene_module_names)
    module_expression = mc_ad.X[:,gene_module_mask].sum(axis=1)
    norm_module_expression = module_expression / mc_ad.X.sum(axis=1)
    log_norm_module_expression = np.log2(norm_module_expression + epsilon_to_add_to_mc_gene_norm_expression)
    
    if 0:
        # this is how MCView calculates enrichment. i don't understand why to add epsilon to the numerator. perhaps a bug?
        log_module_enrichment = np.log2(
            (norm_module_expression + epsilon_to_add_to_mc_gene_norm_expression) / (np.median(norm_module_expression) + epsilon_to_add_to_mc_gene_norm_expression))
    else:
        log_module_enrichment = np.log2(norm_module_expression / (np.median(norm_module_expression) + epsilon_to_add_to_mc_gene_norm_expression))
    
    return log_norm_module_expression, log_module_enrichment

def get_lowest_relative_variance_genes_for_metacell_set(
        mc_ad, metacell_mask, n=20, cell_type_median_log_enrichment_across_mcs_df=None, max_rank_for_cell_type=20,
        median_log_enrichment_threshold_for_exclusively_enriched=None,
):
    # i think i should just use
    # (mc_ad.X.var(axis=0) / (mc_ad.X.mean(axis=0) + 1e-5))
    # i.e., relative variance of binomial vars, which should be 1 (assuming p is very low, which is the case - should be 1e-5 to 1e-2)
    pass

def get_exclusively_enriched_cell_type_median_log_enrichment_across_mcs_df(
        cell_type_median_log_enrichment_across_mcs_df, median_log_enrichment_threshold_for_exclusively_enriched):
    df = cell_type_median_log_enrichment_across_mcs_df
    filtered_df = df.loc[df['median_log_enrichment'] >= median_log_enrichment_threshold_for_exclusively_enriched]
    filtered_df = filtered_df[filtered_df['cell_type'] != 'Unknown']
    non_exclusive_genes_series = filtered_df.loc[filtered_df['gene'].duplicated(), 'gene'].drop_duplicates()
    return df.loc[df.merge(non_exclusive_genes_series, how='left', indicator=True)['_merge'] == 'left_only']

def get_mcs_top_enriched_genes(
        mc_ad, metacell_mask, n_per_metacell=5,
):
    raise RuntimeError('i think you want get_top_genes_dfs')
    dfs = []
    for metacell in mc_ad.obs_names[metacell_mask]:
        # print(f'metacell: {metacell}', type(metacell))
        curr_df = get_metacell_top_enriched_genes(mc_ad, metacell, n=n_per_metacell)
        curr_df['metacell'] = metacell
        dfs.append(curr_df)
    return pd.concat(dfs, ignore_index=True)


def get_metacell_i(metacell_name_or_i):
    if isinstance(metacell_name_or_i, str):
        assert re.fullmatch('MC[0-9]+', metacell_name_or_i)
        return int(metacell_name_or_i.partition('MC')[-1])
    else:
        return metacell_name_or_i

def get_neighbor_mcs_df(mc_ad, metacell_name_or_i):
    mc_i = get_metacell_i(metacell_name_or_i)
    outgoing_weights = mc.ut.to_numpy_vector(mc_ad.obsp['obs_outgoing_weights'][mc_i, :])
    nearest_neighbors_mask = outgoing_weights > 0
    return pd.DataFrame({
        'neighbor_mc_name': mc_ad.obs_names[nearest_neighbors_mask], 'outgoing_weight': outgoing_weights[nearest_neighbors_mask], 
        'neighbor_projected_type': mc_ad.obs.loc[nearest_neighbors_mask, 'projected_type'],
        'neighbor_projected_secondary_type': mc_ad.obs.loc[nearest_neighbors_mask, 'projected_secondary_type'],
    }).sort_values('outgoing_weight', ascending=False)


def get_metacell_zm_score_df(
        mc_ad, metacell_name_or_i, zm_score_layer_name='projected_type_zm_score', min_log_norm_expression_for_positive_zm_score=None,
):
    mc_i = get_metacell_i(metacell_name_or_i)
    df = pd.DataFrame({
        zm_score_layer_name: mc_ad.layers[zm_score_layer_name][mc_i,:], 
        'log_enrichment': mc_ad.layers['log_enrichment'][mc_i,:], 
        'expr': mc_ad.layers['expr'][mc_i,:], 
        'gene': mc_ad.var_names,
    }).sort_values(zm_score_layer_name, ascending=False)
    if min_log_norm_expression_for_positive_zm_score is not None:
        df = df[
            (df['expr'] >= min_log_norm_expression_for_positive_zm_score) |
            (df[zm_score_layer_name] < 0)
        ]
    return df

def get_metacell_top_enriched_genes(
        mc_ad, metacell_name_or_i, n=20, cell_type_median_log_enrichment_across_mcs_df=None, max_rank_for_cell_type=20,
        median_log_enrichment_threshold_for_exclusively_enriched=None,
):
    mc_i = get_metacell_i(metacell_name_or_i)
    df = pd.DataFrame({'log_enrichment': mc_ad.layers['log_enrichment'][mc_i,:], 'expr': mc_ad.layers['expr'][mc_i,:], 'gene': mc_ad.var_names}).sort_values(
        'log_enrichment', ascending=False).head(n)
    if cell_type_median_log_enrichment_across_mcs_df is not None:
        if median_log_enrichment_threshold_for_exclusively_enriched is not None:
            cell_type_median_log_enrichment_across_mcs_df = get_exclusively_enriched_cell_type_median_log_enrichment_across_mcs_df(
                cell_type_median_log_enrichment_across_mcs_df, median_log_enrichment_threshold_for_exclusively_enriched)
        df = df.merge(cell_type_median_log_enrichment_across_mcs_df)
        df = df.loc[df['rank_for_cell_type'] <= max_rank_for_cell_type]
    return df

def get_metacell_top_zm_score_genes(
        mc_ad, metacell_name_or_i, n_per_metacell=None, min_norm_expression_std=None, min_log_norm_expression_for_pos_zm_score=None,
):
    mc_i = get_metacell_i(metacell_name_or_i)
    df = pd.DataFrame({'zm_score': mc_ad.layers['zm_score'][mc_i,:], 'expr': mc_ad.layers['expr'][mc_i,:], 'gene': mc_ad.var_names, 'norm_expression_std': mc_ad.var['norm_expression_std']})
    
    df['abs_zm_score'] = df['zm_score'].abs()
    df.sort_values('abs_zm_score', ascending=False, inplace=True)
    df.drop('abs_zm_score', axis=1, inplace=True)

    if min_norm_expression_std is not None:
        df = df[df['norm_expression_std'] >= min_norm_expression_std]
    
    if min_log_norm_expression_for_pos_zm_score is not None:
        df = df[(df['zm_score'] < 0) | (df['expr'] >= min_log_norm_expression_for_pos_zm_score)]

    if n_per_metacell is not None:
        df = df.head(n_per_metacell)
    return df


def get_mcs_top_zm_score_genes(
        mc_ad, metacell_mask, n_per_metacell=5, min_norm_expression_std=None, min_log_norm_expression_for_pos_zm_score=None, min_count=None,
):
    dfs = []
    for metacell in mc_ad.obs_names[metacell_mask]:
        # print(f'metacell: {metacell}', type(metacell))
        curr_df = get_metacell_top_zm_score_genes(
            mc_ad, metacell, n_per_metacell=n_per_metacell, min_norm_expression_std=min_norm_expression_std,
            min_log_norm_expression_for_pos_zm_score=min_log_norm_expression_for_pos_zm_score)
        curr_df['metacell'] = metacell
        dfs.append(curr_df)
    
    all_df = pd.concat(dfs, ignore_index=True)
    if min_count is not None:
        gene_counts = all_df['gene'].value_counts()
        all_df = all_df.merge(pd.Series(np.array(gene_counts[gene_counts >= min_count].index), name='gene'))
    
    return all_df

def get_palette_from_cell_type_color_file(cell_type_colors_csv_file_path):
    df = pd.read_csv(cell_type_colors_csv_file_path)
    return generic_utils.get_dict_mapping_one_df_column_to_other(df, key_column_name='cell_type', val_column_name='color')

def get_palette(
        mc_ad, color_by=None, color_by_name=None, color_to_set_nan_color_to_anyway=None, 
        numeric_colormap_name='bwr', numeric_set_over_color='magenta', numeric_set_under_color='cyan'):
    if color_by_name is None:
        assert isinstance(color_by, str)
        color_by_name = color_by
        if color_by not in mc_ad.obs.columns:
            print(f'{color_by} not in data_ad.obs.columns, so get_palette() returns None')
            return None
        color_by = mc_ad.obs[color_by_name].copy()
    if color_by_name in {'projected_type', 'cell_type', 'state', 'type', 'c_state'}:
        color_column_name = f'{color_by_name}_color'
        if color_column_name in mc_ad.obs.columns:
            # print(color_by.dtype.name)
            # assert isinstance(color_by, pd.Series) and (color_by.dtype.name == object)
            assert isinstance(color_by, pd.Series)
            palette = generic_utils.get_dict_mapping_one_df_column_to_other(
                mc_ad.obs, key_column_name=color_by_name, val_column_name=color_column_name)
            if color_to_set_nan_color_to_anyway is not None:
                palette[np.nan] = color_to_set_nan_color_to_anyway
                palette['nan'] = color_to_set_nan_color_to_anyway
            return palette
    # print(type(color_by))
    # print(color_by.dtype.name)
    # print(isinstance(color_by, (pd.Series, np.ndarray)))
    # if isinstance(color_by, (pd.Series, np.ndarray)) and (color_by.dtype.name in {'float32', 'float64', 'int', 'int64', 'int32'}): # NOTE: changed on 240501 to use dtype.kind
    if isinstance(color_by, (pd.Series, np.ndarray)) and (color_by.dtype.kind in 'iuf'):
        # RdBu_r and vlag seem the best to me.
        
        
        # my_palette = matplotlib.cm.get_cmap('RdBu_r')
        my_palette = matplotlib.colormaps.get_cmap(numeric_colormap_name)
        # my_palette.set_bad("magenta")
        my_palette.set_over(numeric_set_over_color)
        my_palette.set_under(numeric_set_under_color)
        return my_palette


        # matplotlib.colormaps.register(my_palette, name="my_palette", force=True)
        # my_mpl_palette = sb.color_palette(
        #     "my_palette", 
        #     # n_colors=1024, 
        #     # desat=0.2,
        # )

        # return my_mpl_palette


        # return 'RdBu_r'
        # return 'vlag'
        
        # return 'seismic'
        # return 'bwr'
    return None

def get_state_palette(mc_ad, state_column_name='state'):
    palette = get_palette(mc_ad, color_by=state_column_name)
    if 'Outliers' not in palette:
        palette['Outliers'] = 'black'
    return palette

def get_gene_desc(data_ad, gene_names, max_num_of_chars_to_keep=20, skip_missing_genes=False):
    gene_names = get_gene_names(data_ad, gene_names, skip_missing_genes=skip_missing_genes)
    gene_desc = '+'.join(gene_names)
    if len(gene_desc) > max_num_of_chars_to_keep:
        gene_desc = gene_desc[:max_num_of_chars_to_keep] + '[...]'
    return gene_desc

# previously named plot_gene_umi_counts_hist
def plot_gene_hist(
        c_or_mc_ad, gene_name_or_names, c_or_mc_mask=None, hue_named_series=None, figsize=None, ax=None, 
        skip_missing_genes=False, layer_name=None, return_umi_frac=False, 
        mc_epsilon_to_add_to_fraction_before_log=MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG, 
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        cum_if_c=True,
        **hist_kwargs,
):
    if c_or_mc_mask is None:
        c_or_mc_mask = np.full(c_or_mc_ad.n_obs, True)
    elif (hue_named_series is not None) and (len(hue_named_series) != c_or_mc_ad.n_obs):
        # ugly hack. sorry.
        assert c_or_mc_mask.sum() == len(hue_named_series)
        orig_hue_named_series = hue_named_series
        hue_named_series = pd.Series(np.full(c_or_mc_ad.n_obs, np.nan), name=orig_hue_named_series.name)
        hue_named_series[c_or_mc_mask] = orig_hue_named_series.to_numpy() # without this stuff are bad. ugh.
    
    ad_is_mc = is_mc_ad(c_or_mc_ad)
    if layer_name is None:
        if ad_is_mc:
            layer_name = 'expr' if 'expr' in c_or_mc_ad.layers else 'X'
        else:
            layer_name = 'downsampled'
            

    if layer_name == 'downsampled':
        write_downsampled_layer_if_not_exists(c_or_mc_ad)
        c_or_mc_mask &= get_downsampled_c_mask(c_or_mc_ad)

    gene_vals, gene_desc = get_gene_or_genes_log_norm_or_downsampled_expression_and_desc(
        c_or_mc_ad, gene_name_or_names=gene_name_or_names, skip_missing_genes=skip_missing_genes, layer_name=layer_name,
        mc_or_cell_mask=c_or_mc_mask, return_umi_frac=return_umi_frac, 
        mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log, 
        c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log, 
    )
    
    # print(len(gene_vals), len(c_or_mc_ad.obs_names[mc_or_cell_mask]))
    df = pd.DataFrame({
        'name': c_or_mc_ad.obs_names[c_or_mc_mask],
        gene_desc: gene_vals,
    })
    if hue_named_series is not None:
        hue_name = hue_named_series.name
        df[hue_name] = hue_named_series[c_or_mc_mask].to_numpy()

        print(f'df[hue_name].value_counts(): {df[hue_name].value_counts()}')
    
    # df = df[mc_or_cell_mask]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    final_hist_kwargs = {
        **dict(
            data=df,
            x=gene_desc,
            ax=ax,
        ),
        **hist_kwargs,
    }

    cum_hist = (not ad_is_mc) and cum_if_c
    if cum_hist:
        final_hist_kwargs['cumulative'] = True
        final_hist_kwargs['fill'] = False
        final_hist_kwargs['element'] = 'step'
        final_hist_kwargs['bins'] = max(final_hist_kwargs.get('bins', 0), 1000)
        final_hist_kwargs['stat'] = 'proportion'
        final_hist_kwargs['common_norm'] = False

    if hue_named_series is not None:
        final_hist_kwargs['hue'] = hue_name
        if ('palette' not in hist_kwargs) and (hue_named_series.dtype.kind in 'iuf'):
            # final_hist_kwargs['palette'] = 'RdBu_r'
            final_hist_kwargs['palette'] = 'coolwarm'
    sb.histplot(
        **final_hist_kwargs,
    )

    # if cum_hist:
    #     ax.set_ylim(ax.get_ylim()[0], 1.01)
    #     # ax.set_ylim(-0.01, 1.01)
    
    # if you are looking for generic code, then:
    # if ax.get_legend() is not None:
    #     ax.get_legend().get_frame().set_facecolor('none')
    if ax.get_legend() is not None:
        ax.get_legend().get_frame().set_facecolor('none')
    return fig, ax

def is_mc_ad(c_or_mc_ad):
    return ('__name__' in c_or_mc_ad.uns) and (c_or_mc_ad.uns['__name__'] in {'metacells', 'atlas'})

def expr_df_to_expr_vec(c_or_mc_ad, expr_df):
    return generic_utils.merge_preserving_df1_index_and_row_order(c_or_mc_ad.var_names.to_frame(name='gene'), expr_df, how='left')['expr'].to_numpy()

def plot_gene_gene_scatter(
        c_or_mc_ad, gene_x_name, gene_y_name, 
        color_by=None,
        mask_by_state_names=None,
        max_obs_count=int(50e3),
        sample_max_obs_count=False,
        color_by_name=None, legend=False, mc_to_mark_name_to_color=None, 
        c_or_mc_mask=None,
        jitter_sd=0,
        # s=30,
        min_x_val=None,
        min_y_val=None,
        color_by_gene=None,
        mc_epsilon_to_add_to_fraction_before_log=MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        x_is_metadata=False, y_is_metadata=False,
        ax=None,
        expr_vec_infos=None,
        layer_name=None,
        print_gene_expr_quantile_df=False,
        name_to_c_mask=None,
        c_ad=None,
        allow_common_genes=False,
        figsize=None,
        hist2d=False,
        hist2d_bins=(10,10),
        hist2d_log_count=True,
        **scatter_kwargs,
):
    if (color_by is None) and (color_by_gene is None) and (color_by_name is None):
        color_by, color_by_gene = get_default_color_by_and_color_by_gene(c_or_mc_ad)
    ad_is_mc = is_mc_ad(c_or_mc_ad)
    if ad_is_mc:
        epsilon_to_add_to_fraction_before_log = mc_epsilon_to_add_to_fraction_before_log
    else:
        epsilon_to_add_to_fraction_before_log = c_epsilon_to_add_to_fraction_before_log

    if layer_name is None:
        layer_name = 'X' if ad_is_mc else 'downsampled'

    if mask_by_state_names is not None:
        assert c_or_mc_mask is None
        if isinstance(mask_by_state_names, str):
            mask_by_state_names = [mask_by_state_names]
        c_or_mc_mask = mc.ut.to_numpy_vector(c_or_mc_ad.obs['state'].isin(mask_by_state_names))
    elif c_or_mc_mask is None:
        c_or_mc_mask = np.full(c_or_mc_ad.n_obs, True)
    else:
        c_or_mc_mask = mc.ut.to_numpy_vector(c_or_mc_mask)
    
    if layer_name == 'downsampled':
        write_downsampled_layer_if_not_exists(c_or_mc_ad)
        c_or_mc_mask &= get_downsampled_c_mask(c_or_mc_ad)

    obs_count = c_or_mc_mask.sum()
    if obs_count > max_obs_count:
        if sample_max_obs_count:
            print('NOTE: sampling c_or_mc_mask to max_obs_count')
            c_or_mc_mask = generic_utils.sample_mask(c_or_mc_mask, max_obs_count)
        else:
            raise RuntimeError(f'obs_count ({obs_count}) > max_obs_count ({max_obs_count})')

    if color_by_gene:
        color_by, color_by_name = get_gene_or_genes_log_norm_or_downsampled_expression_and_desc(
            c_or_mc_ad, gene_name_or_names=color_by_gene, layer_name=layer_name,
            mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log,
            c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log,
            mc_or_cell_mask=c_or_mc_mask,
        )
        
    if hasattr(color_by, 'name'):
        color_by_name = color_by.name
    if color_by is None:
        palette = None
    elif color_by_name is None:
        if isinstance(color_by, str):
            color_by_name = color_by
            color_by = c_or_mc_ad.obs.loc[c_or_mc_mask, color_by_name].copy()
        else:
            color_by_name = 'dummy'
    
    if expr_vec_infos is None:
        expr_vec_infos = []

    if name_to_c_mask:
        assert (c_ad.var_names == c_or_mc_ad.var_names).all()
        for name, c_mask in name_to_c_mask.items():
            expr_vec_infos.append({
                'name': name,
                'expr': mc.ut.to_numpy_vector(get_mean_expr(
                    c_ad, c_mask, layer_name=layer_name, use_geomean=True, 
                    epsilon_for_log=mc_epsilon_to_add_to_fraction_before_log, 
                    epsilon_for_c_log=1/16,
                    scale_final_geomean_to_sum_to_1=True,
                )),
            })

    if (expr_vec_infos is not None) and expr_vec_infos:
        if 'expr' in expr_vec_infos[0]:
            expr_vec_infos = [{
                **x,
                'norm_expr': np.power(2, x['expr']) - epsilon_to_add_to_fraction_before_log,
            } for x in expr_vec_infos]
        elif 'umi_count' in expr_vec_infos[0]:
            expr_vec_infos = [{
                **x,
                'norm_expr': x['umi_count'] / x['umi_count'].sum(),
            } for x in expr_vec_infos]
        else:
            assert 'norm_expr' in expr_vec_infos[0], 'expr_vec_infos must have either expr, umi_count or norm_expr'
            expr_vec_infos = [x.copy() for x in expr_vec_infos]
        
        for i, expr_vec_info in enumerate(expr_vec_infos):
            norm_expr = expr_vec_info['norm_expr']
            round_up_to_zero_mask = np.isclose(norm_expr, 0) & (norm_expr < 0)
            if round_up_to_zero_mask.any():
                expr_vec_infos[i]['norm_expr'][round_up_to_zero_mask] = 0

        min_norm_exprs = [x['norm_expr'].min() for x in expr_vec_infos]
        max_norm_exprs = [x['norm_expr'].max() for x in expr_vec_infos]
        # print(min_norm_exprs)
        assert all(x >= 0 for x in min_norm_exprs), 'min norm_expr < 0'
        assert all(x <= 1 for x in max_norm_exprs), 'max norm_expr > 1'
        # assert all(((x['norm_expr'] >= 0) & (x['norm_expr'] <= 1)).all() for x in expr_vec_infos), 'ugh'

    gene_names = []
    if x_is_metadata:
        gene_x_vals = c_or_mc_ad.obs.loc[c_or_mc_mask, gene_x_name].to_numpy()
        gene_x_desc = gene_x_name
    else:
        gene_x_vals, gene_x_desc = get_gene_or_genes_log_norm_or_downsampled_expression_and_desc(
            c_or_mc_ad, gene_name_or_names=gene_x_name, 
            jitter_sd=jitter_sd, 
            mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log,
            c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log,
            layer_name=layer_name,
            mc_or_cell_mask=c_or_mc_mask,
        )
        if gene_x_desc.endswith('[...]'):
            gene_x_desc = gene_x_desc.replace('[...]', '[..]') # this is to prevent gene_x_desc == gene_y_desc
        gene_x_names = get_gene_names(c_or_mc_ad, gene_x_name)
        gene_names.extend(gene_x_names)
        if expr_vec_infos is not None:
            gene_x_mask = get_gene_mask(c_or_mc_ad, gene_x_name)
            expr_vec_infos = [{
                **x,
                gene_x_desc: np.log2(x['norm_expr'][gene_x_mask].sum() + epsilon_to_add_to_fraction_before_log),
            } for x in expr_vec_infos]
    if y_is_metadata:
        gene_y_vals = c_or_mc_ad.obs.loc[c_or_mc_mask, gene_y_name].to_numpy()
        gene_y_desc = gene_y_name
    else:
        gene_y_vals, gene_y_desc = get_gene_or_genes_log_norm_or_downsampled_expression_and_desc(
            c_or_mc_ad, gene_name_or_names=gene_y_name, 
            jitter_sd=jitter_sd, 
            mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log,
            c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log,
            layer_name=layer_name,
            mc_or_cell_mask=c_or_mc_mask,
        )

        gene_y_names = get_gene_names(c_or_mc_ad, gene_y_name)
        gene_names.extend(gene_y_names)
        if expr_vec_infos is not None:
            gene_y_mask = get_gene_mask(c_or_mc_ad, gene_y_name)
            expr_vec_infos = [{
                **x,
                gene_y_desc: np.log2(x['norm_expr'][gene_y_mask].sum() + epsilon_to_add_to_fraction_before_log),
            } for x in expr_vec_infos]

    assert gene_x_desc != gene_y_desc

    if (not x_is_metadata) and (not y_is_metadata):
        common_genes = set(gene_x_names) & set(gene_y_names)
        if common_genes:
            common_genes = sorted(common_genes)
            if allow_common_genes:
                print(f'WARNING: common genes between gene_x and gene_y: {common_genes}')
            else:
                raise RuntimeError(f'common genes between gene_x and gene_y: {common_genes}')

    dict_for_scatter = {
        'name': c_or_mc_ad.obs_names[c_or_mc_mask],
        gene_x_desc: gene_x_vals,
        gene_y_desc: gene_y_vals,
    }
    if color_by is not None:
        dict_for_scatter[color_by_name] = color_by
    

    df_for_scatter = pd.DataFrame(dict_for_scatter)
    
    second_mc_or_cell_mask = np.full(len(df_for_scatter), True)
    if min_x_val is not None:
        second_mc_or_cell_mask &= gene_x_vals >= min_x_val
    if min_y_val is not None:
        second_mc_or_cell_mask &= gene_x_vals >= min_y_val

    # print(df_for_scatter)
    if gene_names and print_gene_expr_quantile_df:
        gene_expr_df = pd.DataFrame({x: get_genes_expr(c_or_mc_ad, x, c_or_mc_mask=c_or_mc_mask) for x in gene_names})
        gene_expr_quantile_df = gene_expr_df.quantile([0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1], axis=0).T.sort_values(by=0.5)
        gene_expr_quantile_df['0.9-0.1'] = gene_expr_quantile_df[0.9] - gene_expr_quantile_df[0.1]
        # if 1:
        #     gene_expr_quantile_df.sort_values('0.9-0.1', ascending=False, inplace=True)
        print('gene_expr_quantile_df')
        print(gene_expr_quantile_df)

    df_for_scatter = df_for_scatter[second_mc_or_cell_mask]
    
    num_of_rows_with_nan = df_for_scatter.isna().any(axis=1).sum()
    if num_of_rows_with_nan:
        print(f'Warning: {num_of_rows_with_nan} rows with NaNs are silently not shown in the scatter plot (if plotting cells, then maybe these are outliers with nan state)')
        # df_for_scatter = df_for_scatter[~df_for_scatter.isna().any(axis=1)]

    palette = None if (color_by is None) else get_palette(c_or_mc_ad, color_by, color_by_name)
    # print(df_for_scatter.sort_values(gene_y_desc))
    # raise
    if ax is None:
        if figsize is None:
            figsize = (10, 6) if legend else (6, 6)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = None

    if hist2d:
        print('NOTE: consider plotting a 1d hist stratified by the other var... (would probably (?) be better than this?)')
        hist2d_kwargs = {}
        if hist2d_log_count:
            hist2d_kwargs['norm'] = matplotlib.colors.LogNorm()

        hist2d_obj = ax.hist2d(
            df_for_scatter[gene_x_desc],
            df_for_scatter[gene_y_desc],
            hist2d_bins,
            **hist2d_kwargs,
        )
        ax.set_xlabel(gene_x_desc)
        ax.set_ylabel(gene_y_desc)
        fig.colorbar(hist2d_obj[3], ax=ax)

        # sb.histplot(
        #     data=df_for_scatter,
        #     x=gene_x_desc,
        #     y=gene_y_desc,
        #     bins=hist2d_bins,
        #     # cmap=plt.cm.jet,
        #     # legend='auto',
        #     cbar=True,
        #     log_scale=hist2d_log_scale,
        #     legend=True,
        # )
    else:
        sb.scatterplot(
            data=df_for_scatter,
            x=gene_x_desc,
            y=gene_y_desc,
            hue=color_by_name,
            # **({} if (color_by_name is None) else dict(hue=color_by_name)),
            palette=palette,
            # hue_norm=(0,0.25),
            ax=ax,
            edgecolor='black',
            legend=legend,
            **scatter_kwargs,
        )
        if mc_to_mark_name_to_color is not None:
            for mc_to_mark_name, color in mc_to_mark_name_to_color.items():
                row = df_for_scatter.loc[df_for_scatter['name'] == mc_to_mark_name]
                assert len(row) == 1
                row = row.iloc[0]
                ax.scatter(
                    row[gene_x_desc],
                    row[gene_y_desc],
                    marker='x',
                    # color='black',
                    color=color,
                    s=150,
                )
        if expr_vec_infos is not None:
            for expr_vec_info in expr_vec_infos:
                ax.scatter(
                    expr_vec_info[gene_x_desc],
                    expr_vec_info[gene_y_desc],
                    marker='x',
                    # color='black',
                    # color=color,
                    s=50,
                )
                if 'text' in expr_vec_info:
                    ax.text(expr_vec_info[gene_x_desc], expr_vec_info[gene_y_desc], expr_vec_info['text'], ha='center', va='bottom')
        
        if fig is not None:
            if legend:
                sb.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
            fig.tight_layout()
    return fig, ax

def get_cells_of_interest_fraction_per_mc(mc_ad, c_ad, c_mask):
    df = generic_utils.merge_preserving_df1_index_and_row_order(
        mc_ad.obs[['grouped', 'metacell_name']], 
        c_ad.obs.loc[c_mask, 'metacell_name'].astype(str).value_counts().reset_index(name='num_of_cells_of_interest'),
        how='left',
    )
    df['num_of_cells_of_interest'].fillna(0, inplace=True)
    return (df['num_of_cells_of_interest'] / df['grouped']).to_numpy()

def get_default_color_by_and_color_by_gene(mc_ad):
    color_by = None
    color_by_gene = None
    for x in ('state', 'projected_type', 'c_state'):
        if (x in mc_ad.obs.columns) and (f'{x}_color' in mc_ad.obs.columns):
            if mc_ad.obs[x].nunique() == 1:
                print(f'NOTE: not choosing {x} for default color as nunique()==1')
                continue

            color_by = x
            break
    if color_by is None:
        color_by_gene = 'AVP'
    return color_by, color_by_gene

def plot_manifold_umap(
        mc_ad, 
        color_by=None, 
        color_by_gene=None,
        color_by_name=None, 
        mc_mask=None, legend=False, func_for_specified_mcs='dont_show', 
        palette=None, hue_norm=None,
        c_ad=None,
        color_by_c_mask=None,
        ax=None,
        # vmin=None, vmax=None, # NOTE: don't uncomment - use hue_norm instead.
        move_legend=False,
        remove_axes_etc=False,
        legend_labels_ordered=None,
        get_palette_kwargs={},
        plot_colorbar=False,
        fig_colorbar_kwargs={},
        **scatter_kwargs,
):
    if (color_by is None) and (color_by_gene is None):
        color_by, color_by_gene = get_default_color_by_and_color_by_gene(mc_ad)
    if color_by_gene:
        assert color_by is None
        color_by, color_by_name = get_gene_or_genes_log_norm_or_downsampled_expression_and_desc(
            mc_ad, gene_name_or_names=color_by_gene, layer_name=None)


    if hasattr(color_by, 'name'):
        color_by_name = color_by.name
    if color_by_name is None:
        assert isinstance(color_by, str)
        color_by_name = color_by
        color_by = mc_ad.obs[color_by_name].copy()
    if color_by_c_mask is not None:
        assert c_ad is not None
        color_by = get_cells_of_interest_fraction_per_mc(mc_ad, c_ad, color_by_c_mask)
        
        # print(pd.Series(color_by).describe())

        num_of_different_vals_to_color_by = pd.Series(color_by).nunique()
        assert num_of_different_vals_to_color_by >= 2
        if pd.Series(color_by).nunique() == 2:
            hue_norm = (color_by.min(), color_by.max())

        color_by_name = 'cells_of_interest_fraction'

    if palette is None:
        palette = get_palette(mc_ad, color_by, color_by_name, **get_palette_kwargs)
    
    # if isinstance(color_by, (pd.Series, np.ndarray)) and (color_by.dtype.name in {'float32', 'float64'}) and (hue_norm is not None):
    #     palette.set_over('cyan')
    #     palette.set_under('magenta')
    #     # out_of_hue_norm_mask = (color_by < hue_norm[0]) | (color_by > hue_norm[1])
    #     # # color_by.loc[out_of_hue_norm_mask] = np.nan
    #     # color_by.loc[out_of_hue_norm_mask] = np.inf
    #     # # color_by.loc[out_of_hue_norm_mask] = 2

    df_for_scatter = pd.DataFrame({
        # 'x': mc.ut.get_o_numpy(mc_ad, 'umap_x'),
        'x': mc.ut.get_o_numpy(mc_ad, 'x'),
        # 'y': mc.ut.get_o_numpy(mc_ad, 'umap_y'),
        'y': mc.ut.get_o_numpy(mc_ad, 'y'),
        color_by_name: color_by,
    })

    if mc_mask is not None:
        mc_mask = np.array(mc_mask)

    if (mc_mask is not None) and (func_for_specified_mcs == 'dont_show'):
        df_for_scatter = df_for_scatter.loc[mc_mask]

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.get_figure()
    
    sb.scatterplot(
        data=df_for_scatter,
        x='x',
        y='y',
        hue=color_by_name,
        palette=palette,
        hue_norm=hue_norm,
        ax=ax,
        legend=legend,
        edgecolor='black',
        # s=10,
        plotnonfinite=True,
        **scatter_kwargs,
    )
    if plot_colorbar:
        if hue_norm is None:
            hue_norm = (df_for_scatter[color_by_name].min(), df_for_scatter[color_by_name].max())
        norm = plt.Normalize(*hue_norm)
        sm = plt.cm.ScalarMappable(cmap=palette.name, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, **fig_colorbar_kwargs)
        
    ax.set_xlabel(None)
    ax.set_ylabel(None)

    # if you are looking for generic code, then:
    # if ax.get_legend() is not None:
    #     ax.get_legend().get_frame().set_facecolor('none')

    if legend != False:
        if legend_labels_ordered:
            handles, labels = ax.get_legend_handles_labels()
            sorted_handles_and_labels = sorted(list(zip(handles, labels)), key=lambda x: legend_labels_ordered.index(x[1]))
            ax.legend([x[0] for x in sorted_handles_and_labels], [x[1] for x in sorted_handles_and_labels])
        ax.get_legend().get_frame().set_facecolor('none')
        if move_legend:
            sb.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

    
    if (mc_mask is not None) and (func_for_specified_mcs == 'show_name'):
        for _, row in mc_ad.obs.loc[mc_mask].iterrows():
            ax.text(row['umap_x'], row['umap_y'], row['metacell'])

    if remove_axes_etc:
        generic_utils.make_all_spines_and_x_and_y_axes_invisible(ax)
    fig.tight_layout()
    return fig, ax

def get_median_log_enrichment_across_mcs_df(mc_ad, metacell_mask, genes_mask=None):
    return pd.DataFrame({
        'median_log_enrichment': np.median(
            mc_ad.layers['log_enrichment'][metacell_mask,:] if (genes_mask is None)
            else mc_ad.layers['log_enrichment'][np.ix_(metacell_mask,genes_mask)], axis=0), 
        'gene': mc_ad.var_names if (genes_mask is None) else mc_ad.var_names[genes_mask],
    }).sort_values('median_log_enrichment', ascending=False)

def get_cell_type_to_median_log_enrichment_across_mcs_df(mc_ad, cell_state_column_name, return_a_single_df=True):
    cell_type_to_median_log_enrichment_across_mcs_df = {}
    for cell_type in mc_ad.obs[cell_state_column_name].unique():
        cell_type_metacell_mask = mc_ad.obs[cell_state_column_name] == cell_type
        median_log_enrichment_across_mcs_df = get_median_log_enrichment_across_mcs_df(mc_ad, cell_type_metacell_mask)
        median_log_enrichment_across_mcs_df = median_log_enrichment_across_mcs_df.loc[
            median_log_enrichment_across_mcs_df['median_log_enrichment'] > 0]
        median_log_enrichment_across_mcs_df['rank_for_cell_type'] = list(range(len(median_log_enrichment_across_mcs_df)))
        cell_type_to_median_log_enrichment_across_mcs_df[cell_type] = median_log_enrichment_across_mcs_df
    
    if return_a_single_df:
        return pd.concat([curr_df.assign(**{cell_state_column_name: cell_type}) for cell_type, curr_df in cell_type_to_median_log_enrichment_across_mcs_df.items()], ignore_index=True)
    else:
        return cell_type_to_median_log_enrichment_across_mcs_df

def write_mcs_non_rare_gene_similarity_and_linkage(
        mc_ad,
        output_file_path_mcs_with_gene_rho,
        epsilon_to_add_to_fraction_before_log=1e-5,
):
    # comment from compute_var_var_similarity:
    # Compute the gene-gene (variable-variable) similarity matrix. Note by default this will use
    #    {compute_var_var_similarity} which aren't the normal defaults for ``compute_var_var_similarity``, in order to
    #    keep just the top correlated genes and bottom (anti-)correlated genes for each gene. Otherwise you will get a
    #    dense matrix of ~X0K by ~X0K entries, which typically isn't what you want.


    metacell_gene_umis = mc.ut.to_numpy_matrix(mc_ad.X[:,mc_ad.var['non_rare_gene']], copy=True)
    total_umis_of_mcs = mc.ut.get_o_numpy(mc_ad, name='__x__', sum=True)
    # TODO: is it really how I should normalize?
    metacell_gene_umis_normalized = metacell_gene_umis / total_umis_of_mcs[:, np.newaxis]
    # assert np.allclose(metacell_gene_umis_normalized.sum(axis=1), 1) # Not true because i take only non_rare_genes. not sure that's the normalization i should do, though.

    log_metacell_gene_umis_normalized = np.log2(metacell_gene_umis_normalized + epsilon_to_add_to_fraction_before_log)
    # maybe using mc.ut.to_layout() to switch the layout to columns would have been better here??
    gene_rho = mc.ut.corrcoef(log_metacell_gene_umis_normalized.transpose().copy(), per='row', reproducible=True)
    mc.ut.set_m_data(mc_ad, "non_rare_gene_similarity", gene_rho)
    # gene_rho_pdist = scipy.spatial.distance.pdist(gene_rho)
    gene_rho_linkage = scipy.cluster.hierarchy.ward(gene_rho)
    mc.ut.set_m_data(mc_ad, "non_rare_gene_similarity_linkage", gene_rho_linkage)

    mc_ad.write(output_file_path_mcs_with_gene_rho)

def get_num_of_cells_per_metacell(mc_ad, c_ad, c_mask=None):
    if c_mask is None:
        c_mask = np.full(c_ad.n_obs, True)

    return generic_utils.merge_preserving_df1_index_and_row_order(
        mc_ad.obs[['metacell_i']], 
        c_ad.obs.loc[c_mask, 'metacell_i'].value_counts().reset_index(name='num_of_cells').rename(columns={'index': 'metacell_i'}),
        how='left',
    )['num_of_cells'].fillna(0).to_numpy()

def get_total_num_of_non_excluded_umis_per_metacell(mc_ad, c_ad):
    return generic_utils.merge_preserving_df1_index_and_row_order(
        mc_ad.obs[['metacell_i']], 
        c_ad.obs.groupby('metacell_i')['num_of_non_excluded_umis'].sum().reset_index(name='total_num_of_non_excluded_umis').rename(columns={'index': 'metacell_i'}),
    )['total_num_of_non_excluded_umis'].to_numpy()

def add_mc_str_metadata_mode(mc_ad, c_ad, col_name, mode_col_name, c_mask=None):
    if c_mask is None:
        c_mask = np.full(c_ad.n_obs, True)
    
    c_mask &= c_ad.obs['metacell_name'] != 'Outliers'

    if mode_col_name in mc_ad.obs.columns:
        print(f'Warning: {mode_col_name} already in mc_ad.obs.columns, dropping it.')
        mc_ad.obs.drop(columns=mode_col_name, inplace=True)

    mc_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(
        mc_ad.obs,
        c_ad.obs.loc[c_mask, ['metacell_name', col_name]].value_counts().reset_index(name='count').drop_duplicates(
            subset=['metacell_name'], keep='first').rename(columns={col_name: mode_col_name})[['metacell_name', mode_col_name]],
        how='left',
    )
    mc_ad.obs[mode_col_name] = mc_ad.obs[mode_col_name].astype(str)

@contextlib.contextmanager
def protect_from_convey_obs_to_group_context_manager(c_ad, col):
    # ugh. convey_obs_to_group changes the column someway, such that c_ad.obs[col_name] == 'a' raises `ValueError: buffer source array is read-only`. maybe by the following chain? ut.get_o_numpy -> _get_o_data -> _get_shaped_data -> _fix_data or utt.freeze (freeze sounds more likely?)
    orig_col_series = c_ad.obs[col].copy() 
    try:
        yield
    finally:
        c_ad.obs[col] = orig_col_series

def add_mc_metadata(
        mc_ad,
        c_ad,
        statistic_type_to_column_names_in_cells_ad_obj,
        epsilon_for_log_enrichment=1/16,
        min_dominant_donor_fraction=None,
        add_total_num_of_non_excluded_umis=True,
        allow_overwriting=False,
):
    orig_metacell_obs_column_names = set(mc_ad.obs)

    assert mc_ad_and_c_ad_mc_names_match(mc_ad, c_ad)
    # mc_ad.obs['num_of_cells'] = get_num_of_cells_per_metacell(mc_ad, c_ad) # NOTE: use "grouped" instead
    if add_total_num_of_non_excluded_umis:
        mc_ad.obs['total_num_of_non_excluded_umis'] = get_total_num_of_non_excluded_umis_per_metacell(mc_ad, c_ad)

    # assert mc.__version__ == '0.9.0-dev.1'
    # no need for discarding cells with metacell=-1 (and it would also require creating a new AnnData) that because convey_obs_to_group uses range(gdata.n_obs) as indices...
    # if mc.__version__ == '0.9.0-dev.1':
    #     cells_df = cells_df[cells_df['metacell'] != -1].copy() # discard cells without mcs
    #     cells_df['metacell_name'] = mc_obs_names.to_numpy()[cells_df['metacell']]


    # # removes columns whose values are identical for all rows.
    # if not statistic_type_to_column_names_in_cells_ad_obj:
    #     # TODO: do we really want this?
    #     raise RuntimeError('do we really want this?')
    #     cells_df = c_ad.obs
    #     for column_name in list(cells_df):
    #         # print(f'column_name: {column_name}')
    #         num_of_unique_vals = cells_df[column_name].nunique()
    #         assert num_of_unique_vals >= 0
    #         if num_of_unique_vals <= 1:
    #             cells_df.drop(column_name, axis=1, inplace=True)
    #     print(f'relevant column names in c_ad.obs: {list(cells_df)}')


    assert {x for x in statistic_type_to_column_names_in_cells_ad_obj.keys()} <= {
        'bool_fraction', 'nan_mean', 'category_enrichment_and_mode_and_more', 'nan_mean_and_median_and_std', 'category_fraction_diff', 'non_nan_count',
        'str_mode_mode2_and_fracs_and_nunique',
    }
    print('starting to work on category_fraction_and_mode_and_nunique')
    
    column_names_for_max_and_suspicious_fractions = []
    if mc.__version__ in ('0.9.0-dev.1', '0.9.0', '0.9.1-dev', '0.9.1', '0.9.4dev', '0.9.4'):
        # at least in '0.9.0-dev.1', the 'metacell' column in c_ad.obs specifies the index of the metacell, which is what 
        # convey_obs_to_group expects.
    
        if 'bool_fraction' in statistic_type_to_column_names_in_cells_ad_obj:
            for col in statistic_type_to_column_names_in_cells_ad_obj['bool_fraction']:
                print(f'column_name: {col}')
                assert not c_ad.obs[col].isna().any()
                with protect_from_convey_obs_to_group_context_manager(c_ad, col):
                    to_property_name = f'{col}_frac'
                    if not allow_overwriting:
                        assert to_property_name not in orig_metacell_obs_column_names, to_property_name
                    mc.tl.convey_obs_to_group(
                        adata=c_ad,
                        gdata=mc_ad,
                        group='metacell', # change this in other mcs version?
                        property_name=col,
                        to_property_name=to_property_name,
                        method=np.mean,
                    )
        
        for statistic_type in ('str_mode_mode2_and_fracs_and_nunique', 'category_enrichment_and_mode_and_more'):
            if statistic_type in statistic_type_to_column_names_in_cells_ad_obj:
                for col in statistic_type_to_column_names_in_cells_ad_obj[statistic_type]:
                    print(f'column_name: {col}')
                    assert not c_ad.obs[col].isna().any()

                    with protect_from_convey_obs_to_group_context_manager(c_ad, col):
                        to_property_name = f'{col}_nunique'
                        if not allow_overwriting:
                            assert to_property_name not in orig_metacell_obs_column_names, to_property_name
                        mc.tl.convey_obs_to_group(
                            adata=c_ad,
                            gdata=mc_ad,
                            group='metacell', # change this in other mcs version?
                            property_name=col,
                            to_property_name=to_property_name,
                            method=(lambda vec: len(np.unique(vec))),
                        )
        
        if 'str_mode_mode2_and_fracs_and_nunique' in statistic_type_to_column_names_in_cells_ad_obj:
            for col in statistic_type_to_column_names_in_cells_ad_obj['str_mode_mode2_and_fracs_and_nunique']:
                print(f'column_name: {col}')
                assert not c_ad.obs[col].isna().any()
                
                mode_col_name = f'{col}_mode'
                mode_frac_col_name = f'{col}_mode_frac'
                mode2_col_name = f'{col}_mode2'
                mode2_frac_col_name = f'{col}_mode2_frac'
                
                if not allow_overwriting:
                    assert mode_col_name not in orig_metacell_obs_column_names, mode_col_name
                    assert mode_frac_col_name not in orig_metacell_obs_column_names, mode_frac_col_name
                    assert mode2_col_name not in orig_metacell_obs_column_names, mode2_col_name
                    assert mode2_frac_col_name not in orig_metacell_obs_column_names, mode2_frac_col_name
                mc_ad.obs.drop(columns=[mode_col_name, mode_frac_col_name, mode2_col_name, mode2_frac_col_name], inplace=True, errors='ignore')

                # with protect_from_convey_obs_to_group_context_manager(c_ad, col):
                # mc.tl.convey_obs_to_group(
                #     adata=c_ad,
                #     gdata=mc_ad,
                #     group='metacell', # change this in other mcs version?
                #     property_name=col_name,
                #     to_property_name=mode_col_name,
                #     method=mc.ut.most_frequent,
                # )
                add_mc_str_metadata_mode(mc_ad, c_ad, col, mode_col_name)

                # didn't manage to quickly make it work using existing functions, such as scipy.stats.mode or lambda x: pd.Series(x).mode(),
                obs_with_mcs = generic_utils.merge_preserving_df1_index_and_row_order(
                    c_ad.obs[['metacell_name', col]], mc_ad.obs[['metacell_name', mode_col_name]], how='left')
                # c_ad.obs['metacell_name'] != 'Outliers'
                # c_ad_without_mode = mc.ut.slice(c_ad, obs=(obs_with_mcs[col_name] != obs_with_mcs[mode_col_name]))
                # mc.tl.convey_obs_to_group(
                #     adata=c_ad_without_mode,
                #     gdata=mc_ad,
                #     group='metacell', # change this in other mcs version?
                #     property_name=col_name,
                #     to_property_name=mode2_col_name,
                #     method=lambda vec: mc.ut.most_frequent(vec) if len(vec) > 0 else 'nan',
                # )
                add_mc_str_metadata_mode(mc_ad, c_ad, col, mode2_col_name, c_mask=(obs_with_mcs[col] != obs_with_mcs[mode_col_name]))
                obs_with_mcs = generic_utils.merge_preserving_df1_index_and_row_order(obs_with_mcs, mc_ad.obs[['metacell_name', mode2_col_name]], how='left')
                mc_state_mode_df = generic_utils.merge_preserving_df1_index_and_row_order(
                    mc_ad.obs[['metacell_name', 'grouped']], [
                        obs_with_mcs[obs_with_mcs[col] == obs_with_mcs[mode_col_name]].groupby('metacell_name').size().reset_index(name='mode_count'),
                        obs_with_mcs[obs_with_mcs[col] == obs_with_mcs[mode2_col_name]].groupby('metacell_name').size().reset_index(name='mode2_count'),
                    ], how='left',
                )
                assert mc_state_mode_df['mode_count'].notna().all()
                mc_state_mode_df['mode2_count'].fillna(0, inplace=True)
                mc_state_mode_df[mode_frac_col_name] = mc_state_mode_df['mode_count'] / mc_state_mode_df['grouped']
                mc_state_mode_df[mode2_frac_col_name] = mc_state_mode_df['mode2_count'] / mc_state_mode_df['grouped']
                mc_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(mc_ad.obs, mc_state_mode_df[[
                    'metacell_name', mode_frac_col_name, mode2_frac_col_name]])

                
        if 'category_enrichment_and_mode_and_more' in statistic_type_to_column_names_in_cells_ad_obj:
            for col in statistic_type_to_column_names_in_cells_ad_obj['category_enrichment_and_mode_and_more']:
                print(f'column_name: {col}')
                assert not c_ad.obs[col].isna().any()
                
                with protect_from_convey_obs_to_group_context_manager(c_ad, col):
                    to_property_name = f'{col}_mode'
                    if not allow_overwriting:
                        assert to_property_name not in orig_metacell_obs_column_names, to_property_name
                    mc.tl.convey_obs_to_group(
                        adata=c_ad,
                        gdata=mc_ad,
                        group='metacell', # change this in other mcs version?
                        property_name=col,
                        to_property_name=to_property_name,
                        method=mc.ut.most_frequent,
                    )
                    
                    column_vals = c_ad.obs[col].unique()
                    print(f'len(column_vals): {len(column_vals)}')

                    fraction_column_names = []
                    for curr_val in column_vals:
                        # print(f'curr_val: {curr_val}')
                        statistic_column_name = f'{col}_{curr_val}_frac'
                        if not allow_overwriting:
                            assert statistic_column_name not in orig_metacell_obs_column_names
                        mc.tl.convey_obs_to_group(
                            adata=c_ad,
                            gdata=mc_ad,
                            group='metacell', # change this in other mcs version?
                            property_name=col,
                            to_property_name=statistic_column_name,
                            method=mc.ut.fraction_of_grouped(curr_val),
                        )
                        # mc_ad.obs[statistic_column_name] = mc_ad.obs[statistic_column_name].to_numpy()
                        
                        fraction_column_names.append(statistic_column_name)

                    if fraction_column_names:
                        curr_max_fraction_column_name = f'max_{col}_frac'
                        if not allow_overwriting:
                            assert curr_max_fraction_column_name not in orig_metacell_obs_column_names
                        mc_ad.obs[curr_max_fraction_column_name] = mc_ad.obs[fraction_column_names].max(axis=1)
                    
                    log_enrichment_column_names = []
                    for curr_val, expected_fraction in c_ad.obs.loc[c_ad.obs['metacell_i'] >= 0, col].value_counts(
                        normalize=True, 
                        # dropna=False, # if we can't assert nothing is nan, we can use this and everything would be fine, i think
                    ).items():
                        log_enrichment_column_name = f'{col}_{curr_val}_log_enrich'
                        with generic_utils.allow_numpy_division_by_zero_context_manager():
                            if not allow_overwriting:
                                assert log_enrichment_column_name not in orig_metacell_obs_column_names
                            mc_ad.obs[log_enrichment_column_name] = np.log2(
                                epsilon_for_log_enrichment + (mc_ad.obs[f'{col}_{curr_val}_frac'] / expected_fraction))
                        
                        log_enrichment_column_names.append(log_enrichment_column_name)

                    if log_enrichment_column_names:
                        curr_max_log_enrichment_column_name = f'max_{col}_log_enrich'
                        if not allow_overwriting:
                            assert curr_max_log_enrichment_column_name not in orig_metacell_obs_column_names
                        mc_ad.obs[curr_max_log_enrichment_column_name] = mc_ad.obs[log_enrichment_column_names].max(axis=1)
                        mc_ad.obs[f'{col}_with_{curr_max_log_enrichment_column_name}'] = mc_ad.obs[log_enrichment_column_names].idxmax(axis=1).str.slice(
                            len(col) + 1, -len('_log_enrich'))

                

        if 'category_fraction_diff' in statistic_type_to_column_names_in_cells_ad_obj:
            for column1_name, column2_name in statistic_type_to_column_names_in_cells_ad_obj['category_fraction_diff']:
                fraction_diff_column_name = f'max_{column1_name}_and_{column2_name}_fraction_diff'
                if not allow_overwriting:
                    assert fraction_diff_column_name not in orig_metacell_obs_column_names
                mc_ad.obs[fraction_diff_column_name] = mc_ad.obs[
                    f'max_{column1_name}_frac'] - mc_ad.obs[f'max_{column2_name}_frac']
        
        if 'non_nan_count' in statistic_type_to_column_names_in_cells_ad_obj:
            for col in statistic_type_to_column_names_in_cells_ad_obj['non_nan_count']:
                to_property_name = f'{col}_non_nan_count'
                if not allow_overwriting:
                    assert to_property_name not in orig_metacell_obs_column_names, to_property_name
                with protect_from_convey_obs_to_group_context_manager(c_ad, col):
                    mc.tl.convey_obs_to_group(
                        adata=c_ad,
                        gdata=mc_ad,
                        group='metacell', # change this in other mcs version?
                        property_name=col,
                        to_property_name=to_property_name,
                        method=(lambda vec: (~np.isnan(vec)).sum()),
                    )

        for statistic_type in ('nan_mean_and_median_and_std', 'nan_mean'):
            if statistic_type in statistic_type_to_column_names_in_cells_ad_obj:
                print(f'starting to work on {statistic_type}')
                for col in statistic_type_to_column_names_in_cells_ad_obj[statistic_type]:
                    print(f'column_name: {col}')
                    
                    # print(c_ad.obs[column_name].describe())
                    # assert not c_ad.obs[column_name].isna().any()
                    # assert np.isfinite(c_ad.obs[column_name]).all()
                    assert (
                        np.isfinite(c_ad.obs[col])
                        | c_ad.obs[col].isna()
                    ).all()

                    with protect_from_convey_obs_to_group_context_manager(c_ad, col):

                        to_property_name = f'{col}_mean'
                        if not allow_overwriting:
                            assert to_property_name not in orig_metacell_obs_column_names, to_property_name
                        mc.tl.convey_obs_to_group(
                            adata=c_ad,
                            gdata=mc_ad,
                            group='metacell', # change this in other mcs version?
                            property_name=col,
                            to_property_name=to_property_name,
                            method=np.nanmean,
                        )

                        if statistic_type == 'nan_mean_and_median_and_std':
                            to_property_name = f'{col}_median'
                            if not allow_overwriting:
                                assert to_property_name not in orig_metacell_obs_column_names, to_property_name
                            mc.tl.convey_obs_to_group(
                                adata=c_ad,
                                gdata=mc_ad,
                                group='metacell', # change this in other mcs version?
                                property_name=col,
                                to_property_name=to_property_name,
                                method=np.nanmedian,
                            )
                            
                            to_property_name = f'{col}_std'
                            if not allow_overwriting:
                                assert to_property_name not in orig_metacell_obs_column_names, to_property_name
                            mc.tl.convey_obs_to_group(
                                adata=c_ad,
                                gdata=mc_ad,
                                group='metacell', # change this in other mcs version?
                                property_name=col,
                                to_property_name=to_property_name,
                                method=np.nanstd,
                            )
    else:
        raise NotImplementedError('didnt check whether the same implementation would work for the current mcs version.')

    if ({'max_donor_id_frac', 'donor_id_mode'} <= set(list(mc_ad.obs))) and (min_dominant_donor_fraction is not None):
        assert 'any_dominant_donor' not in orig_metacell_obs_column_names
        mc_ad.obs['any_dominant_donor'] = mc_ad.obs['max_donor_id_frac'] >= min_dominant_donor_fraction
        assert 'any_dominant_donor' not in orig_metacell_obs_column_names
        mc_ad.obs.loc[mc_ad.obs['any_dominant_donor'], 'dominant_donor_id'] = mc_ad.obs.loc[
            mc_ad.obs['any_dominant_donor'], 'donor_id_mode']
        mc_ad.obs.loc[~mc_ad.obs['any_dominant_donor'], 'dominant_donor_id'] = 'nan'
        if mc_ad.obs['any_dominant_donor'].any():
            assert mc_ad.obs[mc_ad.obs['any_dominant_donor']].apply(
                lambda row: row[f'donor_id_{row["dominant_donor_id"]}_frac'] >= min_dominant_donor_fraction, axis=1).all()

def write_mcs_metadata(
        input_file_path_mc_ad,
        input_file_path_cells_ad, # aka clean_data
        **kwargs,
):
    print(f'reading {input_file_path_mc_ad}')
    mc_ad = ad.read_h5ad(input_file_path_mc_ad)
    print(f'reading {input_file_path_cells_ad}')
    c_ad = ad.read_h5ad(input_file_path_cells_ad)
    print(f'done reading {input_file_path_cells_ad}')
    add_mc_metadata(
        mc_ad=mc_ad,
        c_ad=c_ad,
        **kwargs,
    )

    

def write_gene_umis_per_cell_and_metacell(
        input_file_path_mc_ad,
        input_file_path_cells_ad, # aka clean_data
        output_file_path_gene_umis_per_cell_and_metacell_csv,
):    
    mc_ad = ad.read_h5ad(input_file_path_mc_ad)
    c_ad = ad.read_h5ad(input_file_path_cells_ad)

    assert mc_ad.n_vars == c_ad.n_vars
    assert (mc_ad.var_names == c_ad.var_names).all()

    total_umis_of_cells = mc.ut.get_o_numpy(c_ad, name='__x__', sum=True)
    cell_normalized_umis = np.asarray(c_ad.X / total_umis_of_cells[:, np.newaxis])
    cell_normalized_umis_masked = np.ma.masked_equal(cell_normalized_umis, 0)
    zero_ignored_mean_cell_normalized_umis_of_genes = cell_normalized_umis_masked.mean(axis=0).data
    zero_ignored_median_cell_normalized_umis_of_genes = np.ma.median(cell_normalized_umis_masked, axis=0).data
    zero_ignored_std_cell_normalized_umis_of_genes = cell_normalized_umis_masked.std(axis=0).data


    total_umis_of_mcs = mc.ut.get_o_numpy(mc_ad, name='__x__', sum=True)
    metacell_normalized_umis = mc_ad.X / total_umis_of_mcs[:, np.newaxis]
    metacell_normalized_umis_masked = np.ma.masked_equal(metacell_normalized_umis, 0)
    zero_ignored_mean_metacell_normalized_umis_of_genes = metacell_normalized_umis_masked.mean(axis=0).data
    zero_ignored_median_metacell_normalized_umis_of_genes = np.ma.median(metacell_normalized_umis_masked, axis=0).data
    zero_ignored_std_metacell_normalized_umis_of_genes = metacell_normalized_umis_masked.std(axis=0).data

    umis_per_cell_and_metacell_df = pd.DataFrame({
        'gene': mc_ad.var_names,
        'zero_ignored_mean_cell_normalized_umis': zero_ignored_mean_cell_normalized_umis_of_genes,
        'zero_ignored_median_cell_normalized_umis': zero_ignored_median_cell_normalized_umis_of_genes,
        'zero_ignored_std_cell_normalized_umis': zero_ignored_std_cell_normalized_umis_of_genes,
        'zero_ignored_mean_metacell_normalized_umis': zero_ignored_mean_metacell_normalized_umis_of_genes,
        'zero_ignored_median_metacell_normalized_umis': zero_ignored_median_metacell_normalized_umis_of_genes,
        'zero_ignored_std_metacell_normalized_umis': zero_ignored_std_metacell_normalized_umis_of_genes,
    })

    umis_per_cell_and_metacell_df.sort_values('zero_ignored_median_cell_normalized_umis', ascending=False, inplace=True)
    umis_per_cell_and_metacell_df.to_csv(output_file_path_gene_umis_per_cell_and_metacell_csv, index=False)
    


def write_gene_module_csv(
        module_name_to_gene_names,
        output_file_path_gene_modules_csv,
        modules_to_ignore_names=set(),
):
    records = []
    for module_name, gene_names in module_name_to_gene_names.items():
        module_size = len(gene_names)
        assert module_size >= 1
        if module_size == 1:
            print(f'skipping module {module_name} because it contains a single gene (which will cause an error in mcview (sometimes?))')
        elif module_name not in modules_to_ignore_names:
            records.extend((x, module_name) for x in gene_names)
    df = pd.DataFrame(records, columns=['gene', 'module'])
    df.to_csv(output_file_path_gene_modules_csv, index=False)

def add_mask_cols(
        data_ad, mask_and_info_list, 
        verbose=True, 
        mc_epsilon_to_add_to_fraction_before_log=MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG, 
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        dist_with_threshs_out_dir_path=None,
):
    for mask_name, mask_info in mask_and_info_list:
        if verbose:
            print(f'\n{mask_name}')
        curr_mc_mask = get_matching_rules_mask(
            data_ad, mask_info['rules'],
            mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log, 
            c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log, 
            dist_with_threshs_out_dir_path=dist_with_threshs_out_dir_path, mask_name=mask_name)
        data_ad.obs[mask_name] = curr_mc_mask
        if verbose:
            print(f'{mask_name}: {curr_mc_mask.sum()}')

def get_matching_rules_mask(
        data_ad, list_of_match_all_rules, mask_name, 
        mc_epsilon_to_add_to_fraction_before_log=MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG, 
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        allow_expr_for_c_ad=False, dist_with_threshs_out_dir_path=None,
):
    if dist_with_threshs_out_dir_path:
        print("running plt.close('all') after each plot. sorry")
        pathlib.Path(dist_with_threshs_out_dir_path).mkdir(parents=True, exist_ok=True)
        fig_num = 0
    
    curr_mask = np.full(data_ad.n_obs, True)
    for list_of_match_any_rules in list_of_match_all_rules:
        any_match_vec = np.full(data_ad.n_obs, False)
        
        for match_rule in list_of_match_any_rules:
            assert len(match_rule) in (3,4), str(match_rule)
            match_vec = None
            is_val_list = len(match_rule) == 3
            if is_val_list:
                attr_type, names, vals = match_rule
            else:
                attr_type, names, min_val, max_val = match_rule
            if attr_type == 'expr':
                if not allow_expr_for_c_ad:
                    assert is_mc_ad(data_ad), 'set allow_expr_for_c_ad=True if you are sure you want this for a c_ad'
                curr_val_vec = get_genes_expr(
                    data_ad, names, 
                    mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log, 
                    c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log, 
                )
                curr_val_desc = 'expr: ' + get_gene_desc(data_ad, names)
            else:
                assert attr_type == 'metadata', f'{attr_type} {match_rule}'
                assert isinstance(names, str), str(match_rule) # i.e., a single name
                curr_val_desc = obs_column_name = names
                if obs_column_name not in data_ad.obs.columns:
                    print(f'assuming {obs_column_name} doesnt satisfy its requirement because it is missing from data_ad')
                    curr_val_vec = None
                else:
                    curr_val_vec = data_ad.obs[names].to_numpy()
            
            if curr_val_vec is None:
                match_vec = np.full(data_ad.n_obs, False)
            else:
                if is_val_list:
                    assert isinstance(vals, list), str(match_rule)
                    match_vec = np.isin(curr_val_vec, vals)
                else:
                    match_vec = np.array(
                        (curr_val_vec >= min_val) &
                        (curr_val_vec <= max_val)
                    )
                    if dist_with_threshs_out_dir_path:
                        fig, ax = plt.subplots()
                        sb.histplot(curr_val_vec[curr_mask], ax=ax)
                        ax.set_xlabel(curr_val_desc)
                        for x in (min_val, max_val):
                            if np.isfinite(x):
                                ax.axvline(x, color='red', alpha=0.5)
                        fig.savefig(os.path.join(dist_with_threshs_out_dir_path, f'{mask_name}_{str(fig_num).zfill(2)}.png'))
                        fig_num += 1
                        plt.close('all')
            any_match_vec |= match_vec
        curr_mask &= any_match_vec
    return curr_mask

def assign_cell_state(
        mc_ad, cell_state_and_info_list, plot_expression_scatters_if_possible=False,
        allow_multiple_types=False, verbose=True, 
        mc_epsilon_to_add_to_fraction_before_log=MC_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG, 
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
        state_column_name='state',
        assign_dissimilar_according_to_similar_field=True,
        dist_with_threshs_out_dir_path=None,
):
    mc_ad.obs[state_column_name] = 'state_unassigned'

    if assign_dissimilar_according_to_similar_field:
        if 'similar' in mc_ad.obs.columns:
            mc_ad.obs.loc[~mc_ad.obs['similar'], state_column_name] = 'Dissimilar' # TODO: maybe we want to assign Dissimilar only if also projected_correlation is lower than 0.85 etc?
        else:
            print('mc_ad.obs["similar"] is missing, so skipping assigning Dissimilar.')
    else:
        print('WARNING: assign_dissimilar_according_to_similar_field=False')
    already_assigned_cell_type_mask = np.full(mc_ad.n_obs, False)
    already_assigned_overriding_cell_type_mask = np.full(mc_ad.n_obs, False)

    prev_cell_type = None
    for cell_state, cell_state_info in cell_state_and_info_list:
        if verbose:
            print(f'\n{cell_state}')
        unexpected_attr_types = set(cell_state_info) - {
            'rules', 'override_others', 
            'still_unassigned',
        }
        if unexpected_attr_types:
            print(f'unexpected_attr_types: {unexpected_attr_types}')
            raise RuntimeError(f'unexpected_attr_types: {unexpected_attr_types}')

        curr_mc_mask = get_matching_rules_mask(
            mc_ad, cell_state_info['rules'], 
            mc_epsilon_to_add_to_fraction_before_log=mc_epsilon_to_add_to_fraction_before_log,
            c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log,
            mask_name=cell_state, dist_with_threshs_out_dir_path=dist_with_threshs_out_dir_path,
        )

        if 'still_unassigned' in cell_state_info:
            assert cell_state_info['still_unassigned']
            curr_mc_mask &= ~already_assigned_cell_type_mask
        
        num_of_mcs_passing_all_thresholds = curr_mc_mask.sum()
        if verbose:
            print(f'num_of_mcs_passing_all_thresholds: {num_of_mcs_passing_all_thresholds}')

        mcs_that_were_already_assigned_this_type = mc_ad.obs[state_column_name] == cell_state
        if cell_state != prev_cell_type:
            assert not mcs_that_were_already_assigned_this_type.any(), f'{cell_state} has different info dicts which are not consecutive (so probably a mistake)'

        if verbose:
            print(f'mcs_that_were_already_assigned_this_type.sum(): {mcs_that_were_already_assigned_this_type.sum()}')
        mcs_to_assign_type_that_were_already_assigned_an_overriding_type_mask = (
            curr_mc_mask & already_assigned_overriding_cell_type_mask
            # it's ok to have different info dicts for the same type, and that's ok if their MCs overlap
            & (~mcs_that_were_already_assigned_this_type)
        )

        if 'override_others' in cell_state_info:
            assert cell_state_info['override_others']
            if mcs_to_assign_type_that_were_already_assigned_an_overriding_type_mask.any():
                print(
                    f'mcs with overriding type to which we tried to assign another overriding type and their already assigned types: {mc_ad.obs.loc[mcs_to_assign_type_that_were_already_assigned_an_overriding_type_mask, state_column_name]}')
                raise RuntimeError(f'tried to assign two overriding types to the same metacell (current assignment: {cell_state})')

            mc_ad.obs.loc[curr_mc_mask, state_column_name] = cell_state
            already_assigned_overriding_cell_type_mask[curr_mc_mask] = True
        else:
            if mcs_to_assign_type_that_were_already_assigned_an_overriding_type_mask.any() and verbose:
                print(
                        f'mcs for which we tried to assign type for the second time and their already assigned types which are overriding types, and so we did not override them but skipped instead:\n{mc_ad.obs.loc[mcs_to_assign_type_that_were_already_assigned_an_overriding_type_mask, state_column_name]}')
            curr_mc_mask = curr_mc_mask & (~already_assigned_overriding_cell_type_mask)

            assigned_type_2nd_time_mask = already_assigned_cell_type_mask & curr_mc_mask & (~mcs_that_were_already_assigned_this_type)
            if assigned_type_2nd_time_mask.any():
                if verbose:
                    print(
                        f'mcs for which we tried to assign type for the second time and their already assigned types:\n{mc_ad.obs.loc[assigned_type_2nd_time_mask, state_column_name]}')
                
                if allow_multiple_types:
                    mc_ad.obs.loc[assigned_type_2nd_time_mask, state_column_name] = (
                        mc_ad.obs.loc[assigned_type_2nd_time_mask, state_column_name] + ';' + cell_state)
                    mc_ad.obs.loc[(curr_mc_mask & (~already_assigned_cell_type_mask)), state_column_name] = cell_state
                    print(f'ATTENTION: assigned {state_column_name} for the second time. the assigned (multiple) types: {mc_ad.obs.loc[assigned_type_2nd_time_mask, state_column_name].unique()}')
                else:
                    raise RuntimeError('assigned_type_2nd_time_mask.any()')
            else:
                mc_ad.obs.loc[curr_mc_mask, state_column_name] = cell_state
                # print(mc_ad.obs[curr_metacell_mask])
        
        already_assigned_cell_type_mask[curr_mc_mask] = True
        prev_cell_type = cell_state
    if verbose:
        print(f'all assigned {state_column_name} types: {set(mc_ad.obs[state_column_name].unique())}')

def get_correlations_with_gene_norm_expressions(mc_ad, vec, metacell_mask=None):
    raise RuntimeError('i guess you want get_corr_df')
    if metacell_mask is None:
        metacell_mask = np.full(mc_ad.n_obs, True)
    assert len(vec) == metacell_mask.sum()
    norm_expression = get_norm_expression(mc_ad)
    corr_df = pd.DataFrame([{'gene': gene_name, 'corr': scipy.stats.pearsonr(vec, norm_expression[:,gene_i])[0]} for gene_i, gene_name in enumerate(mc_ad.var_names)]).sort_values(
        'corr', ascending=False)
    return corr_df

def get_gene_log_norm_expression_correlations(mc_ad, gene_name, metacell_mask=None, other_gene_mask=None, other_gene_min_max_log_norm_expression=-15):
    raise RuntimeError('i guess you want get_corr_df')
    if metacell_mask is None:
        metacell_mask = np.full(mc_ad.n_obs, True)
    if other_gene_mask is None:
        other_gene_mask = np.full(mc_ad.n_vars, True)
    

    gene_i = get_gene_index(mc_ad, gene_name)
    
    max_log_norm_expression_of_genes = mc_ad.layers['expr'][metacell_mask].max(axis=0)
    other_gene_indices = [x for x in np.where(
        other_gene_mask & 
        (max_log_norm_expression_of_genes >= other_gene_min_max_log_norm_expression)
    )[0] if x != gene_i]

    corr_df = pd.DataFrame([
        {
            'other_gene_name': mc_ad.var_names[other_gene_i], 
            'corr': scipy.stats.pearsonr(
                mc_ad.layers['expr'][metacell_mask, gene_i], 
                mc_ad.layers['expr'][metacell_mask, other_gene_i])[0],
            'other_gene_max_log_norm_expression': max_log_norm_expression_of_genes[other_gene_i],
        } for other_gene_i in other_gene_indices
    ])
    if corr_df.empty:
        return corr_df
    corr_df['abs_corr'] = corr_df['corr'].abs()
    corr_df.sort_values('abs_corr', ascending=False, inplace=True)
    return corr_df

def get_genes_log_norm_expression_correlations(
        mc_ad, gene_names, metacell_mask=None, other_gene_mask=None, other_gene_min_max_log_norm_expression=-15, min_abs_corr=0.75,
):
    dfs = []
    for i, gene_name in enumerate(sorted(set(gene_names))):
        if i % 50 == 0:
            print(f'i: {i}')
        # print(f'i: {i}')
        curr_corr_df = get_gene_log_norm_expression_correlations(
            mc_ad, 
            gene_name=gene_name,
            metacell_mask=metacell_mask, 
            other_gene_mask=other_gene_mask,
        )
        dfs.append(curr_corr_df[curr_corr_df['abs_corr'] >= min_abs_corr].assign(gene_name=gene_name))
    corr_df = pd.concat(dfs, ignore_index=True)
    return corr_df

def get_gene_log_norm_expression_range(mc_ad, metacell_mask, low_quantile=0, high_quantile=1):
    low_log_norm_expression_of_genes = np.quantile(mc_ad.layers['expr'][metacell_mask, :], low_quantile, axis=0)
    high_log_norm_expression_of_genes = np.quantile(mc_ad.layers['expr'][metacell_mask, :], high_quantile, axis=0)
    return high_log_norm_expression_of_genes - low_log_norm_expression_of_genes

def add_gene_proj_fold_median(mc_ad):
    if 'proj_fold_median' in mc_ad.var:
        print('ad.var already has proj_fold_median. doing nothing')
        return
    mc_ad.var['proj_fold_median'] = np.median(mc_ad.layers['proj_fold'], axis=0)


def get_proj_fold_df(mc_ad, mc_name, epsilon_for_fixed_expr=1e-5, sort_by_abs_norm_proj_fold=False):
    add_proj_fold_layer(mc_ad)
    add_gene_proj_fold_median(mc_ad)
    write_expr_and_expr_enrich(mc_ad)

    mc_mask = mc.ut.to_numpy_vector(mc_ad.obs_names == mc_name)
    df = pd.DataFrame({
        'gene': mc_ad.var_names, 
        'proj_fold': mc.ut.to_numpy_vector(mc_ad.layers['proj_fold'][mc_mask, :]),
        'fixed_expr': np.log2(mc.ut.to_numpy_vector(mc_ad.layers['corrected_fraction'][mc_mask, :]) + epsilon_for_fixed_expr),
        'proj_expr': np.log2(mc.ut.to_numpy_vector(mc_ad.layers['projected_fraction'][mc_mask, :]) + epsilon_for_fixed_expr),
        'expr': mc.ut.to_numpy_vector(mc_ad.layers['expr'][mc_mask, :]),
        'proj_fold_median': mc_ad.var['proj_fold_median'],
    })
    df['abs_proj_fold'] = df['proj_fold'].abs()
    df['norm_proj_fold'] = df['proj_fold'] - df['proj_fold_median']
    df['abs_norm_proj_fold'] = df['norm_proj_fold'].abs()

    sort_by = 'abs_norm_proj_fold' if sort_by_abs_norm_proj_fold else 'abs_proj_fold'
    df.sort_values(sort_by, ascending=False, inplace=True)
    return df

def get_genes_with_high_max_abs_projected_fold_mask(
        mc_ad, 
        metacell_mask, 
        projected_fold=None,
        low_quantile=0, 
        high_quantile=0, 
        min_max_abs_projected_fold=None,
):
    raise RuntimeError('i guess you want get_proj_fold_df()')
    if projected_fold is None:
        projected_fold = get_metacell_projected_fold(mc_ad)
    
    low_projected_fold_of_genes = np.quantile(projected_fold[metacell_mask, :], low_quantile, axis=0)
    high_projected_fold_of_genes = np.quantile(projected_fold[metacell_mask, :], high_quantile, axis=0)

    max_abs_projected_fold_of_genes = np.maximum(np.abs(low_projected_fold_of_genes), np.abs(high_projected_fold_of_genes))
    

    plt.close('all')
    fig, ax = plt.subplots()
    sb.histplot(
        max_abs_projected_fold_of_genes, bins=np.arange(0, max_abs_projected_fold_of_genes.max() + 0.01, 0.05), ax=ax)
    ax.set_xlabel('max_abs_projected_fold')
    if min_max_abs_projected_fold is None:
        raise RuntimeError('specify min_max_abs_projected_fold according to the distribution...')
    ax.axvline(min_max_abs_projected_fold, linestyle='--', color='grey', alpha=0.7)

    return max_abs_projected_fold_of_genes >= min_max_abs_projected_fold

def get_genes_with_high_log_norm_expression_range_mask(mc_ad, metacell_mask, low_quantile=0, high_quantile=1, min_expr_range=None):
    log_norm_expression_range_of_genes = get_gene_log_norm_expression_range(
        mc_ad, metacell_mask, low_quantile=low_quantile, high_quantile=high_quantile)

    plt.close('all')
    fig, ax = plt.subplots()
    sb.histplot(log_norm_expression_range_of_genes, bins=np.arange(0, log_norm_expression_range_of_genes.max() + 0.01, 0.2), ax=ax)
    ax.set_xlabel(f'log_norm_expression_range ({high_quantile} quantile - {low_quantile} quantile)')
    if min_expr_range is None:
        raise RuntimeError('specify min_log_norm_expression_range according to the distribution...')
    ax.axvline(min_expr_range, linestyle='--', color='grey', alpha=0.7)

    return log_norm_expression_range_of_genes >= min_expr_range
    

def get_best_separating_gene_log_norm_expression_given_gene_i(mc_ad, metacell_mask1, metacell_mask2, gene_i):
    assert not (metacell_mask1 & metacell_mask2).any()

    mcs1_log_norm_expression = mc_ad.layers['expr'][metacell_mask1, gene_i]
    mcs2_log_norm_expression = mc_ad.layers['expr'][metacell_mask2, gene_i]
    # print(f'mcs1_log_norm_expression: {mcs1_log_norm_expression}')
    # print(f'mcs2_log_norm_expression: {mcs2_log_norm_expression}')
    return generic_utils.get_best_separator(mcs1_log_norm_expression, mcs2_log_norm_expression)
    

def get_best_separating_gene_log_norm_expression(mc_ad, metacell_mask1, metacell_mask2, gene_name):
    return get_best_separating_gene_log_norm_expression_given_gene_i(mc_ad, metacell_mask1, metacell_mask2, get_gene_index(mc_ad, gene_name))

def get_best_separating_genes_log_norm_expression(
        mc_ad, 
        metacell_mask1, # test
        metacell_mask2, # control
        genes_to_go_over_mask=None,
        min_separation_score=0.8, debug_max_i_to_go_over=None,
):
    num_of_mask1_mcs = metacell_mask1.sum()
    num_of_mask2_mcs = metacell_mask2.sum()
    assert num_of_mask1_mcs > 0
    assert num_of_mask2_mcs >= 2, 'must have at least two mcs in mask2, which is the control'
    # if (metacell_mask1.sum() == 1) or (metacell_mask2.sum() == 1):
    #     print('one of the masks has only one metacell. this would probably cause an exception in scipy.stats.ttest_ind. i think we should just assume the other set of vals are normally distributed, and calculate a z-score instead of t-statistic, and the pval of the z-score. (TODO)')

    all_metacell_mask = metacell_mask1 | metacell_mask2
    
    if genes_to_go_over_mask is None:
        genes_to_go_over_mask = np.full(mc_ad.n_vars, True)
    genes_to_go_over_mask &= (
        mc_ad.layers['expr'][all_metacell_mask,:].max(axis=0) - 
        mc_ad.layers['expr'][all_metacell_mask,:].min(axis=0)
    ) > 0
    
    num_of_genes_to_go_over = genes_to_go_over_mask.sum()
    print(f'num_of_genes_to_go_over: {num_of_genes_to_go_over}')
    flat_dicts = []
    for i, gene_i in enumerate(np.where(genes_to_go_over_mask)[0]):
        if debug_max_i_to_go_over is not None:
            if i > debug_max_i_to_go_over:
                break
        # if i <= 27719: 
        #     continue
        # print(f'i: {i}')
        if i % 1000 == 999:
            print(f'i: {i}')
            # break
        best_separator, separation_score = get_best_separating_gene_log_norm_expression_given_gene_i(mc_ad, metacell_mask1, metacell_mask2, gene_i)
        if separation_score >= min_separation_score:
            vals1 = mc_ad.layers['expr'][metacell_mask1, gene_i]
            vals2 = mc_ad.layers['expr'][metacell_mask2, gene_i]

            if num_of_mask1_mcs > 1:
                t_statistic, welch_t_test_pval = scipy.stats.ttest_ind(
                    vals1,
                    vals2,
                    equal_var=False,
                    alternative='two-sided',
                )
                flat_dict_extras = {
                    'welch_t_test_pval': welch_t_test_pval,
                    't_statistic': t_statistic,
                }
            else:
                assert num_of_mask1_mcs == 1
                
                vals2_std = np.std(vals2)
                if vals2_std == 0:
                    ztest_pval = None
                    z_score = None
                else:
                    z_score, ztest_pval = generic_utils.get_z_score_and_ztest_pval(vals1[0], vals2)
                    flat_dict_extras = {
                        'ztest_pval': ztest_pval,
                        'z_score': z_score,
                    }
            
            curr_flat_dict = {
                'gene': mc_ad.var_names[gene_i],
                'best_separator': best_separator,
                'separation_score': separation_score,
                **flat_dict_extras,
            }
            flat_dicts.append(curr_flat_dict)
    
    if not flat_dicts:
        print(f'weirdly, found no gene with separation_score >= {min_separation_score}')
        return pd.DataFrame(columns=['gene', 'best_separator', 'separation_score'])
    
    separators_df = pd.DataFrame(flat_dicts)
    separators_df['unified_pval'] = separators_df[[x for x in ['welch_t_test_pval', 'ztest_pval'] if x in list(separators_df)]].fillna(1).min(axis=1)
    return separators_df.sort_values(['unified_pval', 'separation_score'], ascending=[True, False]).drop('unified_pval', axis=1)

def filter_cells(
        data_ad, filtering_params, random_seed, given_init_c_mask=None, 
        cells_discarded_due_to_donor_attributes_obs_names_file_path=None,
        cells_discarded_due_to_many_cells_per_donor_obs_names_file_path=None,
        allow_keeping_unassigned_cells=False,
        exp_name_to_min_umi_count=None,
        exp_name_to_is_ultima=None,
        exp_name_to_soup_vireo_mip_triplets_to_discard=None,
        mds_and_mds_mpn_and_cytopenia_diagnosis_classes=None,
        diagnosis_to_diagnosis_class=None,
):
    if given_init_c_mask is not None:
        orig_num_of_cells = data_ad.n_obs
        assert not given_init_c_mask.all()
        data_ad = mc.ut.slice(data_ad, obs=given_init_c_mask)
        discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
        print(f'discarded_cell_fraction (due to given_init_c_mask): {discarded_cell_fraction}')

    if 'only_bm' in filtering_params:
        assert filtering_params['only_bm'] == True
        cells_to_keep_mask = data_ad.obs['is_bm']
        orig_num_of_cells = data_ad.n_obs
        if cells_to_keep_mask.all():
            print('not filtering only_bm because nothing to discard')
        else:
            data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
            discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
            print(f'discarded_cell_fraction (due to not being BM): {discarded_cell_fraction}')
    
    if 'only_pb' in filtering_params:
        assert filtering_params['only_pb'] == True
        cells_to_keep_mask = ~data_ad.obs['is_bm']
        orig_num_of_cells = data_ad.n_obs
        if cells_to_keep_mask.all():
            print('not filtering only_pb because nothing to discard')
        else:
            data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
            discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
            print(f'discarded_cell_fraction (due to not being PB): {discarded_cell_fraction}')
    
    if 'only_ultima' in filtering_params:
        assert filtering_params['only_ultima'] == True
        assert set(data_ad.obs['exp_name'].unique()) <= set(exp_name_to_is_ultima)
        ultima_exp_names = [x for x, is_ultima in exp_name_to_is_ultima.items() if is_ultima]
        cells_to_keep_mask = data_ad.obs['exp_name'].isin(ultima_exp_names)
        orig_num_of_cells = data_ad.n_obs
        if cells_to_keep_mask.all():
            print('not filtering only_ultima because nothing to discard')
        else:
            data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
            discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
            print(f'discarded_cell_fraction (due to not being ultima): {discarded_cell_fraction}')
    if 'only_non_ultima' in filtering_params:
        assert filtering_params['only_non_ultima'] == True
        assert set(data_ad.obs['exp_name'].unique()) <= set(exp_name_to_is_ultima)
        non_uli_exp_names = [x for x, is_ultima in exp_name_to_is_ultima.items() if not is_ultima]
        cells_to_keep_mask = data_ad.obs['exp_name'].isin(non_uli_exp_names)
        orig_num_of_cells = data_ad.n_obs
        if cells_to_keep_mask.all():
            print('not filtering only_non_ultima because nothing to discard')
        else:
            data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
            discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
            print(f'discarded_cell_fraction (due to not being ultima): {discarded_cell_fraction}')
    
    cells_to_keep_mask = np.full(data_ad.n_obs, True)
    if 'cells_to_discard_paths' in filtering_params:
        for cells_to_discard_df_csv_file_path in filtering_params['cells_to_discard_paths']:
            cells_to_discard_df = pd.read_csv(cells_to_discard_df_csv_file_path)
            curr_cells_to_discard_mask = data_ad.obs.index.isin(cells_to_discard_df['cell'])
            if not curr_cells_to_discard_mask.any():
                print(f'WARNING: no cells to discard from {cells_to_discard_df_csv_file_path} (specified in cells_to_discard_paths). make sure this makes sense!')
            else:
                cells_discarded_again_mask = (~cells_to_keep_mask) & curr_cells_to_discard_mask
                assert not cells_discarded_again_mask.any(), f'{cells_to_discard_df_csv_file_path} contained ({cells_discarded_again_mask.sum()}) cells that were already discarded (following are the first 5): {data_ad.obs_names[cells_discarded_again_mask][:5]}' # this means we don't discard the same cell twice, so we don't have to worry about counting the same cell twice in other places, e.g., when we estimate the number of yet unidentified same donor doublets.
                
                cells_to_keep_mask &= ~curr_cells_to_discard_mask
        
    assert not data_ad.obs['donor_id'].isna().any()
    cells_assigned_to_donors_mask = data_ad.obs['donor_id'] != 'nan'

    if 'genotype_doublet' in data_ad.obs.columns:
        assert set(data_ad.obs['genotype_doublet'].unique()) <= {True, False}
        doublet_mask = data_ad.obs['genotype_doublet']
    else:
        doublet_mask = np.full(data_ad.n_obs, False)
    unassigned_mask = (~doublet_mask) & (~cells_assigned_to_donors_mask)
    unassigned_or_doublet_mask = doublet_mask | unassigned_mask
    assert (unassigned_or_doublet_mask == (~cells_assigned_to_donors_mask)).all()

    all_exp_names = set(data_ad.obs['exp_name'].unique())
    all_donor_ids = set(data_ad.obs['donor_id'].unique())

    diagnosis_normal_mask = data_ad.obs['diagnosis'].isin(['normal'])

    if (diagnosis_to_diagnosis_class is not None):
        diagnosis_class_of_cells = data_ad.obs['diagnosis'].map(diagnosis_to_diagnosis_class)
        assert not diagnosis_class_of_cells.isna().any()
    # if (mds_and_mds_mpn_and_cytopenia_diagnosis_classes is not None):
    #     cytopenic_donor_mask = diagnosis_class_of_cells.isin(mds_and_mds_mpn_and_cytopenia_diagnosis_classes)
    
    # assert set(data_ad.obs['diagnosis']) <= {'normal', *unhealthy_diagnosis_values}
    orig_cells_to_keep_mask = cells_to_keep_mask.copy()
    if 'donors_to_keep_attributes' in filtering_params:
        for donor_attribute in filtering_params['donors_to_keep_attributes']:
            if donor_attribute == 'assigned':
                cells_to_keep_mask &= cells_assigned_to_donors_mask
            # if donor_attribute == 'cytopenic':
            #     cells_to_keep_mask &= cytopenic_donor_mask
            if donor_attribute == 'except_other_exp_condition':
                cells_to_keep_mask &= (diagnosis_class_of_cells != 'other_exp_condition')
            if donor_attribute == 'mpn_and_normal':
                cells_to_keep_mask &= (diagnosis_class_of_cells == 'MPN') | diagnosis_normal_mask
            if donor_attribute == 'assigned_and_doublets':
                cells_to_keep_mask &= cells_assigned_to_donors_mask | doublet_mask
            if donor_attribute == 'unhealthy':
                cells_of_unhealthy_donors_mask = (
                    cells_assigned_to_donors_mask &
                    (~diagnosis_normal_mask)
                # ) # TODO: uncomment this line and remove the next one
                )
                cells_to_keep_mask &= cells_of_unhealthy_donors_mask
            if donor_attribute == 'normal':
                cells_of_healthy_donors_mask = cells_assigned_to_donors_mask & diagnosis_normal_mask
                cells_to_keep_mask &= cells_of_healthy_donors_mask
    
    assert not (('specific_donors_to_exclude' in filtering_params) and ('specific_donors_to_include' in filtering_params))
    if 'specific_donors_to_exclude' in filtering_params:
        unexpected_specific_donors_to_exclude = set(filtering_params['specific_donors_to_exclude']) - all_donor_ids
        assert not unexpected_specific_donors_to_exclude, f'unexpected_specific_donors_to_exclude: {unexpected_specific_donors_to_exclude}'
        cells_to_keep_mask &= ~data_ad.obs['donor_id'].isin(filtering_params['specific_donors_to_exclude'])
    if 'specific_donors_to_include' in filtering_params:
        unexpected_specific_donors_to_include = set(filtering_params['specific_donors_to_include']) - all_donor_ids
        assert not unexpected_specific_donors_to_include, f'unexpected_specific_donors_to_include: {unexpected_specific_donors_to_include}'
        cells_to_keep_mask &= data_ad.obs['donor_id'].isin(filtering_params['specific_donors_to_include'])
        assert cells_to_keep_mask.any(), f'no cells belong to {filtering_params["specific_donors_to_include"]}'
    
    assert not (('specific_exps_to_exclude' in filtering_params) and ('specific_exps_to_include' in filtering_params))
    if 'specific_exps_to_exclude' in filtering_params:
        unexpected_specific_exps_to_exclude = set(filtering_params['specific_exps_to_exclude']) - all_exp_names
        if unexpected_specific_exps_to_exclude:
            print(f'unexpected_specific_exps_to_exclude: {unexpected_specific_exps_to_exclude}. make sure this makes sense.')
        cells_to_keep_mask &= ~data_ad.obs['exp_name'].isin(filtering_params['specific_exps_to_exclude'])
    if 'specific_exps_to_include' in filtering_params:
        unexpected_specific_exps_to_include = set(filtering_params['specific_exps_to_include']) - all_exp_names
        assert not unexpected_specific_exps_to_include, f'unexpected_specific_exps_to_include: {unexpected_specific_exps_to_include}'
        cells_to_keep_mask &= data_ad.obs['exp_name'].isin(filtering_params['specific_exps_to_include'])
        assert cells_to_keep_mask.any(), f'no cells of {filtering_params["specific_donors_to_include"]} appeared in {filtering_params["specific_exps_to_include"]}'
    
    if 'specific_exp_donors_to_exclude' in filtering_params:
        for exp_name, donor_id in filtering_params['specific_exp_donors_to_exclude']:
            if (exp_name not in all_exp_names):
                print(f'{exp_name} doesnt exist in data_ad. make sure this makes sense.')
            elif (donor_id not in all_donor_ids):
                print(f'{donor_id} doesnt exist in data_ad. make sure this makes sense.')
            else:
                cells_to_keep_mask &= ~((data_ad.obs['exp_name'] == exp_name) & (data_ad.obs['donor_id'] == donor_id))

    if 'specific_exp_barcodes_to_include' in filtering_params:
        for exp_name, barcode_file_path in filtering_params['specific_exp_barcodes_to_include']:
            # NOTE: maybe better to simply use cells_to_discard_paths?
            raise RuntimeError('implement this')
    
    if cells_discarded_due_to_donor_attributes_obs_names_file_path is not None:
        cells_discarded_due_to_donor_attribute_mask = orig_cells_to_keep_mask & (~cells_to_keep_mask)
        cells_discarded_due_to_donor_attribute_obs_names = list(data_ad.obs_names[cells_discarded_due_to_donor_attribute_mask])
        generic_utils.write_text_file(
            cells_discarded_due_to_donor_attributes_obs_names_file_path, '\n'.join(cells_discarded_due_to_donor_attribute_obs_names))
    
    if (cells_to_keep_mask & unassigned_mask).any() and (not allow_keeping_unassigned_cells):
        raise RuntimeError('after filtering, the data still contains unassigned. this is probably wrong. even if you are calculating mcs with doublets to then discard all cells in mcs that have higher than expected number of doublets, then you probably want to calculate these mcs without unassigned, because you dont use them in any way, and they could cluster well with doublets, so just interfering with your goal of identifying doublets.')

    if cells_to_keep_mask.all():
        print('not filtering because nothing to discard')
    else:
        orig_num_of_cells = data_ad.n_obs
        data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
        discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
        print(f'discarded_cell_fraction (according to cells_to_discard_paths and donors_to_keep_attributes): {discarded_cell_fraction}')

    if exp_name_to_min_umi_count is not None:
        cells_to_keep_mask = np.full(data_ad.n_obs, True)
        for exp_name, min_umi_count in exp_name_to_min_umi_count.items():
            cells_to_keep_mask &= ~((data_ad.obs['exp_name'] == exp_name) & (data_ad.obs['cell_ranger_umi_count'] < min_umi_count))
        
        orig_num_of_cells = data_ad.n_obs
        data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
        discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
        print(f'discarded_cell_fraction (according to exp_name_to_min_umi_count): {discarded_cell_fraction}')
    
    if exp_name_to_soup_vireo_mip_triplets_to_discard is not None:
        cells_to_keep_mask = np.full(data_ad.n_obs, True)
        for exp_name, soup_vireo_mip_triplets_to_discard in exp_name_to_soup_vireo_mip_triplets_to_discard.items():
            for (soup_donor_name, vireo_donor_name, donor_id) in soup_vireo_mip_triplets_to_discard:
                cells_to_keep_mask &= ~(
                    (data_ad.obs['exp_name'] == exp_name)
                    & (data_ad.obs['soup_donor_name'] == soup_donor_name)
                    & (data_ad.obs['vireo_donor_name'] == vireo_donor_name)
                    & (data_ad.obs['donor_id'] == donor_id)
                )
        
        orig_num_of_cells = data_ad.n_obs
        data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
        discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
        print(f'discarded_cell_fraction (according to exp_name_to_soup_vireo_mip_triplets_to_discard): {discarded_cell_fraction}')
            
    if 'max_num_of_cells_to_keep_as_a_factor_of_median' in filtering_params:
        max_num_of_cells_to_keep_as_a_factor_of_median = filtering_params['max_num_of_cells_to_keep_as_a_factor_of_median']
        median_num_of_cells_per_donor = int(data_ad.obs['donor_id'].value_counts().median())
        max_num_of_cells_per_donor = int(median_num_of_cells_per_donor * max_num_of_cells_to_keep_as_a_factor_of_median)

        cells_to_keep_mask = np.full(data_ad.n_obs, True)
        random.seed(random_seed)
        for donor_id, donor_num_of_cells in data_ad.obs['donor_id'].value_counts().items():
            if donor_num_of_cells > max_num_of_cells_per_donor:
                print(f'{donor_id} has a very high number of cells ({donor_num_of_cells}), so randomly keeping only {max_num_of_cells_per_donor} (median * {max_num_of_cells_to_keep_as_a_factor_of_median}) of them, to avoid {donor_id} dominating mcs')
                donor_cell_indices = list(np.where(data_ad.obs['donor_id'] == donor_id)[0])
                random_cells_to_discard_indices = np.array(random.sample(donor_cell_indices, donor_num_of_cells - max_num_of_cells_per_donor))
                cells_to_keep_mask[random_cells_to_discard_indices] = False

        orig_num_of_cells = data_ad.n_obs
        cells_discarded_due_to_many_cells_per_donor_obs_names = list(data_ad.obs_names[~cells_to_keep_mask])
        if cells_discarded_due_to_many_cells_per_donor_obs_names_file_path is not None:
            generic_utils.write_text_file(
                cells_discarded_due_to_many_cells_per_donor_obs_names_file_path, '\n'.join(cells_discarded_due_to_many_cells_per_donor_obs_names))
        data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
        discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
        print(f'discarded_cell_fraction (due to very high number of cells per donor): {discarded_cell_fraction}')
        assert data_ad.obs['donor_id'].value_counts().max() <= max_num_of_cells_per_donor
    
        # draft for doing max per donor bleeding rather than donor. copy paste to notebook to keep going.
        # curr_unfiltered_cells_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(
        #     curr_unfiltered_cells_ad.obs, 
        #     sc_rna_seq_preprocessing_params.get_nili_experiment_date_and_exp_name_df().rename(columns={'experiment date': 'bleeding_date'}), 
        # )
        # max_num_of_cells_to_keep_as_a_factor_of_median = 3
        # data_ad = curr_unfiltered_cells_ad
        # random_seed = 1234

        # median_num_of_cells_per_donor_bleeding = int(data_ad.obs[['donor_id', 'bleeding_date']].value_counts().median())
        # max_num_of_cells_per_donor_bleeding = int(median_num_of_cells_per_donor_bleeding * max_num_of_cells_to_keep_as_a_factor_of_median)

        # cells_to_keep_mask = np.full(data_ad.n_obs, True)
        # random.seed(random_seed)
        # for (donor_id, bleeding_date), donor_bleeding_num_of_cells in data_ad.obs[['donor_id', 'bleeding_date']].value_counts().items():
        #     if donor_bleeding_num_of_cells > max_num_of_cells_per_donor_bleeding:
        #         print(f'{donor_id},{bleeding_date} has a very high number of cells ({donor_bleeding_num_of_cells}), so randomly keeping only {max_num_of_cells_per_donor_bleeding} (median * {max_num_of_cells_to_keep_as_a_factor_of_median}) of them, to avoid {donor_id} dominating mcs')
        #         donor_cell_indices = list(np.where(data_ad.obs['donor_id'] == donor_id)[0])
        #         random_cells_to_discard_indices = np.array(random.sample(donor_cell_indices, donor_bleeding_num_of_cells - max_num_of_cells_per_donor_bleeding))
        #         cells_to_keep_mask[random_cells_to_discard_indices] = False

        # orig_num_of_cells = data_ad.n_obs
        # cells_discarded_due_to_many_cells_per_donor_obs_names = list(data_ad.obs_names[~cells_to_keep_mask])
        # if cells_discarded_due_to_many_cells_per_donor_obs_names_file_path is not None:
        #     generic_utils.write_text_file(
        #         cells_discarded_due_to_many_cells_per_donor_obs_names_file_path, '\n'.join(cells_discarded_due_to_many_cells_per_donor_obs_names))
        # data_ad = mc.ut.slice(data_ad, obs=cells_to_keep_mask)
        # discarded_cell_fraction = 1 - (data_ad.n_obs / orig_num_of_cells)
        # print(f'discarded_cell_fraction (due to very high number of cells per donor): {discarded_cell_fraction}')
        # assert data_ad.obs['donor_id'].value_counts().max() <= max_num_of_cells_per_donor_bleeding

    return data_ad

def remove_metacell_etc_columns_from_cells_ad_obs(c_ad):
    c_ad.obs.drop([x for x in CELLS_OBS_COLUMN_NAMES_ADDED_BY_METACELL_AND_DOWNSTREAM_ANALYSIS if x in list(c_ad.obs)], axis=1, inplace=True)

def mc_ad_and_c_ad_mc_names_match(mc_ad, c_ad):
    mc_ad_names = set(mc_ad.obs['metacell_name'].astype(str))
    c_ad_names = set(c_ad.obs['metacell_name'].astype(str))
    assert 'Outliers' not in mc_ad_names
    return (c_ad_names == mc_ad_names) if 'Outliers' not in c_ad_names else (c_ad_names == (mc_ad_names | {'Outliers'}))

def is_metacell_col_as_expected(c_ad):
    mc_is = set(c_ad.obs['metacell'])
    # could there also be -2 in some versions??
    non_outlier_mc_is = mc_is - {-1}
    expected_non_outlier_mc_is = set(range(max(non_outlier_mc_is) + 1))
    unexpected_non_outlier_mc_is = non_outlier_mc_is - expected_non_outlier_mc_is
    missing_non_outlier_mc_is = expected_non_outlier_mc_is - non_outlier_mc_is
    if unexpected_non_outlier_mc_is:
        print(f'unexpected_non_outlier_mc_is: {unexpected_non_outlier_mc_is}')
    if missing_non_outlier_mc_is:
        print(f'missing_non_outlier_mc_is: {missing_non_outlier_mc_is}')
    return non_outlier_mc_is == expected_non_outlier_mc_is

def add_metacell_i_and_metacell_name_to_cells_ad(c_ad, mc_ad):
    assert 'metacell_i' not in list(c_ad.obs)
    assert np.issubdtype(c_ad.obs['metacell'].dtype, np.integer)
    assert set(c_ad.obs['metacell'].unique()) - set(range(len(mc_ad))) <= {-2, -1}
    
    cells_with_mcs_mask = c_ad.obs['metacell'] >= 0
    c_ad.obs.loc[cells_with_mcs_mask, 'metacell_i'] = c_ad.obs.loc[cells_with_mcs_mask, 'metacell']

    if 'metacell_name' in c_ad.obs.columns:
        expected_mc_name_mask = c_ad.obs.loc[
            cells_with_mcs_mask, 'metacell_name'
        ].to_numpy() == mc_ad.obs['metacell_name'].iloc[c_ad.obs.loc[cells_with_mcs_mask, 'metacell_i']].to_numpy()
        assert expected_mc_name_mask.all()
        assert (c_ad.obs.loc[~cells_with_mcs_mask, 'metacell_name'] == 'Outliers').all()

    else:
        # c_ad.obs.loc[cells_with_mcs_mask, 'metacell_name'] = 'MC' + c_ad.obs.loc[cells_with_mcs_mask, 'metacell'].astype(str) # replaced it with the next line on 231005
        c_ad.obs.loc[cells_with_mcs_mask, 'metacell_name'] = mc_ad.obs['metacell_name'].iloc[c_ad.obs.loc[cells_with_mcs_mask, 'metacell_i']].astype(str).to_numpy()
        
        c_ad.obs.loc[~cells_with_mcs_mask, 'metacell_name'] = 'Outliers'

    assert generic_utils.set_is_empty_or_contains_only_np_nan(set(c_ad.obs['metacell_i'].unique()) - set(range(len(mc_ad))))
    assert generic_utils.set_is_empty_or_contains_only_np_nan(set(c_ad.obs['metacell_i'].unique()) - set(mc_ad.obs['metacell_i'].unique()))
    # unexpected_mc_names = 
    assert not (set(c_ad.obs['metacell_name'].unique()) - set(mc_ad.obs['metacell_name'].unique()) - {'Outliers'})

def add_metacell_i(mc_ad):
    mc_ad.obs['metacell_i'] = list(range(mc_ad.n_obs))
    assert mc_ad.obs['metacell_i'].is_monotonic_increasing

def add_metacell_i_and_metacell_name(mc_ad, use_index_as_str_as_metacell_name=True):
    add_metacell_i(mc_ad)

    if use_index_as_str_as_metacell_name:
        metacell_names = mc_ad.obs.index.astype(str)
    else:
        metacell_names = [f'MC{i}' for i in range(mc_ad.n_obs)]
    mc_ad.obs.index = metacell_names
    mc_ad.obs['metacell_name'] = metacell_names # for convenience and safety

def fix_mc_indices_after_manual_mc_filtering(c_ad):
    # we need this because, IIUC, we must not have gaps in metacell nums (which makes sense, as they are used as indices)
    old_mc_is = sorted(c_ad.obs['metacell_i'].drop_duplicates().dropna())
    mc_i_map = {old: new for new, old in enumerate(old_mc_is)}
    c_ad.obs['metacell'] = c_ad.obs['metacell'].map(mc_i_map)
    c_ad.obs['metacell'].fillna(-1, inplace=True)
    c_ad.obs['metacell'] = c_ad.obs['metacell'].astype(int)
    assert ((c_ad.obs['metacell'] == -1) == (c_ad.obs['metacell_i'].isna())).all()
    return mc_i_map

def collect_metacells_wrapper(c_ad, add_expr_and_expr_enrich=True, orig_mc_ad=None, ad_attr_and_columns_to_copy_names=None, **kwargs):
    assert 'metacell_i' not in list(c_ad.obs)
    mc_ad = mc.pl.collect_metacells(c_ad, **kwargs)
    add_metacell_i_and_metacell_name(mc_ad)
    add_metacell_i_and_metacell_name_to_cells_ad(c_ad, mc_ad)
    if add_expr_and_expr_enrich:
        write_expr_and_expr_enrich(mc_ad)
    if 'gene_id' not in mc_ad.var.columns:
        add_gene_id_to_mc_ad(mc_ad, c_ad)

    if orig_mc_ad:
        for ad_attr, columns_to_copy_names in ad_attr_and_columns_to_copy_names.items():
            old_column_names_already_in_new_mc_ad = set(columns_to_copy_names) & set(mc_ad.obs.columns if ad_attr == 'obs' else mc_ad.var.columns)
            assert not old_column_names_already_in_new_mc_ad, f'old_column_names_already_in_new_mc_ad: {old_column_names_already_in_new_mc_ad}'
            if ad_attr == 'obs':
                orig_obs = orig_mc_ad.obs[['metacell_i', *columns_to_copy_names]]
                mc_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(mc_ad.obs, orig_obs, how='left')
            else:
                assert ad_attr == 'var'
                orig_var = orig_mc_ad.var[['gnames', *columns_to_copy_names]]
                mc_ad.var = generic_utils.merge_preserving_df1_index_and_row_order(mc_ad.var, orig_var)

    mc.tl.find_metacells_marker_genes(mc_ad) 
    # mc.tl.find_metacells_significant_genes(mc_ad)
    # NOTE: IIUC, what was previously 'significant_gene', i.e., candidates for feature genes, is now 'selected_gene'. see /net/mraid14/export/tgdata/users/orenmil/anaconda3/envs/minimal/lib/python3.11/site-packages/metacells/pipeline/select.py in extract_selected_data. 
    
    # on 230702 it was:
    # if min_gene_top3 is not None:
    #     var_masks.append("&high_top3_gene")
    #     tl.find_high_topN_genes(adata, "downsampled", topN=3, min_gene_topN=min_gene_top3)
    # if min_gene_total is not None:
    #     var_masks.append("&high_total_gene")
    #     tl.find_high_total_genes(adata, "downsampled", min_gene_total=min_gene_total)
    # if min_gene_relative_variance is not None:
    #     var_masks.append("&high_relative_variance_gene")
    #     tl.find_high_relative_variance_genes(
    #         adata, "downsampled", min_gene_relative_variance=min_gene_relative_variance
    #     )
    # results = tl.filter_data(
    #     adata,
    #     name=name,
    #     top_level=top_level,
    #     track_var="full_gene_index",
    #     var_masks=var_masks + additional_gene_masks,
    #     mask_var="selected_gene",
    # )

    return mc_ad

def add_type_to_atlas(atlas_ad, atlas_cell_types_csv_file_path):
    atlas_cell_types_df = pd.read_csv(atlas_cell_types_csv_file_path, sep='\t')
    assert list(atlas_ad.obs_names) == [str(x) for x in list(atlas_cell_types_df['metacell'])]
    atlas_ad.obs['type'] = np.array(atlas_cell_types_df['cell_type'])

def get_atlas_ad(atlas_params, silently_overwrite_state_and_state_color_columns=True, skip_existing_state_col_handling=False, **add_state_kwargs):
    atlas_ad = ad.read_h5ad(atlas_params['mc_ad_file_path'])
    add_metacell_i_and_metacell_name(atlas_ad)
    if 'atlas_cell_types_csv_file_path' in atlas_params:
        add_type_to_atlas(atlas_ad, atlas_params['atlas_cell_types_csv_file_path'])
    else:
        state_column_name = atlas_params['state_column_name']
        if skip_existing_state_col_handling:
            orig_state_column_name = None
        else:
            assert (state_column_name in atlas_ad.obs.columns)
            orig_state_column_name = atlas_params.get('orig_state_column_name')
            if orig_state_column_name is not None:
                assert state_column_name != orig_state_column_name
                if orig_state_column_name not in atlas_ad.obs.columns:
                    atlas_ad.obs[orig_state_column_name] = atlas_ad.obs[state_column_name]
                print('orig_state.value_counts():')
                print(atlas_ad.obs[orig_state_column_name].value_counts())

    write_expr_and_expr_enrich(atlas_ad)
    if ('cell_state_and_info_list' in atlas_params) and (
        ('state' not in atlas_ad.obs.columns) or silently_overwrite_state_and_state_color_columns
    ):
        add_state(
            atlas_ad, 
            cell_state_and_info_list=atlas_params['cell_state_and_info_list'], 
            cell_type_colors_csv_file_path=atlas_params['atlas_cell_type_colors_csv_file_path'],
            default_state_column_name=orig_state_column_name,
            silently_overwrite_state_and_state_color_columns=silently_overwrite_state_and_state_color_columns,
            **add_state_kwargs,
        )
        
        ### atlas_ad.obs['type'] = atlas_ad.obs['state'] # removed on 231123. hopefully won't cause trouble.
        # print(atlas_ad.obs.columns)
        if orig_state_column_name is not None:
            same_state_mask = atlas_ad.obs[orig_state_column_name].astype(str) == atlas_ad.obs['state'].astype(str)
            print('replaced_states:')
            print(atlas_ad.obs[~same_state_mask][[orig_state_column_name, 'state']].astype(str).value_counts())

    add_selected_gene_best_equivalent_for_old_mc_ad(atlas_ad)
    return atlas_ad

def get_atlas_cells_ad(atlas_params):
    c_ad = ad.read_h5ad(atlas_params['c_ad_file_path'])
    if 'donor_id' not in list(c_ad.obs):
        # cells_with_donor_id_ad = ad.read_h5ad(atlas_params['atlas_cells_with_donor_id_ad_file_path']) # NOTE: 230703: i think we want atlas_cells_ad_file_path instead
        cells_with_donor_id_ad = ad.read_h5ad(atlas_params['c_ad_file_path'])
        c_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(
            c_ad.obs, cells_with_donor_id_ad.obs['indiv_id'].rename('donor_id'), left_index=True, right_index=True)
    return c_ad

def projection_pipeline_wrapper(
    # NOTE: CAUTION with project_correction! if the query is a single donor or only a few donors, you probably want project_correction=False. Why? especially if you look at aged and potentially sick people, they might have deletions/duplications etc in all of their sequenced cells. in such a case, the correction might correct the unexpected expression (e.g., lower (0.5x) due to a deletion) if it was in all query mcs.
        mc_ad, atlas_params, project_correction, 
        # epsilon_to_add_to_mc_gene_norm_expression=1e-5,
        epsilon_to_add_to_mc_gene_norm_expression=None,
        projection_weights_df_csv_file_path=None, out_fig_path=None, **kwargs,
):
    # TODO: what if we want to project mcs created from both ultima and illumina? (actually, how do we analyze this in the first place? need to correct as early as possible??)

    atlas_ad = get_atlas_ad(atlas_params)
    
    # mc.tl.find_metacells_significant_genes(atlas_ad)

    # if these exist and we use defaults argument for mc.pl.projection_pipeline(), it fails in renormalize_query_by_atlas().
    for x in ('expr', 'expr_enrich'):
        if x in mc_ad.layers:
            del mc_ad.layers[x]

    projection_weight_mat = mc.pl.projection_pipeline(
        adata=atlas_ad,
        qdata=mc_ad,
        # project_min_similar_essential_genes_fraction=None,
        project_corrections=project_correction,
        reproducible=True,
        atlas_type_property_name=atlas_params['state_column_name'],
        # ignore_atlas_lateral_genes=True, # already the default
        # ignore_query_lateral_genes=True, # already the default
        # consider_atlas_noisy_genes=False, # consider this?
        # consider_query_noisy_genes=False, # consider this?
        **kwargs,
    )

    mc_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(
        mc_ad.obs, pd.read_csv(atlas_params['atlas_cell_type_colors_csv_file_path']).rename(
            columns={'cell_type': 'projected_type', 'color': 'projected_type_color'}))
    mc_ad.obs['projected_type_color'] = mc_ad.obs['projected_type_color'].str.lower()

    for x in [x for x in list(mc_ad.var) if (
        ('misfit' in x)
        # or x.startswith('ignored_gene_') # seems to not be informative at all.
    )]:
        relevant_gene_names = list(mc_ad.var_names[mc_ad.var[x]])
        if relevant_gene_names:
            print(x)
            print(relevant_gene_names)

    if projection_weights_df_csv_file_path is not None:
        mc.pl.write_projection_weights(
            projection_weights_df_csv_file_path,
            adata=atlas_ad,
            qdata=mc_ad,
            weights=projection_weight_mat,
        )
    else:
        mc.ut.set_m_data(mc_ad, name='projection_mat', data=projection_weight_mat)

    # this is here because projection fails if it exists when we project (applies for both write_log_norm_expression_and_log_enrichment and write_gene_relative_variance).
    if epsilon_to_add_to_mc_gene_norm_expression is not None:
        write_log_norm_expression_and_log_enrichment_and_zm_score_and_norm_expression_std(mc_ad, epsilon_to_add_to_mc_gene_norm_expression)

    if out_fig_path is not None:
        pathlib.Path(out_fig_path).mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(figsize=(10,5))
        sb.histplot(mc_ad.obs['projected_correlation'], ax=ax)
        ax.set_xlim(0,1)
        fig.savefig(os.path.join(out_fig_path, f'obs_projected_correlation_hist.png'))

        fig, ax = plt.subplots(figsize=(10,5))
        plot_manifold_umap(
            atlas_ad,
            color_by=mc.ut.to_numpy_vector(projection_weight_mat.mean(axis=0)),
            # color_by=np.log2(mc.ut.to_numpy_vector(mc_ad.uns['projection_mat'].mean(axis=0)) + 1 / mc_ad.n_obs),
            color_by_name='mean_weight',
            legend='auto',
            ax=ax,
        )
        sb.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        fig.tight_layout()
        fig.savefig(os.path.join(out_fig_path, f'atlas_umap_colored_by_mean_projection_weight.png'))

    return mc_ad, atlas_ad


def plot_and_save_basic_metacell_model_plots(c_ad, mc_ad, metacell_model_figs_dir_path):
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    hist_plot_obj = sb.histplot(
        c_ad.obs['metacell'].value_counts().rename('num_of_cells_per_metacell'),
        bins=list(range(0, 200, 2)),
        ax=ax1,
    )
    hist_plot_obj.figure.savefig(os.path.join(metacell_model_figs_dir_path, 'num_of_cells_per_metacell_hist.png'))

    # fig1, ax1 = plt.subplots(figsize=(8, 5))
    # hist_plot_obj = sb.histplot(
    #     pd.Series(mc_ad.X.sum(axis=1), name='num_of_umis_per_metacell'),
    #     # bins=list(range(0, 200, 2)),
    #     ax=ax1,
    # )
    # hist_plot_obj.figure.savefig(os.path.join(metacell_model_figs_dir_path, 'num_of_umis_per_metacell_hist.png'))

# def copy_metacell_attributes_to_cell_attributes(mc_ad, c_ad, mc_ad_obs_column_names):

def write_projection_weight_df_csv(mc_ad, projection_weight_df_csv_file_path):
    flat_dicts = []
    for query_mc_i, atlas_mc_i in zip(*mc_ad.uns['projection_mat'].nonzero()):
        weight = mc_ad.uns['projection_mat'][query_mc_i, atlas_mc_i]
        flat_dicts.append({
            'query': query_mc_i,
            'atlas': atlas_mc_i,
            'weight': weight,
        })

    pd.DataFrame(flat_dicts).to_csv(projection_weight_df_csv_file_path)

def add_selected_gene_best_equivalent_for_old_mc_ad(mc_ad):
    if 'selected_gene' not in mc_ad.var.columns:
        mc_ad.var['selected_gene'] = (mc_ad.var['feature_gene'] > 0) | (mc_ad.var['marker_gene'] & (~mc_ad.var['lateral_gene']))
    else:
        print('WARNING: selected_gene already exists. make sure this makes sense.')


def get_feature_gene_mask(mc_ad):
    return mc_ad.var['selected_gene'].to_numpy()

def get_feature_gene_names(mc_ad):
    return mc_ad.var_names[get_feature_gene_mask(mc_ad)]

def get_mc_clusters(mc_ad, genes_to_cluster_by_mask=None, num_of_clusters=20, mcs_to_cluster_mask=None):
    if genes_to_cluster_by_mask is None:
        genes_to_cluster_by_mask = get_feature_gene_mask(mc_ad)

    if mcs_to_cluster_mask is None:
        mcs_to_cluster_mask = np.full(mc_ad.n_obs, True)

    kmeans = sklearn.cluster.KMeans(n_clusters=num_of_clusters, random_state=0, init='k-means++').fit(
        mc_ad.layers['expr'][np.ix_(mcs_to_cluster_mask, genes_to_cluster_by_mask)]
    )

    cluster_i_of_mcs = np.full(mc_ad.n_obs, -1)
    cluster_i_of_mcs[mcs_to_cluster_mask] = kmeans.labels_.astype(int)
    return cluster_i_of_mcs


def plot_cumulative_downsampled_umi_count_hist_per_exp_donor_heatmap(
        c_ad,
        gene_names,
        exp_names,
        cell_mask=None,
):
    attr_val_of_cells = c_ad.obs['exp_name'].astype(str) + ',' + c_ad.obs['donor_id'].astype(str)
    gene_indices = get_gene_indices(c_ad, gene_names)
    total_gene_downsampled_umis_of_cells = mc.ut.to_numpy_vector(c_ad.layers['downsampled'][:, gene_indices].sum(axis=1))
    exp_name_to_color = {
        k: v
        for k,v in zip(exp_names, generic_utils.get_n_colors(len(exp_names)))
    }
    if cell_mask is None:
        cell_mask = np.full(c_ad.n_obs, True)


    clustermap_obj = plot_heatmap_hist_per_cell_attr(
        attr_val_of_cells=attr_val_of_cells,
        val_of_cells=total_gene_downsampled_umis_of_cells,
        list_of_attr_val_to_color=[
            {x: exp_name_to_color[x.split(',')[0]] for x in attr_val_of_cells},
        ],
        ordered_attr_vals=sorted(
            set(attr_val_of_cells),
            key=lambda x: (
                exp_names.index(x.split(',')[0]),
                x.split(',')[1],
            ),
        ),
        normalize_bin_counts=True,
        # bins=list(range(60)),
        bins=list(range(int(total_gene_downsampled_umis_of_cells.max() + 1))),
        cell_mask=cell_mask,
        dilute_y_labels=False,
        left_y_ticks=False,
        cmap_name='gnuplot_r',
        cumulative=True,
        x_axis_in_log=True,
        add_num_of_cells_row_color=True,
        obj_desc='barcode',
        cbar_kwargs=dict(fraction=0.6, shrink=0.5),
        # log_x_axis_ticks=sorted(set(range(7)) | set(generic_utils.get_int_powers_of_2_in_range(6, 255))),
        # clustermap_kwargs=dict(dendrogram_ratio=(0.07, 0.15)),
    )

    clustermap_obj.fig.suptitle('+'.join(gene_names))
    return clustermap_obj

def get_z_score_df(
        mc_ad, metacell_mask, control_metacell_mask, cluster_i_of_mcs=None,
        only_feature_genes=False, only_feature_and_lateral_genes=False, min_log_norm_expression=-16, min_num_of_control_mcs=10,
):
    if cluster_i_of_mcs is None:
        cluster_i_of_mcs = np.full(mc_ad.n_obs, 0)

    if only_feature_genes:
        assert not only_feature_and_lateral_genes
        gene_mask = get_feature_gene_mask(mc_ad)
    elif only_feature_and_lateral_genes:
        gene_mask = get_feature_gene_mask(mc_ad) | mc_ad.var['forbidden_gene']
    else:
        gene_mask = np.full(mc_ad.n_vars, True)

    cluster_i_to_control_mean_and_sd_log_norm_expression = {}
    for cluster_i in sorted(set(cluster_i_of_mcs)):
        curr_control_mask = control_metacell_mask & (cluster_i_of_mcs == cluster_i)
        num_of_control_mcs = curr_control_mask.sum()
        # print(num_of_control_mcs)
        if num_of_control_mcs >= min_num_of_control_mcs:
            cluster_i_to_control_mean_and_sd_log_norm_expression[cluster_i] = (
                mc.ut.to_numpy_vector(mc_ad.layers['expr'][np.ix_(curr_control_mask, gene_mask)].mean(axis=0)),
                mc.ut.to_numpy_vector(mc_ad.layers['expr'][np.ix_(curr_control_mask, gene_mask)].std(axis=0)),
                num_of_control_mcs,
            )
    assert not (metacell_mask & control_metacell_mask).any()
    metacell_indices = np.where(metacell_mask)[0]
    dfs = []
    for metacell_i in metacell_indices:
        cluster_i = cluster_i_of_mcs[metacell_i]
        if cluster_i in cluster_i_to_control_mean_and_sd_log_norm_expression:
            dfs.append(pd.DataFrame({
                'gene': mc_ad.var_names[gene_mask],
                'metacell_i': metacell_i,
                'cluster_i': cluster_i,
                'expr': mc.ut.to_numpy_vector(mc_ad.layers['expr'][metacell_i, gene_mask]),
                'control_mean_log_norm_expression': cluster_i_to_control_mean_and_sd_log_norm_expression[cluster_i][0],
                'control_sd_log_norm_expression': cluster_i_to_control_mean_and_sd_log_norm_expression[cluster_i][1],
                'num_of_control_mcs': cluster_i_to_control_mean_and_sd_log_norm_expression[cluster_i][2],
            }))
    df = pd.concat(dfs, ignore_index=False)

    df['z_score'] = (df['expr'] - df['control_mean_log_norm_expression']) / df['control_sd_log_norm_expression']
    df['abs_z_score'] = df['z_score'].abs()
    
    # print((df['metacell_i'] == 3446).value_counts())
    # raise
    
    df = df[df['expr'] >= min_log_norm_expression]
    df.sort_values('abs_z_score', ascending=False, inplace=True)
    return df


def z_score_scatter(z_score_df, alpha_max_abs_z_score=2, plot_gene_name_min_abs_z_score=4):
    fig, ax = plt.subplots()
    x_column_name = 'control_mean_log_norm_expression'
    y_column_name = 'expr'
    x_err_column_name = 'control_sd_log_norm_expression'

    alpha_mask = z_score_df['abs_z_score'] <= alpha_max_abs_z_score

    sb.scatterplot(
        data=z_score_df,
        x=x_column_name,
        y=y_column_name,
        ax=ax,
        s=6,
    )

    ax.errorbar(
        x=z_score_df.loc[~alpha_mask, x_column_name],
        y=z_score_df.loc[~alpha_mask, y_column_name],
        xerr=z_score_df.loc[~alpha_mask, x_err_column_name],
        fmt='none',
        alpha=1,
    )
    ax.errorbar(
        x=z_score_df.loc[alpha_mask, x_column_name],
        y=z_score_df.loc[alpha_mask, y_column_name],
        xerr=z_score_df.loc[alpha_mask, x_err_column_name],
        fmt='none',
        alpha=0.2,
    )
    generic_utils.plot_y_equals_x_line_on_ax(ax=ax, linestyle='--', color='grey')

    generic_utils.plot_texts(
        ax, z_score_df, x_column_name, y_column_name, 'gene', plot_mask=(z_score_df['z_score'].abs() >= plot_gene_name_min_abs_z_score))

    return fig, ax

def plot_gene_gene_scatter_for_cells_of_metacell(c_ad, metacell_i, gene1_names, gene2_names):
    cell_mask = c_ad.obs['metacell_i'] == metacell_i
    num_of_cells = cell_mask.sum()
    print(f'num_of_cells: {num_of_cells}')

    gene1_mask = c_ad.var_names.isin(gene1_names)
    gene2_mask = c_ad.var_names.isin(gene2_names)

    jitter_sd = 0.07
    df_for_scatter = pd.DataFrame({
        'cell_name': c_ad.obs_names[cell_mask],
        'total_downsampled_genes1_umis_with_jitter': mc.ut.to_numpy_vector(
            c_ad.layers['downsampled'][np.ix_(cell_mask, gene1_mask)].sum(axis=1)) + np.random.normal(scale=jitter_sd, size=num_of_cells),
        'total_downsampled_genes2_umis_with_jitter': mc.ut.to_numpy_vector(
            c_ad.layers['downsampled'][np.ix_(cell_mask, gene2_mask)].sum(axis=1)) + np.random.normal(scale=jitter_sd, size=num_of_cells),
    })

    fig, ax = plt.subplots()
    sb.scatterplot(
        data=df_for_scatter,
        x='total_downsampled_genes1_umis_with_jitter',
        y='total_downsampled_genes2_umis_with_jitter',
        s=12,
        # alpha=0.3,
    )
    return fig, ax

def get_donor_id_or_doublet_or_unassigned_column(c_ad=None, c_ad_obs=None):
    if c_ad is not None:
        assert c_ad_obs is None
        # cells_ad_obs = c_ad.obs[['exp_name', 'donor_id', 'genotype_doublet']].copy()
        c_ad_obs = c_ad.obs[['donor_id', 'genotype_doublet']].copy()
    else:
        assert c_ad_obs is not None
        c_ad_obs = c_ad_obs.copy()

    assert set(c_ad_obs['genotype_doublet'].unique()) <= {True, False}
    assert (c_ad_obs.loc[c_ad_obs['genotype_doublet'], 'donor_id'].astype(str).isin(['nan', 'genotype_doublet'])).all()
    column = c_ad_obs['donor_id'].astype(str)
    column[c_ad_obs['genotype_doublet']] = 'doublet'
    column[column == 'nan'] = 'unassigned'
    return column

def get_assigned_or_doublet_or_unassigned_column(c_ad):
    column = get_donor_id_or_doublet_or_unassigned_column(c_ad)
    column[(column != 'doublet') & (column != 'unassigned')] = 'assigned'
    return column

def add_cell_type_color_column(
        cell_type_colors_csv_file_path, cell_state_column_name, data_ad=None, data_ad_obs=None, 
        cell_type_column_name_in_color_file='cell_type', color_column_name_in_color_file='color',
        add_outlier_black_color=False,
):
    if data_ad_obs is None:
        data_ad_obs = data_ad.obs
    else:
        assert data_ad is None
    
    color_column_name = f'{cell_state_column_name}_color'

    assert color_column_name not in data_ad_obs.columns
    
    color_df = pd.read_csv(cell_type_colors_csv_file_path)
    if add_outlier_black_color:
        color_df = pd.concat(
            [color_df, pd.DataFrame({cell_type_column_name_in_color_file: ['Outliers'], color_column_name_in_color_file: ['black']})], ignore_index=True)

    # print(data_ad_obs.head(2), print(pd.read_csv(cell_type_colors_csv_file_path).head(2)))
    data_ad_obs = generic_utils.merge_preserving_df1_index_and_row_order(
        data_ad_obs, color_df.rename(
            columns={cell_type_column_name_in_color_file: cell_state_column_name, color_column_name_in_color_file: color_column_name}),
        verbose=False,
        # print_rows_that_were_lost_in_the_merge=True,
    )
    data_ad_obs[color_column_name] = data_ad_obs[color_column_name].str.lower()
    
    if data_ad is not None:
        data_ad.obs = data_ad_obs

    return data_ad_obs

def add_state(
        mc_ad, cell_state_and_info_list, cell_type_colors_csv_file_path, 
        default_state_column_name=None, 
        state_column_name='state',
        assign_dissimilar_according_to_similar_field=True,
        verbose=True, silently_overwrite_state_and_state_color_columns=False,
        dist_with_threshs_out_dir_path=None,
):
    state_color_column_name = f'{state_column_name}_color'
    if state_column_name in list(mc_ad.obs):
        if silently_overwrite_state_and_state_color_columns:
            mc_ad.obs.drop(columns=state_column_name, inplace=True)
        else:
            raise RuntimeError('state column already exists in mc_ad.obs')
    if state_color_column_name in list(mc_ad.obs):
        if not silently_overwrite_state_and_state_color_columns:
            raise RuntimeError(f'{state_color_column_name} (state_color_column_name) column already exists in mc_ad.obs')
    mc_ad.obs[state_color_column_name] = 'black' # to avoid the annoying case of an error later, then trying to plot, with a state column but no state_color column.

    assign_cell_state(
        mc_ad, cell_state_and_info_list, 
        assign_dissimilar_according_to_similar_field=assign_dissimilar_according_to_similar_field, verbose=verbose,
        state_column_name=state_column_name,
        dist_with_threshs_out_dir_path=dist_with_threshs_out_dir_path,
    )

    if default_state_column_name:
        unassigned_mask = mc_ad.obs[state_column_name] == 'state_unassigned'
        mc_ad.obs.loc[unassigned_mask, state_column_name] = mc_ad.obs.loc[unassigned_mask, default_state_column_name].astype(str).copy()
    assert not mc_ad.obs[state_column_name].isna().any()

    if verbose:
        print(f'state_color_column_name: {state_color_column_name}')
    mc_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(
        mc_ad.obs.drop(columns=state_color_column_name), pd.read_csv(cell_type_colors_csv_file_path).rename(
            columns={'cell_type': state_column_name, 'color': state_color_column_name}))
    mc_ad.obs[state_color_column_name] = mc_ad.obs[state_color_column_name].str.lower()
    assert not mc_ad.obs[state_color_column_name].isna().any()
    state_unassigned_mask = mc_ad.obs[state_column_name] == 'state_unassigned'
    if verbose:
        print(f'state_unassigned_mask.sum(): {state_unassigned_mask.sum()}')
    
    

def get_presumably_healthy_metacell_mask(mc_ad, max_max_donor_id_fraction=0.3):
    disease_disease_cell_types = [
        'N205-MDS', 
        'N186-MDS', 
        'N191-AML', 
        # 'N224-high-ribo',
        
        'Doublet', 
        'Mixture', 
    ]
    return np.array(
        np.array(mc_ad.obs['max_donor_id_fraction'] <= max_max_donor_id_fraction)
        # & (mc_ad.obs['projected_type'] != 'Dissimilar')
        & (~mc_ad.obs[state_column_name].isin(disease_disease_cell_types))
        # & (mc_ad.obs['projected_type'].astype(str) == mc_ad.obs[state_column_name].astype(str))
    )

def get_obs_names_of_cells_of_mc_ad_obs_excluding_cells_of_most_common_donor_in_each_metacell(c_ad, min_mc_size_after_filtering=None):
    cells_ad_obs = c_ad.obs.loc[c_ad.obs['metacell_name'] != 'Outliers', ['metacell_i', 'donor_id']]
    cells_ad_obs['metacell_i'] = cells_ad_obs['metacell_i'].astype(int)
    cells_ad_obs['donor_id'] = cells_ad_obs['donor_id'].astype(str)
    metacell_most_common_donor_df = cells_ad_obs[['metacell_i', 'donor_id']].value_counts().reset_index().drop_duplicates(
        subset='metacell_i', keep='first').drop('count', axis=1)

    extended_cells_ad_obs = generic_utils.merge_preserving_df1_index_and_row_order(
        cells_ad_obs, metacell_most_common_donor_df, indicator=True, how='left', on=['metacell_i', 'donor_id'])
    extended_cells_ad_obs = extended_cells_ad_obs[extended_cells_ad_obs['_merge'] == 'left_only']
    if min_mc_size_after_filtering:
        # print(extended_cells_ad_obs['metacell_i'].value_counts())
        # print(len(extended_cells_ad_obs))
        big_metacell_i_values = (extended_cells_ad_obs['metacell_i'].value_counts() >= min_mc_size_after_filtering).loc[lambda x: x].index
        extended_cells_ad_obs = extended_cells_ad_obs[extended_cells_ad_obs['metacell_i'].isin(big_metacell_i_values)]
        # print(extended_cells_ad_obs)
    return extended_cells_ad_obs.index

def plot_exp_fixed_fractions_of_metacell_cells_excluding_exp_max_donor_cells(
        c_ad,
        mc_ad,
        ordered_exp_names,
        min_metacell_size_after_exclusion=10,
        min_max_exp_fixed_fraction_of_metacell_cells_after_exclusion_to_plot_umap=0.3,
):
    cells_ad_obs = c_ad.obs.loc[
        get_obs_names_of_cells_of_mc_ad_obs_excluding_cells_of_most_common_donor_in_each_metacell(c_ad), ['metacell_i', 'exp_name']].copy()

    cells_ad_obs['metacell_i'] = cells_ad_obs['metacell_i'].astype(int)
    cells_ad_obs['exp_name'] = cells_ad_obs['exp_name'].astype(str)

    metacell_size_df = cells_ad_obs['metacell_i'].value_counts().reset_index(name='metacell_size_after_exclusion').rename(columns={'index': 'metacell_i'})

    metacell_exp_df = cells_ad_obs.value_counts().reset_index(name='exp_num_of_metacell_cells_after_exclusion')
    metacell_exp_df = generic_utils.merge_preserving_df1_index_and_row_order(metacell_exp_df, metacell_size_df)
    metacell_exp_df['exp_fraction_of_metacell_cells_after_exclusion'] = (
        metacell_exp_df['exp_num_of_metacell_cells_after_exclusion'] / metacell_exp_df['metacell_size_after_exclusion'])

    metacell_exp_df['exp_fixed_fraction_of_metacell_cells_after_exclusion'] = metacell_exp_df['exp_fraction_of_metacell_cells_after_exclusion']
    metacell_exp_df.loc[metacell_exp_df['metacell_size_after_exclusion'] < min_metacell_size_after_exclusion, 'exp_fixed_fraction_of_metacell_cells_after_exclusion'] = -0.1
    relevant_exp_names = [
        x for x in ordered_exp_names
        if x in set((
            metacell_exp_df.groupby('exp_name')['exp_fixed_fraction_of_metacell_cells_after_exclusion'].max() >= 
            min_max_exp_fixed_fraction_of_metacell_cells_after_exclusion_to_plot_umap
        ).loc[lambda s: s].index)
    ]

    required_num_of_axes = len(relevant_exp_names)
    nrows = 4
    ncols = 3
    assert nrows * ncols >= required_num_of_axes, f'required_num_of_axes: {required_num_of_axes}'
    plt.close('all')
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, 
        figsize=(30,35),
    )
    for ax, exp_name in zip(axes.flatten(), relevant_exp_names):
        curr_df = generic_utils.merge_preserving_df1_index_and_row_order(
            mc_ad.obs[['metacell_i']], 
            metacell_exp_df.loc[metacell_exp_df['exp_name'] == exp_name, ['metacell_i', 'exp_fixed_fraction_of_metacell_cells_after_exclusion']],
            how='left',
        )
        curr_df['exp_fixed_fraction_of_metacell_cells_after_exclusion'].fillna(-0.1, inplace=True)
        
        plot_manifold_umap(
            mc_ad, 
            color_by=curr_df['exp_fixed_fraction_of_metacell_cells_after_exclusion'].rename('fixed_fraction'),
            hue_norm=(0,1),
            legend='auto',
            # alpha=0.7,
            palette='rocket_r',
            ax=ax,
        )
        ax.set_title(exp_name)

    for ax in axes.flatten()[required_num_of_axes:]:
        generic_utils.make_all_spines_and_x_and_y_axes_invisible(ax)
    
    return fig, axes

def plot_sex_fixed_fractions_of_metacell_cells_excluding_sex_max_donor_cells(
        c_ad,
        mc_ad,
        min_metacell_size_after_exclusion=10,
):
    cells_ad_obs = c_ad.obs.loc[
        get_obs_names_of_cells_of_mc_ad_obs_excluding_cells_of_most_common_donor_in_each_metacell(c_ad), ['metacell_i', 'donor_sex']].copy()

    cells_ad_obs['metacell_i'] = cells_ad_obs['metacell_i'].astype(int)
    cells_ad_obs['donor_sex'] = cells_ad_obs['donor_sex'].astype(str)

    metacell_size_df = cells_ad_obs['metacell_i'].value_counts().reset_index(name='metacell_size_after_exclusion').rename(columns={'index': 'metacell_i'})

    metacell_sex_df = cells_ad_obs.value_counts().reset_index(name='sex_num_of_metacell_cells_after_exclusion')
    metacell_sex_df = generic_utils.merge_preserving_df1_index_and_row_order(metacell_sex_df, metacell_size_df)
    metacell_sex_df['sex_fraction_of_metacell_cells_after_exclusion'] = (
        metacell_sex_df['sex_num_of_metacell_cells_after_exclusion'] / metacell_sex_df['metacell_size_after_exclusion'])

    metacell_sex_df['sex_fixed_fraction_of_metacell_cells_after_exclusion'] = metacell_sex_df['sex_fraction_of_metacell_cells_after_exclusion']
    metacell_sex_df.loc[metacell_sex_df['metacell_size_after_exclusion'] < min_metacell_size_after_exclusion, 'sex_fixed_fraction_of_metacell_cells_after_exclusion'] = -0.1
    plt.close('all')
    fig, axes = plt.subplots(
        nrows=2, ncols=1, 
        figsize=(6,10),
    )
    for ax, sex in zip(axes.flatten(), ('female', 'male')):
        curr_df = generic_utils.merge_preserving_df1_index_and_row_order(
            mc_ad.obs[['metacell_i']], 
            metacell_sex_df.loc[metacell_sex_df['donor_sex'] == sex, ['metacell_i', 'sex_fixed_fraction_of_metacell_cells_after_exclusion']],
            how='left',
        )
        curr_df['sex_fixed_fraction_of_metacell_cells_after_exclusion'].fillna(-0.1, inplace=True)
        
        plot_manifold_umap(
            mc_ad, 
            color_by=curr_df['sex_fixed_fraction_of_metacell_cells_after_exclusion'].rename('fixed_fraction'),
            hue_norm=(0,1),
            legend='auto',
            # alpha=0.7,
            palette='rocket_r',
            ax=ax,
        )
        ax.set_title(sex)

    return fig, axes

def plot_donor_fractions_of_metacell_cells(
        c_ad,
        mc_ad,
):
    relevant_donor_ids = [
        x for x in sorted(c_ad.obs['donor_id'].unique())
        if (f'donor_id_{x}_fraction' in list(mc_ad.obs)) and (mc_ad.obs[f'donor_id_{x}_fraction'].max() > 0.5).any()
    ]

    required_num_of_axes = len(relevant_donor_ids)
    # raise RuntimeError(required_num_of_axes)
    nrows = 5
    ncols = 3
    assert nrows * ncols >= required_num_of_axes
    plt.close('all')
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, 
        figsize=(24,25),
    )
    for ax, donor_id in zip(axes.flatten(), relevant_donor_ids):
        plot_manifold_umap(
            mc_ad, 
            color_by=mc_ad.obs[f'donor_id_{donor_id}_fraction'],
            hue_norm=(0,1),
            legend='auto',
            # alpha=0.7,
            palette='rocket_r',
            ax=ax,
        )
        ax.set_title(donor_id)

    for ax in axes.flatten()[required_num_of_axes:]:
        generic_utils.make_all_spines_and_x_and_y_axes_invisible(ax)

    return fig, axes

def plot_exp_fractions_of_metacell_cells_after_excluding_specific_donors(
        c_ad,
        mc_ad,
        fig_out_dir_path,
        list_of_exp_name_and_excluded_donor_ids,
        min_num_of_metacell_cells_after_exclusion=5,
):
    cells_ad_obs = c_ad.obs.loc[c_ad.obs['metacell_name'] != 'Outliers', ['metacell_i', 'exp_name', 'donor_id']]
    cells_ad_obs['metacell_i'] = cells_ad_obs['metacell_i'].astype(int)
    cells_ad_obs[['exp_name', 'donor_id']] = cells_ad_obs[['exp_name', 'donor_id']].astype(str)

    for exp_name, excluded_donor_ids in list_of_exp_name_and_excluded_donor_ids:
        fig, ax = plt.subplots()

        curr_c_ad_obs = cells_ad_obs.copy()
        curr_c_ad_obs = curr_c_ad_obs[~curr_c_ad_obs['donor_id'].isin(excluded_donor_ids)]
        
        curr_metacell_df = generic_utils.merge_preserving_df1_index_and_row_order(
            mc_ad.obs[['metacell_i']], 
            curr_c_ad_obs['metacell_i'].value_counts().reset_index(name='num_of_metacell_cells_after_exclusion').rename(columns={'index': 'metacell_i'}),
            how='left',
        )
        curr_metacell_df['num_of_metacell_cells_after_exclusion'].fillna(0, inplace=True)

        curr_metacell_df = generic_utils.merge_preserving_df1_index_and_row_order(
            curr_metacell_df, 
            curr_c_ad_obs.loc[curr_c_ad_obs['exp_name'] == exp_name, 'metacell_i'].value_counts().reset_index(
                name='exp_num_of_metacell_cells_after_exclusion').rename(columns={'index': 'metacell_i'}),
            how='left',
        )
        curr_metacell_df['exp_num_of_metacell_cells_after_exclusion'].fillna(0, inplace=True)
        
        div_by_zero_mask = curr_metacell_df['num_of_metacell_cells_after_exclusion'] == 0
        curr_metacell_df.loc[~div_by_zero_mask, 'exp_fraction_of_metacell_cells_after_exclusion'] = (
            curr_metacell_df.loc[~div_by_zero_mask, 'exp_num_of_metacell_cells_after_exclusion'] / 
            curr_metacell_df.loc[~div_by_zero_mask, 'num_of_metacell_cells_after_exclusion']
        )
        curr_metacell_df.loc[div_by_zero_mask, 'exp_fraction_of_metacell_cells_after_exclusion'] = 0
        assert not curr_metacell_df['exp_fraction_of_metacell_cells_after_exclusion'].isna().any()
        assert (curr_metacell_df['exp_fraction_of_metacell_cells_after_exclusion'] >= 0).all()
        assert (curr_metacell_df['exp_fraction_of_metacell_cells_after_exclusion'] <= 1).all()

        curr_metacell_df['exp_fixed_fraction_of_metacell_cells_after_exclusion'] = curr_metacell_df['exp_fraction_of_metacell_cells_after_exclusion']
        curr_metacell_df.loc[
            curr_metacell_df['num_of_metacell_cells_after_exclusion'] < 
            min_num_of_metacell_cells_after_exclusion, 'exp_fixed_fraction_of_metacell_cells_after_exclusion'] = -0.1

        plot_manifold_umap(
            mc_ad, 
            color_by=curr_metacell_df['exp_fixed_fraction_of_metacell_cells_after_exclusion'].rename('fixed fraction'),
            hue_norm=(0,1),
            legend='auto',
            # alpha=0.7,
            palette='rocket_r',
            ax=ax,
        )
        excluded_donor_ids_repr = '_'.join(excluded_donor_ids)
        ax_title = f'{exp_name}_fraction_excluding_{excluded_donor_ids_repr}'
        ax.set_title(ax_title)

        fig.savefig(os.path.join(fig_out_dir_path, f'{ax_title}_on_umap.png'))

def get_gene_score_per_attrs(
        c_ad, gene_names, attr_names, cell_mask, min_num_of_cells, epsilon_for_log_norm_expression=1e-5, 
        geomean_across_pooled_cells=True,
):
    
    if isinstance(gene_names, str):
        gene_names = [gene_names]

    gene_indices = get_gene_indices(c_ad, gene_names)
    gene_mask = c_ad.var_names.isin(gene_names)

    num_of_non_excluded_umis_of_cells = mc.ut.to_numpy_vector(c_ad.X.sum(axis=1))
    if 'num_of_non_excluded_umis' in list(c_ad.obs):
        # return num_of_non_excluded_umis_of_cells, c_ad.obs['num_of_non_excluded_umis']
        assert np.isclose(num_of_non_excluded_umis_of_cells, c_ad.obs['num_of_non_excluded_umis']).all()

    gene_norm_expression = mc.ut.to_numpy_matrix(
        c_ad.X[:, gene_mask]) / num_of_non_excluded_umis_of_cells[:, np.newaxis]

    if isinstance(attr_names, str):
        attr_names = [attr_names]

    flat_dicts = []
    for _, attr_val_row in c_ad.obs.loc[cell_mask, attr_names].drop_duplicates().iterrows():
        print(attr_val_row)
        # print(attr_name, attr_val)
        # print(attr_val_row.to_frame().T.reset_index(drop=True))
        # print()
        # print(c_ad.obs[attr_names].head())
        # curr_cell_mask = cell_mask & (c_ad.obs[attr_names] == attr_val_row.to_frame().T.reset_index(drop=True)).all(axis=1)
        curr_cell_mask = cell_mask & (c_ad.obs[attr_names] == attr_val_row).all(axis=1)
        num_of_cells_used = curr_cell_mask.sum()
        # print(f'num_of_cells_used: {num_of_cells_used}')
        # print(f'attr_val_row: {attr_val_row}')
        # raise
        if num_of_cells_used >= min_num_of_cells:
            if geomean_across_pooled_cells:
                expr = np.log2(gene_norm_expression[curr_cell_mask, :] + epsilon_for_log_norm_expression).mean(axis=0)
            else:
                expr = np.log2(
                    mc.ut.to_numpy_vector(c_ad.X[np.ix_(curr_cell_mask, gene_mask)].sum(axis=0)) / 
                    num_of_non_excluded_umis_of_cells[curr_cell_mask].sum()
                    + epsilon_for_log_norm_expression
                )
            assert len(expr) == len(gene_indices)
            log_geomean_norm_expression = expr.mean()
        else:
            log_geomean_norm_expression = np.nan
        flat_dicts.append({
            # **{k: v for k,v in zip(attr_names: attr_val},
            **attr_val_row.to_dict(),
            'sig_score': log_geomean_norm_expression,
            'num_of_cells_used': num_of_cells_used,
        })
    score_df = pd.DataFrame(flat_dicts)
    # score_df['log2_num_of_cells_used'] = np.log2(score_df['num_of_cells_used'])
    score_df.sort_values(attr_names, inplace=True)
    return score_df

def add_bin_assignment_to_mc_ad(
        mc_ad, bin_column_name, num_of_bins, gene_name, mcs_to_assign_bins_mask=None, ascending=True,
):
    if mcs_to_assign_bins_mask is None:
        mcs_to_assign_bins_mask = np.full(mc_ad.n_obs, True)

    sort_key_column_name = f'log_norm_{gene_name}'
    metacell_size_and_bin_df = mc_ad.obs.loc[mcs_to_assign_bins_mask, 'grouped'].to_frame()
    metacell_size_and_bin_df[sort_key_column_name] = mc_ad.layers['expr'][
        mcs_to_assign_bins_mask, get_gene_indices(mc_ad, gene_name)]
    metacell_size_and_bin_df.sort_values(sort_key_column_name, ascending=ascending, inplace=True)
    
    metacell_size_and_bin_df[bin_column_name] = generic_utils.split_sorted_counts_to_roughly_same_size_bins(
        metacell_size_and_bin_df['grouped'], num_of_bins=num_of_bins)

    min_column_name = f'min_{sort_key_column_name}'
    max_column_name = f'max_{sort_key_column_name}'
    metacell_size_and_bin_df = generic_utils.merge_preserving_df1_index_and_row_order(
        metacell_size_and_bin_df, metacell_size_and_bin_df.groupby(
        bin_column_name)[sort_key_column_name].min().reset_index(name=min_column_name))
    metacell_size_and_bin_df = generic_utils.merge_preserving_df1_index_and_row_order(
        metacell_size_and_bin_df, metacell_size_and_bin_df.groupby(
        bin_column_name)[sort_key_column_name].max().reset_index(name=max_column_name))
    bin_df = metacell_size_and_bin_df[[bin_column_name, min_column_name, max_column_name]].drop_duplicates().reset_index(drop=True)

    print(f'bin_df:\n{bin_df}')

    if bin_column_name in list(mc_ad.obs):
        print(f'{bin_column_name} is already a column in mc_ad.obs. dropping it for the upcoming merge.')
        mc_ad.obs.drop(bin_column_name, axis=1, inplace=True)
    mc_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(
        mc_ad.obs, metacell_size_and_bin_df[bin_column_name], left_index=True, right_index=True, how='left')

    plot_manifold_umap(mc_ad)
    plot_manifold_umap(mc_ad, color_by=bin_column_name)

    return bin_df

def add_bin_assignment_to_cells_ad_according_to_mc_expr(
        mc_ad, c_ad, bin_column_name, num_of_bins, gene_names, mcs_to_assign_bins_mask=None, ascending=True,
):
    if mcs_to_assign_bins_mask is None:
        mcs_to_assign_bins_mask = np.full(mc_ad.n_obs, True)

    if isinstance(gene_names, str):
        gene_names = [gene_names]
    
    gene_names_repr = '_'.join(gene_names)
    sort_key_column_name = f'norm_{gene_names_repr}'
    mc_obs = mc_ad.obs.loc[mcs_to_assign_bins_mask, ['metacell_i']].copy()
    mc_obs[sort_key_column_name] = (
        mc_ad.X[np.ix_(mcs_to_assign_bins_mask, mc_ad.var_names.isin(gene_names))].sum(axis=1) / 
        mc_ad.X[mcs_to_assign_bins_mask, :].sum(axis=1)
    )
    if sort_key_column_name in list(c_ad.obs):
        print(f'{sort_key_column_name} is already a column in c_ad.obs. dropping it for the upcoming merge.')
        c_ad.obs.drop(sort_key_column_name, axis=1, inplace=True)
    c_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(c_ad.obs, mc_obs, how='left')

    cells_mask = ~c_ad.obs[sort_key_column_name].isna()
    c_ad.obs.loc[cells_mask, bin_column_name] = generic_utils.equal_size_bin_digitize_by_single_nums_vec(
        c_ad.obs.loc[cells_mask, sort_key_column_name], num_of_bins)

    if not ascending:
        c_ad.obs.loc[cells_mask, bin_column_name] = num_of_bins - 1 - c_ad.obs.loc[cells_mask, bin_column_name]

    # TODO: add mode bin to each metacell, and plot this in the umap
    add_mc_metadata(
        mc_ad,
        c_ad,
        statistic_type_to_column_names_in_cells_ad_obj={'nan_mean_and_median_and_std': [bin_column_name]},
        add_num_of_cells_in_metacell=False,
    )

    plot_manifold_umap(mc_ad)
    plot_manifold_umap(mc_ad, color_by=f'{bin_column_name}_median')


def get_cna_gene_mask(mc_ad, var_df, minimal_filtered_var_df_sorted_by_unified_end_pos, cna_info, num_of_genes_per_margin=0):
    cna_gene_mask = var_df['chrom_name'] == cna_info['chrom']
    # sorted_genes_with_unified_end_pos = var_df.loc[cna_gene_mask, 'unified_end_pos'].sort_values()
    if 'first_gene' in cna_info:
        first_gene_i = list(minimal_filtered_var_df_sorted_by_unified_end_pos.index).index(cna_info['first_gene'])
        cna_gene_mask &= mc_ad.var_names.isin(minimal_filtered_var_df_sorted_by_unified_end_pos.iloc[
            max(0, first_gene_i - num_of_genes_per_margin):].index)
    if 'last_gene' in cna_info:
        last_gene_i = list(minimal_filtered_var_df_sorted_by_unified_end_pos.index).index(cna_info['last_gene'])
        cna_gene_mask &= mc_ad.var_names.isin(minimal_filtered_var_df_sorted_by_unified_end_pos.iloc[:(last_gene_i + num_of_genes_per_margin + 1)].index)
    
    if ('first_gene' in cna_info) and ('last_gene' in cna_info):
        assert first_gene_i < last_gene_i
    elif ('first_gene' not in cna_info) and ('last_gene' not in cna_info):
        cna_gene_mask &= mc_ad.var_names.isin(minimal_filtered_var_df_sorted_by_unified_end_pos.index)
    else:
        print('either first_gene or last_gene is missing, thus assuming the CNA continues to the end of the chromosome.')

    # print(first_gene_i, last_gene_i)
    if 'presumably_unaffected_genes' in cna_info:
        assert np.isin(cna_info['presumably_unaffected_genes'], mc_ad.var_names[cna_gene_mask]).all()
        unaffected_gene_mask = mc_ad.var_names.isin(cna_info['presumably_unaffected_genes'])
        cna_gene_mask &= ~unaffected_gene_mask
    

    return cna_gene_mask


def add_proj_fold_layer(mc_ad, epsilon_for_expr=1e-5, ignore_correction=False):
    if 'proj_fold' in mc_ad.layers:
        print('ad.layers already has proj_fold. doing nothing')
        return
    
    if ignore_correction:
        expr_fracs = mc_ad.X
    else:
        expr_fracs = mc_ad.layers['corrected_fraction']

    assert mc_ad.n_obs * mc_ad.n_vars < 1e9, 'using mc.ut.to_numpy_matrix on a sparse matrix with more than 1e9 elements seems like a bad idea'
    expr = np.log2(mc.ut.to_numpy_matrix(expr_fracs) + epsilon_for_expr)
    proj_expr = np.log2(mc.ut.to_numpy_matrix(mc_ad.layers['projected_fraction']) + epsilon_for_expr)

    mc.ut.set_vo_data(mc_ad, name='proj_fold', data=(expr - proj_expr))

def get_correlations_with_mc_obs_columns():
    raise NotImplementedError('implement this. maybe get_test_vs_control_mc_obs_columns_df could be useful?')
    if column_names is None:
        # column_names = [x for x in list(mc_ad.obs) if x.endswith('_fraction')]
        column_names = [x for x in list(mc_ad.obs) if (mc_ad.obs[x].dtype.kind in 'biuf')]


def get_test_vs_control_categorical_str_c_obs_columns_df(
        c_ad, test_c_mask=None, control_c_mask=None, col_names=['donor_id', 'exp_name'],
        test_mc_mask=None, control_mc_mask=None, mc_ad=None, skip_numerical_cols=True, skip_barcode_col=True,
):
    if test_c_mask is None:
        assert control_c_mask is None
        test_c_mask = mc_mask_to_c_mask(mc_ad, c_ad, test_mc_mask)
        if control_mc_mask is None:
            control_c_mask = ~test_c_mask
        else:
            control_c_mask = mc_mask_to_c_mask(mc_ad, c_ad, control_mc_mask)
    else:
        assert test_mc_mask is None
        assert control_mc_mask is None
        if control_c_mask is None:
            control_c_mask = ~test_c_mask
    
    assert test_c_mask.any()
    assert control_c_mask.any()
    
    num_of_test_cells = test_c_mask.sum()
    num_of_control_cells = control_c_mask.sum()
    print(f'num_of_test_cells: {num_of_test_cells}')
    print(f'num_of_control_cells: {num_of_control_cells}')
    
    # all_log_ratio = np.log2(num_of_test_cells / num_of_control_cells)

    flat_dicts = []
    for col_name in col_names:
        if skip_barcode_col and (col_name == 'barcode'):
            print(f'skipping barcode column')
            continue
        if skip_numerical_cols and (c_ad.obs[col_name].dtype.kind in 'iuf'):
            print(f'skipping {col_name} because it is numerical')
            continue

        print(col_name)
        for val in sorted(c_ad.obs[col_name].astype(str).unique()):
            val_mask = c_ad.obs[col_name] == val
            num_of_val_test_cells = (test_c_mask & val_mask).sum()
            num_of_non_val_test_cells = num_of_test_cells - num_of_val_test_cells
            num_of_val_control_cells = (control_c_mask & val_mask).sum()
            num_of_non_val_control_cells = num_of_control_cells - num_of_val_control_cells
            if (num_of_val_test_cells == 0) & (num_of_val_control_cells == 0):
                print(f'{val}: (num_of_val_test_cells == 0) & (num_of_val_control_cells == 0)')
                continue
            val_log_ratio = (
                np.log2((num_of_val_test_cells + 1) / (num_of_val_control_cells + 1)) - 
                np.log2((num_of_non_val_test_cells + 1) / (num_of_non_val_control_cells + 1))
            )
            if (val_log_ratio > 0) and (num_of_val_test_cells == 0):
                print(f'{val}: (val_log_ratio > 0) and (num_of_val_test_cells == 0)')
                continue
            flat_dicts.append({
                'column_name': col_name,
                'column_val': val,
                'num_of_val_test_cells': num_of_val_test_cells,
                'num_of_val_control_cells': num_of_val_control_cells,
                'log_ratio': val_log_ratio,
            })

    df = pd.DataFrame(flat_dicts)
    df.sort_values('log_ratio', ascending=False, inplace=True)            
    return df

def mc_mask_to_c_mask(mc_ad, c_ad, mc_mask):
    # NOTE: kind of the opposite is implemented by add_c_mask_fraction_of_mc_column
    return c_ad.obs['metacell_name'].isin(mc_ad.obs.loc[mc_mask, 'metacell_name']).to_numpy()

def add_c_mask_fraction_of_mc_column(mc_ad, c_ad, c_mask, column_name, background_c_mask=None):
    # stuff for searching: c_mask_to_mc_frac, c_mask_to_mc
    mc_ad.obs[column_name] = get_c_mask_fraction_of_mcs(mc_ad, c_ad, c_mask=c_mask, background_c_mask=background_c_mask)

def get_test_vs_control_state_enrichment_df(
        c_or_mc_ad, test_mask, control_mask, state_column_name='state',
):
    assert test_mask.any()
    assert control_mask.any()
    
    num_of_test_cs_or_mcs = test_mask.sum()
    num_of_control_cs_or_mcs = control_mask.sum()
    print(f'num_of_test_cs_or_mcs: {num_of_test_cs_or_mcs}')
    print(f'num_of_control_cs_or_mcs: {num_of_control_cs_or_mcs}')
    
    states = sorted(c_or_mc_ad.obs.loc[test_mask | control_mask, state_column_name].astype(str).unique())
    # all_log_ratio = np.log2(num_of_test_cs_or_mcs / num_of_control_cs_or_mcs)

    flat_dicts = []
    for state in states:
        state_mask = c_or_mc_ad.obs[state_column_name] == state
        num_of_state_test_cells = (test_mask & state_mask).sum()
        num_of_other_state_test_cells = num_of_test_cs_or_mcs - num_of_state_test_cells
        num_of_state_control_cells = (control_mask & state_mask).sum()
        num_of_other_state_control_cells = num_of_control_cs_or_mcs - num_of_state_control_cells
        if (num_of_state_test_cells == 0) & (num_of_state_control_cells == 0):
            continue
        state_log_ratio = (
            np.log2((num_of_state_test_cells + 1) / (num_of_state_control_cells + 1)) - 
            np.log2((num_of_other_state_test_cells + 1) / (num_of_other_state_control_cells + 1))
        )
        if (state_log_ratio > 0) and (num_of_state_test_cells == 0):
            continue
        flat_dicts.append({
            'state': state,
            'num_of_state_test_cells': num_of_state_test_cells,
            'num_of_other_state_test_cells': num_of_other_state_test_cells,
            'num_of_state_control_cells': num_of_state_control_cells,
            'num_of_other_state_control_cells': num_of_other_state_control_cells,
            'log_ratio': state_log_ratio,
        })

    df = pd.DataFrame(flat_dicts)
    df.sort_values('log_ratio', ascending=False, inplace=True)            
    return df
    
def get_test_vs_control_mc_obs_columns_df(
        mc_ad, test_mask, control_mask, column_names=None, use_welch_t_test=True, regularization_epsilon=1e-4,
):
    if column_names is None:
        # column_names = [x for x in list(mc_ad.obs) if x.endswith('_fraction')]
        column_names = [x for x in list(mc_ad.obs) if (mc_ad.obs[x].dtype.kind in 'biuf')]
        # column_names = list(mc_ad.obs)
    
    # column_names = np.array(column_names)

    # for x in column_names:
    #     print(x)
    #     np.log2(
    #         (mc_ad.obs.loc[test_mask, x].mean() + regularization_epsilon) / 
    #         (mc_ad.obs.loc[control_mask, x].mean() + regularization_epsilon)
    #     )

    if use_welch_t_test:
        flat_dicts = []
        for column_name in column_names:
            column = mc_ad.obs[column_name]
            if column.dtype.kind == 'b':
                column = column.astype(int)
            test_vals = column[test_mask]
            control_vals = column[control_mask]
            nunique_test_vals = test_vals.nunique()
            nunique_control_vals = control_vals.nunique()
            if (nunique_test_vals == 1) or (nunique_control_vals == 1):
                pval_column_name = 'ztest_pval'
                if (nunique_test_vals == 1) and (nunique_control_vals == 1):
                    if ~np.isclose(test_vals.iloc[0], control_vals.iloc[0]):
                        ztest_pval = 0
                        z_score = np.inf * np.sign(test_vals.iloc[0] - control_vals.iloc[0])
                    else:
                        ztest_pval = 1
                        z_score = 0
                elif nunique_test_vals == 1:
                    z_score, ztest_pval = generic_utils.get_z_score_and_ztest_pval(test_vals.iloc[0], control_vals)
                else:
                    assert nunique_control_vals == 1
                    z_score, ztest_pval = generic_utils.get_z_score_and_ztest_pval(control_vals.iloc[0], test_vals)
                    z_score = -z_score
                
                flat_dicts.append({
                    'column_name': column_name,
                    'z_score': z_score,
                    'ztest_pval': ztest_pval,
                })
            else:
                pval_column_name = 'welch_t_test_pval'
                t_statistic, welch_t_test_pval = scipy.stats.ttest_ind(
                    test_vals,
                    control_vals,
                    equal_var=False,
                    alternative='two-sided',
                )
                flat_dicts.append({
                    'column_name': column_name,
                    'welch_t_test_pval': welch_t_test_pval,
                    't_statistic': t_statistic,
                })

        df = pd.DataFrame(flat_dicts)
        df.sort_values(pval_column_name, inplace=True)
    else:
        non_negative_mask = (mc_ad.obs[column_names] >= 0).all(axis=0)
        raise NotImplementedError('can use log ratio for non_negative ones, and t test for the rest. sounds less convenient to me')
        # df = pd.concat(
        #     pd.DataFrame({
        #         'column_name': column_names[non_negative_mask], 
        #         'log_ratio': [np.log2(
        #             (mc_ad.obs.loc[test_mask, x].mean() + regularization_epsilon) / 
        #             (mc_ad.obs.loc[control_mask, x].mean() + regularization_epsilon)
        #         ) for x in column_names[non_negative_mask]]
        #     }),
        #     pd.DataFrame({
        #         'column_name': column_names[~non_negative_mask], 
        #         't_test_pval': [np.log2(
        #             (mc_ad.obs.loc[test_mask, x].mean() + regularization_epsilon) / 
        #             (mc_ad.obs.loc[control_mask, x].mean() + regularization_epsilon)
        #         ) for x in column_names[non_negative_mask]]
        #     }),
        # )
        # df.sort_values('log_ratio', ascending=False, inplace=True)

    return df

def get_var_projected_correlation(mc_ad):
    observed = mc_ad.X
    projected = mc_ad.layers['projected']
    return mc.ut.pairs_corrcoef_rows(
        mc.ut.to_layout((observed / observed.sum(axis=1)[:, np.newaxis]).transpose(), layout='row_major'),
        mc.ut.to_layout((projected / projected.sum(axis=1)[:, np.newaxis]).transpose(), layout='row_major'),
        reproducible=True,
    )

def get_donor_id_to_exp_sharing_donor_ids(c_ad):
    donor_exp_df = c_ad.obs[['exp_name', 'donor_id']].drop_duplicates().reset_index(drop=True)

    donor_id_to_exp_sharing_donor_ids = {}
    for donor_id in sorted(donor_exp_df['donor_id'].unique()):
        exp_names = donor_exp_df.loc[donor_exp_df['donor_id'] == donor_id, 'exp_name'].unique().tolist()
        donor_id_to_exp_sharing_donor_ids[donor_id] = donor_exp_df.loc[donor_exp_df['exp_name'].isin(exp_names), 'donor_id'].unique().tolist()
        
    return donor_id_to_exp_sharing_donor_ids


def get_donor_id_cell_fraction_of_mc_mask_df(mc_ad, mc_mask, donor_ids):
    mc_obs = mc_ad.obs[mc_mask]
    donor_fraction_column_names = [f'donor_id_{x}_fraction' for x in donor_ids]
    all_donor_all_mc_total_num_of_cells = (mc_obs[donor_fraction_column_names].sum(axis=1) * mc_obs['grouped']).sum()

    flat_dicts = []
    for donor_id in donor_ids:
        num_of_donor_cells = (mc_obs[f'donor_id_{donor_id}_fraction'] * mc_obs['grouped']).sum()
        frac_of_mcs_of_interest = num_of_donor_cells / all_donor_all_mc_total_num_of_cells
        
        flat_dicts.append({
            'donor_id': donor_id,
            'num_of_donor_cells': num_of_donor_cells,
            'frac_of_mcs_of_interest': frac_of_mcs_of_interest,
        })
    df = pd.DataFrame(flat_dicts)
    df.sort_values('frac_of_mcs_of_interest', ascending=False, inplace=True)
    return df

def get_donor_id_to_multiple_exp_names_and_all_donors_with_single_exp_ids(list_of_exp_name_and_donor_id):
    all_donor_exp_df = pd.DataFrame(data=list_of_exp_name_and_donor_id, columns=['exp_name', 'donor_id'])
    
    all_donor_ids = sorted(set(all_donor_exp_df['donor_id'].unique()) - {'nan'})
    donor_id_to_exp_names = {
        donor_id: all_donor_exp_df[all_donor_exp_df['donor_id'] == donor_id]['exp_name'].tolist()
        for donor_id in all_donor_ids
        if donor_id != 'nan'
    }
    donor_id_to_multiple_exp_names = {k: v for k,v in donor_id_to_exp_names.items() if len(v) >= 2}
    all_donors_with_single_exp_ids = sorted(x for x in donor_id_to_exp_names if len(donor_id_to_exp_names[x]) == 1)
    assert not (set(donor_id_to_multiple_exp_names) & set(all_donors_with_single_exp_ids))
    assert (set(donor_id_to_multiple_exp_names) | set(all_donors_with_single_exp_ids)) == set(all_donor_ids)
    return donor_id_to_multiple_exp_names, all_donors_with_single_exp_ids


def draw_mc_model_projected_correlation_plots(
        mc_ad, atlas_ad, cmap, 
        out_dir_path=None, out_file_name_suffix='', differentiated_cell_types=None,
        my_projected_correlation_column_name='my_projected_correlation',
):
    fig, ax = plt.subplots()
    non_nan_mask = ~np.isnan(mc_ad.obs[my_projected_correlation_column_name])
    rho = scipy.stats.pearsonr(
        mc_ad.obs.loc[non_nan_mask, my_projected_correlation_column_name], 
        mc_ad.obs.loc[non_nan_mask, 'projected_correlation'],
    )[0]
    sb.scatterplot(
        data=mc_ad.obs,
        x=my_projected_correlation_column_name, 
        y='projected_correlation',
        ax=ax,
        hue='state',
        palette=get_palette(mc_ad, 'state'),
        legend=False,
    )
    ax.set_title(f'rho={rho:.2f}')
    generic_utils.plot_y_equals_x_line_on_ax(ax)
    if out_dir_path is not None:
        fig.savefig(os.path.join(out_dir_path, f'my_vs_orenbk_projected_correlation_scatter{out_file_name_suffix}.png'))


    # plt.close('all')
    fig, ax = plt.subplots(figsize=(6, 5))
    sb.scatterplot(x=atlas_ad.obs['umap_x'], y=atlas_ad.obs['umap_y'], color='#cccccc', s=10, ax=ax)

    hue_norm = (0.8, 1)
    norm = plt.Normalize(*hue_norm)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    sb.scatterplot(
        data=mc_ad.obs,
        x='atlas_umap_x',
        y='atlas_umap_y',
        # hue='projected_correlation', # NOTE: i think this is not a good idea when we compare different mc models that were projected, because the projection algorithm might have used different genes when calculating projected_correlation
        hue=my_projected_correlation_column_name,
        hue_norm=hue_norm,
        palette=cmap,
        # vmax=1,
        ax=ax,
        s=25,
        legend=False,
        # xticklabels=[],
        # yticklabels=[],
    )
    fig.colorbar(sm)
    generic_utils.make_all_spines_and_x_and_y_axes_invisible(ax)
    fig.tight_layout()
    if out_dir_path is not None:
        fig.savefig(os.path.join(out_dir_path, f'my_projected_correlation_on_atlas_umap{out_file_name_suffix}.png'))

    metacell_mask_and_descs = [
        (np.full(mc_ad.n_obs, True), 'all'),
    ]
    if differentiated_cell_types is not None:
        metacell_mask_and_descs.append(
            ((~mc_ad.obs['state'].isin(differentiated_cell_types)), 'non_differentiated'))

    for mc_mask, mask_desc in metacell_mask_and_descs:
        fig, ax = plt.subplots(figsize=(4, 3))
        curr_mc_obs = mc_ad.obs[mc_mask]
        xs = curr_mc_obs[my_projected_correlation_column_name]
        median_x = xs.median()
        sb.histplot(
            xs,
            bins=np.linspace(0, 1, 51),
            stat='proportion',
            # common_norm=False,
            ax=ax,
        )
        ax.set_title(f'{mask_desc}, median={median_x:.2f}')
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.01, 1.01)
        fig.tight_layout()
        if out_dir_path is not None:
            fig.savefig(os.path.join(out_dir_path, f'my_projected_correlation_hist_{mask_desc}{out_file_name_suffix}.png'))

def get_list_of_exp_name_and_donor_id(c_ad, c_mask=None):
    if c_mask is None:
        c_mask = np.full(c_ad.n_obs, True)
    return c_ad.obs.loc[c_mask, ['exp_name', 'donor_id']].drop_duplicates().to_records(index=False).tolist()

def get_gene_gc_correction_df(
        mc_umi_count_mat, mc_proj, gene_names, 
        gene_read_gc_df, num_of_gc_bins, epsilon_to_add_for_log, before_or_after_correction, 
        out_fig_dir_path, use_median_for_gc_bin_correction=True,
        gene_mask=None,
):
    # mc_x_sum = mc_umi_count_mat.sum()
    # mc_proj_sum = mc_proj.sum()
    # print(f'mc_x_sum, mc_proj_sum: {mc_x_sum, mc_proj_sum}')
    # raise
    if gene_mask is None:
        gene_mask = np.full(mc_umi_count_mat.shape[1], True)
    # print(epsilon_to_add_for_log)
    # print(mc_umi_count_mat)
    # raise
    df = pd.DataFrame({
        'gene': gene_names[gene_mask],
        'total_observed': mc.ut.to_numpy_vector(mc_umi_count_mat[:, gene_mask].sum(axis=0)),
        'total_projected': mc.ut.to_numpy_vector(mc_proj[:, gene_mask].sum(axis=0)),
    })
    total_observed = df['total_observed'].sum()
    total_projected = df['total_projected'].sum()
    # print(f'total_observed: {total_observed}')
    # print(f'total_projected: {total_projected}')
    # print(f'total_projected/total_observed: {total_projected/total_observed}')
    df['norm_total_observed'] = df['total_observed'] / total_observed
    df['norm_total_projected'] = df['total_projected'] / total_projected

    df['mean_gc_fraction'] = generic_utils.merge_preserving_df1_index_and_row_order(
        pd.DataFrame({'gene': gene_names[gene_mask]}), 
        gene_read_gc_df[['gene', 'median_illumina_mean_read2_gc_fraction']],
        how='left',
    )['median_illumina_mean_read2_gc_fraction'].to_numpy()

    na_mean_gc_fraction_mask = df['mean_gc_fraction'].isna() 
    # assert not na_mean_gc_fraction_mask.any() # NOTE: for some reason, this is always true...

    df['log_norm_total_observed'] = np.log2(df['norm_total_observed'] + epsilon_to_add_for_log)
    df['log_norm_total_projected'] = np.log2(df['norm_total_projected'] + epsilon_to_add_for_log)

    # print(df)
    # raise

    fig, ax = plt.subplots()
    sb.scatterplot(
        data=df,
        x='log_norm_total_observed',
        y='log_norm_total_projected',
        hue='mean_gc_fraction',
        ax=ax,
        s=10,
    )
    generic_utils.plot_y_equals_x_line_on_ax(ax)
    fig.savefig(os.path.join(out_fig_dir_path, f'log_norm_total_observed_vs_projected_{before_or_after_correction}_correction.png'))

    df['log_norm_total_diff'] = df['log_norm_total_observed'] - df['log_norm_total_projected']
    fig, ax = plt.subplots()
    sb.histplot(df['log_norm_total_diff'], ax=ax, bins=50)
    ax.set_yscale('log')
    fig.savefig(os.path.join(out_fig_dir_path, f'log_norm_total_diff_hist_{before_or_after_correction}_correction.png'))

    gc_bin_i_of_genes = generic_utils.equal_size_bin_digitize_by_single_nums_vec(
        df.loc[~na_mean_gc_fraction_mask, 'mean_gc_fraction'], num_of_gc_bins)
    df.loc[~na_mean_gc_fraction_mask, 'gc_bin_i'] = gc_bin_i_of_genes
    df['gc_bin_i_as_str'] = df['gc_bin_i'].astype(str)
    fig, ax = plt.subplots()
    sb.scatterplot(
        data=df,
        x='mean_gc_fraction',
        y='log_norm_total_diff',
        # hue='gc_bin_i',
        hue='gc_bin_i_as_str',
        ax=ax,
        s=5,
        legend=False,
    )
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    fig.savefig(os.path.join(out_fig_dir_path, f'log_norm_total_diff_vs_gc_{before_or_after_correction}_correction.png'))


    df_grouped_by_gc_bin = df.groupby('gc_bin_i')
    if use_median_for_gc_bin_correction:
        gc_bin_log_correction_factor_df = df_grouped_by_gc_bin['log_norm_total_diff'].median().reset_index(name='gc_bin_log_correction_factor').merge(
            df_grouped_by_gc_bin['mean_gc_fraction'].median().reset_index(name='gc_bin_avg_mean_gc_fraction'))
    else:
        gc_bin_log_correction_factor_df = df_grouped_by_gc_bin['log_norm_total_diff'].mean().reset_index(name='gc_bin_log_correction_factor').merge(
            df_grouped_by_gc_bin['mean_gc_fraction'].mean().reset_index(name='gc_bin_avg_mean_gc_fraction'))

    fig, ax = plt.subplots()
    sb.scatterplot(
        data=gc_bin_log_correction_factor_df,
        x='gc_bin_avg_mean_gc_fraction',
        y='gc_bin_log_correction_factor',
        ax=ax,
    )
    fig.savefig(os.path.join(out_fig_dir_path, f'gc_vs_correction_per_bin_{before_or_after_correction}_correction.png'))
    # print(gc_bin_log_correction_factor_df)
    gc_bin_log_correction_factor_df.to_csv(os.path.join(out_fig_dir_path, f'gc_bin_log_correction_factor_df_{before_or_after_correction}_correction.csv'), index=False)


    gc_bin_log_correction_factor_df['gc_bin_correction_factor'] = 2 ** gc_bin_log_correction_factor_df['gc_bin_log_correction_factor']
    df = generic_utils.merge_preserving_df1_index_and_row_order(df, gc_bin_log_correction_factor_df, how='left')
    df['gc_bin_correction_factor'].fillna(1, inplace=True)
    # assert not df['gc_bin_correction_factor'].isna().any()
    return df



def plot_stripplot_per_obs_column(
        c_or_mc_ad, df=None, column_name=None, gene_names=None, 
        group_obs_col='state', 
        hue_obs_col='state', 
        palette=None, 
        ordered_group_obs_cols=None, c_or_mc_mask=None, 
        min_group_size=1,
        min_group_size_for_krusk=3,
        group_obs_col_vals_excluded_from_krusk=None,
        **stripplot_kwargs,
):
    fig, ax = plt.subplots(figsize=(10, 5))

    if df is None:
        if c_or_mc_mask is None:
            c_or_mc_mask = np.full(c_or_mc_ad.n_obs, True)
        df = c_or_mc_ad.obs[c_or_mc_mask]

    nan_mask = df[column_name].isna()
    if nan_mask.any():
        print(f'dropping {nan_mask.sum()} rows with nan {column_name}')
        df = df[~nan_mask]

    if column_name is None:
        assert gene_names is not None
        column_name = get_gene_desc(c_or_mc_ad, gene_names)
        df = df[[group_obs_col]].copy()
        df[column_name] = get_genes_expr(c_or_mc_ad, gene_names)

    df = df[df[group_obs_col].isin((df[group_obs_col].value_counts() >= min_group_size).loc[lambda x: x].index)]


    common_kwargs = dict(
        data=df, 
        x=group_obs_col, 
        y=column_name, 
        dodge=False,
        ax=ax,
    )
    # print(get_palette(c_or_mc_ad, color_by=hue_obs_column_name))
    # df = df.copy()
    # df[hue_obs_column_name] = df[hue_obs_column_name].astype(str)
    # print(df[hue_obs_column_name].unique())
    if ordered_group_obs_cols:
        common_kwargs['order'] = [x for x in ordered_group_obs_cols if x in df[group_obs_col].unique()]
    else:
        common_kwargs['order'] = df.groupby(group_obs_col)[column_name].median().sort_values().index.tolist()
    
    df_for_krusk = df if (group_obs_col_vals_excluded_from_krusk is None) else df[~df[group_obs_col].isin(group_obs_col_vals_excluded_from_krusk)]
    df_for_krusk = df_for_krusk[df_for_krusk[group_obs_col].isin(
        (df_for_krusk[group_obs_col].value_counts() >= min_group_size_for_krusk).loc[lambda x: x].index)]
    if df_for_krusk[group_obs_col].nunique() >= 2:
        krusk_stat, krusk_pval = scipy.stats.kruskal(*df_for_krusk.groupby(group_obs_col)[column_name].apply(list))
        ax.set_title(f'krusk_pval={krusk_pval:.2e}')
    
    sb.stripplot(
        jitter=0.25,
        linewidth=0.5,
        hue=hue_obs_col, 
        palette=palette if palette else get_palette(c_or_mc_ad, color_by=hue_obs_col),
        **common_kwargs,
        **stripplot_kwargs,
    )
    sb.boxplot(
        # fliersize=3,
        # flierprops=dict(alpha=0.5),
        showfliers=False,
        showbox=False,
        showcaps=False,
        # meanprops={'visible': False},
        medianprops={'visible': True},
        # medianprops={'color': 'k', 'ls': '-', 'lw': 1},
        whiskerprops={'visible': False},
        # notch=True,
        **common_kwargs,
    )
    ax.get_legend().remove()
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45, ha='right', 
        rotation_mode='anchor', 
        # fontsize='small',
    )
    ax.set_xlabel(None)
    fig.tight_layout()

    return fig, ax

def add_donor_bool_attr_to_cells_and_mc_ad(c_ad, mc_ad, donor_id_and_sample_date_df, bool_attr_name, sample_date_column_name='bleeding_date'):
    assert set(donor_id_and_sample_date_df.columns) == {'donor_id', sample_date_column_name}
    assert bool_attr_name not in c_ad.obs.columns
    donor_id_and_sample_date_df = donor_id_and_sample_date_df.copy()
    donor_id_and_sample_date_df[bool_attr_name] = True
    c_ad.obs = generic_utils.merge_preserving_df1_index_and_row_order(
        c_ad.obs, donor_id_and_sample_date_df, on=['donor_id', sample_date_column_name], how='left')
    c_ad.obs[bool_attr_name].fillna(False, inplace=True)
    add_mc_metadata(mc_ad, c_ad, {'bool_fraction': [bool_attr_name]}, add_num_of_cells_in_metacell=False)

def get_gene_expr_range_df(mc_ad, mc_mask=None):
    if mc_mask is None:
        mc_mask = np.full(mc_ad.n_obs, True)
    df = pd.DataFrame({
        'gene': mc_ad.var_names,
        'max_exp': mc_ad.layers['expr'][mc_mask, :].max(axis=0),
        'min_exp': mc_ad.layers['expr'][mc_mask, :].min(axis=0),
    })
    df['exp_range'] = df['max_exp'] - df['min_exp']
    df.sort_values('exp_range', ascending=False, inplace=True)
    return df

def get_cell_state_freq_df(
        c_ad=None,
        c_ad_obs=None,
        cell_mask=None,
        # group_by_column_names=['donor_id', 'exp_name'], 
        group_by_column_names=['numbered_donor_id'], 
        cell_state_column_name='state', 
        cell_states_to_ignore=None,
        add_count_column=False,
):
    if c_ad is not None:
        assert c_ad_obs is None
        c_ad_obs = c_ad.obs
    
    if cell_mask is None:
        cell_mask = np.full(len(c_ad_obs), True)
    df = c_ad_obs.loc[cell_mask, group_by_column_names + [cell_state_column_name]].copy()
    
    if 0:
        # 240229: i don't think this should be here... what if we asked for c_state? could have an argument saying whether to set state to Outliers...
        outlier_mask = c_ad_obs['metacell_name'] == 'Outliers'
        assert not df.loc[~outlier_mask, cell_state_column_name].isna().any()
        df.loc[outlier_mask, cell_state_column_name] = 'Outliers'
    df[cell_state_column_name] = df[cell_state_column_name].astype(str)
    
    if cell_states_to_ignore:
        df = df[~df[cell_state_column_name].isin(cell_states_to_ignore)]
    df = df.astype(str)

    if group_by_column_names:
        grouped_df = df.groupby(group_by_column_names)
    else:
        grouped_df = df

    df = grouped_df[cell_state_column_name].value_counts(normalize=True).reset_index(name='freq').rename(columns={'index': 'state'})
    if add_count_column:
        df = generic_utils.merge_preserving_df1_index_and_row_order(
            df, grouped_df[cell_state_column_name].value_counts().reset_index(name='count').rename(columns={'index': 'state'}))

    return df

def get_umi_count_mat(mc_ad=None, mc_x=None, mc_mask=None, gene_mask=None):
    if mc_x is None:
        mc_x = mc_ad.X
    if mc_mask is None:
        mc_mask = np.full(mc_x.shape[0], True)
    if gene_mask is None:
        gene_mask = np.full(mc_x.shape[1], True)
    if scipy.sparse.issparse(mc_x) and (mc_x.getformat() == 'coo'):
        mc_x = mc_x.tocsr()
    # 240725: umi_count_mat is currently not an integer matrix, as current metacell version performs geometric mean between its cells to obtain fraction per gene.
    umi_count_mat = mc.ut.to_numpy_matrix(mc_x[np.ix_(mc_mask, gene_mask)]) * mc_ad.obs.loc[mc_mask, 'total_umis'].to_numpy()[:, np.newaxis]
    return umi_count_mat



def get_top_genes_dfs(mc_ad, mc_mask=None, gene_mask=None, num_of_top_enriched_genes_per_mc=2, mc_mask_after_expr_enrich_calc=None, min_expr_enrich=1):
    if mc_mask is None:
        mc_mask = np.full(mc_ad.n_obs, True)

    if gene_mask is None:
        gene_mask = np.full(mc_ad.n_vars, True)
    gene_names = mc_ad.var_names[gene_mask]

    expr = mc_ad.layers['expr'][np.ix_(mc_mask, gene_mask)]
    expr_enrich = expr - np.median(expr, axis=0)
    expr = None

    if mc_mask_after_expr_enrich_calc is None:
        mc_mask_after_expr_enrich_calc = mc_mask
    else:
        assert not (mc_mask_after_expr_enrich_calc & (~mc_mask)).any()
        expr_enrich = expr_enrich[mc_mask_after_expr_enrich_calc[mc_mask], :]

    partition_indices = np.argpartition(-expr_enrich, num_of_top_enriched_genes_per_mc, axis=1)

    flat_dicts = []
    for mc_i, mc_name in enumerate(mc_ad.obs_names[mc_mask_after_expr_enrich_calc]):
        curr_expr_enrich = expr_enrich[mc_i,:]
        curr_partition_indices = partition_indices[mc_i,:][:num_of_top_enriched_genes_per_mc]
        
        # because partition doesnt sort the top k elements, we need to sort them ourselves
        curr_top_gene_indices = curr_partition_indices[np.argsort(-curr_expr_enrich[curr_partition_indices])]

        flat_dict = {
            'metacell_name': mc_name,
        }
        for i, top_enriched_gene_i in enumerate(curr_top_gene_indices):
            top_enriched_gene_name = gene_names[top_enriched_gene_i]
            top_enriched_gene_expr_enrich = curr_expr_enrich[top_enriched_gene_i]
            if top_enriched_gene_expr_enrich < min_expr_enrich:
                continue
            flat_dict = {
                **flat_dict,
                f'top{i + 1}_enriched_gene_name': top_enriched_gene_name,
                f'top{i + 1}_enriched_gene_expr_enrich': top_enriched_gene_expr_enrich,
            }
        if len(flat_dict) > 1:
            flat_dicts.append(flat_dict)
    mc_top_genes_df = pd.DataFrame(flat_dicts)
    top_genes_df = pd.concat(
        [
            mc_top_genes_df[f'top{i + 1}_enriched_gene_name']
            for i in range(num_of_top_enriched_genes_per_mc)
        ],
        ignore_index=True,
    ).value_counts().reset_index(name='num_of_top_enriched_mcs')
    column_to_replace_name = 'top1_enriched_gene_name' if num_of_top_enriched_genes_per_mc == 1 else 'index'
    top_genes_df.rename(columns={column_to_replace_name: 'gene'}, inplace=True)
    return mc_top_genes_df, top_genes_df

def add_top_gene_of_any_mc(mc_ad, min_expr_enrich=1, gene_mask=None):
    mc_top_genes_df, top_genes_df = get_top_genes_dfs(mc_ad, min_expr_enrich=min_expr_enrich, gene_mask=gene_mask)
    mc_ad.var['top_gene_of_any_mc'] = mc_ad.var_names.isin(top_genes_df['gene'])

def plot_obs_column_heatmap(
        c_or_mc_ad, col_names, c_or_mc_mask=None,
        row_color_obs_column_names=('c_state_color', 'state_color', 'projected_type_color', 'projected_correlation'),
        extra_row_color_obs_column_names=('log_num_of_non_excluded_umis', 'log_norm_num_of_cell_ranger_reported_umis_mean',),
        max_obs_count=int(10e3),
        metric='correlation',
        vmin=-3,
        vmax=3,
        **clustermap_kwargs,
):
    if c_or_mc_mask is None:
        c_or_mc_mask = np.full(c_or_mc_ad.n_obs, True)
    
    obs_count = c_or_mc_mask.sum()
    if obs_count > max_obs_count:
        raise RuntimeError(f'obs_count ({obs_count}) > max_obs_count ({max_obs_count})')
    assert obs_count > 1, f'obs_count: {obs_count}'
    
    df = c_or_mc_ad.obs.loc[c_or_mc_mask, col_names]
    df = df - df.median(axis=0)

    num_of_extra_cols = 0
    col_color_vecs = []
    for col_name in row_color_obs_column_names + extra_row_color_obs_column_names:
        if col_name not in c_or_mc_ad.obs.columns:
            print(f'\nNOTE: {col_name} is missing from c_or_mc_ad.obs\n')
            continue
        vals = c_or_mc_ad.obs.loc[c_or_mc_mask, col_name]
        # if vals.dtype.name in {'float32', 'float64', 'int32', 'int64'}: # NOTE: changed on 240501 to use dtype.kind
        if vals.dtype.kind in 'iuf':
            color_vec = generic_utils.get_num_colors(list(vals)) + ['white'] * num_of_extra_cols
        else:
            color_vec = list(vals) + ['white'] * num_of_extra_cols
        col_color_vecs.append(color_vec)

    clustermap_obj = sb.clustermap(
        data=df.transpose(),
        xticklabels=False,
        # yticklabels=gene_names,
        yticklabels=col_names,
        col_colors=col_color_vecs if col_color_vecs else None,
        # cmap='RdBu_r',
        cmap='bwr',
        vmin=vmin,
        vmax=vmax,
        metric=metric,
        **clustermap_kwargs,
    )
    return clustermap_obj



def plot_marker_heatmap(
        mc_ad, mc_mask=None,
        num_of_top_enriched_genes_per_mc=2,
        row_color_obs_column_names=(
            'state_color', 'projected_type_color', 
            'projected_correlation', 
            # 'log_norm_cell_ranger_umi_count_mean',
        ),
        extra_row_color_obs_column_names=(),
        use_these_genes_instead_of_top_genes=None,
        genes_to_include_anyway=('CD34',),
        name_to_c_mask=None,
        name_to_mc_mask_to_agg=None,
        mask_c_mask_size=2000,
        c_mask_to_norm_by=None,
        c_ad=None,
        min_expr_enrich=1,
        mcs_to_mark_mask=None,
        add_prefixes_to_gene_names=True,
        show_mc_names=True,
        show_gene_names=True,
        # clustering_metric='correlation',
        clustering_metric='euclidean', # scipy says "Method 'ward' requires the distance metric to be Euclidean", so i listen to scipy...
        clustering_method='ward', # added only on 240407
        write_ordered_genes_and_mc_names_to_temp_files=True,
        add_binary_color_row_for_added_pooled_c_masks=True,
        agg_by_val_of_col=None,
        equal_contrib_for_each_val_of_col=None,
        equal_contrib_for_each_val_of_col_min_mc_count=3,
        equal_contrib_for_each_val_of_col_max_mc_count=20,
        vmin=-3,
        vmax=3,
        **clustermap_kwargs,
):
    if mc_mask is None:
        mc_mask = np.full(mc_ad.n_obs, True)
    
    mc_count = mc_mask.sum()
    assert mc_count > 1, f'mc_count: {mc_count}'
    
    if equal_contrib_for_each_val_of_col:
        mc_ad_obs_col_as_str = mc_ad.obs[equal_contrib_for_each_val_of_col].astype(str)
        count_series = mc_ad_obs_col_as_str[mc_mask].value_counts()
        rare_val_mask = count_series <= equal_contrib_for_each_val_of_col_min_mc_count
        rare_vals = count_series[rare_val_mask].index

        if rare_val_mask.all():
            mc_mask &= mc_ad_obs_col_as_str.isin(rare_vals)
        else:
            non_rare_vals = count_series[~rare_val_mask].index
            min_common_val_mc_count = count_series[~rare_val_mask].iloc[-1]

            mc_mask &= (
                mc_ad_obs_col_as_str.isin(rare_vals)
                | generic_utils.get_equal_contrib_mask(
                    [(mc_ad_obs_col_as_str == x) & mc_mask for x in non_rare_vals], 
                    min(min_common_val_mc_count, equal_contrib_for_each_val_of_col_max_mc_count)
                )
            )
    
    if use_these_genes_instead_of_top_genes is not None:
        top_gene_names = use_these_genes_instead_of_top_genes
    else:
        mc_top_genes_df, top_genes_df = get_top_genes_dfs(
            mc_ad, mc_mask=mc_mask, num_of_top_enriched_genes_per_mc=num_of_top_enriched_genes_per_mc, min_expr_enrich=min_expr_enrich)
        top_gene_names = list(top_genes_df['gene'])
        
    plt.close('all')
    genes_to_include_anyway_also_in_top_genes = list(set(top_gene_names) & set(genes_to_include_anyway))
    gene_mask = mc_ad.var_names.isin(list(set(top_gene_names) | set(genes_to_include_anyway)))
    expr_mat = mc_ad.layers['expr'][np.ix_(mc_mask, gene_mask)]

    if agg_by_val_of_col:
        assert name_to_mc_mask_to_agg is None
        mc_ad_obs_col_as_str = mc_ad.obs[agg_by_val_of_col].astype(str)
        name_to_mc_mask_to_agg = {x: mc_ad_obs_col_as_str == x for x in set(mc_ad_obs_col_as_str)}

    agg_mc_names = []
    if name_to_mc_mask_to_agg:
        agg_expr_mat_vecs = []
        all_agg_mc_mask = np.full(mc_ad.n_obs, False)
        for name, mc_mask_to_agg in name_to_mc_mask_to_agg.items():
            # print(name, mc_mask_to_agg.sum())
            assert not (all_agg_mc_mask & mc_mask_to_agg).any()
            all_agg_mc_mask |= mc_mask_to_agg

            masked_mc_mask_to_agg = mc_mask_to_agg[mc_mask]
            if masked_mc_mask_to_agg.any():
                agg_mc_names.append(name)
                agg_expr_mat_vecs.append(np.median(expr_mat[masked_mc_mask_to_agg], axis=0))
        
        expr_mat = np.vstack([expr_mat[~all_agg_mc_mask,:], *agg_expr_mat_vecs])
        mc_mask &= ~all_agg_mc_mask
    
    median_expr = np.median(expr_mat, axis=0)
    enriched_expr_mat = expr_mat - median_expr
    

    gene_names = mc_ad.var_names[gene_mask]
    mc_names = list(mc_ad.obs_names[mc_mask]) + agg_mc_names
    enriched_expr_df = pd.DataFrame(enriched_expr_mat, columns=gene_names, index=mc_names).T

    if name_to_c_mask:
        assert c_ad is not None, 'c_ad argument is missing'
        assert (c_ad.var_names == mc_ad.var_names).all()
        for name, mask in name_to_c_mask.items():
            enriched_expr_df[name] = mc.ut.to_numpy_vector(get_mean_expr(
                c_ad, generic_utils.sample_mask(mask, mask_c_mask_size, allow_not_sampling_due_to_too_few_pos=True), layer_name='X',
                use_geomean=True, epsilon_for_log=1e-5, epsilon_for_c_log=1/16,
            )[gene_mask] - median_expr)
        named_c_mask_count = len(name_to_c_mask)
    else:
        named_c_mask_count = 0
    
    
    if show_gene_names:
        if add_prefixes_to_gene_names:
            selected_gene_names = sorted(mc_ad.var_names[mc_ad.var['selected_gene']])
            non_selected_prefixes = ['' if x in selected_gene_names else '**' for x in gene_names]
            genes_to_include_anyway_prefixes = [
                '#' if x in genes_to_include_anyway_also_in_top_genes else ('##' if x in genes_to_include_anyway else '') for x in gene_names]
            y_tick_labels = [x[0] + x[1] + x[2] for x in zip(non_selected_prefixes, genes_to_include_anyway_prefixes, gene_names)]
        else:
            y_tick_labels = gene_names
    else:
        y_tick_labels = False


    col_color_vecs = []
    col_colors_dict = {}
    for col_name in row_color_obs_column_names + extra_row_color_obs_column_names:
        if isinstance(col_name, tuple):
            col_name, curr_vmin, curr_vmax = col_name
        else:
            curr_vmin = curr_vmax = None
        if col_name not in mc_ad.obs.columns:
            print(f'\nNOTE: {col_name} is missing from mc_ad.obs\n')
            continue
        
        vals = mc_ad.obs.loc[mc_mask, col_name]
        vals_kind = vals.dtype.kind
        vals = list(vals)
        # if vals.dtype.name in {'float32', 'float64', 'int32', 'int64'}: # NOTE: changed on 240501 to use dtype.kind
        if vals_kind in 'iuf':
            if agg_mc_names:
                vals += [mc_ad.obs.loc[name_to_mc_mask_to_agg[name], col_name].median() for name in agg_mc_names]
            color_vec = generic_utils.get_num_colors(vals, vmin=curr_vmin, vmax=curr_vmax) + ['white'] * named_c_mask_count
        else:
            if agg_mc_names:
                vals += [mc_ad.obs.loc[name_to_mc_mask_to_agg[name], col_name].mode().iloc[0] for name in agg_mc_names]
            color_vec = vals + ['white'] * named_c_mask_count
        col_color_vecs.append(color_vec)
        col_colors_dict[col_name] = color_vec
    # raise
    if add_binary_color_row_for_added_pooled_c_masks and name_to_c_mask and (not show_mc_names):
        added_pooled_vec = ['white'] * len(mc_names) + ['black'] * named_c_mask_count
        col_color_vecs.append(added_pooled_vec)
        col_colors_dict['custom_pooled_cells'] = added_pooled_vec
    
    if col_colors_dict:
        col_colors_df = pd.DataFrame(col_colors_dict)
        col_colors_df.index = enriched_expr_df.columns
    else:
        col_colors_df = None


    if show_mc_names:
        if mcs_to_mark_mask is not None:
            xticklabels = [f'****{x}' if x in mc_ad.obs_names[mcs_to_mark_mask] else x for x in enriched_expr_df.columns]
        else:
            xticklabels = enriched_expr_df.columns
    else:
        xticklabels = False

    const_mask = enriched_expr_df.var(axis=1) == 0

    if const_mask.any():
        const_genes = list(enriched_expr_df.index[const_mask])
        print(f'const_genes: {const_genes}. not showing them')
        enriched_expr_df = enriched_expr_df[~const_mask]

    # return enriched_expr_df
    if 'figsize' not in clustermap_kwargs:
        y_size = 2.8 + 0.17 * len(enriched_expr_df) if show_gene_names else 15
        figsize = (10, y_size)
        print(f'clustermap figsize: {figsize}')
        clustermap_kwargs['figsize'] = figsize
    clustermap_obj = sb.clustermap(
        data=enriched_expr_df,
        xticklabels=xticklabels,
        # yticklabels=gene_names,
        yticklabels=y_tick_labels,
        # col_colors=col_color_vecs if col_color_vecs else None,
        col_colors=col_colors_df,
        # cmap='RdBu_r',
        cmap='bwr',
        vmin=vmin,
        vmax=vmax,
        metric=clustering_metric,
        method=clustering_method,
        # metric=metric,
        **clustermap_kwargs,
    )

    ordered_genes, ordered_mc_names = generic_utils.get_clustermap_row_and_col_names_ordered(enriched_expr_df, clustermap_obj)
    if write_ordered_genes_and_mc_names_to_temp_files:
        pathlib.Path('temp').mkdir(parents=True, exist_ok=True)
        generic_utils.write_text_file('temp/marker_heatmaps/ordered_genes.txt', "'" + "','".join(ordered_genes) + "'")
        generic_utils.write_text_file('temp/marker_heatmaps/ordered_mc_names.txt', "'" + "','".join(ordered_mc_names) + "'")
        

    return clustermap_obj

def get_norm_umi_count_geomean(
        x_mat, epsilon_for_log=None, umi_count_of_cells=None, x_already_normed=False, 
        norm_epsilon_for_log=1e-5, scale_final_geomean_to_sum_to_1=True,
        min_gene_count_to_scale_final_geomean_to_sum_to_1=500,
):
    # see metacells/pipeline/collect.py
    if umi_count_of_cells is not None:
        assert (umi_count_of_cells > 5).all(), f'umi_count_of_cells (expected ints): {umi_count_of_cells}'
    elif (not x_already_normed) and (umi_count_of_cells is None):
        umi_count_of_cells = mc.ut.to_numpy_vector(x_mat.sum(axis=1))
    
    if x_already_normed:
        # NOTE: for x_already_normed=True, the geomean will not be weighted.
        norm_umi_count = x_mat + norm_epsilon_for_log
    else:
        norm_umi_count = (x_mat + epsilon_for_log) / umi_count_of_cells[:, np.newaxis]
    weights = np.full(x_mat.shape[0], 1) if x_already_normed else np.log2(umi_count_of_cells)
    geomean_norm_umi_count = scipy.stats.mstats.gmean(norm_umi_count, weights=weights[:, np.newaxis], axis=0)
    if x_already_normed:
        geomean_norm_zero_umi_count = norm_epsilon_for_log
    else:
        norm_zero_umi_count = epsilon_for_log / umi_count_of_cells
        geomean_norm_zero_umi_count = scipy.stats.mstats.gmean(norm_zero_umi_count, weights=weights) # was wrong here on 230919. OrenBK confirmed and said he would fix.
    
    corrected_geomean_norm_umi_count = geomean_norm_umi_count - geomean_norm_zero_umi_count
    
    all_zeros_mask = mc.ut.to_numpy_vector((x_mat == 0).all(axis=0))
    assert np.isclose(corrected_geomean_norm_umi_count[all_zeros_mask], 0, atol=1e-6).all(), f'the following should be close to zero: {corrected_geomean_norm_umi_count[all_zeros_mask]}'
    corrected_geomean_norm_umi_count[all_zeros_mask] = 0
    
    under_zero_mask = corrected_geomean_norm_umi_count < 0
    assert np.isclose(corrected_geomean_norm_umi_count[under_zero_mask], 0, atol=1e-6).all(), f'the following should be close to zero: {corrected_geomean_norm_umi_count[under_zero_mask]}'
    corrected_geomean_norm_umi_count[under_zero_mask] = 0
    if scale_final_geomean_to_sum_to_1:
        gene_count = x_mat.shape[1]
        if gene_count < min_gene_count_to_scale_final_geomean_to_sum_to_1:
            raise RuntimeError(f'gene_count={gene_count} is lower than min_gene_count_to_scale_final_geomean_to_sum_to_1 ({min_gene_count_to_scale_final_geomean_to_sum_to_1})')
        corrected_geomean_norm_umi_count /= corrected_geomean_norm_umi_count.sum()
    return mc.ut.to_numpy_vector(corrected_geomean_norm_umi_count)

def get_cell_and_proj_mean_expr_df(c_ad, mc_ad, c_mask, epsilon_for_c_log=1/16, epsilon_for_log=1e-5):
    write_proj_expr(mc_ad, silent=True)

    orig_num_of_cells = c_mask.sum()
    c_mask &= c_ad.obs['metacell_name'] != 'Outliers'
    num_of_cells = c_mask.sum()
    if num_of_cells < orig_num_of_cells:
        print(f'ignoring {orig_num_of_cells - num_of_cells} outlier cells')
    
    metacell_name_counts = c_ad.obs.loc[c_mask, 'metacell_name'].astype(str).value_counts().to_frame(name='num_of_donor_cells_in_mc')
    assert metacell_name_counts.index.isin(mc_ad.obs_names).all()

    filtered_x = c_ad.X[c_mask, :]
    # filtered_x = c_ad.layers['downsampled'][c_mask, :] # this makes the effect much worse.
        
    filtered_x = mc.ut.to_numpy_matrix(filtered_x)

    corrected_geomean_norm_umi_count = get_norm_umi_count_geomean(filtered_x, epsilon_for_c_log)

    # if scipy.sparse.issparse(norm_umi_count):
    #     norm_umi_count = mc.ut.to_numpy_matrix(norm_umi_count)
    donor_c_mean_expr = mc.ut.to_numpy_vector(np.log2(corrected_geomean_norm_umi_count + epsilon_for_log))

    ordered_metacell_name_counts = generic_utils.merge_preserving_df1_index_and_row_order(mc_ad.obs[[]], metacell_name_counts, left_index=True, right_index=True, how='left').fillna(0)
    donor_mc_mean_expr = mc.ut.to_numpy_vector(
        (ordered_metacell_name_counts.T.to_numpy() @ mc_ad.layers['expr']) / ordered_metacell_name_counts.to_numpy().sum())
    donor_mc_mean_proj_expr = mc.ut.to_numpy_vector(
        (ordered_metacell_name_counts.T.to_numpy() @ mc_ad.layers['proj_expr']) / ordered_metacell_name_counts.to_numpy().sum())

    df = pd.DataFrame({
        'gene': mc_ad.var_names,
        'c_expr': donor_c_mean_expr,
        'mc_expr': donor_mc_mean_expr,
        'proj_expr': donor_mc_mean_proj_expr,
    })
    df['c_proj_log_ratio'] = df['c_expr'] - df['proj_expr']
    df['c_mc_log_ratio'] = df['c_expr'] - df['mc_expr']
    df.sort_values('c_proj_log_ratio', ascending=False, inplace=True)

    return df

def get_c_proj_log_ratio_df(
        c_ad, mc_ad, name_to_c_mask, 
        log_ratio_df_csv_file_path=None, 
        cell_count_df_csv_file_path=None,
        c_mc_log_ratio_df_csv_file_path=None,
        proj_expr_df_csv_file_path=None,
        c_expr_df_csv_file_path=None,
        general_c_mask=None, min_num_of_cells=300, 
        # gene_mask=None, 
        # min_max_c_proj_log_ratio=1, 
        # min_max_c_proj_log_ratio_range=1,
        epsilon_for_c_log=1/16, epsilon_for_log=1e-5, 
        out_dir_path=None, use_existing_df_csv_file=False, cache_dir_path=None,
):
    if out_dir_path is not None:
        pathlib.Path(out_dir_path).mkdir(parents=True, exist_ok=True)
    if cache_dir_path is not None:
        pathlib.Path(cache_dir_path).mkdir(parents=True, exist_ok=True)
    
    if general_c_mask is None:
        general_c_mask = np.full(c_ad.n_obs, True)

    # gene_passing_thresh_names = set()
    name_to_df = {}
    name_to_num_of_cells = {}
    names_analyzed = []
    num_of_names = len(name_to_c_mask)
    for i, (name, c_mask) in enumerate(name_to_c_mask.items()):
        # print(name)
        curr_c_mask = c_mask & general_c_mask
        num_of_cells = curr_c_mask.sum()
        name_to_num_of_cells[name] = num_of_cells
        if num_of_cells < min_num_of_cells:
            print(f'skipping {name} due to too low #cells ({num_of_cells})')
            continue
        names_analyzed.append(name)
        if i % 20 == 15:
            print(f'processing {name} ({num_of_cells} cells) ({i+1}/{num_of_names})')
            # break
        df_csv_file_path = None if cache_dir_path is None else os.path.join(cache_dir_path, f'{name}_cell_and_proj_mean_expr_df.csv')
        if (not use_existing_df_csv_file) or (cache_dir_path is None) or (not os.path.isfile(df_csv_file_path)):
            curr_df = get_cell_and_proj_mean_expr_df(
                c_ad, mc_ad, curr_c_mask, epsilon_for_c_log=epsilon_for_c_log, epsilon_for_log=epsilon_for_log)
            if df_csv_file_path:
                curr_df.to_csv(df_csv_file_path, index=False)
        else:
            curr_df = pd.read_csv(df_csv_file_path)
        
        # gene_passing_thresh_names |= (
        #     set(curr_df.loc[curr_df['c_proj_log_ratio'] >= min_max_c_proj_log_ratio, 'gene'].unique())
        #     | set(curr_df.loc[np.abs(curr_df['c_proj_log_ratio'] - curr_df['c_proj_log_ratio']) >= min_max_c_proj_log_ratio_range, 'gene'].unique())
        # )

        name_to_df[name] = curr_df
        # if name != 'N227':
        #     raise
    # gene_passing_thresh_names = sorted(gene_passing_thresh_names)

    # if gene_mask is None:
    #     gene_names = gene_passing_thresh_names
    # else:
    #     gene_names = sorted(mc_ad.var_names[gene_mask])

    gene_names = sorted(mc_ad.var_names)
    atlas_gene_names = mc_ad.var_names[mc_ad.var['atlas_gene']]
    gene_names = [x for x in gene_names if x in atlas_gene_names]

    log_ratio_df = pd.DataFrame({'gene': gene_names})
    c_mc_log_ratio_df = pd.DataFrame({'gene': gene_names})
    c_expr_df = pd.DataFrame({'gene': gene_names})
    proj_expr_df = pd.DataFrame({'gene': gene_names})
    for name, df in name_to_df.items():
        log_ratio_df = generic_utils.merge_preserving_df1_index_and_row_order(log_ratio_df, df[['gene', 'c_proj_log_ratio']].rename(
            columns={'c_proj_log_ratio': name}))
        c_mc_log_ratio_df = generic_utils.merge_preserving_df1_index_and_row_order(c_mc_log_ratio_df, df[['gene', 'c_mc_log_ratio']].rename(
            columns={'c_mc_log_ratio': name}))
        c_expr_df = generic_utils.merge_preserving_df1_index_and_row_order(c_expr_df, df[['gene', 'c_expr']].rename(
            columns={'c_expr': name}))
        proj_expr_df = generic_utils.merge_preserving_df1_index_and_row_order(proj_expr_df, df[['gene', 'proj_expr']].rename(
            columns={'proj_expr': name}))

    for df in (log_ratio_df, c_mc_log_ratio_df, c_expr_df, proj_expr_df):
        df.set_index('gene', inplace=True)
    
    if cell_count_df_csv_file_path:
        cell_count_df = pd.Series(name_to_num_of_cells).reset_index(name='cell_count').rename(columns={'index': 'name'})
        cell_count_df['analyzed'] = cell_count_df['name'].isin(names_analyzed)
        cell_count_df.to_csv(cell_count_df_csv_file_path, index=False)
    if log_ratio_df_csv_file_path:
        log_ratio_df.to_csv(log_ratio_df_csv_file_path)
    if c_mc_log_ratio_df_csv_file_path:
        c_mc_log_ratio_df.to_csv(c_mc_log_ratio_df_csv_file_path)
    if proj_expr_df_csv_file_path:
        proj_expr_df.to_csv(proj_expr_df_csv_file_path)
    if c_expr_df_csv_file_path:
        c_expr_df.to_csv(c_expr_df_csv_file_path)
    
    return log_ratio_df, cell_count_df

def single_cell_proj(
    c_ad,
    atlas_ad,
    c_mask,
    mc_ad=None,
    gene_names=None,
    num_of_top_atlas_mcs=5,
):
    if gene_names is None:
        if mc_ad is not None:
            assert (c_ad.var_names == mc_ad.var_names).all()
            # init_gene_mask = mc_ad.var['top_gene_of_any_mc'].to_numpy()
            init_gene_mask = (mc_ad.var['top_gene_of_any_mc'] & mc_ad.var['selected_gene']).to_numpy()
        else:
            init_gene_mask = c_ad.var['selected_gene'].to_numpy()

        gene_names = sorted(set(c_ad.var_names[init_gene_mask]) & set(atlas_ad.var_names))

    gene_mask = c_ad.var_names.isin(gene_names)
    ordered_gene_names = c_ad.var_names[gene_mask]
    atlas_var_names = list(atlas_ad.var_names)
    atlas_gene_indices = [atlas_var_names.index(x) for x in ordered_gene_names]
    atlas_ad_expr_ordered = atlas_ad.layers['expr'][:, atlas_gene_indices]

    c_log_norm_expr = get_log_norm_expression(ad_x=c_ad.X[c_mask, :], gene_mask=gene_mask, mask_genes_after_normalizing=True)
    # c_log_norm_expr = mc.ut.to_numpy_matrix(np.array([mebemp_m_big_expr[gene_mask]]))

    corr_mat = mc.ut.cross_corrcoef_rows(
        mc.ut.to_layout(c_log_norm_expr, layout='row_major'),
        mc.ut.to_layout(atlas_ad_expr_ordered, layout='row_major'),
        reproducible=True,
    )

    unordered_top_corrs = np.partition(corr_mat, -num_of_top_atlas_mcs, axis=1)[:,-num_of_top_atlas_mcs:]
    unordered_top_corr_indices = np.argpartition(corr_mat, -num_of_top_atlas_mcs, axis=1)[:,-num_of_top_atlas_mcs:]
    top_corr_indices = unordered_top_corr_indices[np.arange(unordered_top_corr_indices.shape[0])[:,np.newaxis], np.argsort(unordered_top_corrs, axis=1)][:, ::-1]
    top_corrs = corr_mat[np.arange(corr_mat.shape[0])[:,np.newaxis], top_corr_indices]

    assert (np.diff(top_corrs, axis=1) <= 0).all()

    c_proj_df = pd.DataFrame(
        {
            **{
                f'top{i+1}_atlas_mc_i': top_corr_indices[:, i]
                for i in range(num_of_top_atlas_mcs)
            },
            **{
                f'top{i+1}_atlas_mc_corr': top_corrs[:, i]
                for i in range(num_of_top_atlas_mcs)
            },
            'donor_id': c_ad.obs.loc[c_mask, 'donor_id'],
        },
        index=c_ad.obs_names[c_mask],
    )
    for i in range(num_of_top_atlas_mcs):
        c_proj_df[f'top{i+1}_atlas_mc_state'] = atlas_ad.obs['type'][c_proj_df[f'top{i+1}_atlas_mc_i']].to_numpy()

    c_proj_df['top_atlas_mc_state_mode'] = c_proj_df[[f'top{i+1}_atlas_mc_state' for i in range(num_of_top_atlas_mcs)]].mode(axis=1)[0]
    c_proj_df['num_of_top_atlas_mcs_with_state_mode'] = (
        c_proj_df[[f'top{i+1}_atlas_mc_state' for i in range(num_of_top_atlas_mcs)]] == c_proj_df['top_atlas_mc_state_mode'].to_numpy()[:,np.newaxis]).sum(axis=1)
    
    majority_count = num_of_top_atlas_mcs // 2 + 1
    
    print(f'num_of_top_atlas_mcs_with_state_mode >= majority_count ({majority_count})')
    print((c_proj_df['num_of_top_atlas_mcs_with_state_mode'] >= majority_count).value_counts())

    return c_proj_df
    if 0:
        plt.close('all')
        sb.histplot(c_proj_df['top_atlas_mc_corr'])
        sb.histplot(c_proj_df['num_of_top_atlas_mcs_with_state_mode'])

        c_proj_df.loc[c_proj_df['num_of_top_atlas_mcs_with_state_mode'] >= 3, 'top_atlas_mc_state_mode'].value_counts(normalize=True)

def get_exp_donor_count(c_ad):
    exp_donor_count = c_ad.obs.loc[(~c_ad.obs['donor_id'].isna()) & (c_ad.obs['donor_id'] != 'nan'), ['donor_id', 'exp_name']].astype(
        str).drop_duplicates()['exp_name'].value_counts()
    assert (exp_donor_count >= 1).all()
    return exp_donor_count

def get_single_donor_exp_names(c_ad):
    exp_donor_count = get_exp_donor_count(c_ad)
    return sorted((exp_donor_count == 1).loc[lambda x: x].index)

def get_exp_names_with_at_least_x_donors(c_ad, min_num_of_donors_identified_in_exp):
    exp_donor_count = get_exp_donor_count(c_ad)
    return sorted((exp_donor_count >= min_num_of_donors_identified_in_exp).loc[lambda x: x].index)

def get_log_ugly_norm_cell_ranger_umi_count_of_cells(c_ad, min_num_of_mebemp_m_and_mpp=20):
    if 'log_cell_ranger_umi_count' not in c_ad.obs.columns:
        c_ad.obs['log_cell_ranger_umi_count'] = np.log2(c_ad.obs['cell_ranger_umi_count'])
    assert not c_ad.obs['donor_id'].isna().any()
    c_ad.obs['donor_id'] = c_ad.obs['donor_id'].astype(str)
    c_ad.obs['exp_name'] = c_ad.obs['exp_name'].astype(str)
    assert (c_ad.obs['state'] == 'MEBEMP-E').any()
    assert (c_ad.obs['state'] == 'MPP').any()
    mebemp_m_and_mpp_mask = c_ad.obs['state'].isin(['MEBEMP-E', 'MPP']) & (c_ad.obs['donor_id'] != 'nan')
    donor_exps_with_sufficient_mebemp_m_and_mpp_cells = (c_ad.obs.loc[mebemp_m_and_mpp_mask, ['donor_id', 'exp_name']].value_counts() >= min_num_of_mebemp_m_and_mpp).loc[
        lambda x: x].index.to_frame().reset_index(drop=True)
    sufficient_mebemp_m_and_mpp_mask = c_ad.obs[['donor_id', 'exp_name']].isin(donor_exps_with_sufficient_mebemp_m_and_mpp_cells).all(axis=1) & mebemp_m_and_mpp_mask
    sufficient_mebemp_m_and_mpp_mask = generic_utils.get_mask_of_df1_rows_that_inner_join_will_keep(
        c_ad.obs[['donor_id', 'exp_name']], donor_exps_with_sufficient_mebemp_m_and_mpp_cells) & mebemp_m_and_mpp_mask

    grouped_df = c_ad.obs[sufficient_mebemp_m_and_mpp_mask].groupby(['donor_id', 'exp_name'])
    mebemp_m_and_mpp_median_log_umi_count_df = grouped_df['log_cell_ranger_umi_count'].median().reset_index(name='mebemp_m_and_mpp_median_log_cell_ranger_umi_count')

    log_norm_umi_count_df = generic_utils.merge_preserving_df1_index_and_row_order(
        c_ad.obs[['donor_id', 'exp_name', 'log_cell_ranger_umi_count']], mebemp_m_and_mpp_median_log_umi_count_df, on=['donor_id', 'exp_name'], how='left')
    log_norm_umi_count_df['log_norm_cell_ranger_umi_count'] = log_norm_umi_count_df['log_cell_ranger_umi_count'] - log_norm_umi_count_df['mebemp_m_and_mpp_median_log_cell_ranger_umi_count']
    return log_norm_umi_count_df['log_norm_cell_ranger_umi_count'].to_numpy()

def get_c_mask_freq_per_obs_column(c_ad, c_mask, obs_col_names, background_c_mask=None):
    if background_c_mask is None:
        background_c_mask = np.full(c_ad.n_obs, True)

    if (c_mask & (~background_c_mask)).any():
        print('WARNING (get_c_mask_freq_per_obs_column): c_mask is not a subset of background_c_mask')
        c_mask &= background_c_mask

    if isinstance(obs_col_names, str):
        obs_col_names = [obs_col_names]
    else:
        assert isinstance(obs_col_names, list), f'obs_column_names must be str or list, got {type(obs_col_names)}'

    if c_ad.obs.loc[background_c_mask, obs_col_names].isna().any().any():
        for col in obs_col_names:
            if c_ad.obs.loc[background_c_mask, col].isna().any():
                raise RuntimeError(f'c_ad.obs[{col}] has NaNs in background_c_mask')
    df = c_ad.obs.loc[background_c_mask, obs_col_names].value_counts().reset_index(name='background_count')
    df = generic_utils.merge_preserving_df1_index_and_row_order(
        df,
        c_ad.obs.loc[c_mask, obs_col_names].value_counts().reset_index(name='count'),
        on=obs_col_names,
        how='left',
    )
    df['count'].fillna(0, inplace=True)
    df['freq'] = df['count'] / df['background_count']
    df.sort_values('freq', ascending=False, inplace=True)
    if (len(obs_col_names) == 1) and (len(df) == 2) and (set(df[obs_col_names[0]]) == {True, False}):
        col_name = obs_col_names[0]
        true_freq = df.loc[df[col_name], 'freq'].iloc[0]
        false_freq = df.loc[~df[col_name], 'freq'].iloc[0]
        assert (true_freq > 0) or (false_freq > 0)
        if true_freq == 0:
            log_ratio = -np.inf
        elif false_freq == 0:
            log_ratio = np.inf
        else:
            log_ratio = np.log2(true_freq) - np.log2(false_freq)
        print(f'log_ratio: {log_ratio:.2f}')
    return df


def get_c_mask_fraction_of_mcs(mc_ad, c_ad, c_mask, background_c_mask=None):
    if background_c_mask is None:
        background_c_mask = np.full(c_ad.n_obs, True)
    elif not (c_mask & (~background_c_mask)).any():
        print('WARNING (get_c_mask_fraction_of_mcs): c_mask is not a subset of background_c_mask')
        c_mask &= background_c_mask
    numerator_vec = get_num_of_cells_per_metacell(mc_ad, c_ad, c_mask=c_mask)
    denom_vec = get_num_of_cells_per_metacell(mc_ad, c_ad, c_mask=background_c_mask)
    res_vec = np.zeros(len(numerator_vec))
    non_zero_denom_mask = denom_vec != 0
    res_vec[non_zero_denom_mask] = numerator_vec[non_zero_denom_mask] / denom_vec[non_zero_denom_mask]
    return res_vec
    
def get_same_exp_donor_pairs_df(df_with_exp_name_and_donor_id):
    donor_exp_df = df_with_exp_name_and_donor_id[['exp_name', 'donor_id']].drop_duplicates().astype(str)
    df = donor_exp_df.merge(donor_exp_df, on='exp_name', suffixes=('1', '2'))
    df = df[df['donor_id1'] < df['donor_id2']]
    return df

def get_pooled_expr_per_state(c_ad, gene_names, c_mask=None, min_cell_count=1, layer_name='downsampled', epsilon_for_log=1e-5):
    if c_mask is None:
        c_mask = np.full(c_ad.n_obs, True)
        
    flat_dicts = []
    for state in sorted(c_ad.obs['state'].astype(str).unique()):
        curr_c_mask = c_mask & (c_ad.obs['state'] == state)
        cell_count = curr_c_mask.sum()
        if cell_count < min_cell_count:
            continue
        flat_dicts.append({
            'state': state,
            'cell_count': cell_count,
            'expr': get_genes_pooled_log_total_norm_expression(
                c_ad, curr_c_mask, gene_names, downsample_layer_name=layer_name, epsilon_for_log=epsilon_for_log),
        })
    return pd.DataFrame(flat_dicts).sort_values('expr', ascending=False)

def get_donor_pooled_expr_per_state(c_ad, gene_names, donor_id, exp_name=None, min_cell_count=1):
    # if donor_id is None:
    #     c_mask = np.full(c_ad.n_obs, True)
    donor_c_mask = c_ad.obs['donor_id'] == donor_id
    assert donor_c_mask.any(), 'donor_id missing from c_ad'
    if exp_name is not None:
        exp_c_mask = c_ad.obs['exp_name'] == exp_name
        donor_c_mask &= exp_c_mask
    donor_pooled_expr_per_state = get_pooled_expr_per_state(c_ad, gene_names, c_mask=donor_c_mask, min_cell_count=min_cell_count)
    if exp_name is not None:
        exp_rest_pooled_expr_per_state = get_pooled_expr_per_state(
            c_ad, gene_names, c_mask=(exp_c_mask & (~donor_c_mask)), min_cell_count=min_cell_count)
        donor_pooled_expr_per_state = generic_utils.merge_preserving_df1_index_and_row_order(
            donor_pooled_expr_per_state, exp_rest_pooled_expr_per_state, on='state', how='left', suffixes=('', '_exp_rest'))
        donor_pooled_expr_per_state['expr_log_ratio'] = donor_pooled_expr_per_state['expr'] - donor_pooled_expr_per_state['expr_exp_rest']
        donor_pooled_expr_per_state.sort_values('expr_log_ratio', ascending=False, inplace=True)
    return donor_pooled_expr_per_state

def get_var_mean_df(
        c_or_mc_ad, c_or_mc_mask=None, min_umi_count=None, show_and_return_scatter_fig_and_ax=True, 
        min_var_mean_log_ratio_to_show_gene_name=None, min_mean_expr=1e-5, color_by='selected_gene',
):
    # NOTE: 240110: seems like this function only works for c_ad currently, which sounds fine to me.
    
    # i tried to simulate this in 231210_trying_to_understand_why_var_mean_ratio_higher_in_strong_genes.ipynb. our sequencing is sampling which is approximately binomial (IIUC). but what is the sampling that the cells perform? i tried to simulate their sampling as binomial, but didn't get a higher var/mean for stronger genes. when simulating their sampling as n*U(p-0.2p,p+0.2p) (n is the number of mRNA moledules in the cell) while p is the mean frequency of the gene, i got a higher var/mean for stronger genes.
    # so essentially, the variance should be larger than the mean because poisson/binomial is not a good approximation. there is truly more variance between cells. intuitively, the lower the mean, the higher the noise we get by our poisson/binomial (sequencing) sampling, and so our noise masks the true variation in frequencies between cells. taking this to the extreme, the noise completely masks the true variation, and we are left with a poisson/binomial distribution with mean=var. on the other hand, if the mean is higher, the noise from our sampling is smaller, and we can observe the true variation in frequencies between cells.
    # but why is the true variance higher in stronger genes? simply because of the true distribution, which we don't know, IIUC. that is, this mysterious distribution simply doesn't behave as poisson/binomial. for this mysterious distribution, the variance is larger for stronger genes. if, for example, the true distribution is uniform, e.g., n*U(p-0.2p,p+0.2p) (n is the number of mRNA moledules in the cell), then the variance is 1/12*(b-a)^2 = n**2 * p**2 * 4/300. so var/mean = np * 4/300, which is higher for stronger genes.

    if c_or_mc_mask is not None:
        c_or_mc_ad = mc.ut.slice(c_or_mc_ad, obs=c_or_mc_mask)

    if min_umi_count is not None:
        curr_c_ad_only_large = mc.ut.slice(c_or_mc_ad, obs=c_or_mc_ad.obs['num_of_non_excluded_umis'] >= min_umi_count)
        curr_c_ad_only_large.layers.pop('downsampled', None)
        mc.tl.downsample_cells(
            curr_c_ad_only_large,
            '__x__',
            downsample_min_samples=min_umi_count,
            random_seed=1234,
        )
    else:
        assert 'downsampled' in c_or_mc_ad.layers
        min_umi_count = c_or_mc_ad.uns['downsample_samples']
        curr_c_ad_only_large = mc.ut.slice(c_or_mc_ad, obs=c_or_mc_ad.obs['num_of_non_excluded_umis'] >= min_umi_count)

    ds_mat = mc.ut.to_numpy_matrix(curr_c_ad_only_large.layers['downsampled'])
    # ds_mat /= curr_c_ad_only_large.uns['downsample_samples'] # don't do this because the var is changed as the square of this...



    print(f'calculating var and mean across {len(ds_mat)} cells')
    var_mean_df = pd.DataFrame({
        'gene': curr_c_ad_only_large.var_names,
        'var': ds_mat.var(axis=0),
        'mean': ds_mat.mean(axis=0),
    })
    if color_by is not None:
        if hasattr(color_by, 'name'):
            hue_column_name = color_by.name
            var_mean_df[hue_column_name] = color_by.to_numpy()
        elif isinstance(color_by, str):
            assert color_by in c_or_mc_ad.var.columns
            hue_column_name = color_by
            var_mean_df[hue_column_name] = c_or_mc_ad.var[color_by].to_numpy()
    else:
        hue_column_name = None
    var_mean_df = var_mean_df[var_mean_df['mean'] > min_mean_expr * curr_c_ad_only_large.uns['downsample_samples']]
    var_mean_df['fixed_mean'] = var_mean_df['mean'] / curr_c_ad_only_large.uns['downsample_samples']
    var_mean_df['log_var'] = np.log2(var_mean_df['var'])
    var_mean_df['log_mean'] = np.log2(var_mean_df['mean'])
    var_mean_df['log_fixed_mean'] = np.log2(var_mean_df['fixed_mean'])

    var_mean_df['var_mean_log_ratio'] = var_mean_df['log_var'] - var_mean_df['log_mean']
    var_mean_df.sort_values('var_mean_log_ratio', ascending=False, inplace=True)

    # jitter_sd = 0.05
    # var_mean_df['jit_log_var'] = var_mean_df['log_var'] + np.random.normal(0, jitter_sd, len(var_mean_df))
    # var_mean_df['jit_log_mean'] = var_mean_df['log_mean'] + np.random.normal(0, jitter_sd, len(var_mean_df))

    if show_and_return_scatter_fig_and_ax:
        fig, ax = plt.subplots()
        sb.scatterplot(
            data=var_mean_df, 
            # x='jit_log_mean', y='jit_log_var', 
            x='log_fixed_mean', y='var_mean_log_ratio', 
            s=8,
            hue=hue_column_name,
            ax=ax,
        )
        if min_var_mean_log_ratio_to_show_gene_name is None:
            min_var_mean_log_ratio_to_show_gene_name = max(
                var_mean_df['var_mean_log_ratio'].sort_values(ascending=False).iloc[14],
                var_mean_df['var_mean_log_ratio'].quantile(0.995),
            )
        for _, row in var_mean_df[var_mean_df['var_mean_log_ratio'] > min_var_mean_log_ratio_to_show_gene_name].iterrows():
            ax.annotate(row['gene'], (row['log_fixed_mean'], row['var_mean_log_ratio']))
        # ax.axhline(11.25, color='red', alpha=0.3)
        # ax.axhline(0.19, color='red', alpha=0.3)

        return var_mean_df, fig, ax
    return var_mean_df

def get_quantiles_across_multiple_mc_ads(mc_ads_and_masks, genes, quantiles=(0.25, 0.5, 0.75), ignore_masks=False):
    all_expr = []
    for mc_ad, mc_mask in mc_ads_and_masks:
        if ignore_masks:
            mc_mask = None
        all_expr.extend(get_genes_expr(mc_ad, genes, c_or_mc_mask=mc_mask))
    return np.quantile(all_expr, quantiles)

def get_non_overlapping_expr_range_log_ratio_and_is_1_higher(mc_ad_and_mask1, mc_ad_and_mask2, genes, quantiles=(0.25, 0.75)):
    assert len(quantiles) == 2
    
    expr1 = get_genes_expr(mc_ad_and_mask1[0], genes, c_or_mc_mask=mc_ad_and_mask1[1])
    low_q1, high_q1 = np.quantile(expr1, quantiles)
    expr2 = get_genes_expr(mc_ad_and_mask2[0], genes, c_or_mc_mask=mc_ad_and_mask2[1])
    low_q2, high_q2 = np.quantile(expr2, quantiles)
    
    log_ratio1 = low_q1 - high_q2
    log_ratio2 = low_q2 - high_q1
    if log_ratio1 > log_ratio2:
        return log_ratio1, True
    return log_ratio2, False
    
def plot_comparison_scatters_between_models(mc_ads_and_masks, x_genes, y_genes, truncate_y_lims=True):
    flat_dicts = []
    for gene in x_genes:
        log_ratio, is_1_higher = get_non_overlapping_expr_range_log_ratio_and_is_1_higher(*mc_ads_and_masks, gene)
        flat_dicts.append({
            'gene': gene,
            'log_ratio': log_ratio,
            'is_1_higher': is_1_higher,
        })
    df = pd.DataFrame(flat_dicts).sort_values(
        [
            # 'is_1_higher', 
            'log_ratio',
        ], 
        ascending=False,
    )

    if truncate_y_lims:
        y_quantiles = get_quantiles_across_multiple_mc_ads(mc_ads_and_masks, y_genes, quantiles=(0, 0.25, 0.5, 0.75, 1))
        y_range = y_quantiles[-1] - y_quantiles[0]
        y_lims = (y_quantiles[0] - 0.1 * y_range, y_quantiles[-1] + 0.1 * y_range)
    else:
        y_lims = [*get_quantiles_across_multiple_mc_ads(mc_ads_and_masks, y_genes, quantiles=(0, 1), ignore_masks=True)]

    plt.close('all')
    fig, axes = generic_utils.get_multi_panel_fig_and_axes(len(df) * 2, ncols=2, figsize=(10,50))
    for ax_i, gene in enumerate(df['gene']):
        for ax, mc_ad_and_mask in zip(axes[ax_i], mc_ads_and_masks):
            
            plot_gene_gene_scatter(
                mc_ad_and_mask[0], gene, y_genes, ax=ax, 
                alpha=0.5,
            )
            plot_gene_gene_scatter(mc_ad_and_mask[0], gene, y_genes, c_or_mc_mask=mc_ad_and_mask[1], ax=ax)
            ax.set_ylim(y_lims)

            x_quantiles = get_quantiles_across_multiple_mc_ads(mc_ads_and_masks, gene)

            for q_val in x_quantiles:
                ax.axvline(q_val, color='grey', alpha=0.3)
        # raise
    fig.tight_layout()

    return fig, df

def add_cell_scores(
        c_ad, cell_score_desc_to_genes, overwrite_existing_cols=False, genes_to_ignore=None, layer_name=None, 
        add_as_umi_count=False, 
        # epsilon_for_log=1e-5,
        c_epsilon_to_add_to_fraction_before_log=C_EPSILON_TO_ADD_TO_FRACTION_BEFORE_LOG,
):
    ds_umi_count_of_cells = None
    cell_score_name_to_col_name = {}
    for cell_score_name, genes in cell_score_desc_to_genes.items():
        col_name = cell_score_name
        if layer_name == 'downsampled':
            col_name = f'ds_{col_name}'
        if genes_to_ignore and c_ad.var_names.isin(list(genes_to_ignore)).any():
            col_name = f'filtered_{col_name}'

        print(col_name)
        cell_score_name_to_col_name[cell_score_name] = col_name
        if col_name in c_ad.obs.columns:
            if overwrite_existing_cols:
                print(f'overwriting {col_name}')
            else:
                print(f'skipping {col_name}, as it already exists in c_ad.obs')
                continue
        
        if genes_to_ignore:
            genes = list(set(genes) - set(genes_to_ignore))
        if add_as_umi_count:
            assert layer_name == 'downsampled'
            if ds_umi_count_of_cells is None:
                ds_umi_count_of_cells = c_ad.layers['downsampled'].sum(axis=1)
                assert (ds_umi_count_of_cells <= c_ad.uns['downsample_samples']).all()
                c_mask = mc.ut.to_numpy_vector(ds_umi_count_of_cells == c_ad.uns['downsample_samples'])
            c_ad.obs.loc[c_mask, col_name] = get_genes_expr(c_ad, genes, layer_name=layer_name, return_umi_count=True, c_or_mc_mask=c_mask)
        else:
            c_ad.obs[col_name] = get_genes_expr(c_ad, genes, layer_name=layer_name, c_epsilon_to_add_to_fraction_before_log=c_epsilon_to_add_to_fraction_before_log)
    return cell_score_name_to_col_name

            

def merge_c_ads_after_running_divide_and_conquer_with_same_args(c_ads_to_merge):
    all_obs_names = pd.concat([pd.Series(x.obs_names) for x in c_ads_to_merge])
    assert all_obs_names.is_unique

    curr_highest_mc_i = -1
    var_names = c_ads_to_merge[0].var_names
    rare_gene_col = np.full(len(var_names), False)
    selected_gene_col = np.full(len(var_names), False)
    mc_i_series_list = []
    for curr_c_ad in c_ads_to_merge:
        assert (curr_c_ad.var_names == var_names).all()
        assert is_metacell_col_as_expected(curr_c_ad)
        
        mc_i_series = curr_c_ad.obs['metacell'].copy()
        update_mask = mc_i_series >= 0
        mc_i_series.loc[update_mask] += curr_highest_mc_i + 1
        mc_i_series_list.append(mc_i_series)
        curr_highest_mc_i = mc_i_series.max()

        rare_gene_col |= curr_c_ad.var['rare_gene']
        selected_gene_col |= curr_c_ad.var['selected_gene']

    merged_c_ad = ad.concat(
        c_ads_to_merge, 
        merge='same',
    )
    print(f'merged_c_ad.n_obs: {merged_c_ad.n_obs}')
    merged_c_ad.var['gene'] = merged_c_ad.var_names
    merged_c_ad.var.drop(columns='rare_gene_module', inplace=True, errors='ignore')

    assert (merged_c_ad.var_names == var_names).all()
    merged_c_ad.var['rare_gene'] = rare_gene_col
    merged_c_ad.var['selected_gene'] = selected_gene_col

    # for some reason, ad.concat drops some of the columns, so we add them back.
    concat_obs = pd.concat([x.obs for x in c_ads_to_merge])
    assert (concat_obs.index == merged_c_ad.obs.index).all()
    merged_c_ad.obs = concat_obs

    concat_mc_i_series = pd.concat(mc_i_series_list)
    merged_c_ad.obs.loc[concat_mc_i_series.index, 'metacell'] = concat_mc_i_series
    assert is_metacell_col_as_expected(merged_c_ad)

    return merged_c_ad

def get_agg_expr_per_gene(c_ad, agg_func, c_mask=None, epsilon_for_log=1e-5):
    if c_mask is None:
        c_mask = np.full(c_ad.n_obs, True)
    
    filtered_x = c_ad.X[c_mask,:]
    return np.log2(mc.ut.to_numpy_vector(agg_func(filtered_x / filtered_x.sum(axis=1), axis=0)) + epsilon_for_log)

def get_max_expr_per_gene(c_ad, c_mask=None, epsilon_for_log=1e-5):
    return get_agg_expr_per_gene(c_ad, np.max, c_mask=c_mask, epsilon_for_log=epsilon_for_log)

def get_mean_expr_per_gene_but_arith_for_performance(c_ad, c_mask=None, epsilon_for_log=1e-5):
    return get_agg_expr_per_gene(c_ad, np.mean, c_mask=c_mask, epsilon_for_log=epsilon_for_log)

def get_c_state_state_matching_df(c_ad):
    c_states = set(c_ad.obs['c_state'].unique())
    dfs = []
    for c_state in c_states:
        print(c_state)
        enrich_df = get_c_mask_freq_per_obs_column(
            c_ad, 
            c_ad.obs['c_state'] == c_state,
            ['state'],
        ).head(2)
        enrich_df['c_state'] = c_state
        enrich_df['state2'] = enrich_df['state'].iloc[1]
        enrich_df['freq2'] = enrich_df['freq'].iloc[1]
        enrich_df = enrich_df.head(1)
        dfs.append(enrich_df)
    df = pd.concat(dfs, ignore_index=True)
    df.sort_values('freq', inplace=True)
    df['state_equal_c_state'] = df['state'] == df['c_state']
    # df.sort_values('freq2', inplace=True)
    return df

def get_obs_col_val_to_mask(data_ad, obs_col):
    return {x: data_ad.obs[obs_col] == x for x in sorted(data_ad.obs[obs_col].unique())}

def print_most_enriched_per_mask(mask_names, obs_col_name, c_ad, min_freq=0.1):
    final_reprs = []
    for mask_name in mask_names:
        enrich_df = get_c_mask_freq_per_obs_column(
            c_ad, 
            c_ad.obs[mask_name],
            obs_col_name,
        )
        enrich_df = enrich_df[enrich_df['freq'] >= min_freq]
        freq_reprs = []
        for _, row in enrich_df.iterrows():
            freq_reprs.append(f'{row[obs_col_name]}:{row["freq"]:.2f}')
        freq_repr = ', '.join(freq_reprs)
        
        final_reprs.append(f'{mask_name} -- {freq_repr}')
    print('\n'.join(final_reprs))

def plot_x_across_traj_per_donor(
        c_ad, 
        numbered_donor_df,
        diagnosis_class_to_color, exp_names,
        traj_pos_obs_col_name, traj_pos_min_truncate_val, traj_pos_max_truncate_val, 
        bin_count, 
        numbered_donor_id_to_exp_date_color,
        analysis_name=None,
        # row_cluster=True, 
        ignore_lower_than_truncate_val=False, 
        min_donor_cell_count=10,
        min_bin_cell_count=1, # irrelevant if x_is_cell_freq=True
        x_obs_col_name=None,
        x_obs_col_agg_func=np.median,
        x_is_cell_freq=False,
        x_is_cell_cum_freq=False,
        x_is_max_selected_gene_corr_with_ref_mc=False,
        x_is_cell_freq_or_cum_freq_out_of_background=False,
        c_mask=None, 
        c_mask_col=None, 
        c_mask_col_and_vals=None,
        background_c_mask_for_freq=None,
        background_c_mask_for_freq_col=None, background_c_mask_for_freq_c_states=None,
        show_only_best_tech_rep=False,
        bin_edge_to_tick_label_func=lambda x: generic_utils.remove_redundant_trailing_zeros(f'{x:.1f}'),
        clustermap_kwargs=None,
        norm_each_traj_bin=False,
        norm_each_donor=False,
        df_for_clustermap_index=None,
        ref_mc_ad=None,
):
    
    if x_is_cell_freq or x_is_cell_cum_freq:
        if min_bin_cell_count > 0:
            print('plot_x_across_traj_per_donor: min_bin_cell_count is ignored as (x_is_cell_freq or x_is_cell_cum_freq)==True')
        min_bin_cell_count = None
    
    numbered_donor_id_to_numbered_donor_id_without_tech_rep_suffix = generic_utils.get_dict_mapping_one_df_column_to_other(
        numbered_donor_df, 'numbered_donor_id', 'numbered_donor_id_without_tech_rep_suffix')
    numbered_donor_id_to_exp_name = generic_utils.get_dict_mapping_one_df_column_to_other(numbered_donor_df, 'numbered_donor_id', 'exp_name')
    numbered_donor_id_to_diagnosis_class = generic_utils.get_dict_mapping_one_df_column_to_other(numbered_donor_df, 'numbered_donor_id', 'diagnosis_class')
    numbered_donor_id_to_diagnosis_class_color = {
        k: diagnosis_class_to_color[v] for k,v in numbered_donor_id_to_diagnosis_class.items()}
    exp_name_to_color = {
        k: v
        for k,v in zip(exp_names, generic_utils.get_n_colors(len(exp_names)))
    }
    numbered_donor_id_to_exp_color = {k: exp_name_to_color[v] for k,v in numbered_donor_id_to_exp_name.items()}
    all_numbered_donor_ids = sorted(numbered_donor_df['numbered_donor_id'])

    if c_mask is None:
        c_mask = np.full(c_ad.n_obs, True)
    if c_mask_col:
        c_mask &= c_ad.obs[c_mask_col]
    if c_mask_col_and_vals:
        c_mask &= c_ad.obs[c_mask_col_and_vals[0]].isin(c_mask_col_and_vals[1])

    if x_is_cell_freq_or_cum_freq_out_of_background:
        if background_c_mask_for_freq is None:
            background_c_mask_for_freq = np.full(c_ad.n_obs, True)
        if background_c_mask_for_freq_col:
            background_c_mask_for_freq &= c_ad.obs[background_c_mask_for_freq_col]
        if background_c_mask_for_freq_c_states:
            background_c_mask_for_freq &= c_ad.obs['c_state'].isin(background_c_mask_for_freq_c_states)
        assert not (c_mask & (~background_c_mask_for_freq)).any(), 'c_mask must be a subset of background_c_mask_for_freq'
        background_cell_count_df = c_ad.obs.loc[background_c_mask_for_freq, 'numbered_donor_id'].value_counts().reset_index(name='background_count')

    if show_only_best_tech_rep:
        best_tech_rep_numbered_donor_ids = numbered_donor_df.sort_values('cell_count', ascending=False).drop_duplicates(
            subset=['donor_id', 'bleeding_date', 'is_bm', 'diagnosis'])['numbered_donor_id']

        best_tech_rep_mask = c_ad.obs['numbered_donor_id'].isin(best_tech_rep_numbered_donor_ids)
        c_mask &= best_tech_rep_mask
    if ignore_lower_than_truncate_val:
        c_mask &= (c_ad.obs[traj_pos_obs_col_name] >= traj_pos_min_truncate_val)

    cell_count_df = c_ad.obs.loc[c_mask, 'numbered_donor_id'].value_counts().reset_index(name='count')
    enough_cells_numbered_donor_ids = list(cell_count_df.loc[cell_count_df['count'] >= min_donor_cell_count, 'numbered_donor_id'])
    too_few_cells_numbered_donor_ids = sorted(set(all_numbered_donor_ids) - set(enough_cells_numbered_donor_ids))
    print(f'discarded due to too low total cell count: {too_few_cells_numbered_donor_ids}')
    c_mask &= c_ad.obs['numbered_donor_id'].isin(enough_cells_numbered_donor_ids)

    numbered_donor_id_of_cells = c_ad.obs.loc[c_mask, 'numbered_donor_id'].astype(str)
    traj_pos_of_cells = c_ad.obs.loc[c_mask, traj_pos_obs_col_name]
    if traj_pos_min_truncate_val is not None:
        traj_pos_of_cells = np.maximum(traj_pos_of_cells, traj_pos_min_truncate_val)
    if traj_pos_max_truncate_val is not None:
        traj_pos_of_cells = np.minimum(traj_pos_of_cells, traj_pos_max_truncate_val)
    
    bins = generic_utils.get_bins(bin_count + 1, traj_pos_of_cells)
    traj_bin_i_of_cells = generic_utils.safe_digitize(traj_pos_of_cells, bins=bins)
    df = pd.DataFrame({
        'numbered_donor_id': numbered_donor_id_of_cells,
        'traj_bin_i': traj_bin_i_of_cells,
        'traj_pos': traj_pos_of_cells,
    })
    
    cell_count_pivot_table = df[['numbered_donor_id', 'traj_bin_i']].value_counts().reset_index().pivot(
        index='numbered_donor_id', columns='traj_bin_i', values='count').fillna(0)

    assert (x_obs_col_name is not None) + (x_is_cell_freq != False) + (x_is_cell_cum_freq != False) + (x_is_cell_freq_or_cum_freq_out_of_background != False) == 1, 'more than one x was specified (via x_obs_col_name, x_is_cell_freq, x_is_cell_cum_freq, x_is_cell_freq_or_cum_freq_out_of_background)'
    if x_is_cell_freq or x_is_cell_cum_freq:
        df_for_clustermap = cell_count_pivot_table.divide(cell_count_pivot_table.sum(axis=1), axis=0)
        if x_is_cell_cum_freq:
            assert pd.Series(df.columns).is_monotonic_increasing
            df_for_clustermap = df_for_clustermap.cumsum(axis=1)
    elif x_is_cell_freq_or_cum_freq_out_of_background:
        ordered_filtered_background_cell_count_df = generic_utils.merge_preserving_df1_index_and_row_order(
            pd.DataFrame({'numbered_donor_id': cell_count_pivot_table.index}),
            background_cell_count_df,
        )
        df_for_clustermap = cell_count_pivot_table.divide(ordered_filtered_background_cell_count_df['background_count'].to_numpy(), axis=0)
    elif x_obs_col_name:
        df[x_obs_col_name] = c_ad.obs.loc[c_mask, x_obs_col_name]
        df_for_clustermap = df.groupby(['numbered_donor_id', 'traj_bin_i'])[x_obs_col_name].agg(x_obs_col_agg_func).reset_index().pivot(
            index='numbered_donor_id', columns='traj_bin_i', values=x_obs_col_name)
        
        assert min_bin_cell_count >= 1
        any_bin_with_too_few_cells_mask = (cell_count_pivot_table < min_bin_cell_count).any(axis=1)
        any_bin_with_too_few_cells_numbered_donor_ids = sorted(cell_count_pivot_table.index[any_bin_with_too_few_cells_mask])
        print(f'discarded due to any bin with too few cells: {any_bin_with_too_few_cells_numbered_donor_ids}')
        df_for_clustermap = df_for_clustermap[~any_bin_with_too_few_cells_mask]
    elif x_is_max_selected_gene_corr_with_ref_mc:
        raise NotImplementedError('would be pretty slow, i guess. anyway, i think: generate mean_expr for each bin and donor, then use corrcoef to get corr with each atlas_mc_ad mc')

    if norm_each_traj_bin:
        assert not norm_each_donor
        median_of_cols = df_for_clustermap.median(axis=0)
        # print(f'median_of_cols: {median_of_cols}')
        df_for_clustermap = df_for_clustermap.subtract(median_of_cols, axis=1)
    if norm_each_donor:
        assert not norm_each_traj_bin
        median_of_rows = df_for_clustermap.median(axis=1)
        # print(f'median_of_rows: {median_of_rows}')
        df_for_clustermap = df_for_clustermap.subtract(median_of_rows, axis=0)

    if clustermap_kwargs is None:
        clustermap_kwargs = dict()

    if 'cmap' not in clustermap_kwargs:
        clustermap_kwargs['cmap'] = 'bwr' if (norm_each_traj_bin or x_is_cell_cum_freq or norm_each_donor) else 'rocket_r'
    
    if ('vmin' not in clustermap_kwargs) and ('vmax' not in clustermap_kwargs):
        if norm_each_traj_bin or norm_each_donor:
            abs_val_90_percentile = np.percentile(np.abs(df_for_clustermap.to_numpy().flatten()), 90)
            clustermap_kwargs['vmin'] = -abs_val_90_percentile
            clustermap_kwargs['vmax'] = abs_val_90_percentile
        if x_is_cell_cum_freq:
            clustermap_kwargs['vmin'] = 0
            clustermap_kwargs['vmax'] = 1
    
    if 'method' not in clustermap_kwargs:
        clustermap_kwargs['method'] = 'ward'
        clustermap_kwargs['metric'] = 'euclidean' # scipy says "Method 'ward' requires the distance metric to be Euclidean", so i listen to scipy...
    if 'figsize' not in clustermap_kwargs:
        clustermap_kwargs['figsize'] = (10, 2.8 + 0.15 * len(df_for_clustermap))
    if 'row_cluster' in clustermap_kwargs:
        assert not (clustermap_kwargs['row_cluster'] and x_is_cell_cum_freq), '_Amos_ said that row_cluster=True is a bad idea for x_is_cell_cum_freq=True. instead we should just sort by median etc'
    else:
        if df_for_clustermap_index:
            clustermap_kwargs['row_cluster'] = False
            df_for_clustermap = df_for_clustermap.loc[df_for_clustermap_index]
        elif x_is_cell_cum_freq:
            clustermap_kwargs['row_cluster'] = False
            ordered_numbered_donor_ids = list(df.groupby('numbered_donor_id')['traj_pos'].median().sort_values().index)
            assert len(ordered_numbered_donor_ids) == len(df_for_clustermap)
            assert set(ordered_numbered_donor_ids) == set(df_for_clustermap.index)
            df_for_clustermap = df_for_clustermap.loc[ordered_numbered_donor_ids]
        else:
            clustermap_kwargs['row_cluster'] = True

    clustermap_obj = sb.clustermap(
        df_for_clustermap,
        col_cluster=False,
        row_colors=[
            [numbered_donor_id_to_exp_color[x] for x in df_for_clustermap.index],
            # [['black', 'white'][numbered_donor_id_to_is_ultima[x]] for x in df_for_clustermap.index],
            [numbered_donor_id_to_exp_date_color[x] for x in df_for_clustermap.index],
            [numbered_donor_id_to_diagnosis_class_color[x] for x in df_for_clustermap.index],
        ],
        # cbar_kws=None,
        yticklabels=[numbered_donor_id_to_numbered_donor_id_without_tech_rep_suffix[x] for x in df_for_clustermap.index],
        # yticklabels=df_for_clustermap.index,
        xticklabels=True,
        # cbar_pos=None, # if we do this, then clustermap_obj.ax_cbar gives another ax.
        **clustermap_kwargs,
    )
    heatmap_ax = clustermap_obj.ax_heatmap
    heatmap_ax.set_xlabel(traj_pos_obs_col_name)
    heatmap_ax.set_ylabel(None)
    title = ''
    if analysis_name:
        title += analysis_name
    if norm_each_traj_bin:
        title += '\n(norm_each_bin)'
    if norm_each_donor:
        title += '\n(norm_each_donor)'

    if title:
        heatmap_ax.set_title(title)

    generic_utils.set_heatmap_bin_edge_xticks(
        heatmap_ax, bins,
        bin_edge_to_tick_label_func=bin_edge_to_tick_label_func,
        min_truncate_val=traj_pos_min_truncate_val,
        max_truncate_val=traj_pos_max_truncate_val,
        dilute_x_labels=x_is_cell_cum_freq,
    )
    
    return clustermap_obj

def get_pooled_expr_mat(c_ad, c_mask, obs_col_name, min_cell_count=3):
    obc_vals = list((c_ad.obs.loc[c_mask, obs_col_name].value_counts() >= min_cell_count).loc[lambda x: x].index)
    val_to_expr_vec = {}
    for val in obc_vals:
        curr_c_mask = c_mask & (c_ad.obs[obs_col_name] == val)
        val_to_expr_vec[val] = get_mean_expr(c_ad, curr_c_mask)
    expr_mat = pd.DataFrame(val_to_expr_vec, index=c_ad.var_names)
    return expr_mat

def get_norm_pooled_diff_expr_mat(c_ad, c_mask1, c_mask2, obs_col_name, min_cell_count=3, min_max_abs_norm_diff=None, expr_mats=None):
    if expr_mats is None:
        expr_mat1 = get_pooled_expr_mat(c_ad, c_mask1, obs_col_name=obs_col_name, min_cell_count=min_cell_count)
        expr_mat2 = get_pooled_expr_mat(c_ad, c_mask2, obs_col_name=obs_col_name, min_cell_count=min_cell_count)
        expr_mats = (expr_mat1, expr_mat2)
    else:
        expr_mat1, expr_mat2 = expr_mats
    diff_expr_mat = expr_mat1 - expr_mat2
    diff_expr_mat = diff_expr_mat.loc[:, ~diff_expr_mat.isna().any(axis=0)]
    norm_pooled_diff_expr_mat = diff_expr_mat.subtract(diff_expr_mat.median(axis=1), axis=0)
    abs_norm_pooled_diff_expr_mat = norm_pooled_diff_expr_mat.abs()
    max_abs_norm_diff_of_genes = abs_norm_pooled_diff_expr_mat.max(axis=1)
    
    if min_max_abs_norm_diff is None:
        plt.close('all')
        sb.histplot(max_abs_norm_diff_of_genes, bins=50)
        print('specify min_max_abs_norm_diff according to the hist.')
        return None, expr_mats
    norm_pooled_diff_expr_mat = norm_pooled_diff_expr_mat[max_abs_norm_diff_of_genes >= min_max_abs_norm_diff]
    print(f'len(norm_pooled_diff_expr_mat): {len(norm_pooled_diff_expr_mat)}')

    return norm_pooled_diff_expr_mat, expr_mats

def add_c_state_and_mc_c_state_stats(
        mc_ad, c_ad, mask_and_info_list, cell_state_and_info_list, cell_type_colors_csv_file_path, 
        silently_overwrite_state_and_state_color_columns=True, verbose=False,
        add_all_c_state_fracs=True, only_add_c_state=False,
):
    add_mask_cols(
        c_ad, 
        mask_and_info_list=mask_and_info_list, 
        # dist_with_threshs_out_dir_path='temp/aa',
        verbose=verbose,
    )
    add_state(
        c_ad, 
        cell_state_and_info_list=cell_state_and_info_list, 
        cell_type_colors_csv_file_path=cell_type_colors_csv_file_path,
        silently_overwrite_state_and_state_color_columns=silently_overwrite_state_and_state_color_columns,
        state_column_name='c_state',
        verbose=verbose,
    )

    if only_add_c_state:
        return

    add_c_mask_fraction_of_mc_column(mc_ad, c_ad, c_ad.obs['c_state'] == 'state_unassigned', 'unassigned_c_state_freq')
    add_mc_metadata(mc_ad, c_ad, {'str_mode_mode2_and_fracs_and_nunique': ['c_state']}, allow_overwriting=True)
    if add_all_c_state_fracs:
        add_mc_metadata(mc_ad, c_ad, {'category_enrichment_and_mode_and_more': ['c_state']}, allow_overwriting=True)


    if verbose:
        print(mc_ad.obs[['c_state_mode', 'c_state_mode2']].value_counts())

def get_potential_artifact_mixed_mcs(mc_ad, states=['state_unassigned', 'Doublet']):
    return mc_ad.obs.loc[mc_ad.obs['state'].isin(states), ['c_state_mode', 'c_state_mode2', 'c_state_mode_frac', 'c_state_mode2_frac']].sort_values(
        ['c_state_mode', 'c_state_mode2'])

