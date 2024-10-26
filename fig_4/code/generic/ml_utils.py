import pathlib
import random
import sys
import warnings
import io
import unittest
import hashlib
import math
import os.path
import os
import itertools
import collections
import pandas as pd
import numpy as np
import subprocess
import contextlib
import dis
import string
import types
import shutil
import time
import pickle
import logging
import textwrap
import scipy.stats
import datetime
import matplotlib
import matplotlib.pyplot as plt

import sklearn
import sklearn.pipeline
import sklearn.feature_selection
import seaborn as sb
import metacells as mc
import xgboost as xgb
import shap

from generic import generic_utils

try: import __builtin__ as builtins
except ImportError: import builtins

# 200624: It seems that this is the default (what I get from np.geterr() after importing numpy): {'divide': 'warn', 'over': 'warn', 'under': 'ignore', 'invalid': 'warn'}
# https://numpy.org/doc/stable/reference/generated/numpy.seterr.html says: Underflow: result so close to zero that some precision was lost.
# so I guess it is Ok that I just ignore underflow problems?
# np.seterr(all='raise')
np.seterr(divide='raise', over='raise', invalid='raise')


def get_tpr(ground_truth_vec, prediction_vec):
    if (ground_truth_vec == prediction_vec).all() or (ground_truth_vec == 0).all():
        return 1
    confusion_mat = sklearn.metrics.confusion_matrix(ground_truth_vec, prediction_vec)
    raveled_confusion_mat = confusion_mat.ravel()
    # if len(raveled_confusion_mat) == 1:
    #     print(confusion_mat, raveled_confusion_mat, ground_truth_vec, prediction_vec)
    #     print(ground_truth_vec)
    #     print(prediction_vec)
    tn, fp, fn, tp = raveled_confusion_mat
    tpr = tp / (tp + fn)
    return tpr

def lasso_cv(feature_df, target_vec, cv, return_lasso_cv_obj_on_all=True, alphas=np.logspace(-5, 1, 20), **lasso_cv_kwargs):
    # modified from chatGPT answer
    X = feature_df.values
    y = target_vec
    
    flat_dicts = []
    for alpha in alphas:
        print(f'log(alpha): {np.log2(alpha)}')
        # Train LASSO model with the specific alpha
        
        # # Predict on the test set
        # y_pred = curr_alpha_lasso.predict(X)
        # y_pred_binary = np.where(y_pred > 0.5, 1, 0)
        
        kfold = sklearn.model_selection.KFold(n_splits=cv, shuffle=True, random_state=0)

        # Calculate TPR for each fold
        tpr_values = []
        nonzero_coef_counts = []
        for train_idx, test_idx in kfold.split(X):
            X_cv_train, X_cv_test = X[train_idx], X[test_idx]
            y_cv_train, y_cv_test = y[train_idx], y[test_idx]
            
            curr_alpha_lasso = sklearn.linear_model.Lasso(alpha=alpha, **lasso_cv_kwargs)
            curr_alpha_lasso.fit(X_cv_train, y_cv_train)
            y_cv_pred = curr_alpha_lasso.predict(X_cv_test)
            y_cv_pred_binary = np.where(y_cv_pred > 0.5, 1, 0)
            
            tpr = get_tpr(y_cv_test, y_cv_pred_binary)
            tpr_values.append(tpr)
            nonzero_coef_counts.append(np.sum(curr_alpha_lasso.coef_ != 0))

        
        # Store the mean and std of TPR values
        flat_dicts.append({
            'alpha': alpha,
            'tpr_mean': np.mean(tpr_values),
            'tpr_std': np.std(tpr_values),
            'nonzero_coef_count_mean': np.mean(nonzero_coef_counts).astype(int),
        })
    
    df = pd.DataFrame(flat_dicts)
    df['log_alpha'] = np.log2(df['alpha'])
    # Plot the results
    fig, ax = plt.subplots()
    # plt.figure(figsize=(10, 6))
    ax.errorbar(df['log_alpha'], df['tpr_mean'], yerr=df['tpr_std'], fmt='-o', ecolor='r', capsize=5)
    ax.set_xlabel('log($\lambda$)')
    ax.set_ylabel('TPR')
    for _, row in df.iterrows():
        ax.text(
            row['log_alpha'], 0.96, 
            str(int(row['nonzero_coef_count_mean'])),
            ha='center', 
            transform=matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes),
            fontsize='small',
        )
    # ax.grid(True)
    
    if return_lasso_cv_obj_on_all:
        lasso_cv = sklearn.linear_model.LassoCV(alphas=alphas, **lasso_cv_kwargs)
        lasso_cv.fit(X, y)
        pred_y = lasso_cv.predict(X)
        return lasso_cv, pred_y, fig, ax
    
    return fig, ax
        


def amos_feature_df_and_prediction_barplot(
        feature_df, prediction_df, title_prefix_for_plot='', name_to_hue=None, palette=None, show_prediction_plus1=False, 
        subplots_adjust_top=0.92,
):

    feature_df = feature_df.copy()
    for col in feature_df.columns:
        if not (feature_df[col] > 0).any(): # do it this way because of nans
            feature_df[col] = -feature_df[col]
            feature_df.rename(columns={col: f'-{col}'}, inplace=True)
    feature_cols = list(feature_df.columns)

    df = generic_utils.merge_preserving_df1_index_and_row_order(prediction_df[['name', 'pred_y']], feature_df.reset_index().rename(columns={'index': 'name'}))


    if show_prediction_plus1:
        df['pred_y+1'] = df['pred_y'] + 1
        pred_y_col = 'pred_y+1'
    else:
        pred_y_col = 'pred_y'
    if name_to_hue:
        df['hue'] = df['name'].map(name_to_hue)
        assert df['hue'].notna().all()
    df.sort_values('pred_y', ascending=False, inplace=True)
    
    barplot_kwargs = dict(
        data=df, 
        y='name', 
        dodge=False,
    )
    if 'hue' in df:
        barplot_kwargs['hue'] = 'hue'
        barplot_kwargs['palette'] = palette
    
    n_axes = len(feature_cols) + 2
    fig, axes = plt.subplots(
        # ncols=n_axes, 
        nrows=1, ncols=n_axes,
        gridspec_kw=dict(
            wspace=0.07,
            width_ratios=[1] * (n_axes - 2) + [1.2, 2],
        ),
        figsize=(1.1 * n_axes, len(df) * 0.15),
    )

    fig.subplots_adjust(
        top=subplots_adjust_top,
        left=0.03,
        bottom=0.03,
        right=0.97,
    )
    for ax in axes[:-1]:
        generic_utils.make_all_spines_and_x_and_y_axes_invisible(ax)
    
    for ax, feature_col in zip(axes[:-2], feature_cols):
        sb.barplot(
            x=feature_col, 
            **barplot_kwargs, ax=ax,
        )
        quantiles = df[feature_col].quantile([0.25, 0.5, 0.75])
        ax.axvline(quantiles[0.5], color='black')
        for q in (0.25, 0.75):
            ax.axvline(quantiles[q], color='black', linestyle='--')
        ax.get_legend().remove()
        # ax.set_title(feature_col, fontsize='x-small')
        ax.text(
            0.5, 1, feature_col,
            size="x-small", rotation=15,
            horizontalalignment='center', verticalalignment='bottom',
            transform=ax.transAxes,
            # rotation_mode='anchor',
        )
        generic_utils.set_ax_lim_by_vals(ax, 'x', vals=df[feature_col])
        # break

    ax = axes[-1]
    sb.barplot(
        x=pred_y_col, 
        **barplot_kwargs, ax=ax,
    )
    ax.set_title('prediction', fontsize='x-small')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize='small')
    ax.get_legend().remove()
    generic_utils.set_ax_lim_by_vals(ax, 'x', vals=df[pred_y_col])
    generic_utils.make_all_spines_invisible(ax)
    ax.get_xaxis().set_visible(False)
    ax.set_ylabel(None)

    title = f'{title_prefix_for_plot}'
    fig.suptitle(title)
    
    fig.savefig(f'temp/features_and_preds/{generic_utils.strip_special_chars_on_edges_and_replace_others_with_underscores(title_prefix_for_plot)}_features_and_preds_barplot.png')

def regression_prediction_scatter(
        pred_df, ax_title='', name_to_hue=None, palette=None,
        show_names=False,
):
    if name_to_hue:
        pred_df['hue'] = pred_df['name'].map(name_to_hue)
        assert pred_df['hue'].notna().all()
    
    fig, ax = plt.subplots(figsize=(5, 5))
    scatterplot_kwargs = dict(
        data=pred_df,
        x='pred_y',
        y='y',
        ax=ax
    )
    if 'hue' in pred_df:
        scatterplot_kwargs['hue'] = 'hue'
        scatterplot_kwargs['palette'] = palette
    sb.scatterplot(**scatterplot_kwargs)
    if show_names:
        for i, row in pred_df.iterrows():
            ax.text(row['pred_y'], row['y'], row['name'], fontsize='x-small')
    # sb.scatterplot(**scatterplot_kwargs)
    ax.get_legend().remove()
    ax.set_title(ax_title)
    fig.tight_layout()
    fig.savefig(f'temp/regress_pred/{generic_utils.strip_special_chars_on_edges_and_replace_others_with_underscores(ax_title)}_preds_scatter.png')

def binary_classification_prediction_barplot(
        pred_df, title_prefix_for_plots='', name_to_hue=None, palette=None, show_prediction_plus1=False, 
        threshs_to_also_plot_only_low_and_only_high_preds=(0.05, 0.95),
):
    pred_df = pred_df.copy()

    pred_df['pred_y+1'] = pred_df['pred_y'] + 1
    if name_to_hue:
        pred_df['hue'] = pred_df['name'].map(name_to_hue)
        assert pred_df['hue'].notna().all()
    pred_df.sort_values('pred_y', inplace=True)

    
    pred_df['name_with_#_if_y'] = pred_df['name']
    pred_df.loc[pred_df['y'] == 1, 'name_with_#_if_y'] += '#'


    descs_and_pred_dfs_for_plots = [
        ('all', pred_df),
    ]
    below_low_thresh_mask = pred_df['pred_y'] < threshs_to_also_plot_only_low_and_only_high_preds[0]
    above_high_thresh_mask = pred_df['pred_y'] > threshs_to_also_plot_only_low_and_only_high_preds[1]
    if (below_low_thresh_mask | above_high_thresh_mask).all():
        descs_and_pred_dfs_for_plots += [
            (f'only preds below {threshs_to_also_plot_only_low_and_only_high_preds[0]}', pred_df[below_low_thresh_mask]), 
            (f'only preds above {threshs_to_also_plot_only_low_and_only_high_preds[1]}', pred_df[above_high_thresh_mask]), 
        ]

    for desc, pred_df_for_plot in descs_and_pred_dfs_for_plots:
        pred_df_for_plot = pred_df_for_plot.copy()
        # pred_df_for_plot['']
        fig, ax = plt.subplots(figsize=(13, 5))
        
        y_col = 'pred_y+1' if show_prediction_plus1 else 'pred_y'
        barplot_kwargs = dict(
            data=pred_df_for_plot, 
            # x='name', 
            x='name_with_#_if_y', 
            y=y_col,
            ax=ax, 
            dodge=False,
        )
        if 'hue' in pred_df_for_plot:
            barplot_kwargs['hue'] = 'hue'
            barplot_kwargs['palette'] = palette
        sb.barplot(**barplot_kwargs)
        # sb.scatterplot(**barplot_kwargs)
        ax.get_legend().remove()
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=90, 
            # ha='right', 
            # rotation_mode='anchor', 
            # fontsize='small',
        )
        generic_utils.set_ax_lim_by_vals(ax, 'y', vals=pred_df_for_plot[y_col])
        title = f'{title_prefix_for_plots}{desc}'
        ax.set_title(title)
        ax.set_xlabel(None)
        fig.tight_layout()
        pathlib.Path('temp/binary_classif_pred').mkdir(parents=True, exist_ok=True)
        fig.savefig(f'temp/binary_classif_pred/{generic_utils.strip_special_chars_on_edges_and_replace_others_with_underscores(title_prefix_for_plots)}{desc}_preds_barplot.png')

def plot_auc(pred_df_or_name_to_pred_df, title_suffix=''):
    plt.close('all')
    fig, ax = plt.subplots(figsize=(4, 4))
    auc_reprs = []

    name_to_pred_df = pred_df_or_name_to_pred_df if isinstance(pred_df_or_name_to_pred_df, dict) else {'': pred_df_or_name_to_pred_df}
    
    for name, pred_df in name_to_pred_df.items():
        sample_count = len(pred_df)
        # print(name, pred_df['pred_y'].value_counts())
        fpr, tpr, _ = sklearn.metrics.roc_curve(
            pred_df['y'], pred_df['pred_y'], 
            # drop_intermediate=False,
        )
        # print(tpr, fpr)
        # print(name, pd.DataFrame({'fpr': fpr, 'tpr': tpr}))
        auc = sklearn.metrics.roc_auc_score(pred_df['y'], pred_df['pred_y'])
        print(f'{auc:.3f}')
        if 0:
            roc_df = pd.DataFrame({'fpr': fpr, 'tpr': tpr})
            sb.lineplot(roc_df, x='fpr', y='tpr', ax=ax) # this looks weird for some reason...
        ax.plot(fpr, tpr, label=name)
        auc_reprs.append(f'{name}(n={sample_count}) auc={auc:.3f}')

    ax.set_xlabel('FPR')
    ax.set_ylabel('TPR')
    ax.legend()
    ax.set_aspect('equal')
    title_suffix_repr = f'\n{title_suffix}' if title_suffix else ''
    ax.set_title(f'{", ".join(auc_reprs)}{title_suffix_repr}', fontsize='xx-small')
    fig.tight_layout()
    return fig, ax

def get_pipeline_transformed_x(pipe, X):
    transformed_X = X
    for step_name, transformer in pipe.steps[:-1]:
        transformed_X = transformer.transform(transformed_X)
    return transformed_X
    

def xgboost_cv(
        is_classif,
        train_feature_df, train_target_vec, 
        test_feature_df=None,
        test_target_vec=None,
        prediction_barplot_kwargs=dict(),
        prediction_scatterplot_kwargs=dict(),
        # default xgboost params: https://xgboost.readthedocs.io/en/latest/parameter.html
        # arg_name_and_vals_and_plot_log=('max_depth', range(2, 9), False),
        arg_name_and_vals_and_plot_log=('dummy', [0], False),
        cv=None,
        leave_one_out=True,
        title_prefix_for_plots='',
        get_pipe=None,
        modify_pipe_inside_cv_func=None,
        score_train_subgroup_name_to_mask=dict(),
        score_train_subgroup_name_for_best_cv='all',
        skip_cv=False,
        shap_expected_approx_zero_vals_allclose_atol=1e-5,
        out_dir_path='temp/xgboost_cv',
        **xgboost_kwargs,
):
    
    # which hyperparameters to tune?
    # https://xgboost.readthedocs.io/en/stable/tutorials/param_tuning.html
    # https://www.kaggle.com/code/prashant111/a-guide-on-xgboost-hyperparameters-tuning
    # maybe? https://datascience.stackexchange.com/questions/108233/recommendations-for-tuning-xgboost-hyperparams/108242#108242
    # default params: https://xgboost.readthedocs.io/en/stable/parameter.html

    assert is_classif in {True, False}
    xgb_model_class = xgb.XGBClassifier if is_classif else xgb.XGBRegressor
    score_name = 'roc_auc' if is_classif else 'r'

    if get_pipe is None:
        get_pipe = lambda clf: sklearn.pipeline.Pipeline(steps=[('clf', clf)])

    pathlib.Path(out_dir_path).mkdir(parents=True, exist_ok=True)
    
    train_X = train_feature_df
    train_y = train_target_vec

    if not leave_one_out:
        raise NotImplementedError('only leave_one_out is supported. sorry')

    flat_dicts = []
    assert len(arg_name_and_vals_and_plot_log) == 3
    arg_name = arg_name_and_vals_and_plot_log[0]
    arg_val_to_pred_df = {}
    
    if is_classif:
        score_train_subgroup_name_to_mask = {'all': np.full(len(train_y), True), **score_train_subgroup_name_to_mask}
    
    if skip_cv:
        chosen_arg_val = arg_name_and_vals_and_plot_log[1][0]
        cv_res_df = None
        cv_pipes = None
    else:
        for arg_val in arg_name_and_vals_and_plot_log[1]:
            if arg_name != 'dummy':
                print(f'{arg_name}: {arg_val}')
                xgboost_kwargs[arg_name] = arg_val
            clf = xgb_model_class(**xgboost_kwargs)
            pred_valid_ys = []
            valid_ys = []
            valid_names = []
            cv_pipes = []
            n_splits = cv.get_n_splits(train_X)
            for cv_i, (train, valid) in enumerate(cv.split(train_X, train_y)):
                assert len(valid) == 1 # only leave one out is implemented...
                if cv_i % 20 == 0:
                    print(f'CV split iter: {cv_i+1}/{n_splits}')
                curr_train_X = train_X.iloc[train]
                curr_train_y = train_y[train]
                
                curr_valid_X = train_X.iloc[valid]
                curr_valid_y = train_y[valid][0]
                
                
                pipe = get_pipe(sklearn.base.clone(clf))
                
                if modify_pipe_inside_cv_func is not None:
                    modify_pipe_inside_cv_func(pipe, train)

                # pipe = get_pipe(clf)
                # print(pipe)
                pipe.fit(curr_train_X, curr_train_y)
                
                pred_valid_y = pipe.predict_proba(curr_valid_X)[:, 1][0] if is_classif else pipe.predict(curr_valid_X)[0]

                pred_valid_ys.append(pred_valid_y)
                valid_ys.append(curr_valid_y)
                valid_names.append(train_feature_df.index[valid[0]])

                cv_pipes.append(pipe)

            assert valid_names == list(train_feature_df.index), 'leave one out in unexpected order' # need this for score_train_subgroup_name_to_mask (here and later)
            
            curr_flat_dict = {arg_name: arg_val}
            for subgroup_name, subgroup_mask in score_train_subgroup_name_to_mask.items():
                curr_score_name = f'{score_name}_{subgroup_name}'
                curr_valid_ys = np.array(valid_ys)[subgroup_mask]
                curr_pred_valid_ys = np.array(pred_valid_ys)[subgroup_mask]
                curr_score = (
                    sklearn.metrics.roc_auc_score(curr_valid_ys, curr_pred_valid_ys) if is_classif
                    else scipy.stats.pearsonr(curr_valid_ys, curr_pred_valid_ys)[0]
                )
                curr_flat_dict[curr_score_name] = curr_score


            flat_dicts.append(curr_flat_dict)
            arg_val_to_pred_df[arg_val] = pd.DataFrame({
                'name': valid_names,
                'y': valid_ys,
                'pred_y': pred_valid_ys,
            })
        cv_res_df = pd.DataFrame(flat_dicts)
        
        score_name_for_best_cv = f'{score_name}_{score_train_subgroup_name_for_best_cv}'
        chosen_arg_val = cv_res_df.loc[cv_res_df[score_name_for_best_cv].idxmax(), arg_name]
        best_score = cv_res_df.loc[cv_res_df[score_name_for_best_cv].idxmax(), score_name_for_best_cv]
        print(f'best_{score_name_for_best_cv}: {best_score}')



    final_clf = xgb_model_class(**{arg_name: chosen_arg_val, **xgboost_kwargs})
    final_pipe = get_pipe(final_clf)
    final_pipe.fit(train_X, train_y)

    train_X_transformed_by_final_pipe = get_pipeline_transformed_x(final_pipe, train_X)
    final_pipe_model = final_pipe.steps[-1][1]
    raw_pred = final_pipe_model.predict(train_X_transformed_by_final_pipe, output_margin=True)
    explainer = shap.TreeExplainer(final_pipe_model)
    final_pipe_model_shap_explanation = explainer(train_X_transformed_by_final_pipe)
    shap_values = final_pipe_model_shap_explanation.values
    # make sure the SHAP values add up to marginal predictions
    shap_expected_approx_zero_vals = shap_values.sum(axis=1) + final_pipe_model_shap_explanation.base_values - raw_pred
    assert np.allclose(shap_expected_approx_zero_vals, 0, atol=shap_expected_approx_zero_vals_allclose_atol)
    final_pipe_model_features = list(train_X_transformed_by_final_pipe.columns[shap_values.max(axis=0) > 0])

    relevant_feature_count = len(final_pipe_model_features)
    feature_to_show_count = relevant_feature_count + 1
    plt.close('all')
    ax = shap.plots.beeswarm(
        final_pipe_model_shap_explanation, plot_size=(10, 1 + feature_to_show_count*0.27), show=False, max_display=feature_to_show_count)
    fig = ax.get_figure()
    fig.suptitle(title_prefix_for_plots, fontsize='x-small')
    fig.tight_layout()
    pathlib.Path(f'{out_dir_path}/shap').mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{out_dir_path}/shap/shap_{generic_utils.strip_special_chars_on_edges_and_replace_others_with_underscores(title_prefix_for_plots)}.png')
    plt.close('all')


    if arg_name_and_vals_and_plot_log[2]:
        cv_res_df[f'log_{arg_name}'] = np.log2(cv_res_df[arg_name])
        arg_name = f'log_{arg_name}'
        chosen_arg_val_to_plot = np.log2(chosen_arg_val)
    else:
        chosen_arg_val_to_plot = chosen_arg_val

    if skip_cv:
        best_pred_df = None
    else:
        fig, ax = plt.subplots(figsize=(12, 5))
        sb.lineplot(cv_res_df, x=arg_name, y=score_name_for_best_cv, ax=ax)
        ax.axvline(chosen_arg_val_to_plot, color='red', alpha=0.3, linestyle='--')
        fig.tight_layout()
        fig.savefig(f'{out_dir_path}/xgboost_{generic_utils.strip_special_chars_on_edges_and_replace_others_with_underscores(title_prefix_for_plots)}_{score_name}_vs_{arg_name}.png')

        best_pred_df = arg_val_to_pred_df[chosen_arg_val]
        # print(best_leave_one_out_df)

        # print(f'score_train_subgroup_name_to_mask: {score_train_subgroup_name_to_mask}')
        for subgroup_name, subgroup_mask in score_train_subgroup_name_to_mask.items():
            # print(f'subgroup_name: {subgroup_name}')
            curr_repr = f'{title_prefix_for_plots}{subgroup_name}'
            curr_pred_df = best_pred_df[subgroup_mask]
            if is_classif:
                fig, ax = plot_auc(curr_pred_df, title_suffix=curr_repr)
                fig.tight_layout()
                fig.savefig(f'{out_dir_path}/xgboost_{generic_utils.strip_special_chars_on_edges_and_replace_others_with_underscores(curr_repr)}_best_{score_name}.png')
            else:
                regression_prediction_scatter(
                    curr_pred_df, 
                    ax_title=f'xgboost train ({score_name})\n{curr_repr}',
                    **prediction_scatterplot_kwargs,
                )


        if is_classif:
            binary_classification_prediction_barplot(
                best_pred_df, 
                title_prefix_for_plots=f'xgboost train ({score_name}), {title_prefix_for_plots}',
                **prediction_barplot_kwargs,
            )
        else:
            regression_prediction_scatter(
                best_pred_df, 
                ax_title=f'xgboost train ({score_name})\n{title_prefix_for_plots}',
                **prediction_scatterplot_kwargs,
            )
    res_dict = {
        'cv_res_df': cv_res_df,
        'train_X': train_X,
        'final_pipe': final_pipe,
        'train_X_transformed_by_final_pipe': train_X_transformed_by_final_pipe,
        'title_prefix_for_plots': title_prefix_for_plots,
        'arg_val_to_leave_one_out_df': arg_val_to_pred_df,
        'train_prediction_df': best_pred_df,
        'cv_pipes': cv_pipes,
        'final_pipe_model_shap_explanation': final_pipe_model_shap_explanation,
        'final_pipe_model_features': final_pipe_model_features,
    }

    if test_feature_df is not None:
        assert (test_feature_df.columns == train_feature_df.columns).all()

        test_prediction_df = pd.DataFrame({
            'name': test_feature_df.index,
            'y': test_target_vec,
            'pred_y': final_pipe.predict_proba(test_feature_df)[:, 1] if is_classif else final_pipe.predict(test_feature_df),
        })
        res_dict['test_prediction_df'] = test_prediction_df
        
        # return test_feature_df, res_dict['test_predicted_y']
        if is_classif:
            binary_classification_prediction_barplot(
                test_prediction_df, 
                title_prefix_for_plots=f'xgboost test ({score_name}), {title_prefix_for_plots}',
                **prediction_barplot_kwargs,
            )
        else:
            regression_prediction_scatter(
                test_prediction_df, 
                ax_title=f'xgboost test ({score_name})\n{title_prefix_for_plots}',
                **prediction_scatterplot_kwargs,
            )
    return res_dict


def logistic_regr_cv(
        train_feature_df, train_target_vec, 
        test_feature_df=None,
        test_target_vec=None,
        scoring='accuracy',
        prediction_barplot_kwargs=dict(),
        Cs=np.logspace(-3, 3, 10),
        cv=None,
        score_statistics_to_plot='mean_and_std', # 'quantiles'
        title_prefix_for_plots='',
        refit=True,
        use_manual_refit=None,
        get_pipe=None,
        **logistic_regr_cv_kwargs,
):
    pathlib.Path('temp/logistic_regr_cv').mkdir(parents=True, exist_ok=True)

    if get_pipe is None:
        get_pipe = lambda clf: sklearn.pipeline.Pipeline(steps=[('clf', clf)])

    if ('penalty' not in logistic_regr_cv_kwargs) and (logistic_regr_cv_kwargs['penalty'] != 'l1'):
        raise NotImplementedError("only penalty='l1' is supported")
    if (use_manual_refit is None):
        use_manual_refit = refit # for some reason, with refit=True, seems like always intercept_==0.
    
    train_X = train_feature_df.values
    train_y = train_target_vec
    test_X = None if test_feature_df is None else test_feature_df.values
    lambda_to_leave_one_out_df = {}
    
    if isinstance(cv, sklearn.model_selection.LeaveOneOut) and (scoring == 'roc_auc'):
        roc_auc_cv = True
        print('manual CV because of LeaveOneOut and roc_auc')
        lambda_flat_dicts = []
        for curr_c in Cs:
            curr_lambda = 1 / curr_c
            print(f'curr_lambda: {curr_lambda}')
            clf = sklearn.linear_model.LogisticRegression(
                C=curr_c, 
                **{k: v for k,v in logistic_regr_cv_kwargs.items() if k not in {'cv', 'Cs', 'refit'}},
            )
            pipe = get_pipe(clf)

            pred_valid_ys = []
            valid_ys = []
            valid_names = []
            non_zero_coef_counts = []
            n_splits = cv.get_n_splits(train_X)
            for cv_i, (train, valid) in enumerate(cv.split(train_X, train_y)):
                print(f'CV split iter: {cv_i+1}/{n_splits}')
                curr_train_X = train_X[train]
                curr_train_y = train_y[train]
                
                curr_valid_X = train_X[valid]
                curr_valid_y = train_y[valid][0]
                
                pred_valid_y = pipe.fit(curr_train_X, curr_train_y).predict_proba(curr_valid_X)[:, 1][0]
                
                pred_valid_ys.append(pred_valid_y)
                valid_ys.append(curr_valid_y)
                non_zero_coef_counts.append((pipe['clf'].coef_ != 0).sum())
                assert len(valid) == 1
                valid_names.append(train_feature_df.index[valid[0]])
                # if cv_i == 10:
                #     break
            roc_auc = sklearn.metrics.roc_auc_score(np.array(valid_ys), np.array(pred_valid_ys))
            non_zero_coef_count_mean = np.mean(non_zero_coef_counts)
            lambda_flat_dicts.append({
                'lambd': curr_lambda,
                'roc_auc': roc_auc,
                'non_zero_coef_count_mean': non_zero_coef_count_mean,
            })
            lambda_to_leave_one_out_df[curr_lambda] = pd.DataFrame({
                'name': valid_names,
                'y': valid_ys,
                'pred_y': pred_valid_ys,
            })
        lambda_df = pd.DataFrame(lambda_flat_dicts)
        lambda_df['log_lambda'] = np.log2(lambda_df['lambd'])

        fig, ax = plt.subplots(figsize=(12, 5))
        sb.lineplot(lambda_df, x='log_lambda', y='roc_auc', ax=ax)

        chosen_lambda = lambda_df.loc[lambda_df['roc_auc'].idxmax(), 'lambd']
        chosen_lambda_leave_one_out_df = lambda_to_leave_one_out_df[chosen_lambda]
        clf = sklearn.linear_model.LogisticRegression(
            C=1 / chosen_lambda, 
            **{k: v for k,v in logistic_regr_cv_kwargs.items() if k not in {'cv', 'Cs', 'refit'}},
        )

        final_pipe = get_pipe(clf)
        final_pipe.fit(train_X, train_y)
        final_features = list(train_feature_df.columns[final_pipe['select'].get_support()])
        res_dict = {}

    else:
        roc_auc_cv = False
        cv_obj = sklearn.linear_model.LogisticRegressionCV(scoring=scoring, cv=cv, **logistic_regr_cv_kwargs)
        cv_obj.fit(train_X, train_y)

        chosen_c = cv_obj.C_
        assert len(chosen_c) == 1
        chosen_c = chosen_c[0]
        chosen_lambda = 1 / chosen_c

        manual_refit_cv_obj = sklearn.linear_model.LogisticRegression(
            C=chosen_c, 
            **{k: v for k,v in logistic_regr_cv_kwargs.items() if k not in {'cv', 'Cs', 'refit'}},
        )
        manual_refit_cv_obj.fit(train_X, train_y)
        
        final_pipe = (manual_refit_cv_obj if use_manual_refit else cv_obj)

        lambda_df = pd.DataFrame({
            'log_lambda': np.log2(1 / cv_obj.Cs_),
            'score_median': np.median(cv_obj.scores_[1], axis=0),
            'score_0.9': np.quantile(cv_obj.scores_[1], 0.9, axis=0),
            'score_0.1': np.quantile(cv_obj.scores_[1], 0.1, axis=0),
            'score_std': cv_obj.scores_[1].std(axis=0),
            'score_mean': cv_obj.scores_[1].mean(axis=0),
            'non_zero_coef_count_mean': (cv_obj.coefs_paths_[1] != 0).sum(axis=2).mean(axis=0),
        })

        # Plot the results
        fig, ax = plt.subplots(figsize=(12, 5))
        # plt.figure(figsize=(10, 6))
        
        errorbar_kwargs = dict(fmt='-o')
        if isinstance(cv, sklearn.model_selection.LeaveOneOut):
            y_label = f'{scoring} (mean)'
            median_or_mean = 'mean'
        else:
            errorbar_kwargs['ecolor'] = 'grey'
            errorbar_kwargs['capsize'] = 0
            if score_statistics_to_plot == 'quantiles':
                lambda_df['score_0.9_median_diff'] = lambda_df['score_0.9'] - lambda_df['score_median']
                lambda_df['score_0.1_median_diff'] = lambda_df['score_median'] - lambda_df['score_0.1']
                errorbar_kwargs['yerr'] = lambda_df[['score_0.1_median_diff', 'score_0.9_median_diff']].T.to_numpy()
                y_label = f'{scoring} (0.1, 0.5, 0.9 quantiles)'
                median_or_mean = 'median'
            elif score_statistics_to_plot == 'mean_and_std':
                errorbar_kwargs['yerr'] = lambda_df['score_std']
                y_label = f'{scoring} (mean +- std)'
                median_or_mean = 'mean'
            else:
                raise NotImplementedError(f'score_statistics_to_plot: {score_statistics_to_plot}')
        
        res_dict = {
            'cv_obj': cv_obj,
            'manual_refit_cv_obj': manual_refit_cv_obj,
        }

        ax.errorbar(lambda_df['log_lambda'], lambda_df[f'score_{median_or_mean}'], **errorbar_kwargs)
        ax.set_ylabel(y_label)
    
    ax.set_xlabel('log($\lambda$)')
    for _, row in lambda_df.iterrows():
        ax.text(
            row['log_lambda'], 0.96, 
            str(int(row['non_zero_coef_count_mean'])),
            ha='center', 
            transform=matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes),
            fontsize='small',
        )
    
    coef_vec = mc.ut.to_numpy_vector(final_pipe['clf'].coef_)
    nonzero_coef_count = (coef_vec != 0).sum()

    train_sample_count = train_X.shape[0]
    feature_count = final_pipe['select'].get_support().sum()
    ax.set_title(title_prefix_for_plots + f'train_sample_count: {train_sample_count}, feature_count: {feature_count}, fold_count={cv.get_n_splits(train_X)}, #coefs={nonzero_coef_count}')
    ax.axvline(np.log2(chosen_lambda), color='red', linestyle='--', alpha=0.3)
    fig.tight_layout()
    fig.savefig(f'temp/logistic_regr_cv/logist_regr_{generic_utils.strip_special_chars_on_edges_and_replace_others_with_underscores(title_prefix_for_plots)}_{scoring}_vs_lambda.png')

    # ax.grid(True)
    train_prediction_df = chosen_lambda_leave_one_out_df if roc_auc_cv else pd.DataFrame({
        'name': train_feature_df.index,
        'y': train_y,
        'pred_y': final_pipe.predict_proba(train_X)[:, 1],
    })

    binary_classification_prediction_barplot(
        train_feature_df, 
        train_prediction_df['pred_y'], 
        title_prefix_for_plots=f'logist_regr train ({scoring}), {title_prefix_for_plots}',
        **prediction_barplot_kwargs,
    )
    res_dict = {
        **res_dict,
        'final_pipe': final_pipe,
        'train_prediction_df': train_prediction_df,
        'final_features': final_features,
        'score_lambda_fig_ax': (fig, ax),
    }

    if test_feature_df is not None:
        assert (test_feature_df.columns == train_feature_df.columns).all()

        test_prediction_df = pd.DataFrame({
            'name': test_feature_df.index,
            'y': test_target_vec,
            'pred_y': final_pipe.predict_proba(test_X)[:, 1],
        })
        res_dict['test_prediction_df'] = test_prediction_df
        
        # return test_feature_df, res_dict['test_predicted_y']
        binary_classification_prediction_barplot(
            test_feature_df, 
            test_prediction_df['pred_y'], 
            title_prefix_for_plots=f'logist_regr test ({scoring}), {title_prefix_for_plots}',
            **prediction_barplot_kwargs,
        )
        # if needed to do this from a new cv_obj:
        # logistic_regr_obj = sklearn.linear_model.LogisticRegression()
        # logistic_regr_obj.coef_ = cv_res[0].coef_
        # logistic_regr_obj.intercept_ = cv_res[0].intercept_
        # logistic_regr_obj.classes_ = np.array([0,1])
        # pred_y = logistic_regr_obj.predict_proba(curr_log_ratio_df.T.values)
        # (pred_y == cv_res[1]).all()

    return res_dict
        



def get_bh_fdr_corrected_pvals(pvals, alpha):
    raise RuntimeError('use scipy.stats.false_discovery_control. the following doesnt follow the exact algorithm that should be used, i.e., if the lowest pval is 1e-5 and the second lowest is 1.3e-5, i would give a lower corrected pval for the second, which doesnt make sense. the right correction would be to give the same corrected pval for both. see https://stats.stackexchange.com/questions/238458/whats-the-formula-for-the-benjamini-hochberg-adjusted-p-value/402217#402217')
    df = pd.DataFrame({'pval': pvals, 'orig_i': range(len(pvals))})
    df.sort_values('pval', inplace=True)
    df['pval_rank'] = range(1, len(pvals) + 1)
    df['approx_expected_pval_if_uniform_between_0_and_alpha'] = df['pval_rank'] / len(pvals) * alpha
    df['bh_fdr_corrected_pval'] = df['pval'] / df['approx_expected_pval_if_uniform_between_0_and_alpha']
    return df.sort_values('orig_i')['bh_fdr_corrected_pval'].to_numpy()
    
def mw_separating_y(X, y, only_higher_in=None, ignore_nans=False, sample_mask=None):
    if sample_mask is None:
        sample_mask = np.full(len(y), True)
    # else:
    #     print('mw_separating_y got sample_mask')
    
    x_is_df = isinstance(X, pd.DataFrame)
    if x_is_df:
        x1 = X[(y == 1) & sample_mask]
        x0 = X[(y == 0) & sample_mask]
    else:
        assert False, 'mw_separating_y(): please pass a df'
        x1 = X[(y == 1) & sample_mask, :]
        x0 = X[(y == 0) & sample_mask, :]
    
    if only_higher_in == 1:
        alternative = 'greater'
    elif only_higher_in == 0:
        alternative = 'less'
    else:
        assert only_higher_in is None
        alternative = 'two-sided'

    pvals = []
    u_minus_mean_u_vals = []
    for i in range(X.shape[1]):
        if x_is_df:
            curr_x1 = x1.iloc[:, i]
            curr_x0 = x0.iloc[:, i]
        else:
            curr_x1 = x1[:, i]
            curr_x0 = x0[:, i]
        if ignore_nans:
            # print(f'curr_x0: {type(curr_x0)} {curr_x0}')
            # print(f'curr_x0.dtype: {curr_x0.dtype}')
            curr_x0 = curr_x0[~np.isnan(curr_x0)]
            curr_x1 = curr_x1[~np.isnan(curr_x1)]
            if len(curr_x0) == 0 or len(curr_x1) == 0:
                pvals.append(1)
                u_minus_mean_u_vals.append(0)
                continue
        mw_res = generic_utils.perform_mw_test(curr_x1, curr_x0, alternative=alternative)
        pvals.append(mw_res['pvalue'])
        u_minus_mean_u_vals.append(mw_res['u_minus_mean_u'])
    return u_minus_mean_u_vals, pvals


# class SelectFdrAllowingNans(sklearn.feature_selection.GenericUnivariateSelect, sklearn.base.BaseEstimator):
# class SelectFdrAllowingNans(sklearn.feature_selection.SelectorMixin, sklearn.base.BaseEstimator):
class SelectorAllowingNans(sklearn.base.TransformerMixin, sklearn.base.BaseEstimator):
    # ugh. seems like we cant use sklearn.feature_selection.GenericUnivariateSelect nor sklearn.feature_selection.SelectorMixin, because they don't allow nans.

    # https://xgboost.readthedocs.io/en/stable/faq.html: "XGBoost supports missing values by default. In tree algorithms, branch directions for missing values are learned during training"
    # for a slightly more elaborate explanation: https://datascience.stackexchange.com/questions/15305/how-does-xgboost-learn-what-are-the-inputs-for-missing-values/15306#15306

    # used https://stackoverflow.com/questions/25250654/how-can-i-use-a-custom-feature-selection-function-in-scikit-learns-pipeline/53870180#53870180 to help with understanding how to implement this

    def __init__(self, get_scores_and_support_mask_func, get_scores_and_support_mask_func_kwargs=dict()):
        self.get_scores_and_support_mask_func = get_scores_and_support_mask_func
        self.get_scores_and_support_mask_func_kwargs = get_scores_and_support_mask_func_kwargs

    def fit(self, X, y=None, get_scores_and_support_mask_func_kwargs=dict()):
        # NOTE: ugh. seems like sklearn.compose.ColumnTransformer.set_params() works differently - it actually passes params to fit(), while sklearn.pipeline.Pipeline.set_params() directly sets attributes of the object. IIUC
        if get_scores_and_support_mask_func_kwargs:
            assert not self.get_scores_and_support_mask_func_kwargs
            self.get_scores_and_support_mask_func_kwargs = get_scores_and_support_mask_func_kwargs
        del get_scores_and_support_mask_func_kwargs # shouldn't use it anymore - as it is different for sklearn.pipeline.Pipeline.set_params() and sklearn.compose.ColumnTransformer.set_params()
        
        # print('SelectorAllowingNans.fit(), self.get_scores_and_support_mask_func_kwargs:', self.get_scores_and_support_mask_func_kwargs)

        self.feature_names_in_ = X.columns
        self.n_features_in_ = X.shape[1]
        self._support_mask, self._scores = self.get_scores_and_support_mask_func(X, y, **self.get_scores_and_support_mask_func_kwargs)
        return self

    def get_support(self):
        return self._support_mask
    
    def get_scores(self):
        return self._scores
    
    def transform(self, X):
        return X.loc[:, self._support_mask]

    

def get_y_associated_feature_mask(X, y, fdr_alpha=0.05, also_return_corrected_pvals=False, separating_y_func=mw_separating_y, **separating_y_func_kwargs):
    # print(f'get_y_associated_feature_mask(): {separating_y_func.__name__}, {separating_y_func_kwargs}')
    _, pvals = separating_y_func(X, y, **separating_y_func_kwargs)
    
    bh_fdr_corrected_pvals = scipy.stats.false_discovery_control(pvals)
    
    if 0:
        plt.close('all')
        fig, ax = plt.subplots()
        sb.histplot(np.log2(pvals), bins=20, ax=ax)

        fig, ax = plt.subplots()
        sb.histplot(np.log2(bh_fdr_corrected_pvals), bins=20, ax=ax)
        raise

    feature_mask = bh_fdr_corrected_pvals < fdr_alpha
    print(f'get_y_associated_feature_mask(): feature_mask.sum(): {feature_mask.sum()}')
    if also_return_corrected_pvals:
        return feature_mask, bh_fdr_corrected_pvals
    return feature_mask

def get_y_pos_and_neg_associated_feature_masks(X, y, fdr_alpha=0.05, separating_y_func=mw_separating_y, **separating_y_func_kwargs):
    scores, pvals = separating_y_func(X, y, **separating_y_func_kwargs)
    scores = np.array(scores)
    bh_fdr_corrected_pvals = scipy.stats.false_discovery_control(pvals)
    

    # NOTE: assuming scores are positive for positive associations and negative for negative associations
    
    feature_mask = bh_fdr_corrected_pvals < fdr_alpha

    # print(bh_fdr_corrected_pvals)
    # print(f'feature_mask.sum(): {feature_mask.sum()}')
    pos_feature_mask = feature_mask & (scores > 0)
    neg_feature_mask = feature_mask & (scores < 0)

    # print(f'get_y_pos_and_neg_associated_feature_masks(): pos_feature_mask.sum(): {pos_feature_mask.sum()}')
    # print(f'get_y_pos_and_neg_associated_feature_masks(): neg_feature_mask.sum(): {neg_feature_mask.sum()}')
    return bh_fdr_corrected_pvals, [pos_feature_mask, neg_feature_mask]

def get_pooled_y_associated_features(X, y, **associated_feature_func_kwargs):
    print('get_pooled_y_associated_features()')
    print(f'X.shape: {X.shape}')
    raise
    return X

class poolColsTransformer(sklearn.base.TransformerMixin, sklearn.base.BaseEstimator):
    def __init__(self, feature_names_out, get_col_scores_and_masks_func, agg_func=np.mean, get_col_scores_and_masks_func_kwargs=dict()):
        self.agg_func = agg_func
        self.get_col_scores_and_masks_func = get_col_scores_and_masks_func
        self.get_col_scores_and_masks_func_kwargs = get_col_scores_and_masks_func_kwargs
        self.feature_names_out = feature_names_out
    
    def fit(self, X, y=None, get_col_scores_and_masks_func_kwargs=dict()):
        # print('poolColsTransformer.fit() self: ', self)
        # NOTE: ugh. seems like sklearn.compose.ColumnTransformer set_params() works differently - it actually passes params to fit(), while sklearn.pipeline.Pipeline.set_params() sets attributes of the object. IIUC
        if get_col_scores_and_masks_func_kwargs:
            assert not self.get_col_scores_and_masks_func_kwargs
            self.get_col_scores_and_masks_func_kwargs = get_col_scores_and_masks_func_kwargs
        del get_col_scores_and_masks_func_kwargs # shouldn't use it anymore - as it is different for sklearn.pipeline.Pipeline.set_params() and sklearn.compose.ColumnTransformer.set_params()

        self.feature_names_in_ = X.columns

        # assert X.dtype.kind == 'f', 'maybe X contains booleans? should only contain floats...'
        assert X.values.dtype.kind == 'f', 'maybe X contains booleans? should only contain floats...'

        self._col_scores, self._col_masks = self.get_col_scores_and_masks_func(X, y, **self.get_col_scores_and_masks_func_kwargs)
        assert isinstance(self._col_masks, list)
        assert len(self._col_scores) == X.shape[1]

        pooled_feature_count = len(self._col_masks)
        assert len(self.feature_names_out) == pooled_feature_count
        for i in range(pooled_feature_count):
            self.feature_names_out[i] += f'_{self._col_masks[i].sum()}'
            # print(f'self.feature_names_out[{i}]: {self.feature_names_out[i]}')

        return self

    def get_col_scores(self):
        return self._col_scores
    
    def get_col_masks(self):
        return self._col_masks

    def transform(self, X, *args, **kwargs):
        # print(f'poolColsTransformer.transform()')
        if args:
            # print(f'args: {args}')
            # ugh. i don't understand why more args are passed to transform
            assert len(args) == 1
            assert args[0] is None
        assert not kwargs
        # print(f'kwargs: {kwargs}')

        agg_col_name_to_vals = {}
        assert len(self._col_masks) == len(self.feature_names_out)
        for col_mask, agg_feature_name in zip(self._col_masks, self.feature_names_out):
            if col_mask.any():
                # print(f'poolColsTransformer.transform(): col_mask.sum(): {col_mask.sum()}')
                vals = self.agg_func(X.loc[:, col_mask], axis=1)
            else:
                vals = np.full(X.shape[0], np.nan)
            agg_col_name_to_vals[agg_feature_name] = vals
        
        trans_X = pd.DataFrame(agg_col_name_to_vals, index=X.index)
        assert trans_X.shape[0] == X.shape[0]
        assert trans_X.shape[1] == len(self._col_masks)

        return trans_X

    def fit_transform(self, X, y, sample_weight=None, groups=None):
        return self.fit(X, y, sample_weight).transform(X, groups)
    
    def get_feature_names_out(self, input_features=None):
        # print('get_feature_names_out:', self.feature_names_out)
        return self.feature_names_out


def get_feature_means(X, y, y_val=None):
    sample_mask = np.full(len(y), True) if (y_val is None) else y == y_val
    X = X[sample_mask,:]
    return X.mean(axis=0)

class PassthroughTransformer(sklearn.base.TransformerMixin, sklearn.base.BaseEstimator):
    def __init__(self, feature_names):
        self.feature_names = feature_names

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        assert isinstance(X, pd.DataFrame), 'PassthroughTransformer.transform(): please pass a df'
        return X

    def get_feature_names_out(self, input_features=None):
        # print('get_feature_names_out:', self.feature_names)
        return self.feature_names




