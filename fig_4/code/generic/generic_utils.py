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


try: import __builtin__ as builtins
except ImportError: import builtins

# 200624: It seems that this is the default (what I get from np.geterr() after importing numpy): {'divide': 'warn', 'over': 'warn', 'under': 'ignore', 'invalid': 'warn'}
# https://numpy.org/doc/stable/reference/generated/numpy.seterr.html says: Underflow: result so close to zero that some precision was lost.
# so I guess it is Ok that I just ignore underflow problems?
# np.seterr(all='raise')
np.seterr(divide='raise', over='raise', invalid='raise')

INPUT_FILE_PATH_ARG_PREFIX = 'input_file_path_'
MULTIPLE_INPUT_FILE_PATH_ARG_PREFIX = 'multiple_input_file_paths_'
OUTPUT_FILE_PATH_ARG_PREFIX = 'output_file_path_'
OUTPUT_FILE_PATH_FOR_CACHING_ONLY_ARG = f'{OUTPUT_FILE_PATH_ARG_PREFIX}for_caching_only' # output_file_path_for_caching_only
HASH_FILE_NAME_SUFFIX = '.kind_of_hash_by_oren.txt'
FUNC_INPUT_ARGS_IN_PREV_CALL_FILE_NAME_SUFFIX = '.input_args_in_prev_call.txt'
INTERNAL_FILES_TO_SKIP_REDUNDANT_CALCULATIONS_DIR_PATH = 'internal___hashes_etc_to_skip_redundant_calculations'

ColumnNameAndValuesAndRanges = collections.namedtuple('ColumnNameAndValuesAndRanges', ['column_name', 'values_and_ranges'])
Range = collections.namedtuple('Range', ['start', 'end'])
Range.__str__ = lambda r: f'Range({r.start}, {r.end})'

WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES = ['word', 'count', 'freq']
WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES_SET = set(WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES)

WORD_COUNTS_AND_FREQS_COLUMN_NAME_TO_DTYPE = {
    'word': 'string',
    'count': int,
    'freq': np.float64,
}

prop_cycle = plt.rcParams['axes.prop_cycle']
TEN_DEFAULT_MATPLOTLIB_COLORS = prop_cycle.by_key()['color']
TEN_DEFAULT_MATPLOTLIB_COLORS_AND_SOME_MORE_KIND_OF_DISTINGUISHABLE = TEN_DEFAULT_MATPLOTLIB_COLORS + [
    'gold',
    'teal',
    'olive',
    'lawngreen',
    'olivedrab',
    'cadetblue',
    'mediumseagreen',
    # 'mediumorchid',
    'mediumpurple',
    'deepskyblue',
    'indianred',
    'tan',
    'darkkhaki',
    'steelblue',
    'lightgreen',
    'powderblue',
    'plum',
    'mediumaquamarine',
    'navajowhite',
]

COLORBREWER_SET1_WITHOUT_GREY_AND_YELLOW = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00' ,'#a65628','#f781bf']
# was used as NON_REF_VARIANT_COLORS in the programmed inversions paper code
MANY_KIND_OF_DISTINGUISHABLE_COLORS = COLORBREWER_SET1_WITHOUT_GREY_AND_YELLOW + [
    'gold',
    'grey',
    'teal',
    'olive',
    'olivedrab',
    'lawngreen',
    'cadetblue',
    'mediumseagreen',
    # 'mediumorchid',
    'mediumpurple',
    'deepskyblue',
    'indianred',
    'tan',
    'darkkhaki',
    'steelblue',
    'lightgreen',
    'powderblue',
    'plum',
    'mediumaquamarine',
    'navajowhite',
]

PYTHON_BUILTIN_TYPES = (int, float, complex, bool, str, type(None), list, tuple, dict, set, frozenset)

# pathlib.Path(INTERNAL_FILES_TO_SKIP_REDUNDANT_CALCULATIONS_DIR_PATH).mkdir(parents=True, exist_ok=True)

def merge_mutual_exclusive_dicts(dicts):
    assert not get_common_vals_in_sets([set(x) for x in dicts])

    merged_dict = {}
    for d in dicts:
        merged_dict = {
            **merged_dict,
            **d,
        }
    return merged_dict

def get_max_finite(num_series):
    return num_series[np.isfinite(num_series)].max()

def get_min_finite(num_series):
    return num_series[np.isfinite(num_series)].min()

def get_min_strictly_positive(num_series):
    return num_series[num_series > 0].min()

def get_range_tuple_from_edges(edges):
    assert pd.Series(edges).is_monotonic_increasing
    return tuple(Range(x,y) for x,y in zip(edges[:-1], edges[1:]))

def print_and_write_to_log(msg):
    print(msg)
    logging.debug(msg)

def read_text_file(file_path, if_file_doesnt_exist_return_empty_str=False):
    if if_file_doesnt_exist_return_empty_str and (not os.path.isfile(file_path)):
        return ''
    with open(file_path) as f:
        return f.read()

def get_file_size(file_path):
    return os.stat(file_path).st_size

def is_file_empty(file_path):
    return get_file_size(file_path) == 0

def is_text_file_empty_or_containing_only_whitespaces(file_path):
    return not read_text_file(file_path).strip()

def read_bin_file(file_path):
    with open(file_path, 'rb') as f:
        return f.read()

def write_text_file(file_path, file_contents):
    with open(file_path, 'w') as f:
        return f.write(file_contents)

def append_to_text_file(file_path, str_to_append, append_newline_first=True):
    if str_to_append:
        if append_newline_first:
            str_to_append = '\n' + str_to_append
        with open(file_path, 'a') as f:
            return f.write(str_to_append)
    else:
        print('append_to_text_file received empty str_to_append, so doing nothing')

def write_empty_file(file_path):
    with open(file_path, 'w') as f:
        pass

def add_suffix_to_file_name_while_keeping_extension(file_path, suffix):
    file_path_without_extension, extension = os.path.splitext(file_path)
    return file_path_without_extension + suffix + extension
def add_prefix_to_file_name_while_keeping_extension(file_path, prefix):
    dir_path, file_name = os.path.split(file_path)
    return os.path.join(dir_path, prefix + file_name)

def get_num_of_lines_in_text_file(file_path):
    return int(subprocess.check_output(['wc', '-l', file_path]).split()[0])

def get_len_of_longest_line_in_text_file(file_path):
    return int(subprocess.check_output(['wc', '-L', file_path]).split()[0])

def get_num_of_false_values(values):
    return sum(x == False for x in values)

def does_any_str_in_strs1_contain_any_str_in_strs2(strs1, strs2):
    # print('\nstrs1')
    # print(strs1)
    return any(str2 in str1 for (str1, str2) in itertools.product(strs1, strs2))

def does_any_str_in_strs1_ends_with_any_str_in_strs2(strs1, strs2):
    # print('\nstrs1')
    # print(strs1)
    return any(str1.endswith(str2) for (str1, str2) in itertools.product(strs1, strs2))

def get_public_attribute_names_set(obj):
    return set(attr for attr in dir(obj) if not attr.startswith('_'))


def get_file_sha256_in_hex(file_path):
    return hashlib.sha256(read_bin_file(file_path)).hexdigest()


def get_str_sha256_in_hex(str_to_hash):
    return hashlib.sha256(str_to_hash.encode()).hexdigest()


def get_opposite_key_to_set_dict(key_to_set):
    opposite_key_to_set = collections.defaultdict(set)
    for key, set_ in key_to_set.items():
        for val in set_:
            opposite_key_to_set[val].add(key)
    opposite_key_to_set = dict(opposite_key_to_set) # I don't want a defaultdict moving around.
    return opposite_key_to_set

def print_dict(dicti, file=None, dict_name=None):
    print(file=file)
    if dict_name:
        print(f'{dict_name}:', file=file)
    for k, v in dicti.items():
        print(f'{k}: {v}', file=file)
    print(file=file)


def print_iterable(list_to_print, file=None, iterable_name=None):
    print(file=file)
    if iterable_name:
        print(f'{iterable_name}:', file=file)
    for x in list_to_print:
        print(x, file=file)
    print(file=file)

class KeysStep:
    pass

def get_values_in_nested_obj(obj, path_to_values, predicate_func=None, return_whole_paths=False, curr_path=()):
    assert type(path_to_values) == tuple
    if not path_to_values:
        if predicate_func:
            try:
                predicate_satisfied = predicate_func(obj)
            except Exception as err:
                print(f'predicate_func(obj) raised an error when called by get_values_in_nested_obj().')
                print(f'err:\n{err}')
                print(f'predicate_func(obj):\n{predicate_func.__name__}({obj})')
                raise
            if not predicate_satisfied:
                return []
        if return_whole_paths:
            return [(obj, curr_path)]
        return [obj]

    if obj is None:
        return []

    curr_step = path_to_values[0]
    next_steps = path_to_values[1:]
    if type(curr_step) == str:
        assert type(obj) == dict
        if curr_step not in obj:
            return []
        return get_values_in_nested_obj(obj[curr_step], next_steps, predicate_func=predicate_func, return_whole_paths=return_whole_paths, curr_path=(curr_path + (curr_step,)))
    elif type(curr_step) in {set, frozenset}:
        assert type(obj) == dict
        curr_steps = curr_step
        curr_step = None # just playing it safe.
        values = []
        for curr_step in curr_steps:
            if curr_step in obj:
                values.extend(get_values_in_nested_obj(obj[curr_step], next_steps, predicate_func=predicate_func, return_whole_paths=return_whole_paths,
                                                       curr_path=(curr_path + (curr_step,))))
        return values
    elif curr_step == KeysStep:
        values = []
        assert type(obj) == dict
        for k in obj:
            values.extend(get_values_in_nested_obj(k, next_steps, predicate_func=predicate_func, return_whole_paths=return_whole_paths, curr_path=(curr_path + (KeysStep,))))
        return values
    elif curr_step == Ellipsis:
        values = []
        if type(obj) == dict:
            for k, v in obj.items():
                values.extend(get_values_in_nested_obj(v, next_steps, predicate_func=predicate_func, return_whole_paths=return_whole_paths, curr_path=(curr_path + (k,))))
        elif type(obj) in (tuple, list, set, frozenset):
            for x in obj:
                values.extend(get_values_in_nested_obj(x, next_steps, predicate_func=predicate_func, return_whole_paths=return_whole_paths,
                                                       curr_path=(curr_path + ('element_in_iterable',))))
        else:
            print('obj')
            print(obj)
            print('type(obj)')
            print(type(obj))
            raise NotImplementedError
        return values
    else:
        print('curr_step')
        print(curr_step)
        print('type(curr_step)')
        print(type(curr_step))
        raise NotImplementedError

# def get_paths_to_values_in_nested_obj_that_satisfy_predicate(obj, path_to_values, predicate_func):


def get_substr_indices_including_overlapping(a_str, substr):
    # We trust str.find to do things efficiently, so we just try until it returns -1.
    substr_indices = []
    search_start_index = 0
    while True:
        found_index = a_str.find(substr, search_start_index)
        if found_index == -1:
            return substr_indices
        substr_indices.append(found_index)
        search_start_index = found_index + 1


assert get_substr_indices_including_overlapping('GATATATGCATATACTT', 'ATAT') == [1, 3, 9]
assert get_substr_indices_including_overlapping('.12.456.', '.') == [0, 3, 7]
assert get_substr_indices_including_overlapping('...', '.') == [0, 1, 2]
assert get_substr_indices_including_overlapping('...a', '.') == [0, 1, 2]
assert get_substr_indices_including_overlapping('a...', '.') == [1, 2, 3]
assert get_substr_indices_including_overlapping('...', '..') == [0, 1]
assert get_substr_indices_including_overlapping('...', '...') == [0]


def get_defaultdict_counter_using_numpy_unique___probably_faster_than_Counter_if_np_or_pd_obj(values):

    # Kind of equivalent to the following, but hopefully faster:
    # return collections.Counter(values)

    # When I tried this on a very long list of numbers, Counter was (very roughly) 4 times slower if I gave it a np.array or a Series, and the np.unique based implementation was
    # (very roughly) 1.5 times slower if I gave it a plain python list.
    unique_values, counts = np.unique(values, return_counts=True)
    return collections.defaultdict(int, zip(unique_values, counts))


@contextlib.contextmanager
def chdir_context_manager(new_dir):
    orig_dir = os.getcwd()
    os.chdir(os.path.expanduser(new_dir))
    try:
        yield
    finally:
        os.chdir(orig_dir)

@contextlib.contextmanager
def allow_numpy_division_by_zero_context_manager():
    orig_err = np.geterr()
    np.seterr(divide='ignore', invalid='ignore')
    try:
        yield
    finally:
        np.seterr(divide=orig_err['divide'], invalid=orig_err['invalid'])

@contextlib.contextmanager
def keep_ax_lims_context_manager(ax):
    orig_xlim = ax.get_xlim()
    orig_ylim = ax.get_ylim()
    try:
        yield
    finally:
        ax.set_xlim(orig_xlim)
        ax.set_ylim(orig_ylim)

def get_roughly_equal_chunk_is(total_len, chunk_count, random_state=0):
    assert chunk_count <= total_len
    chunk_sizes = [total_len // chunk_count] * chunk_count
    leftovers = total_len % chunk_count
    random.seed(random_state)
    for i in random.sample(range(chunk_count), leftovers):
        # print(i)
        chunk_sizes[i] += 1
    assert sum(chunk_sizes) == total_len
    return np.repeat(range(chunk_count), chunk_sizes)

def merge_and_verify_no_df1_row_was_lost_but_allow_dups(df1, df2, **kwargs):
    extended_df1 = df1.copy()
    orig_num_of_rows = len(extended_df1)
    
    # # random to avoid collisions. maybe an overkill, but who cares.
    # random.seed(None) # use system time as seed
    # orig_row_i_random_column_name = f'orig_row_i_random_column_name_{get_random_str(10)}'
    # NOTE: formerly randomized this column name each time. cant remember why i did that...
    orig_row_i_random_column_name = 'orig_row_i_random_column_name_15714826758276298527'
    
    orig_row_is = list(range(orig_num_of_rows))
    extended_df1[orig_row_i_random_column_name] = orig_row_is
    
    merged_df = extended_df1.merge(df2, **kwargs)
    missing_row_is = set(orig_row_is) - set(merged_df[orig_row_i_random_column_name].unique())
    assert not missing_row_is, 'some of df1 rows were lost in the merge'
    merged_df.drop(columns=orig_row_i_random_column_name, inplace=True)
    
    return merged_df

def get_mask_of_df1_rows_that_inner_join_will_keep(df1, df2, **kwargs):
    merged_df = merge_preserving_df1_index_and_row_order(df1, df2, how='left', indicator=True, allow_multiple_common_column_names=True, **kwargs)
    return merged_df['_merge'] == 'both'

def inner_join_preserving_df1_row_order(df1, df2, **kwargs):
    return merge_preserving_df1_index_and_row_order(
        df1[get_mask_of_df1_rows_that_inner_join_will_keep(df1, df2, **kwargs)], 
        df2, 
        allow_multiple_common_column_names=True,
        **kwargs,
    )

def get_df1_masked_by_inner_join_with_df2_preserving_df1_row_order(df1, df2, **kwargs):
    return df1[get_mask_of_df1_rows_that_inner_join_will_keep(df1, df2, **kwargs)]

def merge_preserving_df1_row_order_and_verify_no_df1_row_was_duplicated_but_allow_losing_df1_rows(
        df1, df2, 
        allow_multiple_common_column_names=False, 
        **kwargs,
):
    raise RuntimeError('seems like you want inner_join_preserving_df1_row_order')


def merge_preserving_df1_index_and_row_order(
        df1, df2, 
        require_preserving_all_df1_rows=True, allow_multiple_common_column_names=False, 
        print_rows_that_were_lost_in_the_merge=True,
        verbose=True,
        allow_adding_suffixes=False,
        **kwargs,
):
    if isinstance(df2, list):
        merged_df = df1
        for curr_df2 in df2:
            merged_df = merge_preserving_df1_index_and_row_order(
                merged_df, curr_df2, 
                require_preserving_all_df1_rows=require_preserving_all_df1_rows, 
                allow_multiple_common_column_names=allow_multiple_common_column_names, 
                print_rows_that_were_lost_in_the_merge=print_rows_that_were_lost_in_the_merge,
                verbose=verbose,
                **kwargs,
            )
        return merged_df

    common_cols = sorted(set(df1.columns) & set(df2.columns))
    multiple_common_column_names = len(common_cols) > 1
    if verbose and multiple_common_column_names and (not allow_multiple_common_column_names) and ('on' not in kwargs) :
        print(f'NOTE: merging on multiple columns(!!!): {common_cols}')
    if (not allow_multiple_common_column_names) and ('on' not in kwargs) and ('how' in kwargs) and (kwargs['how'] != 'inner') and multiple_common_column_names:
        # thought it would be too annoying and unnecessary for inner join, so only checking in case of non-inner join...
        raise RuntimeError(f'multiple column names ({common_cols}) appear both in df1 and df2. if you really wanted that, specify allow_multiple_common_column_names=True')
    if ('on' in kwargs):
        on_cols = kwargs['on']
        if isinstance(on_cols, str):
            on_cols = [on_cols]
        on_cols_set = set(on_cols)
        common_column_names_set = set(common_cols)
        assert on_cols_set <= common_column_names_set, f'on: {on_cols}, common_cols: {common_cols}'
        if (not allow_adding_suffixes) and ('suffixes' not in kwargs):
            unexpected_common_cols = common_column_names_set - on_cols_set
            assert not unexpected_common_cols, f'unexpected_common_cols: {unexpected_common_cols}'


    extended_df1 = df1.copy()
    orig_num_of_rows = len(extended_df1)
    
    # # random to avoid collisions. maybe an overkill, but who cares.
    # random.seed(None) # use system time as seed
    # orig_row_i_random_column_name = f'orig_row_i_random_column_name_{get_random_str(10)}'
    # NOTE: formerly randomized this column name each time. cant remember why i did that...
    orig_row_i_random_column_name = 'orig_row_i_random_column_name_15714826758276298527'

    assert orig_row_i_random_column_name not in extended_df1.columns
    
    orig_row_is = list(range(orig_num_of_rows))
    extended_df1[orig_row_i_random_column_name] = orig_row_is
    assert extended_df1[orig_row_i_random_column_name].is_monotonic_increasing
    
    orig_indices = extended_df1.index.to_numpy() # use to_numpy() to avoid a weird behavior that uses the strings in some way...
    
    merged_df = extended_df1.merge(df2, **kwargs)
    duplicated_mask = merged_df[orig_row_i_random_column_name].duplicated(keep=False)
    assert not duplicated_mask.any(), f'the merge resulted in duplicates, so the index cant be preserved. duplicated rows:\n{merged_df[duplicated_mask]}'
    if require_preserving_all_df1_rows:
        no_rows_lost = len(merged_df) == orig_num_of_rows
        if print_rows_that_were_lost_in_the_merge and (not no_rows_lost):
            print('rows that were lost in the merge:')
            print(extended_df1[~extended_df1[orig_row_i_random_column_name].isin(merged_df[orig_row_i_random_column_name])])
        assert no_rows_lost, 'some of the rows were lost in the merge. use require_preserving_all_df1_rows=False if that isnt a problem'
    
    merged_df.sort_values(orig_row_i_random_column_name, inplace=True)
    if require_preserving_all_df1_rows:
        # if 'assert len(merged_df) == orig_num_of_rows' passed, this one should also pass, but adding this feels safer
        assert (merged_df[orig_row_i_random_column_name] == orig_row_is).all()

    assert merged_df[orig_row_i_random_column_name].is_monotonic_increasing
    # merged_df.index = orig_indices[merged_df[orig_row_i_random_column_name].isin(orig_row_is)] # 230623: commented out. i think it was a bug that didn't do anything because i never used require_preserving_all_df1_rows=False...
    merged_df.index = orig_indices[np.isin(orig_row_is, merged_df[orig_row_i_random_column_name])] 

    merged_df.drop(orig_row_i_random_column_name, axis=1, inplace=True)
    
    return merged_df

def get_df_sorted_by_column_vals(df, column_name, ordered_column_vals, left_join=False):
    sorted_df = merge_preserving_df1_index_and_row_order(pd.DataFrame({column_name: ordered_column_vals}), df, how='left' if left_join else 'inner')
    sorted_df_with_orig_column_order = sorted_df[list(df)]
    return sorted_df_with_orig_column_order

@contextlib.contextmanager
def timing_context_manager(block_name):
    start_time = time.perf_counter()
    try:
        yield
    finally:
        end_time = time.perf_counter()
        running_time_in_seconds = end_time - start_time
        print(f'block {block_name} running time (in seconds):\n{running_time_in_seconds:.3f}')


def rmtree_silent(dir_path):
    if os.path.exists(dir_path) and os.path.isdir(dir_path):
        shutil.rmtree(dir_path)


def remove_files_silently(file_paths):
    for file_path in file_paths:
        remove_file_silently(file_path)


def remove_file_silently(file_path):
    if os.path.isfile(file_path):
        os.remove(file_path)
    assert not os.path.isdir(file_path)


def silently_make_empty_dir(dir_path):
    rmtree_silent(dir_path)
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=False)

def run_subprocess(
        cmd_str_or_cmd_words, shell=True, verbose=True, 
        raise_exception_if_subproc_returned_non_zero=True, capture_output=True, log_file_path=None, append_to_log_file=True, **subprocess_run_kwargs,
):
    proc = subprocess.run(cmd_str_or_cmd_words, shell=shell, capture_output=capture_output, **subprocess_run_kwargs)
    return_code = proc.returncode
    stdout_str = proc.stdout.decode() if proc.stdout is not None else None
    stderr_str = proc.stderr.decode() if proc.stderr is not None else None

    if isinstance(cmd_str_or_cmd_words, str):
        cmd_str_or_cmd_words_as_iterable_of_str = [cmd_str_or_cmd_words]
    else:
        cmd_str_or_cmd_words_as_iterable_of_str = cmd_str_or_cmd_words

    log_str = (
        f'cmd_str:\n{" ".join(cmd_str_or_cmd_words_as_iterable_of_str)}\n'
        f'return_code:\n{return_code}\n'
    )
    log_str += f'stdout_str:\n{stdout_str}\n' if stdout_str else 'no stdout\n'
    log_str += f'stderr_str:\n{stderr_str}\n' if stderr_str else 'no stderr\n'

    if verbose:
        print(log_str)

    if log_file_path:
        if append_to_log_file:
            append_to_text_file(log_file_path, log_str)
        else:
            write_text_file(log_file_path, log_str)

    if return_code:
        print(log_str)

        if raise_exception_if_subproc_returned_non_zero:
            raise RuntimeError(proc)

    return proc

def run_cmd_and_get_stdout_and_stderr(cmd_line_words, raise_exception_if_subproc_returned_non_zero=True, also_return_return_code=False, verbose=False, **subprocess_run_kwargs):
    assert raise_exception_if_subproc_returned_non_zero or also_return_return_code
    cmd_as_str = ' '.join(cmd_line_words)
    if verbose:
        print('\nrunning cmd: ' + cmd_as_str)

    subproc = subprocess.run(cmd_line_words, check=False, capture_output=True, **subprocess_run_kwargs)
    subproc_stdout = subproc.stdout.decode()
    subproc_stderr = subproc.stderr.decode()
    if raise_exception_if_subproc_returned_non_zero and subproc.returncode:
        print(f'cmd that returned non-zero:\n{cmd_as_str}')
        print(f'cmd result stderr:\n{subproc_stderr}')
        print(f'cmd result stdout:\n{subproc_stdout}')
        print(f'cmd result ret_code: {subproc.returncode}')
        raise subprocess.SubprocessError(f'Process failed with return code {subproc.returncode}')
    if also_return_return_code:
        return subproc_stdout, subproc_stderr, subproc.returncode
    return subproc_stdout, subproc_stderr


def run_cmd_and_write_stdout_and_stderr_to_files(cmd_line_words, stdout_file_path=None, stderr_file_path=None, raise_exception_if_subproc_returned_non_zero=True, verbose=True, **subprocess_run_kwargs):
    if verbose:
        print('\nrunning cmd: ' + ' '.join(cmd_line_words))

    subproc_stdout, subproc_stderr = run_cmd_and_get_stdout_and_stderr(
        cmd_line_words=cmd_line_words,
        raise_exception_if_subproc_returned_non_zero=raise_exception_if_subproc_returned_non_zero,
        **subprocess_run_kwargs,
    )

    if stdout_file_path:
        with open(stdout_file_path, 'w') as f:
            f.write(subproc_stdout)
    if stderr_file_path:
        with open(stderr_file_path, 'w') as f:
            f.write(subproc_stderr)


def run_cmd_and_assert_stdout_and_stderr_are_empty(cmd_line_words, **subprocess_run_kwargs):
    subproc_stdout, subproc_stderr, subproc_ret_code = run_cmd_and_get_stdout_and_stderr(
        cmd_line_words, raise_exception_if_subproc_returned_non_zero=False,
        also_return_return_code=True, verbose=True, **subprocess_run_kwargs)
    if (subproc_ret_code != 0) or subproc_stdout or subproc_stderr:
        if subproc_stderr:
            print(f'cmd result stderr:\n{subproc_stderr}')
        if subproc_stdout:
            print(f'cmd result stdout:\n{subproc_stdout}')
        if subproc_ret_code != 0:
            print(f'cmd result ret_code: {subproc_ret_code}')

        assert False


def run_cmd_and_check_ret_code_and_return_stdout(
        cmd_line_words, stdout_file_path=None, raise_exception_if_stderr_isnt_empty=False, print_stdout_if_verbose=True, verbose=True, **subprocess_run_kwargs,
):
    if verbose:
        print('\nrunning cmd: ' + ' '.join(cmd_line_words))
    # print('\nrunning cmd: ' + ' '.join(cmd_line_words))
    # exit()

    subproc_stdout, subproc_stderr = run_cmd_and_get_stdout_and_stderr(cmd_line_words, **subprocess_run_kwargs)

    if stdout_file_path is not None:
        with open(stdout_file_path, 'w') as f:
            f.write(subproc_stdout)
    else:
        if verbose and print_stdout_if_verbose:
            if subproc_stdout:
                print('cmd result stdout:\n' + subproc_stdout)
            else:
                print('cmd result stdout: <empty str>\n')

    if subproc_stderr:
        if raise_exception_if_stderr_isnt_empty:
            raise RuntimeError('cmd result stderr:\n' + subproc_stderr)
        if verbose:
            print('cmd result stderr:\n' + subproc_stderr)
    else:
        if verbose:
            print('cmd result stderr: <empty str>\n')
    return subproc_stdout

def download_file_from_ftp_server(ftp_url, output_file_path, dont_download_if_file_exists=False, verbose=False):
    if (not dont_download_if_file_exists) or (not os.path.isfile(output_file_path)):
        run_cmd_and_check_ret_code_and_return_stdout(['curl', '-o', output_file_path, ftp_url], verbose=verbose)

def get_file_hash_file_path(file_path):
    # 200522: --------ATTENTION--------
    # Emabarrassingly, until today, this function was returning the same file path for two different files with the same file name (but in different dirs). An epic fail.
    # This comment is here to make sure you don't change the code to what it was earlier, and thus repeat this epic fail. Thank you.
    return file_path + HASH_FILE_NAME_SUFFIX


def get_file_path_with_func_name_and_input_args_hash(func_name, args_hash, suffix):
    func_dir_path = os.path.join(INTERNAL_FILES_TO_SKIP_REDUNDANT_CALCULATIONS_DIR_PATH, func_name)
    assert type(args_hash) == str
    args_hash_start = args_hash[:8]
    # args_hash_end = args_hash[8:]
    curr_hash_dir_path = os.path.join(func_dir_path, *textwrap.wrap(args_hash_start, 2))
    pathlib.Path(curr_hash_dir_path).mkdir(parents=True, exist_ok=True)
    # print(f'curr_hash_dir_path: {curr_hash_dir_path}')
    # exit()
    return os.path.join(curr_hash_dir_path, f'{args_hash}{suffix}')


def get_func_input_args_file_path(func_name, args_hash):
    return get_file_path_with_func_name_and_input_args_hash(func_name, args_hash, FUNC_INPUT_ARGS_IN_PREV_CALL_FILE_NAME_SUFFIX)


def remove_files_and_their_file_hash_files(file_paths):
    for file_path in file_paths:
        file_hash_file_path = get_file_hash_file_path(file_path)
        remove_file_silently(file_hash_file_path)
        remove_file_silently(file_path)

def get_file_last_modification_time(file_path):
    assert os.path.isfile(file_path)
    return os.path.getmtime(file_path)

def get_file_last_modification_time_str_repr(file_path):
    return datetime.datetime.fromtimestamp(get_file_last_modification_time(file_path)).strftime('%d/%m/%Y:%H:%M')

def get_file_or_dir_last_modification_time(file_or_dir_path):
    if os.path.isfile(file_or_dir_path):
        return get_file_last_modification_time(file_or_dir_path)

    assert os.path.isdir(file_or_dir_path)
    # thanks https://stackoverflow.com/questions/29685069/get-the-last-modified-date-of-a-directory-including-subdirectories-using-pytho/29685234#29685234
    return max(os.path.getmtime(root) for root, _, _ in os.walk(file_or_dir_path))


def do_output_files_exist(output_file_paths):
    if not output_file_paths:
        return True
    hashes = set()
    for file_path in output_file_paths:
        file_hash_file_path = get_file_hash_file_path(file_path)
        if (not os.path.exists(file_path)) or (not os.path.isfile(file_hash_file_path)) or (os.path.getmtime(file_hash_file_path) < os.path.getmtime(file_path)):
            return False
        hashes.add(read_text_file(file_hash_file_path))
    return len(hashes) == 1

def get_object_deterministic_str_repr_lines_by_recursive_sort(obj):
    if isinstance(obj, (int, float, complex, bool, str, type(None))):
        return [str(obj)]
    if isinstance(obj, (set, frozenset)):
        return [str(sorted(obj))]

    deterministic_str_repr_lines = []
    if isinstance(obj, (list, tuple)):
        for x in obj:
            deterministic_str_repr_lines += get_object_deterministic_str_repr_lines_by_recursive_sort(x)
    elif isinstance(obj, dict):
        # print('\n\nprint each key')
        # for k in obj:
        #     print(k)
        for k in sorted(obj):
            deterministic_str_repr_lines.append(f'{k}:')
            deterministic_str_repr_lines += get_object_deterministic_str_repr_lines_by_recursive_sort(obj[k])
    else:
        raise NotImplementedError(f'obj is of type that is not supported yet (type(obj): {type(obj)})')
    return deterministic_str_repr_lines

def get_input_file_path_arg_text_repr_lines_and_file_last_modification_time_and_name(input_arg_name, input_arg_val):
    if not isinstance(input_arg_val, str):
        print(f'input_arg_val should be a str. ({input_arg_name}, {input_arg_val})')
    assert isinstance(input_arg_val, str)
    input_file_path = input_arg_val

    # I don't want to receive a directory. Too much added complexity, and it isn't really needed.
    assert not os.path.isdir(input_file_path)
    assert os.path.isfile(input_file_path)


    input_file_hash_file_path = get_file_hash_file_path(input_file_path)
    input_file_path_last_modification_time = get_file_last_modification_time(input_file_path)
    if os.path.isfile(input_file_hash_file_path) and (os.path.getmtime(input_file_hash_file_path) >= input_file_path_last_modification_time):
        input_file_hash = read_text_file(input_file_hash_file_path)
    else:
        # useful_debug_str = (f'input_file_path: {input_file_path}\n'
        #                     f'os.path.getmtime(input_file_hash_file_path): {os.path.getmtime(input_file_hash_file_path)}\n'
        #                     f'input_file_path_last_modification_time:      {input_file_path_last_modification_time}\n')

        # print(useful_debug_str)
        # # if os.path.isfile(input_file_hash_file_path):
        # #     raise RuntimeError(f'seems like someone overwrote a file without updating its hash file. I think this shouldnt happen silently...\n'
        # #                        f'{useful_debug_str}')
        # print('######################## recalculating input_file_hash')
        input_file_hash = get_file_sha256_in_hex(input_file_path)
        write_text_file(input_file_hash_file_path, input_file_hash)

    input_file_name = os.path.basename(input_file_path)
    input_arg_text_repr_lines = [f'{input_arg_name}: {input_file_path}: {input_file_hash}']
    return (input_arg_text_repr_lines, input_file_path_last_modification_time, input_file_name)

def execute_if_output_doesnt_exist_already(func):
    func_name = func.__name__

    def new_func(*args, **kwargs):
        # print(f'skipping execution of ############################################ func_name: {func_name}')

        assert args == ()  # this is just for my convenience (and laziness) - this way i don't need to use reflection to retrieve the names of the arguments...
        output_args_dict = {arg_name: arg_val for arg_name, arg_val in kwargs.items() if arg_name.startswith(OUTPUT_FILE_PATH_ARG_PREFIX)}
        output_file_for_caching_only_path = kwargs.get(OUTPUT_FILE_PATH_FOR_CACHING_ONLY_ARG)
        output_file_paths = [output_path for output_path in output_args_dict.values() if output_path is not None]
        for output_file_path in output_file_paths:
            assert isinstance(output_file_path, str)
            # I am afraid of removing important stuff...
            # assert not os.path.isabs(output_file_path)
            assert not os.path.isdir(output_file_path)
        input_args_dict = {arg_name: arg_val for arg_name, arg_val in kwargs.items() if not arg_name.startswith(OUTPUT_FILE_PATH_ARG_PREFIX)}

        # add the hash of the inputs without their hashes to the file names. nah. you use dirs instead.

        input_files_names = []
        func_curr_input_args_lines = []
        input_files_last_modification_times = []
        for input_arg_name, input_arg_val in sorted(input_args_dict.items()):
            if not isinstance(input_arg_val, PYTHON_BUILTIN_TYPES):
                print(f'input_arg_name: {input_arg_name}')
                print(f'input_arg_val: {input_arg_val}')
                print(f'type(input_arg_val): {type(input_arg_val)}')
                assert isinstance(input_arg_val, PYTHON_BUILTIN_TYPES)
            if input_arg_name.startswith(INPUT_FILE_PATH_ARG_PREFIX):
                input_arg_text_repr_lines, input_file_path_last_modification_time, input_file_name = (
                    get_input_file_path_arg_text_repr_lines_and_file_last_modification_time_and_name(input_arg_name, input_arg_val))

                input_files_names.append(input_file_name)
                input_files_last_modification_times.append(input_file_path_last_modification_time)
            elif input_arg_name.startswith(MULTIPLE_INPUT_FILE_PATH_ARG_PREFIX):
                if (not isinstance(input_arg_val, (list, tuple, set, frozenset))) or (not (len(input_arg_val) == len(set(input_arg_val)))):
                    print(f'input_arg_name: {input_arg_name}')
                    print(f'input_arg_val: {input_arg_val}')
                    print(f'type(input_arg_val): {type(input_arg_val)}')
                    assert isinstance(input_arg_val, (list, tuple, set, frozenset))
                    assert len(input_arg_val) == len(set(input_arg_val))
                input_arg_text_repr_lines = []
                for i, input_file_path in enumerate(sorted(input_arg_val)):
                    curr_input_arg_text_repr_lines, input_file_path_last_modification_time, input_file_name = (
                        get_input_file_path_arg_text_repr_lines_and_file_last_modification_time_and_name(f'{input_arg_name}{i}', input_file_path))

                    input_files_names.append(input_file_name)
                    input_files_last_modification_times.append(input_file_path_last_modification_time)
                    input_arg_text_repr_lines += curr_input_arg_text_repr_lines
            else:
                input_arg_val_text_repr_lines = get_object_deterministic_str_repr_lines_by_recursive_sort(input_arg_val)
                input_arg_text_repr_lines = [f'{input_arg_name}:'] + input_arg_val_text_repr_lines

            func_curr_input_args_lines += input_arg_text_repr_lines

        func_curr_input_args = '\n'.join(func_curr_input_args_lines)
        func_curr_input_args_hash = get_str_sha256_in_hex(func_curr_input_args)
        func_curr_input_args_file_path = get_func_input_args_file_path(func_name, func_curr_input_args_hash)
        need_to_execute_func = (
            # I suspect that I could remove func_curr_input_args_file_path and still make the caching work, ***BUT*** I think that would be a bad idea. Currently, a convenient
            # way to make my code execute the cached function anyway is to remove all func_curr_input_args_file_path in internal___hashes_etc_to_skip_redundant_calculations.
            # Removing func_curr_input_args_file_path would make it harder to make the function execute anyway, I think.
                (not os.path.isfile(func_curr_input_args_file_path)) or
                any(not os.path.isfile(x) for x in output_file_paths)
        )
        if need_to_execute_func:
            # print(f'skipping execution of ############################################ need_to_execute_func1: {need_to_execute_func} {func_name}')
            # print(f'skipping execution of ############################################ (not os.path.isfile(func_curr_input_args_file_path)): {(not os.path.isfile(func_curr_input_args_file_path))}')
            # print(f'skipping execution of ############################################ (not do_output_files_and_potentially_valid_hash_files_exist(output_file_paths)): {(not do_output_files_and_potentially_valid_hash_files_exist(output_file_paths))}')
            pass
        # print(f'skipping execution of ############################################')
        # print(f'skipping execution of ############################################ {func_name}')
        # print(f'skipping execution of ############################################ (not os.path.isfile(func_curr_input_args_file_path)): {(not os.path.isfile(func_curr_input_args_file_path))}')

        if (not need_to_execute_func) and output_file_paths:
            output_files_last_modification_times = [get_file_or_dir_last_modification_time(x) for x in output_file_paths]

            if 1: # 220104: added this only now. so if some bug arises, it might be the culprit.
                if not need_to_execute_func:
                    if max(output_files_last_modification_times) > get_file_or_dir_last_modification_time(func_curr_input_args_file_path):
                        need_to_execute_func = True
                        # print(f'skipping execution of ############################################ need_to_execute_func3: {need_to_execute_func}')

        # print(f'skipping execution of ############################################ need_to_execute_func4: {need_to_execute_func}')

        if need_to_execute_func:
            remove_file_silently(func_curr_input_args_file_path)
            if output_file_for_caching_only_path in output_file_paths:
                remove_file_silently(output_file_for_caching_only_path)


            func_ret = func(*args, **kwargs)
            assert func_ret is None

            write_text_file(func_curr_input_args_file_path, func_curr_input_args)
            for output_file_path in output_file_paths:
                if output_file_path == output_file_for_caching_only_path:
                    # Why is this needed? For things like make_blast_nucleotide_db and execute_blast_nucleotide - we don't know exactly what makeblastdb creates, so we create this
                    # file. then, when we call execute_blast_nucleotide, the arguments include the path of the blast db, and we also give this file's path as an argument.
                    assert not os.path.isfile(output_file_path)
                    # this was problematic (IIUC) when less than a second passed between the two calls! use a random number instead!
                    # write_text_file(output_file_path, time.strftime('%Y_%m_%d_%H_%M_%S') + str(time.time()))
                    write_text_file(output_file_path, str(random.random()))
                elif not os.path.isfile(output_file_path):
                    raise RuntimeError(f'Func {func_name} did not create the output file {output_file_path}. Shame.')

        else:
            comma_separated_input_file_names_as_str_literals = ','.join(f"'{file_name}'" for file_name in input_files_names)
            print(f'skipping execution of {func_name}({comma_separated_input_file_names_as_str_literals})')

    new_func.__name__ = f'execute_if_output_doesnt_exist_already___{func_name}'
    return new_func

def unzip_gz_file(gz_file_path, unzipped_file_path=None, do_nothing_if_unzipped_file_already_exists=False):
    # NOTE: 
    # NOTE: you probably want to use zcat instead. gz files are stored as gz because storing as unzipped is wasteful.
    # NOTE: 
    
    if unzipped_file_path is None:
        unzipped_file_path = gz_file_path.partition('.gz')[0]

    if do_nothing_if_unzipped_file_already_exists and os.path.isfile(unzipped_file_path):
        return

    import gzip
    import shutil
    with gzip.open(gz_file_path, 'rb') as f_in:
        with open(unzipped_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

@execute_if_output_doesnt_exist_already
def cached_unzip_gz_file(
        input_file_path_gz,
        output_file_path_unzipped,
):
    # if gunzip is too slow, consider using pigz?
    # https://unix.stackexchange.com/questions/623881/how-to-gzip-100-gb-files-faster-with-high-compression/623883#623883
    # https://superuser.com/questions/599329/why-is-gzip-slow-despite-cpu-and-hard-drive-performance-not-being-maxed-out
    assert input_file_path_gz.endswith('.gz')
    assert output_file_path_unzipped == input_file_path_gz[:-3]
    run_cmd_and_check_ret_code_and_return_stdout(['gunzip', '-k', input_file_path_gz], verbose=False)

def unzip_gz_file_with_caching(
        input_file_path_gz,
        output_file_path_unzipped,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_unzip_gz_file(
        input_file_path_gz=input_file_path_gz,
        output_file_path_unzipped=output_file_path_unzipped,
    )

@contextlib.contextmanager
def redirect_stdout_and_verify_that_each_line_starts_with_skipping_execution_of():
    orig_stdout = sys.stdout
    captured_output_stringio = io.StringIO()
    sys.stdout = captured_output_stringio
    try:
        yield
    finally:
        sys.stdout = orig_stdout
        captured_output = captured_output_stringio.getvalue()
        captured_output_lines = captured_output.splitlines()
        for line in captured_output_lines:
            if not line.startswith('skipping execution of '):
                raise RuntimeError(f'Expected all output lines to start with "skipping execution of ", but this line does not:\n{line}')

def get_func_ret_val_and_stdout(f, *args, **kwargs):
    orig_stdout = sys.stdout
    captured_output_stringio = io.StringIO()
    sys.stdout = captured_output_stringio
    ret_val = f(*args, **kwargs)
    sys.stdout = orig_stdout
    captured_output = captured_output_stringio.getvalue()
    return ret_val, captured_output

def feature_scale_pd_series(pd_series):
    series_max = pd_series.max()
    series_min = pd_series.min()
    series_range_size = series_max - series_min
    if series_range_size == 0:
        return pd.Series(np.ones(len(pd_series)))
    return pd_series.apply(lambda x: (x - series_min) / series_range_size)


def np_array_contains_only_zeroes(arr):
    return not arr.any()


def get_hamming_dist_between_same_len_strs(str1, str2):
    hamming_dist = 0
    for char1, char2 in zip(str1, str2):
        if char1 != char2:
            hamming_dist += 1
    return hamming_dist


def is_hamming_dist_at_most_x(str1, str2, max_hamming_dist):
    hamming_dist = 0
    for char1, char2 in zip(str1, str2):
        if char1 != char2:
            hamming_dist += 1
            if hamming_dist > max_hamming_dist:
                return False
    return True


def get_indices_of_mismatches_of_same_len_strs(str1, str2):
    assert len(str1) == len(str2)
    indices_of_mismatches = []
    for i, (char1, char2) in enumerate(zip(str1, str2)):
        if char1 != char2:
            indices_of_mismatches.append(i)

    # print(str1)
    # print(str2)
    # print(indices_of_mismatches)
    return indices_of_mismatches


class MyTestCase(unittest.TestCase):
    @contextlib.contextmanager
    def assertPrints(self, expected_output):
        orig_stdout = sys.stdout
        captured_output = io.StringIO()
        sys.stdout = captured_output
        try:
            yield
        finally:
            sys.stdout = orig_stdout
            self.assertEqual(captured_output.getvalue(), expected_output)

def replace_file_name_in_path(file_path, new_file_name):
    file_dir_path, file_name = os.path.split(file_path)
    return os.path.join(file_dir_path, new_file_name)


def does_str_contain_any_whitespace(s):
    splitted_s = s.split()
    num_of_splitted_parts = len(splitted_s)
    assert num_of_splitted_parts >= 1
    return num_of_splitted_parts > 1


def read_object_from_pickle_file(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)

def write_object_to_pickle_file(file_path, obj):
    with open(file_path, 'wb') as f:
        pickle.dump(obj, f, protocol=4)

def strip_special_chars_on_edges_and_replace_others_with_underscores(orig_str):
    str_with_special_replaced_with_underscores = ''.join(char if char.isalnum() else '_' for char in orig_str)
    return str_with_special_replaced_with_underscores.strip('_')

def replace_special_chars_with_underscores(orig_str):
    return ''.join(char if char.isalnum() else '_' for char in orig_str)

def str_contains_only_alnum_and_underscores(s):
    return (s.replace('_', '')).isalnum()

def filter_df_by_column_name_and_value_or_values_or_range(df, column_name, value_or_values_or_range):
    # if df.empty:
    #     return df.copy()

    if isinstance(value_or_values_or_range, Range):
        # if column_name not in list(df):
        #     print('df')
        #     print(df)
        rang = value_or_values_or_range
        return df[(df[column_name] >= rang.start) & (df[column_name] <= rang.end)]
    elif isinstance(value_or_values_or_range, tuple):
        values = value_or_values_or_range
        any_value_filter = pd.Series([False] * len(df), index=df.index)
        for value in values:
            any_value_filter = any_value_filter | (df[column_name] == value)
        return df[any_value_filter]
    elif (value_or_values_or_range is None) or pd.isna(value_or_values_or_range):
        return df[df[column_name].isna()]
    else:
        assert isinstance(value_or_values_or_range, (int, float, complex, bool, str))
        # print(f'df: {df}')
        # print(f'value_or_values_or_range: {value_or_values_or_range}')
        value = value_or_values_or_range
        return df[df[column_name] == value]

@execute_if_output_doesnt_exist_already
def cached_filter_csv_by_column_name_and_values_range(
        input_file_path_csv,
        column_name,
        values_range,
        output_file_path_filtered_csv,
        csv_missing_column_names,
):
    df = pd.read_csv(input_file_path_csv, sep='\t', names=csv_missing_column_names)
    filtered_df = filter_df_by_column_name_and_value_or_values_or_range(df, column_name, Range(*values_range))
    filtered_df.to_csv(output_file_path_filtered_csv, sep='\t', index=False)

def filter_csv_by_column_name_and_values_range(
        input_file_path_csv,
        column_name,
        values_range,
        output_file_path_filtered_csv,
        csv_missing_column_names=None,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_filter_csv_by_column_name_and_values_range(
        input_file_path_csv=input_file_path_csv,
        column_name=column_name,
        values_range=values_range,
        output_file_path_filtered_csv=output_file_path_filtered_csv,
        csv_missing_column_names=csv_missing_column_names,
    )

def filter_df_to_keep_only_rows_with_extreme_val_for_each_group(
        df,
        group_by_column_name,
        extreme_val_column_name,
        max_or_min,
):
    assert max_or_min in ('max', 'min')

    filtered_df = df.groupby(group_by_column_name, as_index=False, sort=False).apply(
        lambda group_df: group_df.loc[(group_df[extreme_val_column_name].idxmax() if max_or_min == 'max' else group_df[extreme_val_column_name].idxmin()),:]
    )
    return filtered_df

@execute_if_output_doesnt_exist_already
def cached_filter_csv_to_keep_only_rows_with_extreme_val_for_each_group(
        input_file_path_csv,
        group_by_column_name,
        extreme_val_column_name,
        max_or_min,
        output_file_path_filtered_csv,
        csv_missing_column_names,
):
    df = pd.read_csv(input_file_path_csv, sep='\t', names=csv_missing_column_names)

    filtered_df = filter_df_to_keep_only_rows_with_extreme_val_for_each_group(
        df=df,
        group_by_column_name=group_by_column_name,
        extreme_val_column_name=extreme_val_column_name,
        max_or_min=max_or_min,
    )

    filtered_df.to_csv(output_file_path_filtered_csv, sep='\t', index=False)

def filter_csv_to_keep_only_rows_with_extreme_val_for_each_group(
        input_file_path_csv,
        group_by_column_name,
        extreme_val_column_name,
        max_or_min,
        output_file_path_filtered_csv,
        csv_missing_column_names=None,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_filter_csv_to_keep_only_rows_with_extreme_val_for_each_group(
        input_file_path_csv=input_file_path_csv,
        group_by_column_name=group_by_column_name,
        extreme_val_column_name=extreme_val_column_name,
        max_or_min=max_or_min,
        output_file_path_filtered_csv=output_file_path_filtered_csv,
        csv_missing_column_names=csv_missing_column_names,
    )

def get_filtered_df_by_fixed_column_names_and_values_or_ranges(df, fixed_column_names_and_values_or_ranges):
    fixed_columns_str_reprs = []
    for column_name, value_or_values_or_range_single_elem_list in fixed_column_names_and_values_or_ranges:
        assert len(value_or_values_or_range_single_elem_list) == 1
        value_or_values_or_range = value_or_values_or_range_single_elem_list[0]
        df = filter_df_by_column_name_and_value_or_values_or_range(df, column_name, value_or_values_or_range)
        fixed_columns_str_reprs.append(f'{column_name}: {value_or_values_or_range}')
    return (df, fixed_columns_str_reprs)


def load_np_array_from_single_array_npz_file(npz_file_path):
    loaded_npz = np.load(npz_file_path)
    assert len(loaded_npz.files) == 1
    return loaded_npz[loaded_npz.files[0]]

def write_word_counts_and_freqs_df_to_csv(word_counts_and_freqs_df, csv_file_path):
    if word_counts_and_freqs_df.empty:
       word_counts_and_freqs_df = pd.DataFrame(columns=WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES)

    assert set(word_counts_and_freqs_df) == WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES_SET
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        word_counts_and_freqs_df.to_csv(csv_file_path, sep='\t', index=False)

def read_word_counts_and_freqs_from_csv(csv_file_path):
    word_counts_and_freqs_df = pd.read_csv(csv_file_path, sep='\t', dtype=WORD_COUNTS_AND_FREQS_COLUMN_NAME_TO_DTYPE)
    assert set(word_counts_and_freqs_df) == WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES_SET
    return word_counts_and_freqs_df

def add_freq_column_and_write_word_counts_and_freqs_to_csv_and_total_word_count_to_text_file(word_counts_and_freqs_df, word_counts_and_freqs_csv_file_path,
                                                                                             total_word_count_txt_file_path):
    total_count = word_counts_and_freqs_df['count'].sum()
    if total_count == 0:
        word_counts_and_freqs_df.loc[:, 'freq'] = np.nan
    else:
        word_counts_and_freqs_df.loc[:, 'freq'] = word_counts_and_freqs_df['count'] / total_count

    write_word_counts_and_freqs_df_to_csv(word_counts_and_freqs_df, word_counts_and_freqs_csv_file_path)
    write_text_file(total_word_count_txt_file_path, str(total_count))

def write_word_counts_and_freqs_to_csv_and_total_word_count_to_text_file(text_series, word_counts_and_freqs_csv_file_path, total_word_count_txt_file_path):
    # print(f'text_series:\n{text_series}')
    if text_series.empty:
        word_counts_and_freqs_df = pd.DataFrame(columns=WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES)
    else:
        # Counting the number of rows for each word would be much less efficient, i think, because i don't know how to do that with the nice split(' ', expand=True).
        # So let's just do this. It seems to me like it won't change much (there is a difference only if a word appears more than once in the same row).
        word_counts_and_freqs_df = text_series.str.lower().str.split(' ', expand=True).stack().value_counts().reset_index(name='count').rename(columns={'index': 'word'})

    add_freq_column_and_write_word_counts_and_freqs_to_csv_and_total_word_count_to_text_file(
        word_counts_and_freqs_df=word_counts_and_freqs_df,
        word_counts_and_freqs_csv_file_path=word_counts_and_freqs_csv_file_path,
        total_word_count_txt_file_path=total_word_count_txt_file_path,
    )

def merge_word_counts_and_freqs_dfs_and_write_to_csv_and_write_total_word_count_to_text_file(word_counts_and_freqs_dfs_list, merged_word_counts_and_freqs_csv_file_path,
                                                                                             merged_total_word_count_txt_file_path):
    if not word_counts_and_freqs_dfs_list:
        merged_word_counts_and_freqs_df = pd.DataFrame(columns=WORD_COUNTS_AND_FREQS_DF_COLUMN_NAMES)
    else:
        merged_word_counts_and_freqs_df = pd.concat(word_counts_and_freqs_dfs_list)
        merged_word_counts_and_freqs_df = merged_word_counts_and_freqs_df.groupby('word', as_index=False, sort=False).agg({'count': np.sum})

    add_freq_column_and_write_word_counts_and_freqs_to_csv_and_total_word_count_to_text_file(
        word_counts_and_freqs_df=merged_word_counts_and_freqs_df,
        word_counts_and_freqs_csv_file_path=merged_word_counts_and_freqs_csv_file_path,
        total_word_count_txt_file_path=merged_total_word_count_txt_file_path,
    )

def get_num_of_strs_in_which_substr_appears(strs, substr):
    return sum((substr in s) for s in strs)

def get_paths_of_files_and_dirs_in_dir(dir_path, assert_only_files=False):
    paths = [os.path.join(dir_path, filename) for filename in os.listdir(dir_path)]
    if assert_only_files:
        assert all(os.path.isfile(path) for path in paths)
    return paths

def get_df_filtered_by_allowed_column_values(df, column_name, allowed_values):
    filter_column = pd.Series(np.zeros(len(df)), dtype='bool', index=df.index)
    for allowed_val in allowed_values:
        filter_column.loc[df[column_name] == allowed_val] = True
    return df[filter_column]

def get_df_filtered_by_unallowed_column_values(df, column_name, unallowed_values):
    filter_column = pd.Series(np.ones(len(df)), dtype='bool', index=df.index)
    for unallowed_val in unallowed_values:
        filter_column.loc[df[column_name] == unallowed_val] = False
    return df[filter_column]


@execute_if_output_doesnt_exist_already
def cached_write_filtered_df_to_csv(
        input_file_path_csv,
        column_name,
        column_value,
        output_file_path_filtered_csv,
        csv_separator,
):
    df = pd.read_csv(input_file_path_csv, sep=csv_separator)
    filtered_df = df[df[column_name] == column_value]
    filtered_df.to_csv(output_file_path_filtered_csv, sep=csv_separator, index=False)
    # with warnings.catch_warnings():
    #     warnings.simplefilter('ignore', DeprecationWarning)
    #     filtered_df.to_csv(output_file_path_filtered_csv, sep=csv_separator, index=False)

def write_filtered_df_to_csv(
        input_file_path_csv,
        column_name,
        column_value,
        output_file_path_filtered_csv,
        csv_separator='\t',
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_filtered_df_to_csv(
        input_file_path_csv=input_file_path_csv,
        column_name=column_name,
        column_value=column_value,
        output_file_path_filtered_csv=output_file_path_filtered_csv,
        csv_separator=csv_separator,
    )

def is_point_touching_interval(point, interval):
    return interval[0] <= point <= interval[1]

@execute_if_output_doesnt_exist_already
def cached_write_index_to_is_covered_by_any_interval_array(
        intervals_csv_file_path,
        array_len,
        output_file_path_index_to_is_covered_by_any_interval_array_npz,
):
    intervals_df = pd.read_csv(intervals_csv_file_path, sep='\t')
    index_to_is_covered_by_any_interval_array = np.zeros(array_len)
    for _, row in intervals_df.iterrows():
        index_to_is_covered_by_any_interval_array[(row['start'] - 1):row['end']] = 1

    np.savez_compressed(output_file_path_index_to_is_covered_by_any_interval_array_npz,
                        index_to_is_covered_by_any_interval_array=index_to_is_covered_by_any_interval_array)

def write_index_to_is_covered_by_any_interval_array(
        intervals_csv_file_path,
        array_len,
        output_file_path_index_to_is_covered_by_any_interval_array_npz,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_index_to_is_covered_by_any_interval_array(
        intervals_csv_file_path=intervals_csv_file_path,
        array_len=array_len,
        output_file_path_index_to_is_covered_by_any_interval_array_npz=output_file_path_index_to_is_covered_by_any_interval_array_npz,
    )

@execute_if_output_doesnt_exist_already
def cached_write_index_to_num_of_intervals_covering_it_array(
        intervals_csv_file_path,
        array_len,
        output_file_path_index_to_num_of_intervals_covering_it_array_npz,
):
    intervals_df = pd.read_csv(intervals_csv_file_path, sep='\t')
    index_to_num_of_intervals_covering_it_array = np.zeros(array_len)
    for _, row in intervals_df.iterrows():
        index_to_num_of_intervals_covering_it_array[(row['start'] - 1):row['end']] += 1

    np.savez_compressed(output_file_path_index_to_num_of_intervals_covering_it_array_npz,
                        index_to_num_of_intervals_covering_it_array=index_to_num_of_intervals_covering_it_array)

def write_index_to_num_of_intervals_covering_it_array(
        intervals_csv_file_path,
        array_len,
        output_file_path_index_to_num_of_intervals_covering_it_array_npz,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_index_to_num_of_intervals_covering_it_array(
        intervals_csv_file_path=intervals_csv_file_path,
        array_len=array_len,
        output_file_path_index_to_num_of_intervals_covering_it_array_npz=output_file_path_index_to_num_of_intervals_covering_it_array_npz,
    )

def get_max_num_of_overlapping_intervals(intervals_df):
    assert (intervals_df.iloc[:,0] <= intervals_df.iloc[:,1]).all()
    start_to_count = intervals_df.iloc[:,0].value_counts().to_dict()
    end_to_count = intervals_df.iloc[:,1].value_counts().to_dict()
    starts_and_ends_sorted = sorted(set(start_to_count) | set(end_to_count))
    curr_num_of_overlapping_intervals = 0
    max_num_of_overlapping_intervals = 0
    for pos in starts_and_ends_sorted:
        if pos in start_to_count:
            curr_num_of_overlapping_intervals += start_to_count[pos]
        if pos in end_to_count:
            max_num_of_overlapping_intervals = max(max_num_of_overlapping_intervals, curr_num_of_overlapping_intervals)
            curr_num_of_overlapping_intervals -= end_to_count[pos]
    return max(max_num_of_overlapping_intervals, curr_num_of_overlapping_intervals)

assert get_max_num_of_overlapping_intervals(pd.DataFrame([(1,1),(1,24),(3,5)])) == 2
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(1,1),(1,24),(3,5),(4,8)])) == 3
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(1,1),(1,4),(3,5),(4,8)])) == 3
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(1,1),(1,3),(3,5),(4,8)])) == 2
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(1,1),(1,3),(3,5),(3,8)])) == 3
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(1,3),(1,3),(3,5),(3,8)])) == 4
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(1,3),(1,2),(3,5),(3,8)])) == 3
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(2,2),(1,2),(3,5),(7,8)])) == 2
assert get_max_num_of_overlapping_intervals(pd.DataFrame([(2,2),(1,1),(3,5),(7,8)])) == 1

def is_interval(interval):
    return (type(interval) == tuple) and (len(interval) == 2) and (interval[0] <= interval[1])

def is_interval1_contained_in_interval2(interval1, interval2):
    assert is_interval(interval1)
    assert is_interval(interval2)
    return (interval2[0] <= interval1[0]) and (interval1[1] <= interval2[1])

def is_iterable_of_intervals(intervals):
    return all(is_interval(interval) for interval in intervals)

def old_slow_very_naive_get_merged_intervals(intervals):
    assert is_iterable_of_intervals(intervals)

    merged_intervals = set(intervals)
    while True:
        merged_something_this_cycle = False

        for interval1, interval2 in itertools.permutations(merged_intervals, 2):
            if not (
                (interval1[1] < interval2[0]) or
                (interval2[1] < interval1[0])
            ):
                both_intervals_edges = interval1 + interval2
                assert len(both_intervals_edges) == 4
                new_interval = (min(both_intervals_edges), max(both_intervals_edges))
                merged_intervals.remove(interval1)
                merged_intervals.remove(interval2)
                merged_intervals.add(new_interval)
                merged_something_this_cycle = True
                break

        if not merged_something_this_cycle:
            return merged_intervals

def get_merged_intervals(intervals):
    assert is_iterable_of_intervals(intervals)

    if not intervals:
        return set()

    sorted_intervals = sorted(intervals)
    merged_intervals = set()
    curr_merged_interval_start, curr_merged_interval_end = sorted_intervals[0]
    for interval in sorted_intervals[1:]:
        if interval[0] > curr_merged_interval_end:
            merged_intervals.add((curr_merged_interval_start, curr_merged_interval_end))
            curr_merged_interval_start, curr_merged_interval_end = interval
        else:
            curr_merged_interval_end = max(curr_merged_interval_end, interval[1])
    merged_intervals.add((curr_merged_interval_start, curr_merged_interval_end))

    return merged_intervals

def naive_get_interval_to_merged_interval(intervals):
    intervals = set(intervals)
    merged_intervals = get_merged_intervals(intervals)
    interval_to_merged_interval = {}
    for interval in intervals:
        for merged_interval in merged_intervals:
            if is_interval1_contained_in_interval2(interval, merged_interval):
                assert interval not in interval_to_merged_interval
                interval_to_merged_interval[interval] = merged_interval
                # break # not breaking makes things slower but more safe - it is unnecessary if we trust get_merged_intervals() completely.
    return interval_to_merged_interval

def naive_get_interval_to_merged_interval_and_merged_interval_to_intervals(intervals):
    interval_to_merged_interval = naive_get_interval_to_merged_interval(intervals)

    merged_interval_to_intervals = collections.defaultdict(set)
    for interval, merged_interval in interval_to_merged_interval.items():
        merged_interval_to_intervals[merged_interval].add(interval)
    merged_interval_to_intervals = dict(merged_interval_to_intervals) # I don't want a defaultdict moving around.

    return interval_to_merged_interval, merged_interval_to_intervals

def naive_get_merged_interval_to_intervals(intervals):
    return naive_get_interval_to_merged_interval_and_merged_interval_to_intervals(intervals)[1]

def do_intervals_overlap(interval1, interval2):
    # an alternative would be: is the distance between the interval middles smaller or equal to the sum of half their lengths?
    return len(get_merged_intervals((interval1, interval2))) == 1

# https://en.wikipedia.org/wiki/Partition_of_an_interval
def find_partition_to_subintervals_and_return_subinterval_to_containing_given_intervals(given_intervals, region_to_partition_start=None, region_to_partition_end=None):
    assert given_intervals
    assert is_iterable_of_intervals(given_intervals)
    if region_to_partition_start is None:
        region_to_partition_start = min(start for start, _ in given_intervals)
    if region_to_partition_end is None:
        region_to_partition_end = max(end for _, end in given_intervals)
    assert min(start for start, _ in given_intervals) >= region_to_partition_start
    assert max(end for _, end in given_intervals) <= region_to_partition_end

    ordered_partition_subintervals = []
    subinterval_to_containing_given_intervals = {}

    given_intervals_edges = set(itertools.chain.from_iterable(given_intervals))
    given_interval_edge_to_edge_type_to_given_intervals = {
        edge: collections.defaultdict(set)
        for edge in given_intervals_edges
    }
    for given_interval in given_intervals:
        start, end = given_interval
        given_interval_edge_to_edge_type_to_given_intervals[start]['start'].add(given_interval)
        given_interval_edge_to_edge_type_to_given_intervals[end]['end'].add(given_interval)
    given_interval_edge_to_edge_type_to_given_intervals = {
        edge: dict(edge_type_to_given_intervals) # I don't want a defaultdict moving around.
        for edge, edge_type_to_given_intervals in given_interval_edge_to_edge_type_to_given_intervals.items()
    }

    # print(f'given_interval_edge_to_edge_type_to_given_intervals: {given_interval_edge_to_edge_type_to_given_intervals}')
    # print()

    curr_containing_given_intervals = set()
    curr_subinterval_start = region_to_partition_start
    for given_interval_edge in sorted(given_interval_edge_to_edge_type_to_given_intervals):
        assert curr_subinterval_start <= given_interval_edge
        given_interval_edge_type_to_given_intervals = given_interval_edge_to_edge_type_to_given_intervals[given_interval_edge]
        if 'start' in given_interval_edge_type_to_given_intervals:
            if curr_subinterval_start < given_interval_edge:
                new_subinterval = (curr_subinterval_start, given_interval_edge - 1)
                ordered_partition_subintervals.append(new_subinterval)
                subinterval_to_containing_given_intervals[new_subinterval] = frozenset(curr_containing_given_intervals) # use frozenset to prevent mutability bugs...

            assert (not ordered_partition_subintervals) or (ordered_partition_subintervals[-1][1] == given_interval_edge - 1)
            curr_subinterval_start = given_interval_edge

            curr_containing_given_intervals = curr_containing_given_intervals | given_interval_edge_type_to_given_intervals['start']
            
        if 'end' in given_interval_edge_type_to_given_intervals:
            assert curr_containing_given_intervals
            assert curr_subinterval_start <= given_interval_edge
            new_subinterval = (curr_subinterval_start, given_interval_edge)
            ordered_partition_subintervals.append(new_subinterval)
            subinterval_to_containing_given_intervals[new_subinterval] = frozenset(curr_containing_given_intervals) # use frozenset to prevent mutability bugs...

            curr_subinterval_start = given_interval_edge + 1

            curr_containing_given_intervals = curr_containing_given_intervals - given_interval_edge_type_to_given_intervals['end']

    assert not curr_containing_given_intervals
    assert curr_subinterval_start <= region_to_partition_end + 1
    if curr_subinterval_start < region_to_partition_end + 1:
        assert ordered_partition_subintervals[-1][1] < region_to_partition_end
        new_subinterval = (curr_subinterval_start, region_to_partition_end)
        ordered_partition_subintervals.append(new_subinterval)
        subinterval_to_containing_given_intervals[new_subinterval] = frozenset() # use frozenset to prevent mutability bugs...

    assert sorted(subinterval_to_containing_given_intervals) == ordered_partition_subintervals
    assert is_iterable_an_ordered_partition_to_subintervals(ordered_partition_subintervals)
    # return (ordered_partition_subintervals, subinterval_to_containing_given_intervals)
    return subinterval_to_containing_given_intervals

def find_partition_to_subintervals_and_return_subinterval_to_num_of_containing_given_intervals(given_intervals, region_to_partition_start=None, region_to_partition_end=None):
    assert is_iterable_of_intervals(given_intervals)

    given_intervals_counter = collections.Counter(given_intervals)

    subinterval_to_containing_given_intervals = find_partition_to_subintervals_and_return_subinterval_to_containing_given_intervals(
        given_intervals=given_intervals,
        region_to_partition_start=region_to_partition_start,
        region_to_partition_end=region_to_partition_end,
    )
    return {
        subinterval: sum(given_intervals_counter[given_interval] for given_interval in containing_given_intervals)
        for subinterval, containing_given_intervals in subinterval_to_containing_given_intervals.items()
    }

def naive_get_n_positions_that_touch_all_intervals(intervals, n):
    assert intervals
    assert is_iterable_of_intervals(intervals)
    assert n > 0

    intervals = set(intervals) # remove duplicates
    intervals_edges = set(itertools.chain.from_iterable(intervals))
    for curr_edges in itertools.permutations(intervals_edges, n):
        curr_edges_touch_all_intervals = True
        for interval in intervals:
            if all(not is_point_touching_interval(edge, interval) for edge in curr_edges):
                curr_edges_touch_all_intervals = False
                break
        if curr_edges_touch_all_intervals:
            return set(curr_edges)

def is_iterable_an_ordered_partition_to_subintervals(intervals):
    if not is_iterable_of_intervals(intervals):
        return False

    for curr_interval, next_interval in zip(intervals[:-1], intervals[1:]):
        if curr_interval[1] + 1 != next_interval[0]:
            return False

    return True

def get_intersection_of_2_intervals(interval1, interval2):
    assert is_interval(interval1)
    assert is_interval(interval2)

    intersection_interval = (max(interval1[0], interval2[0]),
                             min(interval1[1], interval2[1]))
    if intersection_interval[0] > intersection_interval[1]:
        return None

    return intersection_interval

def get_intersections_of_intervals_with_interval(intervals, interval):
    intersections = {get_intersection_of_2_intervals(x, interval) for x in intervals}
    return {x for x in intersections if x is not None}

def replace_chars_with_char(s, chars_to_replace, new_char):
    chars_in_s = set(s)
    new_str = s
    for char_to_replace in chars_in_s & set(chars_to_replace):
        new_str = new_str.replace(char_to_replace, new_char)
    return new_str

def get_count_dict_with_low_counts_merged(count_dict, low_count_max_fraction_of_total_count):
    total_count = sum(count_dict.values())
    low_counts = [x for x in count_dict.values() if (x / total_count <= low_count_max_fraction_of_total_count)]
    total_low_count = sum(low_counts)
    new_dict = {k: count for k, count in count_dict.items() if (count / total_count > low_count_max_fraction_of_total_count)}
    new_dict['other'] = total_low_count
    assert sum(new_dict.values()) == total_count
    return new_dict

def remove_redundant_trailing_zeros(str_repr, unit_str_repr=''):
    if '.' in str_repr:
        while str_repr.endswith(f'0{unit_str_repr}'):
            str_repr = str_repr[:-(len(unit_str_repr) + 1)] + unit_str_repr
        if str_repr.endswith(f'.{unit_str_repr}'):
            str_repr = str_repr[:-(len(unit_str_repr) + 1)] + unit_str_repr
    return str_repr

def get_num_rounded_to_thousands_str_repr(num, num_of_digits_after_decimal_point=0):
    # assert num >= 1e3
    str_repr = str(round(num / 1e3, num_of_digits_after_decimal_point)) + 'K'
    str_repr = remove_redundant_trailing_zeros(str_repr, unit_str_repr='K')
    return str_repr

def get_num_rounded_to_millions_str_repr(num, num_of_digits_after_decimal_point=0):
    # assert num >= 1e6
    str_repr = str(round(num / 1e6, num_of_digits_after_decimal_point)) + 'M'
    str_repr = remove_redundant_trailing_zeros(str_repr, unit_str_repr='M')
    return str_repr

def perform_mw_test(nums1, nums2, alternative, **kwargs):
    assert alternative in {'two-sided', 'less', 'greater'}
    u, pvalue = scipy.stats.mannwhitneyu(
        nums1,
        nums2,
        alternative=alternative,
        **kwargs,
    )
    mean_u = len(nums1) * len(nums2) / 2
    u_minus_mean_u = u - mean_u
    return {
        'pvalue': pvalue,
        'u': u,
        'u_minus_mean_u': u_minus_mean_u,
    }

def perform_mw_test_for_df(df, test_column_name, bool_column_name, alternative, **kwargs):
    nums1 = df[df[bool_column_name]][test_column_name]
    nums2 = df[~df[bool_column_name]][test_column_name]
    return perform_mw_test(nums1, nums2, alternative=alternative, **kwargs)

def perform_g_test_or_fisher_exact_test(df, boolean_column1_name, boolean_column2_name, alternative, return_matrix_in_4_keys=False, use_fisher_exact_anyway=False):
    # Require expected and observed frequencies to be at least 5, according to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.power_divergence.html:
    # "This test is invalid when the observed or expected frequencies in each category are too small. A typical rule is that all of the observed and expected frequencies
    # should be at least 5."
    MIN_EXPECTED_AND_OBSERVED_FREQ_FOR_G_TEST = 5

    # it seems to me that G-test the way I use it is equivalent to fisher exact with two-sided alternative.
    # scipy.stats.chi2_contingency([[10000,10000],[10000,12000]], lambda_="log-likelihood")
    # scipy.stats.fisher_exact([[10000,10000],[10000,12000]], alternative='two-sided')

    assert not df[boolean_column1_name].isna().any()
    assert not df[boolean_column2_name].isna().any()
    counts_df = df[[boolean_column1_name, boolean_column2_name]].value_counts()
    false_false_count = counts_df[(False, False)] if ((False, False) in counts_df) else 0
    false_true_count = counts_df[(False, True)] if ((False, True) in counts_df) else 0
    true_false_count = counts_df[(True, False)] if ((True, False) in counts_df) else 0
    true_true_count = counts_df[(True, True)] if ((True, True) in counts_df) else 0
    matrix_for_test = [[false_false_count, false_true_count],
                       [true_false_count, true_true_count]]
    matrix_total_count = int(np.array(matrix_for_test).sum())

    # non two-sided requires fisher just because I didn't try to figure out how to perform the test with another type of alternative hypothesis. as mentioned above, as far as I
    # can tell, the default result that the G test gives is for the two-sided alternative hypothesis.
    if use_fisher_exact_anyway or (alternative != 'two-sided') or (not (np.array(matrix_for_test) >= MIN_EXPECTED_AND_OBSERVED_FREQ_FOR_G_TEST).all()):
        use_fisher_exact = True
    else:
        use_fisher_exact = not (scipy.stats.contingency.expected_freq(matrix_for_test) >= MIN_EXPECTED_AND_OBSERVED_FREQ_FOR_G_TEST).all()
        # print(scipy.stats.contingency.expected_freq(matrix_for_test))

    if use_fisher_exact:
        test_performed = 'fisher_exact'
        odds_ratio, pvalue = scipy.stats.fisher_exact(matrix_for_test, alternative=alternative)
    else:
        test_performed = 'g'
        g_test_result = scipy.stats.chi2_contingency(matrix_for_test, lambda_="log-likelihood")
        assert true_false_count > 0
        assert false_true_count > 0
        odds_ratio = true_true_count * false_false_count / (true_false_count * false_true_count)
        pvalue = g_test_result[1]

    if return_matrix_in_4_keys:
        test_result = {
            'matrix_for_test_false_false': false_false_count,
            'matrix_for_test_false_true': false_true_count,
            'matrix_for_test_true_false': true_false_count,
            'matrix_for_test_true_true': true_true_count,
            'matrix_total_count': matrix_total_count,
            'odds_ratio': odds_ratio,
            'pvalue': pvalue,
            'test_performed': test_performed,
        }
    else:
        test_result = {
            'matrix_for_test': matrix_for_test,
            'matrix_total_count': matrix_total_count,
            'odds_ratio': odds_ratio,
            'pvalue': pvalue,
            'test_performed': test_performed,
        }
    return test_result

# test_df = pd.DataFrame([
#     *([(False, False)] * 3500),
#     *([(False, True)] * 70),
#     *([(True, False)] * 400),
#     *([(True, True)] * 4),
# ], columns=['a', 'b'])
# # observed is too low
# print(perform_g_test_or_fisher_exact_test(test_df, 'a', 'b', alternative='two-sided'), scipy.stats.fisher_exact([[3500,70],[400, 4]]))
#
# test_df = pd.DataFrame([
#     *([(False, False)] * 3500),
#     *([(False, True)] * 70),
#     *([(True, False)] * 400),
#     *([(True, True)] * 5),
# ], columns=['a', 'b'])
# print(perform_g_test_or_fisher_exact_test(test_df, 'a', 'b', alternative='two-sided'), scipy.stats.chi2_contingency([[3500,70],[400, 5]], lambda_="log-likelihood"))
#
# test_df = pd.DataFrame([
#     *([(False, False)] * 3500),
#     *([(False, True)] * 30),
#     *([(True, False)] * 400),
#     *([(True, True)] * 5),
# ], columns=['a', 'b'])
# # expected is too low
# print(perform_g_test_or_fisher_exact_test(test_df, 'a', 'b', alternative='two-sided'), scipy.stats.fisher_exact([[3500,30],[400, 5]]))

@execute_if_output_doesnt_exist_already
def cached_perform_linear_fit_of_column_histogram(
        input_file_path_df_csv,
        column_name,
        bins,
        fit_log10_of_hist,
        output_file_path_linear_fit_result_pickle,
):
    df = pd.read_csv(input_file_path_df_csv, sep='\t', low_memory=False)
    column = df[column_name]
    assert not column.isna().any()
    hist, _ = np.histogram(column, bins=bins)

    hist_sum = hist.sum()
    hist_median = np.median(hist)

    std_err_of_estimated_gradient = r_value = p_value = intercept = slope = np.nan
    # r_value = p_value = slope_std_err = intercept_std_err = intercept = slope = np.nan
    reason_for_not_performing_fit = None

    if fit_log10_of_hist:
        if (hist > 0).all():
            ys = np.log10(hist)
        else:
            reason_for_not_performing_fit = 'fit_log10_of_hist=True but one of the hist values was zero'
    else:
        ys = hist

    if reason_for_not_performing_fit is None:
        bins = np.array(bins)
        xs = (bins[:-1] + bins[1:]) / 2
        # only numpy 1.7 has the alternative argument. If I understand correctly, the default until then is the same as alternative='two-sided' in newer versions.
        slope, intercept, r_value, p_value, std_err_of_estimated_gradient = scipy.stats.linregress(
            xs, ys,
            # alternative='two-sided',
        )

        # wrote this (but didn't test it) before finding out about scipy.stats.linregress.
        # coefficients, (residual_sum_of_squares, _, _, _) = np.polynomial.Polynomial.fit(xs, ys, deg=1, full=True)
        # intercept, slope = coefficients

    linear_fit_result = {
        'intercept': intercept,
        'slope': slope,
        'r_value': r_value,
        'p_value': p_value,
        # the following two are only available in newer versions, while the version i use instead gives std_err_of_estimated_gradient, if I understand correctly.
        # 'slope_std_err': slope_std_err,
        # 'intercept_std_err': intercept_std_err,
        'std_err_of_estimated_gradient': std_err_of_estimated_gradient,
        'reason_for_not_performing_fit': reason_for_not_performing_fit,
        'hist_sum': hist_sum,
        'hist_median': hist_median,
    }

    with open(output_file_path_linear_fit_result_pickle, 'wb') as f:
        pickle.dump(linear_fit_result, f, protocol=4)

def perform_linear_fit_of_column_histogram(
        input_file_path_df_csv,
        column_name,
        bins,
        fit_log10_of_hist,
        output_file_path_linear_fit_result_pickle,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_perform_linear_fit_of_column_histogram(
        input_file_path_df_csv=input_file_path_df_csv,
        column_name=column_name,
        bins=bins,
        fit_log10_of_hist=fit_log10_of_hist,
        output_file_path_linear_fit_result_pickle=output_file_path_linear_fit_result_pickle,
    )

def get_vertex_to_connected_component_index_for_undirected_graph(vertices, edges):
    assert isinstance(vertices, set)
    assert isinstance(edges, set)

    vertex_index_to_vertex = {i: v for i, v in enumerate(vertices)}
    vertex_to_vertex_index = {v: i for i, v in enumerate(vertices)}
    edges_as_sorted_vertex_indices = {tuple(sorted((vertex_to_vertex_index[v1], vertex_to_vertex_index[v2]))) for v1, v2 in edges}

    vertex_index_to_connected_component = {i: frozenset({i}) for i in vertex_index_to_vertex}
    for i1, i2 in edges_as_sorted_vertex_indices:
        new_connected_component = vertex_index_to_connected_component[i1] | vertex_index_to_connected_component[i2]
        for i in new_connected_component:
            vertex_index_to_connected_component[i] = new_connected_component

    index_connected_components = set(vertex_index_to_connected_component.values())
    index_connected_component_to_connected_component_index = {c: i for i, c in enumerate(index_connected_components)}
    vertex_to_connected_component_index = {vertex_index_to_vertex[i]: index_connected_component_to_connected_component_index[c]
                                           for i, c in vertex_index_to_connected_component.items()}
    return vertex_to_connected_component_index

assert len(set(get_vertex_to_connected_component_index_for_undirected_graph(set(range(9)), set()).values())) == 9
assert len(set(get_vertex_to_connected_component_index_for_undirected_graph(set(range(9)), {(1, 2), (2, 8), (6, 7)}).values())) == 6
assert len(set(get_vertex_to_connected_component_index_for_undirected_graph(
    set(range(9)), {(1, 2), (2, 8), (6, 7), (2, 3), (4, 3), (3, 7), (5, 4)}).values())) == 2
assert len(set(get_vertex_to_connected_component_index_for_undirected_graph(
    set(range(9)), {(1, 2), (2, 8), (6, 7), (2, 3), (4, 3), (3, 7), (5, 4), (6, 3)}).values())) == 2
assert len(set(get_vertex_to_connected_component_index_for_undirected_graph(
    set(range(9)), {(1, 2), (2, 8), (6, 7), (2, 3), (4, 3), (3, 7), (5, 4), (6, 0)}).values())) == 1

def uncomment_python_line(line):
    line_before_hash, _, line_after_hash = line.partition('#')
    return line_before_hash + line_after_hash.lstrip()

def plot_line_on_ax(ax, slope=1, offset=0, color='red', alpha=0.3, **kwargs):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    min_lim = min(xlim[0],ylim[0])
    max_lim = max(xlim[1],ylim[1])
    ax.plot(
        (min_lim, max_lim), (slope * min_lim + offset, slope * max_lim + offset),
        color=color,
        alpha=alpha,
        **kwargs,
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def plot_y_equals_x_line_on_ax(ax, offset=0, linestyle='--', alpha=0.5, color='grey', **kwargs):
    plot_line_on_ax(ax, slope=1, offset=offset, linestyle=linestyle, alpha=alpha, color=color, **kwargs)

def get_num_of_header_lines(file_path, header_line_prefix, allow_only_header_lines=False):
    lines = read_text_file(file_path).splitlines()
    # print('aoeu')
    # print(file_path)
    # print(lines[:6])
    # print('snth')
    for i, line in enumerate(lines):
        if not line.startswith(header_line_prefix):
            return i
    if allow_only_header_lines:
        return len(lines)
    raise RuntimeError('no non-header lines in file. maybe the program creating the file crashed?')
    return i

def get_num_series_with_or_without_zeros(nums, ignore_zeros=False):
    if ignore_zeros:
        return nums[nums != 0]
    return nums

def make_all_spines_invisible(ax):
    for spine_name in ('right', 'left', 'top', 'bottom'):
        ax.spines[spine_name].set_visible(False)

def make_all_spines_visible(ax):
    for spine_name in ('right', 'left', 'top', 'bottom'):
        ax.spines[spine_name].set_visible(True)

def make_all_spines_and_x_and_y_axes_invisible(ax):
    make_all_spines_invisible(ax)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    

def get_dict_mapping_one_df_column_to_other(df, key_column_name, val_column_name):
    return {k: v for k,v in df[[key_column_name, val_column_name]].to_records(index=False).tolist()}

### NOTE: some attempts to produce random colors that seaborn can parse. did not completely succeed...
import colorsys
import binascii
def get_N_HexCol(n=5):
    HSV_tuples = [(x*1.0/n, 0.5, 0.5) for x in range(n)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#' + "".join(map(lambda x: hex(x)[2:].lower(),rgb)))
    return hex_out

def get_2_digit_hex_repr(num):
    assert 0 <= num <= 255
    return hex(num)[2:].zfill(2)

def get_hex_color_from_rgba_color(rgba_color):
    assert len(rgba_color) == 4
    return '#' + ''.join(map(lambda x: get_2_digit_hex_repr(int(x * 255)), rgba_color[:-1]))

def get_n_colors(n, random_seed=0):
    if n <= len(TEN_DEFAULT_MATPLOTLIB_COLORS_AND_SOME_MORE_KIND_OF_DISTINGUISHABLE):
        return TEN_DEFAULT_MATPLOTLIB_COLORS_AND_SOME_MORE_KIND_OF_DISTINGUISHABLE[:n]
    
    rgba_colors = matplotlib.cm.gist_rainbow(np.linspace(0, 1, n))
    np.random.seed(random_seed)
    rgba_colors = np.random.permutation(rgba_colors)
    # print(rgba_colors)
    return [get_hex_color_from_rgba_color(x) for x in rgba_colors]

def add_color_col(df, col, color_col=None):
    if color_col is None:
        color_col = f'{col}_color'
    uniq_col_vals = sorted(df[col].unique())
    val_to_color = {k: v for k, v in zip(uniq_col_vals, get_n_colors(len(uniq_col_vals)))}
    df[color_col] = df[col].map(val_to_color)
    assert not df[color_col].isna().any()
    

def get_num_colors(nums, colormap='binary', vmin=None, vmax=None, also_return_norm=False):
    if vmin is None:
        vmin = min(nums)
    if vmax is None:
        vmax = max(nums)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    rgba_colors = getattr(matplotlib.cm, colormap)(norm(nums))

    num_colors = [get_hex_color_from_rgba_color(x) for x in rgba_colors]
    if also_return_norm:
        return (num_colors, norm)
    return num_colors


def get_dict_str_repr(d, key_val_sep=':', item_sep=';', round_vals=False):
    return item_sep.join(
        f'{k}{key_val_sep}' + (f'{v if np.isnan(v) else round(v)}' if round_vals else f'{v}')
        for k,v in d.items()
    )


def get_best_separator(nums1, nums2):
    # print(set(nums1), set(nums2))
    total_count1 = len(nums1)
    total_count2 = len(nums2)
    total_count = total_count1 + total_count2
    min_val_to_test_smaller_separator = max(nums1.min(), nums2.min())
    max_val_to_test_bigger_separator = min(nums1.max(), nums2.max())

    # print(f'sorted(nums1): {sorted(nums1)}')
    # print(f'sorted(nums2): {sorted(nums2)}')

    sorted_counts_df = pd.Series(nums1).value_counts().reset_index(name='count1').merge(
        pd.Series(nums2).value_counts().reset_index(name='count2'), how='outer',
    ).fillna(0).sort_values('index').rename(columns={'index': 'num'})
    # print(f'sorted_counts_df: {sorted_counts_df}')
    # print(f'sorted_counts_df["num"].max(): {sorted_counts_df["num"].max()}')
    # print(f'sorted_counts_df["num"].min(): {sorted_counts_df["num"].min()}')


    left_vals_to_skip_mask = sorted_counts_df['num'] < min_val_to_test_smaller_separator
    if not left_vals_to_skip_mask.any():
        # nothing to skip, which actually means we can skip the first one.
        left_vals_to_skip_mask = np.full(len(sorted_counts_df), False)
        left_vals_to_skip_mask[0] = True
        
    curr_nums1_to_the_left = sorted_counts_df.loc[left_vals_to_skip_mask, 'count1'].sum()
    curr_nums2_to_the_left = sorted_counts_df.loc[left_vals_to_skip_mask, 'count2'].sum()

    # print(f'curr_nums1_to_the_left: {curr_nums1_to_the_left}')
    # print(f'curr_nums2_to_the_left: {curr_nums2_to_the_left}')

    max_val_to_test_smaller_separator = sorted_counts_df.loc[sorted_counts_df['num'] > max_val_to_test_bigger_separator, 'num'].min()
    if np.isnan(max_val_to_test_smaller_separator):
        max_val_to_test_smaller_separator = sorted_counts_df['num'].max()
    # print(f'max_val_to_test_smaller_separator: {max_val_to_test_smaller_separator}')
    vals_to_go_over_mask = (~left_vals_to_skip_mask) & (sorted_counts_df['num'] <= max_val_to_test_smaller_separator)
    num_of_vals_to_go_over = vals_to_go_over_mask.sum()
    assert num_of_vals_to_go_over > 0
    # print(f'num_of_vals_to_go_over: {num_of_vals_to_go_over}')
    max_prob = 0
    max_prob_num_after_sep_i_in_df_slice = None
    for i, (_, row) in enumerate(sorted_counts_df[vals_to_go_over_mask].iterrows()):
        curr_total_count_to_the_left = curr_nums1_to_the_left + curr_nums2_to_the_left
        curr_total_count_to_the_right = total_count - curr_total_count_to_the_left
        curr_nums1_to_the_right = total_count1 - curr_nums1_to_the_left
        curr_nums2_to_the_right = total_count2 - curr_nums2_to_the_left
        prob_to_get_both_right_if_1_is_left = (curr_nums1_to_the_left / curr_total_count_to_the_left) * (curr_nums2_to_the_right / curr_total_count_to_the_right)
        prob_to_get_both_right_if_2_is_left = (curr_nums2_to_the_left / curr_total_count_to_the_left) * (curr_nums1_to_the_right / curr_total_count_to_the_right)
        curr_max_prob = max(prob_to_get_both_right_if_1_is_left, prob_to_get_both_right_if_2_is_left)
        if curr_max_prob > max_prob:
            max_prob = curr_max_prob
            max_prob_num_after_sep_i_in_df_slice = i

        curr_nums1_to_the_left += row['count1']
        curr_nums2_to_the_left += row['count2']

    max_prob_num_after_sep_i_in_df = max_prob_num_after_sep_i_in_df_slice + left_vals_to_skip_mask.sum()
    best_separator = sorted_counts_df['num'].iloc[(max_prob_num_after_sep_i_in_df - 1):(max_prob_num_after_sep_i_in_df + 1)].mean()
    return best_separator, max_prob


def get_random_str(str_len):
    return ''.join(random.choices(string.ascii_letters + string.digits, k=str_len))


def set_is_empty_or_contains_only_np_nan(s):
    if not s:
        return True
    return (len(s) == 1) and np.isnan(next(iter(s)))

def get_pvalue_repr(pvalue, verbose=False):
    if pvalue == 0:
        return '0'
    pval_exponent = math.floor(np.log10(pvalue))
    if verbose:
        print(f'pvalue: {pvalue}')
    if pval_exponent >= -4:
        pvalue_repr = f'{pvalue:.4f}'
    else:
        pval_exponent_repr = f'\\mathdefault{{10^{{{pval_exponent}}}}}'
        pval_without_exponent = pvalue / (10 ** pval_exponent)
        if pval_without_exponent == 1:
            pvalue_repr = f'${pval_exponent_repr}$'
        else:
            pval_repr_without_exponent = f'{pval_without_exponent:.1f}'
            pvalue_repr = f'${pval_repr_without_exponent}\cdot{pval_exponent_repr}$'
        if verbose:
            print(pval_repr_without_exponent)
            print(f'pval_repr_without_exponent: {pval_repr_without_exponent}')
    
    if verbose:
        print(f'pvalue_repr: {pvalue_repr}')
    
    return pvalue_repr

def dilute_str_list(strs, dilution_factor=None, target_num_of_strs=50, dont_dilute_last=True):
    if target_num_of_strs is not None:
        assert dilution_factor is None
        dilution_factor = len(strs) // target_num_of_strs
    last_i = len(strs) - 1
    def should_dilute(i):
        if dont_dilute_last and (i == last_i):
            return False
        return ((i % dilution_factor) != 0) or (i > last_i - (dilution_factor / 2))
    return ['' if should_dilute(i) else x for i, x in enumerate(strs)]


def get_int_powers_of_2_in_range(min, max):
    exponents = list(range(math.floor(np.log2(min)), math.ceil(np.log2(max) + 1)))
    # print(exponents)
    powers = np.exp2(exponents)
    powers = powers[(powers >= min) & (powers <= max)]
    return powers

def ax_has_a_legend(ax):
    return len(ax.get_legend_handles_labels()[0]) > 0

def plot_texts(ax, df, x_column_name, y_column_name, text_column_name, plot_mask=None):
    if plot_mask is None:
        plot_mask = np.full(len(df), True)
    
    for _, row in df[plot_mask].iterrows():
        ax.text(row[x_column_name], row[y_column_name], row[text_column_name])

def normalize_df_columns(df, allow_nans=False, inplace=False, column_names=None, col_to_min_and_max=dict(), allow_auto_min_max=False):
    assert not df.empty
    if not inplace:
        df = df.copy()

    if column_names is None:
        column_names = list(df)

    if not allow_nans:
        assert not df[column_names].isna().any().any()
    for column_name in column_names:
        if allow_nans and df[column_name].isna().all():
            continue

        non_nan_vals_set = set(df[column_name].dropna().unique())
        if len(non_nan_vals_set) == 1:
            only_val = next(iter(non_nan_vals_set))
            print(f'cant normalize {column_name} such that min is 0 and max is 1, because there is only one unique val ({only_val}). so setting to nan.')
            df[column_name] = np.nan
        
        else:
            if column_name in col_to_min_and_max:
                min_val_for_norm, max_val_for_norm = col_to_min_and_max[column_name]
            else:    
                if allow_auto_min_max or (non_nan_vals_set == {0, 1}):
                    min_val_for_norm = df[column_name].min()
                    max_val_for_norm = df[column_name].max()
                else:
                    raise RuntimeError(
                        f'set allow_auto_min_max=True or provide min and max vals for normalization for {column_name} in column_name_to_min_and_max ({sorted(df[column_name].unique())})')
            # print(column_name, df[column_name].unique())
            df[column_name] -= min_val_for_norm
            df.loc[df[column_name] < 0, column_name] = 0

            df[column_name] /= max_val_for_norm - min_val_for_norm
            df.loc[df[column_name] > 1, column_name] = 1

            assert np.nanmin(df[column_name]) >= 0
            assert np.nanmax(df[column_name]) <= 1

    if not inplace:
        return df

def get_color_df(
        df, allow_nans=False, 
        col_to_cmap_and_min_and_max=dict(),
        col_to_val_to_color=dict(),
):
    assert (set(col_to_cmap_and_min_and_max) | set(col_to_val_to_color)) == set(df.columns)
    assert not (set(col_to_cmap_and_min_and_max) & set(col_to_val_to_color))

    df = normalize_df_columns(
        df, allow_nans=allow_nans, column_names=list(col_to_cmap_and_min_and_max), 
        col_to_min_and_max={k: v[1:] for k, v in col_to_cmap_and_min_and_max.items()},
    )

    # if not in_place_add_color_cols:
    #     assert (set(col_to_cmap_and_min_and_max) | set(col_to_val_to_color)) == set(df.columns)
    # assert not (set(col_to_cmap_and_min_and_max) & set(col_to_val_to_color))

    # if in_place_add_color_cols:
    #     for col in col_to_cmap_and_min_and_max:
    #         color_col = f'{col}_color'
    #         df[f'']
    #     cols = sorted(set(col_to_cmap_and_min_and_max) | set(col_to_val_to_color))
    #     for col in cols:
    #         df[]
    
    for col in df.columns:
        if col in col_to_cmap_and_min_and_max:
            cmap = col_to_cmap_and_min_and_max[col][0]
            # df[col] = cmap(df[col])
            df[col] = df[col].apply(lambda x: matplotlib.colors.to_hex(cmap(x)))
        else:
            df[col] = df[col].map(col_to_val_to_color[col])

    return df

def get_color_arr(
        df, allow_nans=False, 
        col_to_cmap_and_min_and_max=dict(),
        col_to_val_to_color=dict(),
):
    assert (set(col_to_cmap_and_min_and_max) | set(col_to_val_to_color)) == set(df.columns)
    assert not (set(col_to_cmap_and_min_and_max) & set(col_to_val_to_color))

    df = normalize_df_columns(
        df, allow_nans=allow_nans, column_names=list(col_to_cmap_and_min_and_max), 
        col_to_min_and_max={k: v[1:] for k, v in col_to_cmap_and_min_and_max.items()},
    )
    
    color_arr = np.full((*df.shape, 4), np.nan)
    for col_i, col in enumerate(df.columns):
        if col in col_to_cmap_and_min_and_max:
            cmap = col_to_cmap_and_min_and_max[col][0]
            color_arr[:,col_i] = cmap(df[col])
        else:
            col_as_colors = df[col].map(col_to_val_to_color[col])

            assert col_as_colors.notna().all()
            color_arr[:,col_i] = np.array(list(col_as_colors.apply(lambda x: np.array(matplotlib.colors.to_rgba(x)))))

    return color_arr


def split_sorted_counts_to_roughly_same_size_bins(counts, num_of_bins):
    total_count = sum(counts)
    wanted_bin_total_count = total_count / num_of_bins
    print(f'wanted_bin_total_count: {wanted_bin_total_count}')
    cum_counts_series = pd.Series(np.array(counts).cumsum())
    
    wanted_bin_ends_cum_counts = [wanted_bin_total_count * x for x in range(1, num_of_bins + 1)]
    assert np.isclose(wanted_bin_ends_cum_counts[-1], cum_counts_series.iloc[-1])

    bin_indices = np.full(len(counts), np.nan)
    for curr_bin_i, wanted_bin_end_cum_count in enumerate(wanted_bin_ends_cum_counts):
        closest_cum_count_i = (cum_counts_series - wanted_bin_end_cum_count).abs().idxmin()
        curr_assign_bin_mask = np.isnan(bin_indices) & (np.arange(len(bin_indices)) <= closest_cum_count_i)
        assert curr_assign_bin_mask.any()
        bin_indices[curr_assign_bin_mask] = curr_bin_i

    assert not np.isnan(bin_indices).any()
    bin_indices = bin_indices.astype(int)

    print(f'bin_indices_value_counts:\n{pd.Series(bin_indices).value_counts()}')
    print(pd.DataFrame({'count': counts, 'bin_i': bin_indices}).groupby('bin_i')['count'].sum().reset_index(name='bin_total_count'))

    return bin_indices


def equal_size_bin_digitize_by_single_nums_vec(nums, num_of_bins):
    num_of_nums = len(nums)
    bin_size = num_of_nums // num_of_bins
    bin_start_indices = np.arange(num_of_bins) * bin_size
    # print(f'bin_start_indices: {bin_start_indices}')
    sorted_nums = sorted(nums)
    bins = list(np.array(sorted_nums)[bin_start_indices]) + [sorted_nums[-1]]
    # print(f'bins: {bins}')
    assert len(bins) == len(set(bins)), f'bins: {bins} - probably too many duplicated vals (many zeros?), so bins are not unique.'

    assert len(bins) == num_of_bins + 1
    
    bin_indices = safe_digitize(nums, bins=bins)
    
    return bin_indices

def equal_size_bin_digitize(list_of_nums_and_num_of_bins):
    # num_of_dimensions = len(list_of_nums_and_num_of_bins)
    df = pd.DataFrame({
        i: equal_size_bin_digitize_by_single_nums_vec(nums, num_of_bins) 
        for i, (nums, num_of_bins) in enumerate(list_of_nums_and_num_of_bins)
    })
    df = merge_preserving_df1_index_and_row_order(df, df.drop_duplicates().sort_values(list(df)).reset_index(drop=True).reset_index())
    return df['index'].to_numpy()

def get_bin_median_and_repr_vecs(
        nums, bins=None, num_of_bins=None, add_bin_size_if_non_equal_size=False,
        bin_edge_to_repr_func=lambda x: remove_redundant_trailing_zeros(f'{x:.1f}'),
):
    if num_of_bins:
        assert bins is None
        bin_indices = equal_size_bin_digitize_by_single_nums_vec(nums, num_of_bins)
        equal_bin_size = True
    else:
        assert num_of_bins is None
        bin_indices = safe_digitize(nums, bins=bins)
        num_of_bins = len(bins) - 1
        equal_bin_size = False
        # bin_sizes = pd.Series(bin_indices).value_counts()
        # print(f'bin_sizes: {bin_sizes}')

    bin_medians = np.array([np.median(nums[bin_indices == i]) for i in range(num_of_bins)])
    median_vec = bin_medians[bin_indices]
    
    # nums_series = pd.Series(nums) # TODO: discard if not needed
    # median_vec = nums_series.map({i: v for i, v in enumerate(bin_medians)}) # TODO: discard if not needed
    
    bin_reprs = []
    for i in range(num_of_bins):
        curr_nums = nums[bin_indices == i]
        min_val = curr_nums.min()
        max_val = curr_nums.max()
        curr_repr = f'{bin_edge_to_repr_func(min_val)}_{bin_edge_to_repr_func(max_val)}'
        if add_bin_size_if_non_equal_size and (not equal_bin_size):
            curr_repr += f' ({len(curr_nums)})'
        bin_reprs.append(curr_repr)
    bin_reprs = np.array(bin_reprs)
    repr_vec = bin_reprs[bin_indices]

    return median_vec, repr_vec

def verbose_drop_duplicates(df, df_desc='df'):
    num_of_rows_before_drop_duplicates = len(df)
    df = df.drop_duplicates()
    num_of_rows_after_drop_duplicates = len(df)
    num_of_rows_dropped = num_of_rows_before_drop_duplicates - num_of_rows_after_drop_duplicates
    if num_of_rows_dropped > 0:
        print(f'dropped {num_of_rows_dropped} duplicate rows from {df_desc}.')
    return df

def safe_digitize(nums, bins):
    assert len(bins) == len(set(bins))
    
    nums = np.array(nums)
    # print(nums.min())
    assert nums.min() >= bins[0]
    assert nums.max() <= bins[-1]

    num_of_bins = len(bins) - 1
    bin_indices = np.digitize(nums, bins=bins) - 1
    
    equal_to_last_bin_upper_mask = nums == bins[-1]
    assert (equal_to_last_bin_upper_mask == (bin_indices == num_of_bins)).all()
    bin_indices[equal_to_last_bin_upper_mask] = num_of_bins - 1

    return bin_indices

def add_bin_desc_column(df, val_column_name, bins, bin_desc_column_name=None):
    if bin_desc_column_name is None:
        bin_desc_column_name = f'{val_column_name}_range'
    
    assert bin_desc_column_name not in df.columns
    
    bin_indices = safe_digitize(df[val_column_name], bins=bins)
    for bin_i, (low, high) in enumerate(zip(bins[:-1], bins[1:])):
        range_mask = bin_indices == bin_i
        if range_mask.any():
            df.loc[range_mask, bin_desc_column_name] = f'{low}_to_{high}'

def get_2d_array_row_major(arr):
    raise RuntimeError('ugh. my implementation was wrong. use mc.ut.to_layout, which uses np.ravel. ugh.')

def get_2d_array_col_major(arr):
    raise RuntimeError('ugh. my implementation was wrong. use mc.ut.to_layout, which uses np.ravel. ugh.')

def get_z_score_and_ztest_pval(val, vals):
    z_score = scipy.stats.zmap(val, vals)[0]
    ztest_pval = scipy.stats.norm.sf(abs(z_score)) * 2
    return z_score, ztest_pval

def set_pd_display_large_repr_and_row_options():
    pd.set_option('display.large_repr', 'truncate')
    pd.set_option('display.max_rows', 40)
    pd.set_option('display.min_rows', 40)

def get_imported_module_names_and_aliases(globals_dict):
    imported_module_names_and_aliases = set()
    for name, val in globals_dict.items():
        if isinstance(val, types.ModuleType):
            imported_module_names_and_aliases |= {val.__name__, name}
    return sorted(imported_module_names_and_aliases)

def print_dis_output_relevant_for_detecting_global_variables_used_by_function(f, globals_dict):
    # print_dis_output_relevant_for_detecting_global_variables_used_by_function(f, globals())
    _, dis_output = get_func_ret_val_and_stdout(dis.dis, f)
    dis_output_lines = dis_output.splitlines()
    
    # uninteresting_global_names = set(get_imported_module_names_and_aliases(globals_dict) + dir(__builtins__))
    # assert 'str' in dir(__builtins__) # was False for some reason
    uninteresting_global_names = set(get_imported_module_names_and_aliases(globals_dict) + list(builtins.__dict__.keys()))
    
    print('\n'.join([
        dis_output_line for dis_output_line in dis_output_lines 
        if ('LOAD_GLOBAL' in dis_output_line) and (not any(f'({x})' in dis_output_line for x in uninteresting_global_names))
    ]))

def print_dis_output_relevant_for_detecting_global_variables_used_for_each_global_functions(globals_dict):
    # print_dis_output_relevant_for_detecting_global_variables_used_for_each_global_functions(globals())
    for func_name, func in [(k,v) for k,v in globals_dict.items() if isinstance(v, types.FunctionType)]:
        print(func_name)
        print_dis_output_relevant_for_detecting_global_variables_used_by_function(func, globals_dict)
        print()

def get_k_to_vs_as_k_v_pairs(k_to_vs):
    k_v_pairs = []
    for k, vs in k_to_vs.items():
        for v in vs:
            k_v_pairs.append((k, v))
    return k_v_pairs

def df_to_index_and_column_namd_and_val_df(df):
    flat_dicts = []
    for index, row in df.iterrows():
        for column, val in row.items():
            flat_dicts.append({
                'orig_index': index,
                'orig_column_name': column,
                'orig_index_and_column_name': f'{index},{column}',
                'val': val,
            })
    return pd.DataFrame(flat_dicts)

def replace_vals_by_join(
        df1, df2, cols_to_join_on, 
        ignore_cols_without_nans_in_df1=False, allow_df1_rows_missing_from_df2=False,
        ignore_conflicts=False, is_allowed_conflict=None,
):
    common_cols = set(df1) & set(df2)
    
    if isinstance(cols_to_join_on, str):
        cols_to_join_on = [cols_to_join_on]
    else:
        assert isinstance(cols_to_join_on, list)
    
    assert set(cols_to_join_on) <= common_cols

    df1_cols_with_nans = list(df1.isna().any(axis=0).loc[lambda x: x].index)
    cols_to_overwrite = common_cols - set(cols_to_join_on)
    if ignore_cols_without_nans_in_df1:
        print(f'NOTE: ignore_cols_without_nans_in_df1=True sounds like a bad idea. better to make sure there are no conflicts...ignoring {cols_to_overwrite - set(df1_cols_with_nans)}')
        cols_to_overwrite &= set(df1_cols_with_nans)
    col_to_temp_col = {x: f'{x}_new' for x in cols_to_overwrite}
    if ignore_cols_without_nans_in_df1:
        common_columns_to_ignore_names = sorted(common_cols - set(col_to_temp_col) - set(cols_to_join_on))
        if common_columns_to_ignore_names:
            print(f'ignoring {common_columns_to_ignore_names} in df2 because they contain no nan values in df1.')
            df2 = df2.drop(columns=common_columns_to_ignore_names)
    
    kwargs = {'how': 'left'} if allow_df1_rows_missing_from_df2 else {}
    df = merge_preserving_df1_index_and_row_order(df1, df2.rename(columns=col_to_temp_col), on=cols_to_join_on, **kwargs)
    for col, col_with_suffix in col_to_temp_col.items():
        if not ignore_conflicts:
            non_nan_mask = df[[col, col_with_suffix]].notna().all(axis=1)
            conflict_mask = non_nan_mask & (df[col] != df[col_with_suffix])
            conflict_df = df.loc[conflict_mask, [*cols_to_join_on, col, col_with_suffix]]
            if is_allowed_conflict is not None:
                conflict_df = conflict_df[~(is_allowed_conflict(col, conflict_df[col], conflict_df[col_with_suffix]))]
            assert conflict_df.empty, f'conflict_df:\n{conflict_df}'
        df[col].fillna(df[col_with_suffix], inplace=True)
    
    df.drop(columns=list(col_to_temp_col.values()), inplace=True)
    return df

def get_non_exclusive_vals(iters):
    return set((pd.Series(list(itertools.chain.from_iterable(iters))).value_counts() > 1).loc[lambda x: x].index)

def sample_mask(mask, size, rand_seed=0, allow_not_sampling_due_to_too_few_pos=False):
    pos_count = mask.sum()
    if pos_count < size:
        if allow_not_sampling_due_to_too_few_pos:
            return mask
        raise RuntimeError(f'pos_count ({pos_count}) < size ({size})')
    elif pos_count == size:
        return mask

    if rand_seed is not None:
        np.random.seed(rand_seed)
    sampled_mask = np.full(len(mask), False)
    sampled_mask[np.random.choice(np.where(mask)[0], size, replace=False)] = True
    assert sampled_mask.sum() == size
    return sampled_mask

def sample_mask_per_group(mask, size, group_by_vec, rand_seed=0, allow_not_sampling_due_to_too_few_pos=False, drop_if_too_few_pos=False):
    if rand_seed is not None:
        np.random.seed(rand_seed)

    assert (not drop_if_too_few_pos) or (not allow_not_sampling_due_to_too_few_pos), 'cant have both drop_if_too_few_pos and allow_not_sampling_due_to_too_few_pos'

    sampled_mask = np.full(len(mask), False)

    group_vals = sorted(set(group_by_vec[mask]))

    num_of_group_vals = len(group_vals)
    for group_val in group_vals:
        # print(f'group_val: {group_val}')
        group_mask = (group_by_vec == group_val) & mask # 240109: ugh. only now added `& mask` here. so previous results using this code (in some contexts) were wrong.
        if drop_if_too_few_pos and (group_mask.sum() < size):
            continue
        group_sampled_mask = sample_mask(
            group_mask, size, rand_seed=None, allow_not_sampling_due_to_too_few_pos=allow_not_sampling_due_to_too_few_pos)
        sampled_mask |= group_sampled_mask

    if not drop_if_too_few_pos:
        if not allow_not_sampling_due_to_too_few_pos:
            assert sampled_mask.sum() == size * num_of_group_vals
        else:
            assert sampled_mask.sum() <= size * num_of_group_vals
    return sampled_mask
    
def get_equal_contrib_mask(masks, size=None):
    mask_sizes = [x.sum() for x in masks]
    zero_size_mask_indices = [i for i, x in enumerate(mask_sizes) if x == 0]
    if zero_size_mask_indices:
        raise RuntimeError(f'zero_size_mask_indices: {zero_size_mask_indices}')
    if size is None:
        size = min(mask_sizes)
    mask = np.full(len(masks[0]), False)
    for curr_mask in masks:
        mask |= sample_mask(curr_mask, size)
    if not isinstance(mask, np.ndarray):
        mask = mask.to_numpy()
    return mask

def get_multi_panel_fig_and_axes(num_of_fig_panels, nrows=None, ncols=None, figsize=(30,30), **kwargs):
    if nrows is None:
        if ncols is None:
            nrows = int(np.ceil(np.sqrt(num_of_fig_panels)))
            ncols = int(np.ceil(num_of_fig_panels / nrows))
        else:
            nrows = int(np.ceil(num_of_fig_panels / ncols))
    if ncols is None:
        ncols = int(np.ceil(num_of_fig_panels / nrows))
    
    assert nrows * ncols >= num_of_fig_panels, f'num_of_fig_panels: {num_of_fig_panels}, but nrows*ncols={nrows*ncols}'
    
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, 
        figsize=figsize,
        **kwargs,
    )
    flattened_axes = [axes] if num_of_fig_panels == 1 else axes.flatten()
    for ax in flattened_axes[num_of_fig_panels:]:
        make_all_spines_and_x_and_y_axes_invisible(ax)
    return fig, axes

def get_corr_with_num_columns_or_rows_df(vec, df, cols=True):
    import metacells as mc
    vec = mc.ut.to_numpy_vector(vec).astype(np.float32)
    if cols:
        assert len(vec) == len(df)
        names = [x for x in list(df) if ((df[x].dtype.kind in 'biuf') and (df[x].nunique() > 1))]
        mat = df[names].to_numpy().T.astype(np.float32)
    else:
        # we assume all rows are numerical
        mat = df.to_numpy().astype(np.float32)
        names = df.index
    
    corrs = mc.ut.to_numpy_vector(mc.ut.cross_corrcoef_rows(
        mc.ut.to_layout(mc.ut.to_numpy_matrix(vec[np.newaxis,:]), layout='row_major'),
        mc.ut.to_layout(mat, layout='row_major'),
        reproducible=True,
    ))
    
    corr_df = pd.DataFrame({'corr': corrs, 'name': names}).sort_values('corr', ascending=False)
    return corr_df

def get_clustermap_row_and_col_names_ordered(df, clustermap_obj):
    ordered_row_is = clustermap_obj.dendrogram_row.reordered_ind
    ordered_row_names = list(df.index[ordered_row_is])
    ordered_col_is = clustermap_obj.dendrogram_col.reordered_ind
    ordered_col_names = list(df.columns[ordered_col_is])
    return ordered_row_names, ordered_col_names

def norm_cols_and_plot_clustermap(df, metric='correlation', **clustermap_kwargs):
    import seaborn as sb
    
    norm_df = df - df.median(axis=0)
    # return norm_df
    clustermap_obj = sb.clustermap(
        data=norm_df,
        xticklabels=True,
        # yticklabels=True,
        # col_colors=col_color_vecs if col_color_vecs else None,
        # cmap='RdBu_r',
        cmap='bwr',
        vmin=-3,
        vmax=3,
        metric=metric,
        **clustermap_kwargs,
    )
    return clustermap_obj

def jit_scatter(df, x, jit_sd=None, **kwargs):
    import seaborn as sb
    if jit_sd is None:
        print('jit_sd is None. so just passing args to sb.scatterplot.')
        sb.scatterplot(data=df, x=x, **kwargs)
        return
    
    df = df.copy()
    x_jit_col_name = f'{x} (jit)'
    assert x_jit_col_name not in df.columns

    df[x_jit_col_name] = df[x] + np.random.normal(scale=jit_sd, size=len(df))
    
    sb.scatterplot(data=df, x=x_jit_col_name, **kwargs)

def plot_value_counts_heatmap(
        df, x_col_name, y_col_name, norm_by_row_sum=True, ordered_x_vals=None, ordered_y_vals=None, metric='correlation', **clustermap_kwargs):
    import seaborn as sb
    
    df = df[[x_col_name, y_col_name]].value_counts().reset_index().pivot(index=x_col_name, columns=y_col_name, values='count').fillna(0)
    if ordered_x_vals:
        assert set(df.columns) <= set(ordered_x_vals)
        df = df[[x for x in ordered_x_vals if x in df.columns]]
    if ordered_y_vals:
        assert set(df.index) <= set(ordered_y_vals)
        df = df.loc[[x for x in ordered_y_vals if x in df.index]]

    df = df.div(df.sum(axis=(1 if norm_by_row_sum else 0)), axis=('index' if norm_by_row_sum else 'columns'))

    clustermap_obj = sb.clustermap(
        data=df,
        xticklabels=True,
        yticklabels=True,
        row_cluster=(ordered_x_vals is None),
        col_cluster=(ordered_y_vals is None),
        # cmap='RdBu_r',
        # cmap='bwr',
        cmap='rocket_r',
        # vmin=-3,
        # vmax=3,
        metric=metric,
        **clustermap_kwargs,
    )
    return clustermap_obj

def sort_ad_x_to_make_it_canonical(data_ad):
    # metacells/utilities/typing.py, mustbe_canonical() says: "For sparse matrices, it means the data is in COO format, or compressed (CSC or CSR format), with sorted indices and no duplicates"
    import metacells as mc

    orig_row_sums = mc.ut.to_numpy_vector(data_ad.X.sum(axis=1))
    orig_col_sums = mc.ut.to_numpy_vector(data_ad.X.sum(axis=0))
    data_ad.X = data_ad.X.sorted_indices() # IIUC, this is just the indices used in the "internal representation" of the sparse matrix, such that the data in the matrix doesn't really change.
    assert (mc.ut.to_numpy_vector(data_ad.X.sum(axis=1)) == orig_row_sums).all()
    assert (mc.ut.to_numpy_vector(data_ad.X.sum(axis=0)) == orig_col_sums).all()
    assert data_ad.X.has_canonical_format

def get_common_vals_in_sets(sets):
    # i.e., return something in case the sets are not disjoint.
    all_vals = []
    for x in sets:
        all_vals.extend(x)
    count_df = pd.Series(all_vals).value_counts()
    return set(count_df.index[count_df > 1])

def are_sets_disjoint(sets):
    return not get_common_vals_in_sets(sets)

def assert_sets_are_disjoint(sets):
    unexpected_common_vals_in_sets = get_common_vals_in_sets(sets)
    assert not unexpected_common_vals_in_sets, f'unexpected_common_vals_in_sets: {unexpected_common_vals_in_sets}'


def or_sets(sets):
    return set.union(*sets)

def append_if_not_in(l, x):
    if x not in l:
        l.append(x)

def collapse_df_rows(df, consider_nan_str_as_nan=True, allow_conflicts_cols=()):
    # print(f'collapse_df_rows:\n{df}')
    collapsed_row_as_dict = {}
    for col in df.columns:
        nan_val_mask = df[col].isna()
        if consider_nan_str_as_nan:
            nan_val_mask |= df[col] == 'nan'
        if nan_val_mask.all():
            collapsed_row_as_dict[col] = df[col].iloc[0]
        else:
            curr_df = df[~nan_val_mask]
            if curr_df[col].nunique() == 1:
                collapsed_row_as_dict[col] = curr_df[col].iloc[0]
            else:
                if col in allow_conflicts_cols:
                    collapsed_row_as_dict[col] = np.nan
                else:
                    print(f'col: {col}, curr_df[col].unique(): {curr_df[col].unique()}')
                    raise RuntimeError('conflict in rows to collapse')
    collapsed_df = pd.DataFrame([collapsed_row_as_dict])
    # print(f'collapsed_df:\n{collapsed_df}')
    return collapsed_df

def plot_boxplot_per_bin(
        df, x, y, bins, ax=None, xticklabel_round_ndigits=None, 
        showmeans=True,
        boxprops={'facecolor':'none', 'edgecolor':'#BBBBBB'},
        medianprops={'visible': True, 'linewidth': 3},
        whiskerprops={'visible': False},
        meanprops={'marker':'s','markerfacecolor':'red', 'markeredgecolor':'red'},
        boxplot_kwargs=dict(), **stripplot_kwargs,
):
    import seaborn as sb
    
    x_col = x
    y_col = y
    orig_x_col = f'orig_{x_col}'
    df = df.copy().rename(columns={x_col: orig_x_col})
    min_x = df[orig_x_col].min()
    max_x = df[orig_x_col].max()
    if isinstance(bins, int):
        bins = np.linspace(min_x, max_x, bins + 1)
    else:
        bins = np.array(bins)
    # print(bins, min_x, max_x)
    df['bin_i'] = safe_digitize(df[orig_x_col], bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_i_to_center = {i: c for i, c in enumerate(bin_centers)}
    assert set(df['bin_i'].unique()) <= set(bin_i_to_center)
    df[x_col] = df['bin_i'].map(bin_i_to_center)
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    common_kwargs = dict(
        data=df, 
        x=x_col, 
        y=y_col, 
        dodge=False,
        ax=ax,
    )
    # print(get_palette(c_or_mc_ad, color_by=hue_obs_column_name))
    # if ordered_hue_obs_column_names:
    #     common_kwargs['order'] = [x for x in ordered_hue_obs_column_names if x in df[hue_obs_column_name].unique()]
    # else:
    #     common_kwargs['order'] = df.groupby(hue_obs_column_name)[column_name].median().sort_values().index.tolist()
    
    sb.stripplot(
        jitter=0.25,
        linewidth=0.5,
        **common_kwargs,
        **stripplot_kwargs,
    )
    sb.boxplot(
        # fliersize=3,
        # flierprops=dict(alpha=0.5),
        showfliers=False,
        showbox=True,
        showcaps=False,
        showmeans=showmeans,
        medianprops=medianprops,
        meanprops=meanprops,
        whiskerprops=whiskerprops,
        boxprops=boxprops,
        # notch=True,
        **common_kwargs,
        **boxplot_kwargs,
    )

    if xticklabel_round_ndigits is not None:
        print('starting xticklabel_round_ndigits handling')
        sorted_shown_bin_centers = sorted(df[x_col].unique())
        orig_xticklabels = ax.get_xticklabels()
        assert len(sorted_shown_bin_centers) == len(orig_xticklabels)
        orig_labels_as_floats = [float(x.get_text()) for x in orig_xticklabels]
        # print(f'sorted_shown_bin_centers: {sorted_shown_bin_centers}')
        # print(f'orig_labels_as_floats: {orig_labels_as_floats}')
        assert all(np.isclose(x - y, 0) for x, y in zip(sorted_shown_bin_centers, orig_labels_as_floats))
        ax.set_xticklabels(str(round(x, xticklabel_round_ndigits)) for x in sorted_shown_bin_centers)
    
    # ax.get_legend().remove()
    # ax.set_xticklabels(
    #     ax.get_xticklabels(),
    #     rotation=45, ha='right', 
    #     rotation_mode='anchor', 
    #     # fontsize='small',
    # )
    # ax.set_xlabel(None)
    # fig.tight_layout()
    
    return fig, ax

def get_mask_filtered_by_submask(mask, submask):
    filtered_mask = mask.copy()
    false_indices = np.where(mask)[0][~submask]
    filtered_mask[false_indices] = False
    return filtered_mask

def clustermap_zoom(
        clustermap_obj, 
        row_name=None, row_radius=20, first_row_i=None, last_row_i=None,
        col_name=None, col_radius=20, first_col_i=None, last_col_i=None,
        reset_zoom=False,
):
    df = clustermap_obj.data
    ax_heatmap = clustermap_obj.ax_heatmap
    ax_row_dendrogram = clustermap_obj.ax_row_dendrogram
    ax_col_dendrogram = clustermap_obj.ax_col_dendrogram
    
    if hasattr(clustermap_obj, 'orig_min_x'):
        min_x = clustermap_obj.orig_min_x
        max_x = clustermap_obj.orig_max_x
        min_y = clustermap_obj.orig_min_y
        max_y = clustermap_obj.orig_max_y
        row_dendr_min_y = clustermap_obj.orig_row_dendr_min_y
        row_dendr_max_y = clustermap_obj.orig_row_dendr_max_y
        col_dendr_min_x = clustermap_obj.orig_col_dendr_min_x
        col_dendr_max_x = clustermap_obj.orig_col_dendr_max_x
    else:
        min_x, max_x = ax_heatmap.get_xlim()
        max_y, min_y = ax_heatmap.get_ylim()
    
        row_dendr_max_y, row_dendr_min_y = ax_row_dendrogram.get_ylim()
        col_dendr_min_x, col_dendr_max_x = ax_col_dendrogram.get_xlim()

        clustermap_obj.orig_min_x = min_x
        clustermap_obj.orig_max_x = max_x
        clustermap_obj.orig_min_y = min_y
        clustermap_obj.orig_max_y = max_y
        clustermap_obj.orig_row_dendr_min_y = row_dendr_min_y
        clustermap_obj.orig_row_dendr_max_y = row_dendr_max_y
        clustermap_obj.orig_col_dendr_min_x = col_dendr_min_x
        clustermap_obj.orig_col_dendr_max_x = col_dendr_max_x
    
    assert min_x <= max_x
    assert min_y <= max_y
    assert row_dendr_min_y <= row_dendr_max_y
    assert col_dendr_min_x <= col_dendr_max_x
    assert min_x == 0
    assert min_y == 0
    assert row_dendr_min_y == 0
    assert col_dendr_min_x == 0
    print(min_x, max_x)
    print(min_y, max_y)
    print(row_dendr_min_y, row_dendr_max_y)
    print(col_dendr_min_x, col_dendr_max_x)

    row_dendr_y_ratio = row_dendr_max_y / max_y
    col_dendr_x_ratio = col_dendr_max_x / max_x
    print(row_dendr_y_ratio)
    print(col_dendr_x_ratio)
    print(len(df))

    ordered_row_is = clustermap_obj.dendrogram_row.reordered_ind
    ordered_row_names = list(df.index[ordered_row_is])
    ordered_col_is = clustermap_obj.dendrogram_col.reordered_ind
    ordered_col_names = list(df.columns[ordered_col_is])
    
    
    if row_name:
        assert first_row_i is None
        assert last_row_i is None

        row_i = ordered_row_names.index(row_name)
        first_row_i = max(-2, row_i - row_radius)
        last_row_i = min(max_y + 2, row_i + row_radius)
    elif first_row_i and last_row_i:
        assert row_name is None

    if first_row_i and last_row_i:
        ax_heatmap.set_ylim(last_row_i, first_row_i)
        ax_row_dendrogram.set_ylim(last_row_i * row_dendr_y_ratio, first_row_i * row_dendr_y_ratio)
    
    if col_name:
        assert first_col_i is None
        assert last_col_i is None

        col_i = ordered_col_names.index(col_name)
        first_col_i = max(-2, col_i - col_radius)
        last_col_i = min(max_x + 2, col_i + col_radius)
    elif first_col_i and last_col_i:
        assert col_name is None

    if first_col_i and last_col_i:
        ax_heatmap.set_xlim(first_col_i, last_col_i)
        ax_col_dendrogram.set_xlim(first_col_i * col_dendr_x_ratio, last_col_i * col_dendr_x_ratio)
        print('asdf')


    if reset_zoom:
        ax_heatmap.set_ylim(max_y, min_y)
        ax_row_dendrogram.set_ylim(row_dendr_max_y, row_dendr_min_y)
        ax_heatmap.set_xlim(max_x, min_x)
        ax_col_dendrogram.set_xlim(col_dendr_max_x, col_dendr_min_x)

def get_bins(bins, vals, handle_only_int_vals_differently=False):
    if handle_only_int_vals_differently:
        vals_as_series = pd.Series(vals)
        all_vals_are_ints = (vals_as_series == vals_as_series.astype(int)).all()
        print(f'all_vals_are_ints: {all_vals_are_ints}')
    max_val = vals.max()
    min_val = vals.min()

    if bins is None:
        if handle_only_int_vals_differently and all_vals_are_ints:
            bins = np.arange(0, max_val + 2)
        else:
            bins = 20

    if isinstance(bins, int):
        bins = np.linspace(min_val, max_val, bins)
    else:
        bin_size = bins[1] - bins[0]
        assert np.isclose(np.diff(bins), bin_size).all()
        bins = np.array(bins)
    
    if handle_only_int_vals_differently:
        bins_are_consecutive_integers = all_vals_are_ints and np.isclose(1, np.diff(bins)).all()
        if bins_are_consecutive_integers:
            bins = np.array(list(bins - 0.5) + [bins[-1] + 0.5])

    # remove bins higher than max val (no sense in plotting on the right completely empty columns)
    if bins[-1] > max_val:
        first_bin_edge_higher_than_max_val_i = np.where(bins > max_val)[0][0]
        bins = bins[:(first_bin_edge_higher_than_max_val_i+1)]
        assert max_val < bins[-1]

    assert bins[0] < bins[1]
    return bins

def set_heatmap_bin_edge_xticks(
        heatmap_ax, bins,
        bin_edge_to_tick_label_func=lambda x: remove_redundant_trailing_zeros(f'{x:.1f}'),
        min_truncate_val=None,
        max_truncate_val=None,
        dilute_x_labels=False,
        x_axis_target_num_of_strs=10,
):
    heatmap_ax_x_of_ticks = [x.get_position()[0] for x in heatmap_ax.get_xticklabels()]
    half_bin_size_in_heatmap_ax_x = 0.5 * (heatmap_ax_x_of_ticks[1] - heatmap_ax_x_of_ticks[0])
    leftmost_tick_x = heatmap_ax_x_of_ticks[0] - half_bin_size_in_heatmap_ax_x
    rightmost_tick_x = heatmap_ax_x_of_ticks[-1] + half_bin_size_in_heatmap_ax_x
    
    
    new_heatmap_ax_x_of_ticks = np.linspace(leftmost_tick_x, rightmost_tick_x, len(bins))

    new_xtick_labels = [bin_edge_to_tick_label_func(x) for x in bins]
    if min_truncate_val is not None:
        new_xtick_labels[0] = '-inf'
    if max_truncate_val is not None:
        new_xtick_labels[-1] = 'inf'

    if dilute_x_labels and (len(new_xtick_labels) > x_axis_target_num_of_strs):
        new_xtick_labels = dilute_str_list(new_xtick_labels, target_num_of_strs=x_axis_target_num_of_strs)

    heatmap_ax.set_xticks(new_heatmap_ax_x_of_ticks)
    heatmap_ax.get_xaxis().set_major_formatter(matplotlib.ticker.FixedFormatter(new_xtick_labels))

def assert_iter_is_unique(x):
    assert len(x) == len(set(x)), 'the following appear more than once: ' + str(sorted((pd.Series(x).value_counts() > 1).loc[lambda x: x].index))

def set_ax_lim_by_vals(ax, x_or_y, vals, margin=0.05):
    min_val = np.min(vals)
    max_val = np.max(vals)
    val_range = max_val - min_val
    assert x_or_y in {'x', 'y'}
    f = ax.set_xlim if (x_or_y == 'x') else ax.set_ylim
    f(min_val - margin * val_range, max_val + margin * val_range)


def convert_two_dot_two_dot_two_date_col_to_days_since_ref_date(df, col, ref_date_as_y_m_d_tuple, skip_na_vals=False):
    not_na_mask = df[col].notna() & (df[col] != 'nan')
    if not skip_na_vals:
        assert not_na_mask.all()

    result_vec = np.full(len(df), np.nan)
    ref_date = datetime.date(*ref_date_as_y_m_d_tuple)
    splitted_exp_date = df.loc[not_na_mask, col].str.split('.', expand=True).astype(int)
    filtered_result_vec = splitted_exp_date.apply(
        lambda row: (datetime.date(row[2] + (2000 if (row[2] < 40) else 1900), row[1], row[0]) - ref_date).days, axis=1)
    result_vec[not_na_mask] = filtered_result_vec
    return result_vec

def inefficient_downsample_dfs_to_have_same_dist_of_vals_in_col(df1, df2, col, random_state=0):
    assert set(df1.columns) == set(df2.columns)
    
    shuffeled_df1 = df1.sample(frac=1, random_state=random_state)
    shuffeled_df2 = df2.sample(frac=1, random_state=random_state)
    counts1 = df1[col].value_counts()
    counts2 = df2[col].value_counts()

    filtered_df1_indices = []
    filtered_df2_indices = []
    for val in sorted(set(counts1.index) & set(counts2.index)):
        min_count = min(counts1[val], counts2[val])
        filtered_df1_indices.extend(list(shuffeled_df1[shuffeled_df1[col] == val].index[:min_count]))
        filtered_df2_indices.extend(list(shuffeled_df2[shuffeled_df2[col] == val].index[:min_count]))
    
    return df1.loc[filtered_df1_indices], df2.loc[filtered_df2_indices]

def get_py_file_contents_without_regex_removing_comments_with_regex(file_path, compiled_regex, allowed_lines, preprocessing_func=None):
    def assert_regex_not_inside_str_or_end_of_line(suspect_line, line_without_comment):
        assert (line_without_comment.count("'") % 2 == 0) and (line_without_comment.count('"') % 2 == 0), f'# inside str: {suspect_line}'
        assert '\n' not in suspect_line[:-1], f'suspect word at end of line: {suspect_line}'

    
    file_contents = read_text_file(file_path)
    if preprocessing_func is not None:
        file_contents = preprocessing_func(file_contents)

    file_type = file_path.rpartition('.')[-1]
    assert file_type in {'py', 'ipynb'}, f'file_type: {file_type}, but only py and ipynb are supported'
    is_ipynb = file_type == 'ipynb'
    if is_ipynb:
        assert '"outputs": [\n' not in file_contents, f'clear all outputs first: {file_path}'

    new_file_contents = ''
    last_i_in_orig_file = 0
    for m in compiled_regex.finditer(file_contents):
        suspect_line_i = m.start()
        suspect_line = m.group()
        
        if '#' in suspect_line:
            suspect_line_len = len(suspect_line)
            suspect_line_end_i = suspect_line_i + suspect_line_len
            
            if not is_ipynb:
                line_without_comment = suspect_line.partition('#')[0]
                assert_regex_not_inside_str_or_end_of_line(suspect_line, line_without_comment)
                
                new_file_contents += file_contents[last_i_in_orig_file:suspect_line_i] + line_without_comment + '\n'
            
            if is_ipynb:
                suspect_line_stripped = suspect_line.strip()
                assert suspect_line_stripped[0] == '"', f'suspect_line_stripped: {[suspect_line_stripped]}'
                assert suspect_line_stripped.endswith('\\n",'), f'suspect_line_stripped: {[suspect_line_stripped]}'
                bare_suspect_line = suspect_line_stripped[1:-3]
                bare_suspect_line_without_comment = bare_suspect_line.partition('#')[0]
                assert_regex_not_inside_str_or_end_of_line(bare_suspect_line, bare_suspect_line_without_comment)
                
                new_file_contents += file_contents[last_i_in_orig_file:suspect_line_i] + '"' + bare_suspect_line_without_comment + '\\n",' + '\n'
            
            last_i_in_orig_file = suspect_line_end_i
            # raise

    new_file_contents += file_contents[last_i_in_orig_file:]

    unexpected_lines_with_regex = []
    for m in compiled_regex.finditer(new_file_contents):
        suspect_line = m.group()
        if suspect_line in allowed_lines:
            continue
        suspect_line_i = m.start()
        unexpected_lines_with_regex.append(suspect_line)
        # raise

    if unexpected_lines_with_regex:
        print(f'unexpected_lines_with_regex: {unexpected_lines_with_regex}')
        raise
    
    return new_file_contents
