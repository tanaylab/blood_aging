import itertools
import os.path
import pickle
import random
import subprocess
import time
import warnings

import Bio
import numpy as np
import pandas as pd
from Bio import pairwise2

with warnings.catch_warnings():
    warnings.simplefilter('ignore',
                          PendingDeprecationWarning)  # PendingDeprecationWarning: We intend to remove or replace Bio.Alphabet in 2020, ideally avoid using it explicitly in your code. Please get in touch if you will be adversely affected by this. https://github.com/biopython/biopython/issues/2046
    from Bio import Seq
    from Bio import SeqUtils
# from Bio import Seq # raises a warning because it imports Bio.Alphabet
from Bio import SeqIO
from Bio import Entrez



OREN_EMAIL_ADDRESS = None # enter your email here


OREN_NCBI_API_KEY = None # enter your ncbi api key here

Entrez.api_key = OREN_NCBI_API_KEY
Entrez.email = OREN_EMAIL_ADDRESS

class EntrezEsearchPhraseNotFound(Exception):
    pass

from Bio import SeqRecord

from generic import generic_utils

# 200624: It seems that this is the default (what I get from np.geterr() after importing numpy): {'divide': 'warn', 'over': 'warn', 'under': 'ignore', 'invalid': 'warn'}
# https://numpy.org/doc/stable/reference/generated/numpy.seterr.html says: Underflow: result so close to zero that some precision was lost.
# so I guess it is Ok that I just ignore underflow problems?
# np.seterr(all='raise')
np.seterr(divide='raise', over='raise', invalid='raise')
pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option('display.expand_frame_repr', False)

DNA_BASES = 'ACGT'
LOWERCASE_DNA_BASES = 'acgt'
NUM_OF_DNA_BASES = len(DNA_BASES)
DNA_BASE_TO_BASE_INDEX = dict(zip(DNA_BASES, range(NUM_OF_DNA_BASES)))
PYRIMIDINES = 'CUT'
PURINES = 'AG'
DNA_BASES_SET = set(DNA_BASES)
DNA_BASES_SET_INCLUDING_LOWERCASE = set(DNA_BASES) | set(LOWERCASE_DNA_BASES)
PYRIMIDINES_SET = set(PYRIMIDINES)
PURINES_SET = set(PURINES)
CODON_LEN = 3
TAXON_RANKS = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
TAXON_RANKS_SET = set(TAXON_RANKS)

BASE_TO_COMPLEMENT_BASE = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N',
}

# this is kind of copied from https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py#L1051
SANGER_SCORE_ASCII_TO_PHRED_SCORE = {
    chr(letter): letter - SeqIO.QualityIO.SANGER_SCORE_OFFSET
    for letter in range(0, 256)
}

# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
BACTERIA_CODON_TABLE_ID = 11  #### HA! http://biopython.org/DIST/docs/tutorial/Tutorial.html does: gene.translate(table="Bacterial"). it didn't work in
BACTERIA_START_CODONS = set(Bio.Data.CodonTable.generic_by_id[BACTERIA_CODON_TABLE_ID].start_codons)  # 'Bacterial' didn't work here...
BACTERIA_STOP_CODONS = set(Bio.Data.CodonTable.generic_by_id[BACTERIA_CODON_TABLE_ID].stop_codons)  # 'Bacterial' didn't work here...


BATCH_SIZE_FOR_ENTREZ_QUERIES = 5000
DOWNLOAD_NUCCORE_ENTRY_FROM_NCBI_NUM_OF_RETRIES = 3

# according to https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation
STANDARD_BUT_NON_ACTG_CHARS_IN_DNA_SEQ = set('URYKMSWBDHVN')
STANDARD_CHARS_IN_DNA_SEQ = set(DNA_BASES) | STANDARD_BUT_NON_ACTG_CHARS_IN_DNA_SEQ

def str_to_seq_record(str_to_convert):
    return SeqRecord.SeqRecord(Seq.Seq(str_to_convert))

def silent_seq_record_to_str(seq_record):
    if not isinstance(seq_record, str):
        return seq_record_to_str(seq_record)
    else:
        return seq_record

def seq_record_to_str(seq_record):
    return str(seq_record.seq)

BACTERIA_STOP_CODONS_INVERTED = {seq_record_to_str(str_to_seq_record(x).reverse_complement()) for x in BACTERIA_STOP_CODONS}

def get_non_ACGT_bases_in_seq_record(seq_record, also_allow_acgt):
    seq_record_as_str = seq_record_to_str(seq_record)
    return set(seq_record_as_str) - (DNA_BASES_SET_INCLUDING_LOWERCASE if also_allow_acgt else DNA_BASES_SET)

def does_seq_record_contain_non_ACGT_bases(seq_record, also_allow_acgt=True):
    return len(get_non_ACGT_bases_in_seq_record(seq_record, also_allow_acgt=also_allow_acgt)) >= 1



def get_num_of_non_ACGT_bases_in_seq_record(seq_record):
    seq_record_as_str = seq_record_to_str(seq_record)
    base_counter = collections.Counter(seq_record_as_str)
    for dna_base in DNA_BASES:
        base_counter.pop(dna_base)
    return sum(base_counter.values())

@generic_utils.execute_if_output_doesnt_exist_already
def cached_write_whether_fasta_file_contain_non_ACGT_bases(
        input_file_path_fasta,
        output_file_path_does_seq_record_contain_non_ACGT_bases,
):
    seq_record = get_chrom_seq_from_single_chrom_fasta_file(input_file_path_fasta)
    seq_record_contains_non_ACGT_bases = does_seq_record_contain_non_ACGT_bases(seq_record)
    generic_utils.write_text_file(output_file_path_does_seq_record_contain_non_ACGT_bases,
                                  '1' if seq_record_contains_non_ACGT_bases else '0')

def write_whether_fasta_file_contain_non_ACGT_bases(
        input_file_path_fasta,
        output_file_path_does_seq_record_contain_non_ACGT_bases,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_whether_fasta_file_contain_non_ACGT_bases(
        input_file_path_fasta=input_file_path_fasta,
        output_file_path_does_seq_record_contain_non_ACGT_bases=output_file_path_does_seq_record_contain_non_ACGT_bases,
    )

def get_random_dna_str(seq_len, random_seed=None):
    if random_seed is not None:
        random.seed(random_seed)
    return ''.join(random.choices(DNA_BASES, k=seq_len))

def get_chrom_seq_by_name(chrom_seqs, chrom_name):
    for chrom_seq in chrom_seqs:
        if chrom_seq.id == chrom_name:
            return chrom_seq

def get_num_of_bases_in_region(region):
    assert region[0] <= region[1]
    return region[1] - region[0] + 1

def get_num_of_phosphodiester_bonds_in_region(region):
    assert region[0] <= region[1]
    assert region[0] - int(region[0]) == 0.5
    assert region[1] - int(region[1]) == 0.5
    return region[1] - region[0] + 1

def get_chrom_seqs_from_fasta_file(ref_genome_fasta_file_path):
    return list(SeqIO.parse(ref_genome_fasta_file_path, 'fasta'))

def get_chrom_seq_info_df_from_fasta_file(ref_genome_fasta_file_path):
    chrom_seq_name_to_seq_desc = get_chrom_seq_name_to_seq_description_and_len_from_fasta_file(ref_genome_fasta_file_path)
    return pd.DataFrame([{'chrom_name': x, 'chrom_desc': y[0], 'chrom_len': y[1]} for x,y in chrom_seq_name_to_seq_desc.items()])

def get_chrom_seq_name_to_seq_description_and_len_from_fasta_file(ref_genome_fasta_file_path):
    return {
        seq.name: (seq.description, len(seq))
        for seq in list(SeqIO.parse(ref_genome_fasta_file_path, 'fasta'))
    }
def get_chrom_seq_name_to_seq_description_from_fasta_file(ref_genome_fasta_file_path):
    return {
        seq.name: seq.description
        for seq in list(SeqIO.parse(ref_genome_fasta_file_path, 'fasta'))
    }

def get_raw_seqs(fastq_file_path):
    return list(SeqIO.parse(fastq_file_path, 'fastq'))

def get_longest_read_in_fastq_file(fastq_file_path):
    return generic_utils.get_len_of_longest_line_in_text_file(fastq_file_path)

def get_read_name_from_read_title_in_fastq(title):
    return title.partition(' ')[0]

def get_names_of_reads_in_fastq_file(fastq_file_path):
    read_names = []
    with open(fastq_file_path) as f:
        for title, seq, qual in SeqIO.QualityIO.FastqGeneralIterator(f):
            read_name = get_read_name_from_read_title_in_fastq(title)
            read_names.append(read_name)
    return read_names

def get_gc_fraction_of_reads_in_fastq_file(fastq_file_path):
    with open(fastq_file_path) as f:
        # for title, seq, qual in SeqIO.QualityIO.FastqGeneralIterator(f):
        return [Bio.SeqUtils.GC(x[1]) / 100 for x in SeqIO.QualityIO.FastqGeneralIterator(f)]

def get_reads_in_fastq_file_as_str_list(fastq_file_path):
    with open(fastq_file_path) as f:
        # for title, seq, qual in SeqIO.QualityIO.FastqGeneralIterator(f):
        return [x[1] for x in SeqIO.QualityIO.FastqGeneralIterator(f)]

@generic_utils.execute_if_output_doesnt_exist_already
def cached_write_names_of_fasta_file_seqs(
        input_file_path_fasta,
        output_file_path_seq_names_pickle,
):
    seqs = get_chrom_seqs_from_fasta_file(input_file_path_fasta)
    seq_names = [seq.name for seq in seqs]
    with open(output_file_path_seq_names_pickle, 'wb') as f:
        pickle.dump(seq_names, f, protocol=4)

def write_names_of_fasta_file_seqs(
        input_file_path_fasta,
        output_file_path_seq_names_pickle,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_names_of_fasta_file_seqs(
        input_file_path_fasta=input_file_path_fasta,
        output_file_path_seq_names_pickle=output_file_path_seq_names_pickle,
    )

@generic_utils.execute_if_output_doesnt_exist_already
def cached_filter_fastq_file(
        input_file_path_fastq,
        names_of_reads_to_keep,
        output_file_path_fastq_with_only_specified_reads,
):
    with open(output_file_path_fastq_with_only_specified_reads, "w") as out_handle:
        with open(input_file_path_fastq) as in_handle:
            for title, seq, qual in SeqIO.QualityIO.FastqGeneralIterator(in_handle):
                read_name = get_read_name_from_read_title_in_fastq(title)
                # print(read_name)
                if read_name in names_of_reads_to_keep:
                    out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

def filter_fastq_file(
        input_file_path_fastq,
        names_of_reads_to_keep,
        output_file_path_fastq_with_only_specified_reads,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_filter_fastq_file(
        input_file_path_fastq=input_file_path_fastq,
        names_of_reads_to_keep=names_of_reads_to_keep,
        output_file_path_fastq_with_only_specified_reads=output_file_path_fastq_with_only_specified_reads,
    )

@generic_utils.execute_if_output_doesnt_exist_already
def cached_filter_fasta_file(
        input_file_path_fasta,
        name_of_read_to_read_region_to_keep,
        output_file_path_fasta_with_only_specified_reads,
):
    if 0:
        # seems like this is extremely slower, at least for long reads.
        # and now this old implementation is outdated, as i moved to use name_of_read_to_read_region_to_keep...
        seqs = get_chrom_seqs_from_fasta_file(input_file_path_fasta)
        seq_name_to_seq = get_chrom_name_to_chrom_seq_dict(seqs)
        write_records_to_fasta_or_gb_file([seq_name_to_seq[x] for x in names_of_reads_to_keep],
                                          output_file_path_fasta_with_only_specified_reads, file_type='fasta')
    else:
        seqs_to_keep = []
        for seq in SeqIO.parse(input_file_path_fasta, 'fasta'):
            if seq.name in name_of_read_to_read_region_to_keep:
                region_to_keep = name_of_read_to_read_region_to_keep[seq.name]
                if region_to_keep is None:
                    seqs_to_keep.append(seq)
                else:
                    seqs_to_keep.append(get_region_in_chrom_seq(seq, *region_to_keep))




        write_records_to_fasta_or_gb_file(seqs_to_keep, output_file_path_fasta_with_only_specified_reads, file_type='fasta')


def filter_fasta_file(
        input_file_path_fasta,
        name_of_read_to_read_region_to_keep,
        output_file_path_fasta_with_only_specified_reads,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_filter_fasta_file(
        input_file_path_fasta=input_file_path_fasta,
        name_of_read_to_read_region_to_keep=name_of_read_to_read_region_to_keep,
        output_file_path_fasta_with_only_specified_reads=output_file_path_fasta_with_only_specified_reads,
    )

def write_clipped_reads_to_fastq_file(fastq_file_path, read_name_to_clipped_read_start_and_end_positions_in_read, clipped_reads_fastq_file_path):
    with open(clipped_reads_fastq_file_path, "w") as out_handle:
        with open(fastq_file_path) as in_handle:
            for title, seq, qual in SeqIO.QualityIO.FastqGeneralIterator(in_handle):
                read_name = get_read_name_from_read_title_in_fastq(title)
                # print(read_name)
                if read_name in read_name_to_clipped_read_start_and_end_positions_in_read:
                    start_pos, end_pos = read_name_to_clipped_read_start_and_end_positions_in_read[read_name]
                    start_index = start_pos - 1
                    stop_index = end_pos
                    out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq[start_index:stop_index], qual[start_index:stop_index]))

@generic_utils.execute_if_output_doesnt_exist_already
def cached_write_read_name_to_len_dict(
        input_file_path_fastq,
        output_file_path_read_name_to_len_pickled,
):
    read_name_to_len = {}
    with open(input_file_path_fastq) as in_handle:
        for title, seq, _ in SeqIO.QualityIO.FastqGeneralIterator(in_handle):
            read_name = get_read_name_from_read_title_in_fastq(title)
            read_name_to_len[read_name] = len(seq)

    with open(output_file_path_read_name_to_len_pickled, 'wb') as f:
        pickle.dump(read_name_to_len, f, protocol=4)


def write_read_name_to_len_dict(
        input_file_path_fastq,
        output_file_path_read_name_to_len_pickled,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_read_name_to_len_dict(
        input_file_path_fastq=input_file_path_fastq,
        output_file_path_read_name_to_len_pickled=output_file_path_read_name_to_len_pickled,
    )


def write_and_get_read_name_to_len_dict(input_file_path_fastq):
    read_name_to_len_pickled_file_path = input_file_path_fastq + '.read_name_to_len.pickle'
    write_read_name_to_len_dict(
        input_file_path_fastq=input_file_path_fastq,
        output_file_path_read_name_to_len_pickled=read_name_to_len_pickled_file_path,
    )
    with open(read_name_to_len_pickled_file_path, 'rb') as f:
        return pickle.load(f)

def get_chrom_name_to_chrom_seq_dict(chrom_seqs):
    chrom_names = {chrom_seq.id for chrom_seq in chrom_seqs}
    return {chrom_name: get_chrom_seq_by_name(chrom_seqs, chrom_name) for chrom_name in chrom_names}


def get_chrom_seq_from_single_chrom_fasta_file(fasta_file_path):
    genome_chrom_seqs = get_chrom_seqs_from_fasta_file(fasta_file_path)
    num_of_chroms = len(genome_chrom_seqs)
    if num_of_chroms != 1:
        raise RuntimeError(f'The input fasta file must contain exactly one sequence.')
    return genome_chrom_seqs[0]

def get_chrom_len_from_single_chrom_fasta_file(fasta_file_path):
    return len(get_chrom_seq_from_single_chrom_fasta_file(fasta_file_path))


@generic_utils.execute_if_output_doesnt_exist_already
def cached_write_whether_fasta_seqs_identical(
        input_file_path_fasta1,
        input_file_path_fasta2,
        output_file_path_are_identical_txt,
):
    seq1 = get_chrom_seq_from_single_chrom_fasta_file(input_file_path_fasta1)
    seq2 = get_chrom_seq_from_single_chrom_fasta_file(input_file_path_fasta2)
    are_identical = (len(seq1) == len(seq2)) and (seq_record_to_str(seq1) == seq_record_to_str(seq2))
    are_identical_str = '1' if are_identical else '0'
    generic_utils.write_text_file(output_file_path_are_identical_txt, are_identical_str)

def write_whether_fasta_seqs_identical(
        input_file_path_fasta1,
        input_file_path_fasta2,
        output_file_path_are_identical_txt,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_whether_fasta_seqs_identical(
        input_file_path_fasta1=input_file_path_fasta1,
        input_file_path_fasta2=input_file_path_fasta2,
        output_file_path_are_identical_txt=output_file_path_are_identical_txt,
    )

def get_chrom_len_of_single_chrom_fasta_file(fasta_file_path):
    chrom_seqs = get_chrom_seqs_from_fasta_file(fasta_file_path)
    assert len(chrom_seqs) == 1
    return len(chrom_seqs[0])


def get_chrom_name_and_len_from_single_chrom_fasta_file(fasta_file_path):
    chrom_seqs = get_chrom_seqs_from_fasta_file(fasta_file_path)
    assert len(chrom_seqs) == 1
    chrom_seq = chrom_seqs[0]
    return (chrom_seq.id, len(chrom_seq))


def write_records_to_fasta_or_gb_file(chrom_seqs, fasta_file_path, file_type='fasta'):
    if file_type == 'gb' and type(chrom_seqs) == list:
        raise NotImplementedError('not implemented yet. why would you want that anyway?')

    with open(fasta_file_path, 'w') as f:
        SeqIO.write(chrom_seqs, f, file_type)


@generic_utils.execute_if_output_doesnt_exist_already
def cached_convert_gb_file_to_fasta_file(
        input_file_path_gb,
        output_file_path_fasta,
):
    SeqIO.convert(input_file_path_gb, "gb", output_file_path_fasta, "fasta")

def convert_gb_file_to_fasta_file(
        input_file_path_gb,
        output_file_path_fasta,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_convert_gb_file_to_fasta_file(
        input_file_path_gb=input_file_path_gb,
        output_file_path_fasta=output_file_path_fasta,
    )

def convert_fastq_file_to_fasta_file(fastq_file_rel_path, fasta_file_rel_path):
    SeqIO.convert(fastq_file_rel_path, "fastq", fasta_file_rel_path, "fasta")


def get_entrez_databases_names():
    with Entrez.einfo() as handle:
        record = Entrez.read(handle)
    assert isinstance(record, dict)
    assert list(record) == ['DbList']
    return record['DbList']


def get_info_about_entrez_database(db_name):
    with Entrez.einfo(db=db_name) as handle:
        record = Entrez.read(handle)
    assert isinstance(record, dict)
    assert list(record) == ['DbInfo']
    db_info = record['DbInfo']
    assert (list(db_info) == ['DbName', 'MenuName', 'Description', 'DbBuild', 'Count', 'LastUpdate', 'FieldList', 'LinkList']) or (
            list(db_info) == ['DbName', 'MenuName', 'Description', 'DbBuild', 'Count', 'FieldList', 'LinkList'])
    return db_info


def run_entrez_esearch_and_return_uids(db_name, search_term, retmax=int(10e6), sort=None, convert_uids_to_python_str_objects=True, verbose=True):
    if verbose:
        print(f'run_entrez_esearch_and_return_uids({db_name}, {search_term})')
    with Entrez.esearch(db=db_name, term=search_term, retmax=retmax, sort=sort) as handle:
        record = Entrez.read(handle)
    assert isinstance(record, dict)
    # print(f'term: {term}')
    # print(f'record: {record}')
    if 'ErrorList' in record:
        err_list = record['ErrorList']
        print(f'ErrorList: {err_list}')
        phrase_not_found_err = err_list['PhraseNotFound']
        if phrase_not_found_err:
            raise EntrezEsearchPhraseNotFound('\n'.join(phrase_not_found_err))

    if 'WarningList' in record:
        warning_list = record['WarningList']
        print(f'WarningList: {warning_list}')
    expected_attributes = ['Count', 'RetMax', 'RetStart', 'IdList', 'TranslationSet', 'TranslationStack', 'QueryTranslation']
    if sorted(record) != sorted(expected_attributes):
        print('got unexpected attributes in entrez esearch result:')
        print(sorted(record))
        print(f'expected - actual: {set(expected_attributes) - set(record)}')
        print(f'actual - expected: {set(record) - set(expected_attributes)}')
        raise RuntimeError('got unexpected attributes in entrez esearch result')

    num_of_search_results = int(record['Count'])
    print(f'term: {search_term}, num_of_search_results: {num_of_search_results}')
    uids = record['IdList']
    num_of_uids = len(uids)
    assert num_of_uids <= min(num_of_search_results, retmax)
    # print(num_of_search_results, num_of_uids)
    if num_of_uids < min(num_of_search_results, retmax):
        raise RuntimeError(
            'It seems that ESearch didnt return all of the search results. Maybe https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_ has the solution...')

    # for k,v in record.items():
    #     print(f'{k}: {v}')
    if convert_uids_to_python_str_objects:
        uids = [str(uid) for uid in uids]
    return uids

@generic_utils.execute_if_output_doesnt_exist_already
def cached_run_entrez_esearch_and_write_uids_with_naive_caching(
        db_name,
        search_term,
        output_file_path_uids_pickle,
        dummy_arg_to_make_caching_mechanism_not_skip_execution,
):
    uids = run_entrez_esearch_and_return_uids(db_name=db_name, search_term=search_term)
    with open(output_file_path_uids_pickle, 'wb') as f:
        pickle.dump(uids, f, protocol=4)

def run_entrez_esearch_and_write_uids_with_naive_caching(
        db_name,
        search_term,
        output_file_path_uids_pickle,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_run_entrez_esearch_and_write_uids_with_naive_caching(
        db_name=db_name,
        search_term=search_term,
        output_file_path_uids_pickle=output_file_path_uids_pickle,
        dummy_arg_to_make_caching_mechanism_not_skip_execution=1,
    )


def run_entrez_func_in_batches_and_return_results(db_name, uids, entrez_func, convert_uids_to_python_str_objects=True, **kwargs_for_entrez_func):
    assert isinstance(uids, list)
    assert entrez_func in (Entrez.esummary, Entrez.efetch, Entrez.elink)
    # print(uids)
    num_of_uids = len(uids)
    if entrez_func == Entrez.elink:
        input_uid_to_linked_uids = {}
    elif entrez_func == Entrez.esummary:
        input_uid_to_summary_record = {}
    else:
        records_list = []

    for start_uid_i in range(0, num_of_uids, BATCH_SIZE_FOR_ENTREZ_QUERIES):
        stop_uid_i = start_uid_i + BATCH_SIZE_FOR_ENTREZ_QUERIES
        print(f'uid_i_range: ({start_uid_i}, {stop_uid_i})')
        batch_uids = uids[start_uid_i:stop_uid_i]
        actual_batch_size = len(batch_uids)
        batch_uids_as_set = set(batch_uids)
        print(f'actual_batch_size: {actual_batch_size}')

        while True:
            if entrez_func == Entrez.esummary:
                with entrez_func(db=db_name, id=','.join(batch_uids), **kwargs_for_entrez_func) as handle:
                    batch_records_list = Entrez.read(handle)
            else:
                with entrez_func(db=db_name, id=batch_uids, **kwargs_for_entrez_func) as handle:
                    batch_records_list = Entrez.read(handle)
                # batch_records_list = list(Entrez.parse(handle))
            # print('batch_records_list')
            # print(batch_records_list)
            if not isinstance(batch_records_list, list):
                print(f'type(batch_records_list): {type(batch_records_list)}')
                if isinstance(batch_records_list, dict):
                    assert len(batch_records_list) == 1
                    doc_summary_set = batch_records_list['DocumentSummarySet']
                    assert list(doc_summary_set) == ['DocumentSummary', 'DbBuild']
                    doc_summary = doc_summary_set['DocumentSummary']
                    # print("doc_summary_set['DbBuild']")
                    # print(doc_summary_set['DbBuild'])
                    # print(doc_summary)
                    # assert isinstance(doc_summary, list)
                    batch_records_list = doc_summary
                    # print(next(iter(batch_records_list.items())))
                else:
                    print(f'batch_records_list: {batch_records_list}')
            assert isinstance(batch_records_list, list)
            # could be bigger, i think. this code is problematic, i think, but i am worried about ruining the MIS pipeline, so I won't touch it for now.
            if len(batch_records_list) == actual_batch_size:
                break

            print(f'len(batch_records_list): {len(batch_records_list)}')
            print('sleeping for a minute before trying again.')
            time.sleep(60)

        assert len(batch_records_list) == actual_batch_size

        if entrez_func == Entrez.elink:
            # print(f'batch_records_list:\n{batch_records_list}')
            for record in batch_records_list:
                assert isinstance(record, dict)

                elink_error = record['ERROR']
                if elink_error:
                    raise RuntimeError(f'Entrez.elink returned an error: {elink_error}, for this record:\n{record}')

                link_set_db = record['LinkSetDb']
                # print(f'link_set_db:\n{link_set_db}')
                # print(f'link_set_db:\n{record}')
                assert len(link_set_db) <= 1
                if link_set_db:
                    links = link_set_db[0]['Link']
                    linked_uids = [link['Id'] for link in links]
                    if convert_uids_to_python_str_objects:
                        linked_uids = [str(uid) for uid in linked_uids]

                    input_uid = record['IdList']
                    assert len(input_uid) == 1
                    input_uid = input_uid[0]
                    if convert_uids_to_python_str_objects:
                        input_uid = str(input_uid)
                    assert input_uid in batch_uids_as_set

                    input_uid_to_linked_uids[input_uid] = linked_uids

        elif entrez_func == Entrez.esummary:
            for record in batch_records_list:
                assert isinstance(record, dict)

                # print(record)
                if 'Id' in record:
                    input_uid = record['Id']
                else:
                    assert False # added on 210815. not sure this is right.
                    # print(f'record: {record}')
                    # print(f'dir(record): {dir(record)}')
                    input_uid = record.attributes['uid']
                if convert_uids_to_python_str_objects:
                    input_uid = str(input_uid)
                assert input_uid in batch_uids_as_set

                input_uid_to_summary_record[input_uid] = record

        else:
            raise NotImplementedError('201011: my current guess is that joining the batch_uids is probably wrong, and that passing as a list is probably the right way')
            # changed on 201011
            records_list.extend(batch_records_list)

    if entrez_func == Entrez.elink:
        return input_uid_to_linked_uids
    elif entrez_func == Entrez.esummary:
        return input_uid_to_summary_record
    else:
        return records_list


def run_entrez_esummary_and_return_input_uid_to_summary_record(db_name, uids):
    input_uid_to_summary_record = run_entrez_func_in_batches_and_return_results(db_name, uids, Entrez.esummary)
    for input_uid, summary_record in input_uid_to_summary_record.items():
        assert isinstance(input_uid, str)
        assert isinstance(summary_record, dict)
        if db_name == 'nuccore':
            assert sorted(summary_record) == ['AccessionVersion', 'Caption', 'Comment', 'CreateDate', 'Extra', 'Flags', 'Gi', 'Id', 'Item',
                                              'Length', 'ReplacedBy', 'Status', 'TaxId', 'Title', 'UpdateDate']
        elif db_name == 'taxonomy':
            assert sorted(summary_record) == ['AkaTaxId', 'CommonName', 'Division', 'Genus', 'Id', 'Item', 'ModificationDate', 'Rank', 'ScientificName', 'Species', 'Status', 'Subsp',
                                      'TaxId']


    return input_uid_to_summary_record

@generic_utils.execute_if_output_doesnt_exist_already
def cached_write_taxon_lineage(
        taxon_uid,
        output_file_path_taxon_lineage_pickle,
        dummy_arg_to_make_caching_mechanism_not_skip_execution,
):
    # taxon_uid_to_efetch_results = run_entrez_efetch_and_return_results('taxonomy', [taxon_uid])
    # assert len(taxon_uid_to_efetch_results) == 1
    # print(taxon_uid_to_efetch_results)
    with Entrez.efetch(db='taxonomy', id=taxon_uid) as handle:
        record = Entrez.read(handle)
    # [{'TaxId': '2589952', 'ScientificName': 'Mesorhizobium sp. B2-3-12', 'ParentTaxId': '325217', 'Rank': 'species', 'Division': 'Bacteria', 'GeneticCode': {'GCId': '11', 'GCName': 'Bacterial, Archaeal and Plant Plastid'}, 'MitoGeneticCode': {'MGCId': '0', 'MGCName': 'Unspecified'}, 'Lineage': 'cellular organisms; Bacteria; Proteobacteria; Alphaproteobacteria; Hyphomicrobiales; Phyllobacteriaceae; Mesorhizobium; unclassified Mesorhizobium',
    # 'LineageEx': [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'}, {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'}, {'TaxId': '1224', 'ScientificName': 'Proteobacteria', 'Rank': 'phylum'}, {'TaxId': '28211', 'ScientificName': 'Alphaproteobacteria', 'Rank': 'class'}, {'TaxId': '356', 'ScientificName': 'Hyphomicrobiales', 'Rank': 'order'}, {'TaxId': '69277', 'ScientificName': 'Phyllobacteriaceae', 'Rank': 'family'}, {'TaxId': '68287', 'ScientificName': 'Mesorhizobium', 'Rank': 'genus'}, {'TaxId': '325217', 'ScientificName': 'unclassified Mesorhizobium', 'Rank': 'no rank'}],
    # 'CreateDate': '2019/06/18 12:58:33', 'UpdateDate': '2019/06/18 12:58:34', 'PubDate': '2019/07/08 19:00:16'}]

    if len(record) != 1:
        # ugh. on 211014 i had an assertion error (before adding a print here - i had assert len(record == 1). after changing the code, there wasn't an assertion failure
        # when i ran again the massive screening again. i guess that ncbi sends you the data in some protocol that allows for very rare mistakes. pretty disturbing.
        # happened again on 211017. oh well.
        print(f'taxonomy record:\n{record}')
        assert False
    result_dict = record[0]

    raw_taxon_lineage_infos = result_dict['LineageEx']
    raw_taxon_lineage_infos.append({
        k: result_dict[k]
        for k in {'TaxId', 'ScientificName', 'Rank'}
    })
    for taxon_lineage_info in raw_taxon_lineage_infos:
        assert set(taxon_lineage_info) == {'TaxId', 'ScientificName', 'Rank'}
        for k,v in taxon_lineage_info.items():
            # print(k, type(k), type(v))
            assert type(k) == str
            assert isinstance(v, str)

    taxon_lineage_infos = [
        {
            k: str(v)
            for k,v in taxon_lineage_info.items()
        }
        for taxon_lineage_info in raw_taxon_lineage_infos
    ]
    # print(taxon_lineage_infos)

    with open(output_file_path_taxon_lineage_pickle, 'wb') as f:
        pickle.dump(taxon_lineage_infos, f, protocol=4)

def write_taxon_lineage(
        taxon_uid,
        output_file_path_taxon_lineage_pickle,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_taxon_lineage(
        taxon_uid=taxon_uid,
        output_file_path_taxon_lineage_pickle=output_file_path_taxon_lineage_pickle,
        dummy_arg_to_make_caching_mechanism_not_skip_execution=1,
    )

def get_rank_to_scientific_name_from_lineage_info(lineage_info):
    # 'LineageEx': [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'}, {'TaxId': '2', 'ScientificName': 'Bacteria', 'Rank': 'superkingdom'}, {'TaxId': '1224', 'ScientificName': 'Proteobacteria', 'Rank': 'phylum'}, {'TaxId': '28211', 'ScientificName': 'Alphaproteobacteria', 'Rank': 'class'}, {'TaxId': '356', 'ScientificName': 'Hyphomicrobiales', 'Rank': 'order'}, {'TaxId': '69277', 'ScientificName': 'Phyllobacteriaceae', 'Rank': 'family'}, {'TaxId': '68287', 'ScientificName': 'Mesorhizobium', 'Rank': 'genus'}, {'TaxId': '325217', 'ScientificName': 'unclassified Mesorhizobium', 'Rank': 'no rank'}],
    rank_to_scientific_name = {}
    for lineage_taxon in lineage_info:
        rank = lineage_taxon['Rank']
        if rank in TAXON_RANKS_SET:
            rank_to_scientific_name[rank] = lineage_taxon['ScientificName']
    return rank_to_scientific_name



def run_entrez_efetch_and_return_results(db_name, uids, **kwargs_for_entrez_func):
    if db_name == 'genome':
        raise RuntimeError('https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly seems to say efetch wasnt implemented '
                           'for the genome db. also, on 200430 it seemed like indeed it isnt implemented. but of course, you could try again...')
    records_list = run_entrez_func_in_batches_and_return_results(db_name, uids, Entrez.efetch, **kwargs_for_entrez_func)
    # for record in records_list:
    #     # assert isinstance(record, str)
    #     print(type(record), record)
    return records_list


def run_entrez_elink_and_return_input_uid_to_linked_uids(db_name, uids, **kwargs_for_entrez_func):
    # This seems very useful. on 200725 i found here a link that didn't appear in 201011_relevant_ncbi_dbs_fields_and_links.txt. oh my.
    # https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html

    input_uid_to_linked_uids = run_entrez_func_in_batches_and_return_results(db_name, uids, Entrez.elink, **kwargs_for_entrez_func)
    for input_uid, linked_uids in input_uid_to_linked_uids.items():
        assert isinstance(input_uid, str)
        assert all(isinstance(linked_uid, str) for linked_uid in linked_uids)
    return input_uid_to_linked_uids


def get_genomes_uids_in_ncbi_nuccore_database_linked_to_given_uids_in_ncbi_genome_database(uids_in_genome_db):
    constituent_genomes_uids = run_entrez_elink_and_return_input_uid_to_linked_uids('nuccore', uids_in_genome_db, dbfrom='genome', LinkName='genome_nuccore')
    return constituent_genomes_uids

    # on 200430, this printed: len(other_genomes_uids): 0, When i ran it on all of the bacteria genomes on 200430, i.e., 26908 genomes.
    # so my best guess is that this linked isn't really used and i should just ignore it.
    # other_genomes_uids = bio_utils.run_entrez_elink_and_return_uids('nuccore', uids_in_genome_db, dbfrom='genome', LinkName='genome_nuccore_samespecies')
    # print(f'len(other_genomes_uids): {len(other_genomes_uids)}')


def run_entrez_efetch_to_get_feature_table_from_ncbi_nuccore_database(nuccore_uid):
    # ft is short for feature table.
    with Entrez.efetch(db='nuccore', id=nuccore_uid, rettype='ft') as handle:
        feature_table = handle.read().strip()
    return feature_table

@generic_utils.execute_if_output_doesnt_exist_already
def cached_download_wgs_assembly_from_ncbi(
        nuccore_uid,
        output_file_type,
        output_file_path_wgs_nuccore_entry,
):
    with generic_utils.timing_context_manager('cached_download_wgs_assembly_from_ncbi'):
        assert output_file_type in {'fasta', 'gb'}

        wgs_gb_file_path = output_file_path_wgs_nuccore_entry
        if output_file_type == 'fasta':
            # Why do this? see https://github.com/kblin/ncbi-acc-download/issues/21#issuecomment-723572502
            wgs_gb_file_path += '.gb'

        # ncbi-acc-download accession
        cmd_line_words = [
            'ncbi-acc-download',
            '--recursive',
            # '--verbose',
            '-o', wgs_gb_file_path,
            '--api-key', OREN_NCBI_API_KEY,
            nuccore_uid,
        ]
        subproc_stdout, subproc_stderr, subproc_ret_code = generic_utils.run_cmd_and_get_stdout_and_stderr(
            cmd_line_words, raise_exception_if_subproc_returned_non_zero=False,
            also_return_return_code=True, verbose=True)
        if subproc_stdout or subproc_stderr or subproc_ret_code:
            raise subprocess.SubprocessError(f'orenmil: Seems like ncbi-acc-download failed.\n'
                                             f'subproc_stdout:\n{subproc_stdout}\n'
                                             f'subproc_stderr:\n{subproc_stderr}\n'
                                             f'subproc_ret_code:\n{subproc_ret_code}\n')
        if not os.path.isfile(wgs_gb_file_path):
            raise subprocess.SubprocessError('orenmil: Seems like ncbi-acc-download failed silently.')

        if output_file_type == 'fasta':
            SeqIO.convert(wgs_gb_file_path, 'genbank', output_file_path_wgs_nuccore_entry, 'fasta')
            os.remove(wgs_gb_file_path)

def download_wgs_assembly_from_ncbi(
        nuccore_uid,
        output_file_type,
        output_file_path_wgs_nuccore_entry,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_download_wgs_assembly_from_ncbi(
        nuccore_uid=nuccore_uid,
        output_file_type=output_file_type,
        output_file_path_wgs_nuccore_entry=output_file_path_wgs_nuccore_entry,
    )


@generic_utils.execute_if_output_doesnt_exist_already
def cached_download_nuccore_entry_from_ncbi(
        nuccore_uid, # should be nuccore_accession, but changing would interfere with the caching, so i leave it this way.
        output_file_type,
        output_file_path_nuccore_entry,
):
    # in case of a download error, maybe we should write an empty file? or just write the stderr and stdout to the output file? e.g., what we currently print...
    with generic_utils.timing_context_manager('cached_download_nuccore_entry_from_ncbi'):
        for num_of_failed_attempts in range(DOWNLOAD_NUCCORE_ENTRY_FROM_NCBI_NUM_OF_RETRIES):
            cmd_line_words = [
                'ncbi-acc-download',
                # '--verbose',
                '--format', output_file_type,
                '-o', output_file_path_nuccore_entry,
                '--api-key', OREN_NCBI_API_KEY,
                nuccore_uid,
            ]
            subproc_stdout, subproc_stderr, subproc_ret_code = generic_utils.run_cmd_and_get_stdout_and_stderr(
                cmd_line_words, raise_exception_if_subproc_returned_non_zero=False, also_return_return_code=True, verbose=True)
            if subproc_stdout or subproc_stderr or subproc_ret_code:
                my_err_str = (f'orenmil: Seems like ncbi-acc-download failed.\n'
                              f'subproc_stdout:\n{subproc_stdout}\n'
                              f'subproc_stderr:\n{subproc_stderr}\n'
                              f'subproc_ret_code:\n{subproc_ret_code}\n')
                print(my_err_str)
                if num_of_failed_attempts == DOWNLOAD_NUCCORE_ENTRY_FROM_NCBI_NUM_OF_RETRIES - 1:
                    raise subprocess.SubprocessError(my_err_str)
            elif not os.path.isfile(output_file_path_nuccore_entry):
                raise subprocess.SubprocessError('orenmil: Seems like ncbi-acc-download failed silently.')
            else:
                return
            time.sleep(60) # Stop annoying the server for 1 minute.


        assert False

def download_nuccore_entry_from_ncbi(
        nuccore_accession,
        output_file_type,
        output_file_path_nuccore_entry,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_download_nuccore_entry_from_ncbi(
        nuccore_uid=nuccore_accession, # should be nuccore_accession, but changing would interfere with the caching, so i leave it this way.
        output_file_type=output_file_type,
        output_file_path_nuccore_entry=output_file_path_nuccore_entry,
    )

@generic_utils.execute_if_output_doesnt_exist_already
def cached_extract_sra_file_to_fasta_with_predetermined_path(
        input_file_path_sra,
        output_file_path_fasta,
):
    assert output_file_path_fasta == f'{input_file_path_sra}_pass.fasta'

    # fastq-dump --outdir s6_output_of_massive_screening/sra_entries --fasta 80 --skip-technical --readids
    #    --read-filter pass --dumpbase --split-spot --clip s6_output_of_massive_screening/sra_entries/ERR768077.1
    out_dir_path = os.path.dirname(input_file_path_sra)
    generic_utils.run_cmd_and_check_ret_code_and_return_stdout(
        [
            'fastq-dump',
            '--outdir', out_dir_path,
            '--fasta', '80',
            '--skip-technical',
            '--readids',
            '--read-filter', 'pass',
            '--dumpbase',
            '--split-spot',
            '--clip',
            input_file_path_sra,
        ],
        verbose=True,
    )

def extract_sra_file_to_fasta_with_predetermined_path(
        input_file_path_sra,
        output_file_path_fasta,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_extract_sra_file_to_fasta_with_predetermined_path(
        input_file_path_sra=input_file_path_sra,
        output_file_path_fasta=output_file_path_fasta,
    )

@generic_utils.execute_if_output_doesnt_exist_already
def cached_extract_paired_end_sra_file_to_fastq_with_predetermined_paths(
        input_file_path_sra,
        output_file_path_fastq1,
        output_file_path_fastq2,
        # output_file_path_fastq_orphan,
):
    assert output_file_path_fastq1 == f'{input_file_path_sra}_pass_1.fastq'
    assert output_file_path_fastq2 == f'{input_file_path_sra}_pass_2.fastq'
    # assert output_file_path_fastq_orphan == f'{input_file_path_sra}_pass_orphan.fastq'

    # according to recommended usage in https://edwards.flinders.edu.au/fastq-dump/:
    # fastq-dump --outdir fastq --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip $SRR

    out_dir_path = os.path.dirname(input_file_path_sra)
    generic_utils.run_cmd_and_check_ret_code_and_return_stdout(
        [
            'fastq-dump',
            '--outdir', out_dir_path,
            '--skip-technical',
            '--readids',
            '--read-filter', 'pass',
            '--dumpbase',
            '--split-3',
            '--clip',
            input_file_path_sra,
        ],
        verbose=True,
    )

def extract_paired_end_sra_file_to_fastq_with_predetermined_paths(
        input_file_path_sra,
        output_file_path_fastq1,
        output_file_path_fastq2,
        # output_file_path_fastq_orphan,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_extract_paired_end_sra_file_to_fastq_with_predetermined_paths(
        input_file_path_sra=input_file_path_sra,
        output_file_path_fastq1=output_file_path_fastq1,
        output_file_path_fastq2=output_file_path_fastq2,
        # output_file_path_fastq_orphan=output_file_path_fastq_orphan,
    )

def get_gb_record(gb_file_path):
    return SeqIO.read(gb_file_path, 'gb')

def get_gb_records(gb_file_path):
    return SeqIO.parse(gb_file_path, 'gb')

def get_ncbi_assembly_accession_to_assembly_uid_and_summary_record(assembly_accessions):
    search_term = ' OR '.join(f'"{acc}"[Assembly Accession]' for acc in assembly_accessions)
    assembly_uids = run_entrez_esearch_and_return_uids('assembly', search_term)
    assert len(assembly_uids) == len(assembly_accessions)

    assembly_uid_to_summary_record = run_entrez_esummary_and_return_input_uid_to_summary_record(
        db_name='assembly',
        uids=assembly_uids,
    )
    assert set(assembly_uid_to_summary_record) == set(assembly_uids)

    assembly_accession_to_assembly_uid_and_summary_record = {}
    for assembly_uid, summary_record in assembly_uid_to_summary_record.items():
        assembly_accession = summary_record['AssemblyAccession']
        assert assembly_accession not in assembly_accession_to_assembly_uid_and_summary_record
        assembly_accession_to_assembly_uid_and_summary_record[assembly_accession] = {
            'assembly_uid': assembly_uid,
            'summary_record': summary_record,
        }

    assert set(assembly_accession_to_assembly_uid_and_summary_record) == set(assembly_accessions)
    return assembly_accession_to_assembly_uid_and_summary_record

def get_transition_or_transversion(orig_base, new_base):
    assert orig_base != new_base
    both_bases_as_set = {orig_base, new_base}
    assert both_bases_as_set.issubset(DNA_BASES_SET)

    if both_bases_as_set.issubset(PURINES_SET) or both_bases_as_set.issubset(PYRIMIDINES_SET):
        return 'transition'
    return 'transversion'


def get_reverse_complement(dna_str):
    # print(type(dna_str))
    if isinstance(dna_str, SeqIO.SeqRecord):
        dna_str = str(dna_str.seq)
    return str(Bio.Seq.Seq(dna_str).reverse_complement())


def get_read_phred_scores_from_sanger_scores_str(sanger_scores_str):
    return np.array([SANGER_SCORE_ASCII_TO_PHRED_SCORE[c] for c in sanger_scores_str])




def print_seq_record_as_str(seq_record):
    print(seq_record_to_str(seq_record))



def get_alignemnt_score_and_begin_and_end(alignment):
    # from https://biopython.org/DIST/docs/api/Bio.pairwise2-pysrc.html#_recover_alignments
    # 902                  tracebacks.append((ali_seqB[::-1], ali_seqA[::-1], score,
    # 903                                     begin, end))
    alignment_score, alignment_begin, alignment_end = alignment[-3:]
    return alignment_score, alignment_begin, alignment_end


def get_simulated_circular_chrom_seq(chrom_seq):
    simulated_circular_chrom_seq = chrom_seq + chrom_seq
    return simulated_circular_chrom_seq


def get_region_in_chrom_seq(chrom_seq, start_position, end_position, is_chrom_circular=False, region_name=None, truncate_non_circular_chrom_if_needed=False,
                            return_also_start_index_of_returned_region=False):
    chrom_len = len(chrom_seq)
    after_ori_seq = None
    before_ori_seq = None
    start_position_without_before_ori = start_position
    end_position_without_after_ori = end_position
    if end_position > chrom_len:
        if is_chrom_circular:
            num_of_bases_after_ori = end_position - chrom_len
            after_ori_seq = chrom_seq[:num_of_bases_after_ori]
        else:
            if not truncate_non_circular_chrom_if_needed:
                raise RuntimeError(f'received a non-circular chromosome and also an end position ({end_position}) above chrom_len ({chrom_len}).')
        end_position_without_after_ori = chrom_len
    if start_position < 1:
        if is_chrom_circular:
            num_of_bases_before_ori = 1 - start_position
            before_ori_seq = chrom_seq[-num_of_bases_before_ori:]
        else:
            if not truncate_non_circular_chrom_if_needed:
                raise RuntimeError(f'received a non-circular chromosome and also a start position ({start_position}) below 1.')
        start_position_without_before_ori = 1

    start_index_without_before_ori = start_position_without_before_ori - 1
    region_in_chrom_seq = chrom_seq[start_index_without_before_ori:end_position_without_after_ori]
    if before_ori_seq:
        region_in_chrom_seq = before_ori_seq + region_in_chrom_seq
    if after_ori_seq:
        region_in_chrom_seq = region_in_chrom_seq + after_ori_seq

    expected_region_len_if_not_truncated = end_position - start_position + 1
    region_in_chrom_len = len(region_in_chrom_seq)
    if region_in_chrom_len > chrom_len:
        raise RuntimeError(f'The resulted region_in_chrom_len is too long ({region_in_chrom_len} bases). Such region doesnt really exist.')
    assert (region_in_chrom_len == expected_region_len_if_not_truncated) or ((not is_chrom_circular) and truncate_non_circular_chrom_if_needed and
                                                                                  ((start_position < 1) or (end_position > chrom_len)))

    if region_name:
        region_in_chrom_seq.name = region_in_chrom_seq.description = region_in_chrom_seq.id = region_name

    if return_also_start_index_of_returned_region:
        start_index_of_returned_region = (-len(before_ori_seq)) if before_ori_seq else start_index_without_before_ori
        return (region_in_chrom_seq, start_index_of_returned_region)

    return region_in_chrom_seq

@generic_utils.execute_if_output_doesnt_exist_already
def cached_write_region_to_fasta_file(
        input_file_path_fasta,
        region,
        output_file_path_region_fasta,
        write_reverse_complement_of_region,
        seq_name,
        region_seq_name,
        dummy_arg_to_make_caching_mechanism_not_skip_execution,
):
    if seq_name:
        chrom_seqs = get_chrom_seqs_from_fasta_file(input_file_path_fasta)
        seq = get_chrom_seq_by_name(chrom_seqs, seq_name)
    else:
        seq = get_chrom_seq_from_single_chrom_fasta_file(input_file_path_fasta)

    assert region[0] >= 1
    region_seq = seq[(region[0] - 1):region[1]]


    if write_reverse_complement_of_region:
        region_seq = region_seq.reverse_complement()
    if not region_seq_name:
        region_seq_name = f'{seq.id}_region_{region[0]}_{region[1]}' + ('_reverse_complement' if write_reverse_complement_of_region else '')

    region_seq.name = region_seq.description = region_seq.id = region_seq_name

    write_records_to_fasta_or_gb_file([region_seq], output_file_path_region_fasta, file_type='fasta')


def write_region_to_fasta_file(
        input_file_path_fasta,
        region,
        output_file_path_region_fasta,
        write_reverse_complement_of_region=False,
        seq_name=None,
        region_seq_name=None,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_region_to_fasta_file(
        input_file_path_fasta=input_file_path_fasta,
        region=region,
        output_file_path_region_fasta=output_file_path_region_fasta,
        write_reverse_complement_of_region=write_reverse_complement_of_region,
        seq_name=seq_name,
        region_seq_name=region_seq_name,
        dummy_arg_to_make_caching_mechanism_not_skip_execution=1,
    )

def print_region_in_fasta(fasta_file_path, region_range_as_positions=None, do_print=True):
    chrom_seqs = get_chrom_seqs_from_fasta_file(fasta_file_path)
    region_as_seq = None
    for chrom_seq in chrom_seqs:
        chrom_len = len(chrom_seq)
        if region_range_as_positions:
            start_position, end_position = region_range_as_positions
            if end_position > chrom_len:
                if do_print:
                    print(f'not printing region for {chrom_seq.id}, as the end position of the range specified is larger than chrom_len ({chrom_len}).')
                continue
            region_as_seq = chrom_seq[start_position - 1:end_position]
        else:
            region_as_seq = chrom_seq
        if do_print:
            print(f'Region in {chrom_seq.id}:')

            print_seq_record_as_str(region_as_seq)

    return region_as_seq


@generic_utils.execute_if_output_doesnt_exist_already
def cached_extract_taxon_uid_from_genbank_file(
        input_file_path_gb,
        output_file_path_taxon_uid,
):
    gb_record = get_gb_record(input_file_path_gb)
    seq_features = gb_record.features
    source_feature = seq_features[0]
    assert source_feature.type == 'source'
    # print(source_feature)
    db_xref_qualifier = source_feature.qualifiers['db_xref']
    assert len(db_xref_qualifier) == 1
    db_xref_qualifier = db_xref_qualifier[0].strip()
    assert db_xref_qualifier.startswith('taxon:')
    taxon_uid = int(db_xref_qualifier.partition('taxon:')[-1])
    # print(db_xref_qualifier, taxon_uid)

    generic_utils.write_text_file(output_file_path_taxon_uid, str(taxon_uid))

def extract_taxon_uid_from_genbank_file(
        input_file_path_gb,
        output_file_path_taxon_uid,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_extract_taxon_uid_from_genbank_file(
        input_file_path_gb=input_file_path_gb,
        output_file_path_taxon_uid=output_file_path_taxon_uid,
    )

def get_chrom_name_and_is_chrom_circular_from_files_otherwise_from_gb_file(
        is_chrom_circular_file_path,
        chrom_name_file_path,
        gb_file_path,
):
    if os.path.isfile(is_chrom_circular_file_path):
        is_chrom_circular_file_contents = generic_utils.read_text_file(is_chrom_circular_file_path)
        assert is_chrom_circular_file_contents in ('True', 'False')
        is_chrom_circular = bool(is_chrom_circular_file_contents)

        assert os.path.isfile(chrom_name_file_path)
        chrom_name = generic_utils.read_text_file(chrom_name_file_path)
    else:
        gb_record = get_gb_record(gb_file_path)
        chrom_topology = gb_record.annotations['topology']
        is_chrom_circular = chrom_topology == 'circular'
        if not is_chrom_circular:
            assert chrom_topology == 'linear'
        chrom_name = gb_record.id

        generic_utils.write_text_file(is_chrom_circular_file_path, str(is_chrom_circular))
        generic_utils.write_text_file(chrom_name_file_path, chrom_name)

    return chrom_name, is_chrom_circular


def get_all_distances_between_subseq_and_its_reverse_complement(seq_str, subseq_str):
    assert isinstance(subseq_str, str)
    assert isinstance(seq_str, str)
    subseq_len = len(subseq_str)
    subseq_rc_str = get_reverse_complement(subseq_str)
    subseq_indices = generic_utils.get_substr_indices_including_overlapping(seq_str, subseq_str)
    subseq_rc_indices = generic_utils.get_substr_indices_including_overlapping(seq_str, subseq_rc_str)
    return [abs(subseq_i - subseq_rc_i) - subseq_len for subseq_i, subseq_rc_i in itertools.product(subseq_indices, subseq_rc_indices)]


def does_subseq_appear_close_to_its_reverse_complement(seq_str, subseq_str, max_min_dist):
    assert isinstance(subseq_str, str)
    assert isinstance(seq_str, str)

    min_dist = min(get_all_distances_between_subseq_and_its_reverse_complement(seq_str, subseq_str), default=np.inf)
    ret_val = min_dist <= max_min_dist
    # if ret_val:
    #     print(seq_str, subseq_str, min_dist)
    return ret_val


@generic_utils.execute_if_output_doesnt_exist_already
def cached_find_genus_of_each_taxa_and_write_genus_to_taxon_uids_pickle(
        taxon_uids,
        output_file_path_genus_to_taxon_uids_pickle,
):
    with generic_utils.timing_context_manager('cached_find_genus_of_each_taxa_and_write_genus_to_taxon_uids_pickle'):
        assert all((type(taxon_uid) == str) for taxon_uid in taxon_uids)

        taxon_uid_to_summary_record = run_entrez_esummary_and_return_input_uid_to_summary_record(
            db_name='taxonomy',
            uids=taxon_uids,
        )

        genus_to_taxon_uids = collections.defaultdict(set)
        for taxon_uid, summary_record in taxon_uid_to_summary_record.items():
            genus = taxon_uid_to_summary_record[taxon_uid]['Genus']
            genus = str(genus)
            genus_to_taxon_uids[genus].add(taxon_uid)
        genus_to_taxon_uids = dict(genus_to_taxon_uids) # I don't want a defaultdict moving around.


        with open(output_file_path_genus_to_taxon_uids_pickle, 'wb') as f:
            pickle.dump(genus_to_taxon_uids, f, protocol=4)


def find_genus_of_each_taxa_and_write_genus_to_taxon_uids_pickle(
        taxon_uids,
        output_file_path_relevant_genus_to_taxon_uids_pickle,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_find_genus_of_each_taxa_and_write_genus_to_taxon_uids_pickle(
        taxon_uids=taxon_uids,
        output_file_path_genus_to_taxon_uids_pickle=output_file_path_relevant_genus_to_taxon_uids_pickle,
    )


def get_dynamic_programming_global_alignments(
        seq1_as_str,
        seq2_as_str,
        match_score=1,
        mismatch_score=0,
        gap_open_score=-1,
        gap_extension_score=-0.5,
        penalize_gaps_on_edges=True,
        num_of_alignments_to_return=1,
        get_score_only=False,
):
    one_alignment_only = num_of_alignments_to_return == 1
    if get_score_only:
        assert one_alignment_only # seems like pairwise2.align.globalms doesn't support another option. Kinds of makes sense.

    alignment_results = pairwise2.align.globalms(
        seq1_as_str,
        seq2_as_str,
        match_score,
        mismatch_score,
        gap_open_score,
        gap_extension_score,
        penalize_end_gaps=penalize_gaps_on_edges,
        one_alignment_only=one_alignment_only,
        score_only=get_score_only,
    )
    if get_score_only:
        return alignment_results

    alignment_results = alignment_results[:num_of_alignments_to_return]

    for alignment in alignment_results:
        # print(pairwise2.format_alignment(*alignment))
        pass

    return alignment_results

def get_matching_positions_df_for_biopython_pairwise2_alignment(alignment, num_of_aligned_bases_of_seqA_for_verification, num_of_aligned_bases_of_seqB_for_verification):
    matching_positions_df_rows = []
    curr_position_in_seqA = 1
    curr_position_in_seqB = 1
    for seqA_base, seqB_base in zip(alignment.seqA, alignment.seqB):
        if seqA_base == seqB_base:
            assert seqA_base != '-'
            matching_positions_df_rows.append({'seqA_position': curr_position_in_seqA,
                                               'seqB_position': curr_position_in_seqB})
        if seqA_base != '-':
            curr_position_in_seqA += 1

        if seqB_base != '-':
            curr_position_in_seqB += 1

    assert curr_position_in_seqA == num_of_aligned_bases_of_seqA_for_verification + 1
    assert curr_position_in_seqB == num_of_aligned_bases_of_seqB_for_verification + 1

    matching_positions_df = pd.DataFrame(matching_positions_df_rows)
    return matching_positions_df




@generic_utils.execute_if_output_doesnt_exist_already
def cached_write_optimal_dynamic_programming_global_alignment_info(
        input_file_path_seq1_fasta,
        input_file_path_seq2_fasta,
        take_reverse_complement_of_seq2,
        output_file_path_best_global_alignment_info_pickle,
        seq1_region,
        seq2_region,
        penalize_gaps_on_edges,
        match_score,
        mismatch_score,
        gap_open_score,
        gap_extension_score,
):
    seq1_as_str = seq_record_to_str(get_chrom_seq_from_single_chrom_fasta_file(input_file_path_seq1_fasta))
    if seq1_region:
        seq1_as_str = seq1_as_str[(seq1_region[0] - 1):seq1_region[1]]

    seq2_as_str = seq_record_to_str(get_chrom_seq_from_single_chrom_fasta_file(input_file_path_seq2_fasta))
    if seq2_region:
        seq2_as_str = seq2_as_str[(seq2_region[0] - 1):seq2_region[1]]

    best_global_alignment = get_dynamic_programming_global_alignments(
        seq1_as_str=seq1_as_str,
        seq2_as_str=(get_reverse_complement(seq2_as_str) if take_reverse_complement_of_seq2 else seq2_as_str),
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extension_score=gap_extension_score,
        penalize_gaps_on_edges=penalize_gaps_on_edges,
        num_of_alignments_to_return=1,
        get_score_only=False,
    )[0]

    num_of_matches = len(best_global_alignment.seqA) - generic_utils.get_hamming_dist_between_same_len_strs(best_global_alignment.seqA, best_global_alignment.seqB)
    num_of_indels_in_seq1 = best_global_alignment.seqA.count('-')
    num_of_indels_in_seq2 = best_global_alignment.seqB.count('-')

    best_global_alignment_info = {
        'score': best_global_alignment.score,
        'num_of_matches': num_of_matches,
        'num_of_indels_in_seq1': num_of_indels_in_seq1,
        'num_of_indels_in_seq2': num_of_indels_in_seq2,
    }
    print(f'best_global_alignment_info: {best_global_alignment_info}')
    with open(output_file_path_best_global_alignment_info_pickle, 'wb') as f:
        pickle.dump(best_global_alignment_info, f, protocol=4)

def write_optimal_dynamic_programming_global_alignment_info(
        input_file_path_seq1_fasta,
        input_file_path_seq2_fasta,
        take_reverse_complement_of_seq2,
        output_file_path_best_global_alignment_info_pickle,
        seq1_region=None,
        seq2_region=None,
        penalize_gaps_on_edges=True,
        match_score=1,
        mismatch_score=0,
        gap_open_score=-1,
        gap_extension_score=-0.5,
):
    # the extra level is needed so that the cached functions will always have all arguments.
    # otherwise, if some of the optional arguments aren't specified, my caching algorithm would think the arguments are
    # different, and run the function again.
    return cached_write_optimal_dynamic_programming_global_alignment_info(
        input_file_path_seq1_fasta=input_file_path_seq1_fasta,
        input_file_path_seq2_fasta=input_file_path_seq2_fasta,
        take_reverse_complement_of_seq2=take_reverse_complement_of_seq2,
        output_file_path_best_global_alignment_info_pickle=output_file_path_best_global_alignment_info_pickle,
        seq1_region=seq1_region,
        seq2_region=seq2_region,
        penalize_gaps_on_edges=penalize_gaps_on_edges,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extension_score=gap_extension_score,
    )

def get_num_of_bases_between_regions(region1, region2):
    assert generic_utils.is_interval(region1)
    assert generic_utils.is_interval(region2)

    region1_len = region1[1] - region1[0] + 1
    region2_len = region2[1] - region2[0] + 1
    region_union_len = max(*region1, *region2) - min(*region1, *region2) + 1 # used * just to prevent confusion for people used to elementwise vector operations...
    num_of_bases_between_regions = max(0, region_union_len - region1_len - region2_len)

    return num_of_bases_between_regions

def get_gc_content(seq):
    if not isinstance(seq, str):
        seq_as_str = seq_record_to_str(seq)
    else:
        seq_as_str = seq
    return SeqUtils.GC(seq_as_str)

def get_min_window_gc_content(seq, window_len):
    seq_as_str = silent_seq_record_to_str(seq)
    seq_as_str_upper = seq_as_str.upper()
    if not (set(seq_as_str_upper) <= DNA_BASES_SET):
        return np.nan

    seq_as_arr = np.array([1 if x in 'GC' else 0 for x in seq_as_str_upper])
    num_of_gcs_per_window = np.convolve(seq_as_arr ,np.ones(window_len, dtype=int) ,'valid')
    return num_of_gcs_per_window.min() / window_len * 100

MISMATCHED_BASES_TO_NUM_OF_MISMATCHES_TEMPLATE = {tuple(sorted(x)):0 for x in itertools.combinations('ACGT', 2)}
def get_mismatched_bases_to_num_of_mismatches(seq1, seq2):
    seq1_as_str = silent_seq_record_to_str(seq1)
    seq2_as_str = silent_seq_record_to_str(seq2)
    assert len(seq1_as_str) == len(seq2_as_str)

    mismatched_bases_to_num_of_mismatches = dict(MISMATCHED_BASES_TO_NUM_OF_MISMATCHES_TEMPLATE)
    for base1, base2 in zip(seq1_as_str, seq2_as_str):
        if base1 != base2:
            mismatched_bases_to_num_of_mismatches[tuple(sorted((base1, base2)))] += 1
    return mismatched_bases_to_num_of_mismatches

def does_cds_encode_a_peptide(seq, codon_table_id):
    # bio_utils.str_to_seq_record('ATGGCTGA').translate(cds=True)
    try:
        peptide = seq.translate(table=codon_table_id, cds=True)
    except Bio.Data.CodonTable.TranslationError:
        return False

    assert len(seq) % 3 == 0
    assert len(peptide) == len(seq) // 3 - 1
    return True
