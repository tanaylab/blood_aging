import os.path

ROOT_DIR_PATH = '/home/orenmil/raid/mds/240801_standalone_blood_aging'

OUTPUT_ROOT_DIR_PATH = os.path.join(ROOT_DIR_PATH, 'output_and_given_intermediate_output')

HUMAN_GENOME_DIR_PATH = os.path.join(ROOT_DIR_PATH, 'input/human_genome_files')
ARCH_MUTATION_DIR_PATH = os.path.join(ROOT_DIR_PATH, 'input/ARCH_mutation_data')
PREPROCESSING_DIR_PATH = os.path.join(ROOT_DIR_PATH, 'input/preprocessing_out')
EXPERIMENT_METADATA_DIR_PATH = os.path.join(ROOT_DIR_PATH, 'input/all_experiments_metadata')
GC_BIAS_ANALYSIS_DIR_PATH = os.path.join(ROOT_DIR_PATH, 'input/gc_bias_analysis')

MISC_INPUT_DIR_PATH = os.path.join(ROOT_DIR_PATH, 'input/misc')


print('mds_in_out_dir_paths was loaded/reloaded')
