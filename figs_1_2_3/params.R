
# output directories - for figures and tables
SUPP.TABLE.DIR = ''
BASE.FIG.DIR = ''

# paths to directories with all the expression data and metacell models (MODEL.DIR, files stores in AWS),
# and with smaller data fiels (DATA.DIR, files downloaded from github)
DATA.DIR = '/net/mraid20/export/data4/users/nimrodra/proj/blood_aging/data/for_uploading/'
MODEL.DIR = '/net/mraid20/export/data4/users/nimrodra/proj/blood_aging/metacells_workdir/'

BM.MODEL.PATH = file.path(MODEL.DIR, 'bm_hca_ref_metacells.h5ad')
BM.MODEL.ANNOTATION = file.path(DATA.DIR, 'bm_colors.csv')

ORIG.CBC.PATH = file.path(DATA.DIR, 'cbc_original_cohort.csv')
ORIG.SEX.AGE.PATH = file.path(DATA.DIR, 'age_sex_original_cohort.csv')
NEW.SEX.AGE.PATH = file.path(DATA.DIR, 'age_sex_new_cohort.csv')

PREV.CLIN.VALUES.PATH = file.path(DATA.DIR, 'cbc_longitudinal_original_cohort.csv')
NEW.CLIN.VALUES.PATH = file.path(DATA.DIR, 'cbc_longitudinal_new_cohort.csv')

INDIV.ID.MAP.PATH = file.path(DATA.DIR, 'indiv_id_map_df')

N122.GOT.PATH = file.path(DATA.DIR, 'n122_got.Rdata')

TFS.FILE.PATH = file.path(DATA.DIR, 'tfs.txt')

RDW.MUT.DIR = file.path(DATA.DIR, 'rdw_mut')
MUTATION.FILES.DIR = file.path(DATA.DIR, 'mutation_files')

BATCH.EFFECT.PVALS.PATH = file.path(DATA.DIR, 'gene_batch_effect_pvals')
SEX.GENES.PATH = file.path(DATA.DIR, 'sex_genes')

