import os.path 
from pathlib import Path

VERSION_NUM = 1

NCBI_ID = 'NCBI Gene ID'
ENS_ID = 'Ensembl ID'
ZFIN_ID = 'ZFIN ID'
SYMBOL = 'Symbol'
HUMAN_ID = 'Human NCBI Gene ID'

FILE_DIR = Path(os.path.dirname(os.path.realpath(__file__)))
DATABASE_DIR = FILE_DIR / Path('database/')
RAW_DATA_DIR = FILE_DIR / Path('database/raw_data/')
KEGG_DATA_DIR = FILE_DIR / Path('database/KEGG/')
MASTER_MAPPING_FILE_PATH = DATABASE_DIR / Path('master_gene_mapping_file_V' + str(VERSION_NUM) + '.txt')
MASTER_ORTHO_FILE_PATH = DATABASE_DIR / Path('master_ortho_mapping_file_V' + str(VERSION_NUM) + '.txt')

