'''
This script contains all the tools necessary to map zebrafish
gene ids from the varying identification types. The mapping 
data has been downloaded directly from NCBI and ZFIN. Check
out the database README for more information. 
'''

from settings import *

class DatabaseNotFoundError(Exception):
    "Raised when the the database directory is not found"
    pass

class InvalidGeneTypeError(Exception):
    "Raised when the the database directory is not found"
    pass

# GENE MAPPING FUNCTIONS
# ----------------------

def convert_ids(gene_list: any, id_from: str, id_to: str, keep_mapping = False, out_format = None) -> pd.Series or pd.DataFrame:
    '''
    Convert a list of Gene IDs.
    Parameters: 
        gene_list - a list of Gene IDs with supported formats list, pd.Series,
                    pd.DataFrame and np.array
        id_from - the current Gene ID type
        id_to   - the Gene ID type to convert to
    Other:
        Gene ID Type Options: NCBI Gene ID, ZFIN ID, Symbol, Ensembl ID
    '''

    # last updated: July 7, 2023
    # last checked: July 7, 2023.]
    # tutorial? yes

    try: 
        # some error handling
        # -------------------
        _check_valid_zebrafish_gene_id_type([id_from, id_to])
        gene_list = _make_sure_is_pandas_series(gene_list, id_from)
        # just in case the NCBI Gene IDs are not strings
        # FLAG ? 
        gene_list = gene_list.astype(str) if (id_from == NCBI_ID) and (gene_list.dtype == int) else gene_list
        gene_list = gene_list.drop_duplicates()

        # read in the master mapping file
        master_mapping = pd.read_csv(FILE_DIR / MASTER_MAPPING_FILE_PATH, sep = '\t')
        master_mapping[NCBI_ID] = master_mapping[NCBI_ID].to_numpy()
        # master_mapping[NCBI_ID] = master_mapping[NCBI_ID].astype(str)
        # master_mapping[NCBI_ID] 
        df_desired_columns = master_mapping[[id_from, id_to]]

        # get the rows where our ids come from
        filtered_df = df_desired_columns[df_desired_columns[id_from].isin(gene_list)]
        # reset index
        filtered_df = filtered_df.reset_index(drop=True)

        if keep_mapping:
            mapped_genes = filtered_df
            return mapped_genes.drop_duplicates()
        else:
            mapped_genes = filtered_df[id_to]
            mapped_genes = mapped_genes.dropna()
            mapped_genes = mapped_genes.drop_duplicates()

            if out_format == list:
                mapped_genes = mapped_genes.to_list()
            return mapped_genes
        

    
    except InvalidGeneTypeError:
        pass

def add_mapped_column(data: pd.DataFrame, id_from: str, id_to: str, 
                      column_name_with_ids = None, keep_old_ids = True,
                      drop_na = False):
    '''
    Parameters:
        data - a pandas DataFrame containing a column that has Gene IDs of some type
        id_from - the current Gene ID type must be: NCBI Gene ID, ZFIN ID, 
                Ensembl ID, or Symbol
        id_to - the Gene ID type to convert to, must be: NCBI Gene ID, ZFIN ID,
                Ensembl ID, or Symbol 
        column_name_with_ids - the name of the column containing the Gene IDs (if
                               the column name does not match id_from)
        keep_old_ids - if you would like to keep the old Gene ID column
        drop_na - if you would like to drop any rows that have a NA in the resulting
                  mapped column
    '''

    # last updated: July 13, 2023
    # last checked: July 13, 2023
    # tutorial? yes

    try:
        # some error handling
        # -------------------    
        if type(data) != pd.DataFrame:
            raise TypeError
        _check_valid_zebrafish_gene_id_type([id_from, id_to])
        _check_column_name_matches_id_choice(data, id_from, column_name_with_ids)

        if column_name_with_ids:
            # rename data column to the id in options
            data = data.rename(columns={column_name_with_ids:id_from})

        gene_list = data[id_from]
        if id_from == NCBI_ID:
            data[NCBI_ID] = gene_list.astype(str) if (id_from == NCBI_ID) and (gene_list.dtype == int) else gene_list
            data[NCBI_ID] = data[NCBI_ID].to_numpy()
            # data[NCBI_ID] = gene_list.replace('nan', np.nan) if (id_from == NCBI_ID) and (gene_list.dtype == int) else gene_list
        

        # convert the ids in the gene list
        converted_ids = convert_ids(gene_list, id_from, id_to, keep_mapping=True)

        # add the converted genes to the dataset
        data = pd.merge(data, converted_ids, on=id_from, how='outer')

        # move the id_to column next to the id_from column
        id_from_col_location = data.columns.get_loc(id_from)
        column_names = data.columns.to_list()
        column_names.remove(id_to)
        column_names.insert(id_from_col_location+1, id_to)
        data = data.reindex(columns=column_names)

        if keep_old_ids == False:
            data = data.drop(id_from, axis = 1)

        if column_name_with_ids:
            # rename column back to original name
            data = data.rename(columns={id_from:column_name_with_ids})
        
        if drop_na:
            data = data.dropna(subset=[id_to])

        return data.drop_duplicates()
    
    except InvalidGeneTypeError:
        pass

# ORTHOLOGY FUNCTIONS
# -------------------

def convert_to_human(gene_list, zfish_gene_type, 
                     keep_mapping = False, keep_missing_orthos = False):
    id_to = HUMAN_ID
    human_ids = get_ortho_ids(gene_list, zfish_gene_type, id_to, 
                              keep_mapping, keep_missing_orthos)
    return human_ids

def convert_to_zebrafish(gene_list, zfish_gene_type, 
                         keep_mapping = False, keep_missing_orthos = False):
    id_from = HUMAN_ID
    zebrafish_ids = get_ortho_ids(gene_list, id_from, zfish_gene_type, 
                                  keep_mapping, keep_missing_orthos)
    return zebrafish_ids

def get_ortho_ids(gene_list: list, id_from: str, id_to: str, 
                  keep_mapping = False, keep_missing_orthos = False):
    '''
    Parameters:
        data - a pandas DataFrame containing a column that has Gene IDs of some type
        id_from - the current Gene ID type must be: NCBI Gene ID, ZFIN ID, 
                Ensembl ID, Symbol, or Human NCBI Gene ID
        id_to - the Gene ID type to convert to, must be: NCBI Gene ID, ZFIN ID,
                Ensembl ID, Symbol, or Human NCBI Gene ID
        column_name_with_ids - the name of the column containing the Gene IDs (if
                               the column name does not match id_from)
        keep_old_ids - if you would like to keep the old Gene ID column
    '''

    # last updated: July 13, 2023
    # last checked: July 13, 2023
    # tutorial? yes

    try:
        # some error handling
        # -------------------    
        _check_valid_gene_id_type_for_orthology(id_from, id_to)
        gene_list = _make_sure_is_pandas_series(gene_list, id_from)
        gene_list = gene_list.drop_duplicates()
        # gene_list = gene_list.to_numpy()
        gene_list = gene_list.astype(str)

        # if id_from is a zebrafish gene, make sure it is in ZFIN ID format
        original_id_from = id_from
        zfish_gene_id_options = [NCBI_ID, ZFIN_ID, ENS_ID, SYMBOL]
        if (id_from in zfish_gene_id_options) and (id_from != ZFIN_ID):
            gene_list = convert_ids(gene_list, original_id_from, ZFIN_ID)
            id_from = ZFIN_ID
        
        # if id_to is a zebrafish gene and not ZFIN, 
        original_id_to = id_to
        if (id_to in zfish_gene_id_options) and (id_to != ZFIN_ID):
            id_to = ZFIN_ID

        # read in the master orthology mapping file
        master_mapping = pd.read_csv(MASTER_ORTHO_FILE_PATH, sep = '\t', dtype=str)

        # grab the columns we want
        df_desired_columns = master_mapping[[id_from, id_to]]
        # grab the rows we want
        filtered_df = df_desired_columns[df_desired_columns[id_from].isin(gene_list)]
        filtered_df = pd.merge(gene_list.to_frame(), filtered_df, on=id_from, how='left')

        if keep_mapping:
            mapped_genes = filtered_df
            # if id_to is a zebrafish gene type and not ZFIN ID, map to 
            # desired type
            if (original_id_to in zfish_gene_id_options) and (original_id_to != ZFIN_ID):
                mapped_genes = add_mapped_column(mapped_genes, ZFIN_ID, original_id_to, 
                                                 keep_old_ids=False)
            # if id_from is a zebrafish gene type and not ZFIN ID, map to 
            # desired type
            if (original_id_from in zfish_gene_id_options) and (original_id_from != ZFIN_ID):
                mapped_genes = add_mapped_column(mapped_genes, ZFIN_ID, original_id_from, 
                                                 keep_old_ids=False)
        else:
            mapped_genes = filtered_df[id_to].dropna()
            # if id_to is a zebrafish gene type and not ZFIN ID, map to 
            # desired type
            if (original_id_to in zfish_gene_id_options) and (original_id_to != ZFIN_ID):
                gene_list = filtered_df[id_to]
                mapped_genes = convert_ids(mapped_genes, ZFIN_ID, original_id_to)

        if keep_missing_orthos == False:
            mapped_genes = mapped_genes.dropna()
        mapped_genes = mapped_genes.reset_index(drop=True)
        return mapped_genes.drop_duplicates()
    except InvalidGeneTypeError:
        pass

def add_mapped_ortholog_column(data: pd.DataFrame, id_from: str, id_to: str, 
                               column_name_with_ids = None, keep_old_ids = True,
                               drop_na = False):
    '''
    Parameters:
        gene_list - a list of Gene IDs with supported formats list, pd.Series,
                    pd.DataFrame and np.array
        id_from - the current Gene ID type must be: NCBI Gene ID, ZFIN ID, 
                Ensembl ID, Symbol, or Human NCBI Gene ID
        id_to - the Gene ID type to convert to, must be: NCBI Gene ID, ZFIN ID,
                Ensembl ID, Symbol, or Human NCBI Gene ID
    '''

    # last updated: July 13, 2023
    # last checked: July 13, 2023.
    # tutorial? yes
    
    try:
        # some error handling
        # -------------------          
        if type(data) != pd.DataFrame:
            raise TypeError
        _check_valid_gene_id_type_for_orthology(id_from, id_to)
        _check_column_name_matches_id_choice(data, id_from, column_name_with_ids)
        # get the gene list from the given data
        gene_list = data[id_from]
        if id_from == NCBI_ID:
                # data[NCBI_ID] = gene_list.astype(str)
                data[NCBI_ID] = (gene_list.astype(str) if (id_from == NCBI_ID) and (gene_list.dtype == int) else gene_list)
                data[NCBI_ID] = data[NCBI_ID].to_numpy()
        if id_from == HUMAN_ID:
                # data[NCBI_ID] = gene_list.astype(str)
                data[HUMAN_ID] = (gene_list.astype(str) if (id_from == HUMAN_ID) and (gene_list.dtype == int) else gene_list)
                data[HUMAN_ID] = data[HUMAN_ID].to_numpy()

        gene_list_with_orthos = get_ortho_ids(gene_list, id_from, id_to, keep_mapping=True)

        # add the converted orthologs to the dataset
        data = pd.merge(data, gene_list_with_orthos, on=id_from, how = 'outer')

        # move the id_to column next to the id_from column
        id_from_col_location = data.columns.get_loc(id_from)
        column_names = data.columns.to_list()
        column_names.remove(id_to)
        column_names.insert(id_from_col_location+1, id_to)
        data = data.reindex(columns=column_names)

        if keep_old_ids == False:
            data = data.drop(id_from, axis = 1)
        
        if column_name_with_ids:
            # rename column back to original name
            data = data.rename(columns={id_from:column_name_with_ids})
        
        if drop_na:
            data = data.dropna(subset=[id_to])

        return data.drop_duplicates()
    
    except InvalidGeneTypeError:
        pass

# DATABASE BUILDING FUNCTIONS
# ----------------------

def build_gene_mapping():
    '''
    purpose:    the purpose of this function is to build the master gene mapping file 
                named 'master_gene_mapping_file_V<VERSION_NUM>.txt
                - the files used to build the mapping file have been downloaded from
                  ZFIN and should be in the database/raw_data sub directory
    dependencies:   zfin_to_ncbi_V<VERSION_NUM>.txt
                    zfin_to_ensembl_V<VERSION_NUM>.txt
    '''

    # last updated: July 7, 2023
    # last checked: July 7, 2023.

    try:
        # read in all the data files we have
        zfin_to_ncbi = RAW_DATA_DIR / Path('zfin_to_ncbi_V' + str(VERSION_NUM) + '.txt')
        zfin_to_ens = RAW_DATA_DIR / Path('zfin_to_ensembl_V' + str(VERSION_NUM) + '.txt')

        # check if database directory exists
        if not DATABASE_DIR.is_dir():
            raise DatabaseNotFoundError

        # check to see if required files exist 
        if (not zfin_to_ncbi.is_file()) or (not zfin_to_ens.is_file()):
            raise FileNotFoundError

        # clean some stuff up 
        zfin_geneid = pd.read_csv(zfin_to_ncbi, sep='\t')
        zfin_ens = pd.read_csv(zfin_to_ens, sep='\t')    

        zfin_geneid.drop(columns=['Unnamed: 4'], inplace=True)
        zfin_ens.drop(columns=['Unnamed: 4'], inplace=True)

        zfin_geneid.drop(columns=['SO ID'], inplace=True)
        zfin_ens.drop(columns=['Gene SO ID'], inplace=True)

        # merge the datasets
        merged_df = zfin_geneid.merge(zfin_ens, on=['ZFIN ID', 'Symbol'], how='outer' )

        out_file_name = DATABASE_DIR / Path('master_gene_mapping_file_V' + str(VERSION_NUM) + '.txt')
        merged_df.to_csv(out_file_name, sep = '\t', index=False)
        
    except FileNotFoundError:
        print('To build the gene mapping file, we require these two files:')
        print(' - zfin_to_ncbi_V<VERSION_NUM>.txt')
        print(' - zfin_to_ensembl_V<VERSION_NUM>.txt')
        print('One (or both) of these files was not found in the raw_data folder')
        print('of the database directory: ', RAW_DATA_DIR)
        print('Therefore, the build cannot be completed.')
    except DatabaseNotFoundError:
        print('The database should be accessible via the current working directory.')
        print('It cannot be found. Therefore, the build cannot be completed.')

def build_ortho_mapping():
    '''
    purpose:    the purpose of this function is to build the master orthology mapping file 
                named 'master_ortho_mapping_file_V<VERSION_NUM>.txt
                - the files used to build the orthology file have been downloaded from
                  ZFIN and should be in the database/raw_data sub directory
    dependencies:   zfin_to_ncbi_V<VERSION_NUM>.txt
                    zfin_to_ensembl_V<VERSION_NUM>.txt
    '''

    # last updated: July 7, 2023
    # last checked: July 7, 2023.

    try:

        # read in the gene data and human ortholog information
        # data downloaded from ZFIN!
        zfish_human_orthology_path = RAW_DATA_DIR / Path('zfish_human_orthology_V' + str(VERSION_NUM) + '.txt')
        # check to see if required files exist 
        if not zfish_human_orthology_path.is_file():
            raise FileNotFoundError
        zfish_human_orthos = pd.read_csv(zfish_human_orthology_path, sep='\t')

        # drop unneccessary columns
        zfin_orthos = zfish_human_orthos.drop(['OMIM ID', 'HGNC ID', 'Evidence', 'Pub ID', 'Unnamed: 10',
                                        'ZFIN Symbol', 'ZFIN Name', 'Human Symbol', 'Human Name'], axis=1)
        # # just include genes
        # zfin_orthos = zfin_orthos[zfin_orthos['ZFIN ID'].str.startswith('ZDB-GENE-') ]
        # remove duplicates
        zfin_orthos = zfin_orthos.drop_duplicates(keep='first')
        # reset index
        zfin_orthos=zfin_orthos.reset_index(drop=True)
        zfin_orthos.rename(columns={"Gene ID": "Human NCBI Gene ID"}, inplace=True)

        # # WRITE TO FILE
        out_file_name = DATABASE_DIR / Path('master_ortho_mapping_file_V'+str(VERSION_NUM)+'.txt')
        zfin_orthos.to_csv(out_file_name, sep = '\t', index=False)
    
    except FileNotFoundError:
        print('To build the orthology mapping file, we require the files:')
        print(' - zfish_human_orthology_V<VERSION_NUM>.txt')
        print('This file was not found in the raw_data folder')
        print('of the database directory: ', RAW_DATA_DIR)
        print('Therefore, the build cannot be completed.')
    except DatabaseNotFoundError:
        print('The database should be accessible via the current working directory.')
        print('It cannot be found. Therefore, the build cannot be completed.')

# PRIVATE FUNCTIONS
# -----------------
def convert_to_str(value):
    if value == 'nan':
        return np.nan  # Preserve NaN values
    else:
        return str(value)

def _make_sure_is_pandas_series(gene_list, id_from):
    if type(gene_list) != pd.Series:
        if type(gene_list) == pd.DataFrame:
            gene_list = gene_list.squeeze()
        else:
            gene_list = pd.Series(gene_list, name = id_from)
    return gene_list

def _check_valid_zebrafish_gene_id_type(gene_id_types):
    gene_id_options = [NCBI_ID, ZFIN_ID, ENS_ID, SYMBOL]

    if type(gene_id_types) == str:
        gene_id_types = [gene_id_types]

    invalid_gene_id_types = [item for item in gene_id_types if item not in gene_id_options]
    if invalid_gene_id_types:
        print('One or more of the zebrafish Gene ID types you gave is invalid. The valid zebrafish')
        print(f'Gene ID options are: {NCBI_ID}, {ZFIN_ID}, {ENS_ID}, and {SYMBOL}. The ID(s) you gave')
        print('that are invalid are:')
        print('----------')
        for item in invalid_gene_id_types:
            print(item)
        print('----------')
        print('Reminder: Gene ID types are case and spelling sensitive.')
        raise InvalidGeneTypeError
    
def _check_valid_gene_id_type_for_orthology(gene_id_from, gene_id_to):

    zfish_gene_id_options = [NCBI_ID, ZFIN_ID, ENS_ID, SYMBOL]
    human_gene_id_options = [HUMAN_ID]

    all_gene_id_options = zfish_gene_id_options + human_gene_id_options
    gene_ids_chosen = [gene_id_from, gene_id_to]

    invalid_gene_id_types = [item for item in gene_ids_chosen if item not in all_gene_id_options]
    if invalid_gene_id_types:
        print('One or more of the Gene ID types you gave is invalid. The valid zebrafish Gene ID')
        print(f'options are: {NCBI_ID}, {ZFIN_ID}, {ENS_ID}, and {SYMBOL}. The valid human Gene ID')
        print(f'type is: {HUMAN_ID}. The ID(s) you gave that are invalid are:')
        print('----------')
        for item in invalid_gene_id_types:
            print(item)
        print('----------')
        print('Reminder: Gene ID types are case and spelling sensitive.')
        raise InvalidGeneTypeError
    if not invalid_gene_id_types:
        if (gene_id_from in zfish_gene_id_options and gene_id_to in zfish_gene_id_options) or \
        (gene_id_from in human_gene_id_options and gene_id_to in human_gene_id_options):
            print('The Gene ID choices you gave are from the same organism. To perform othology mapping,')
            print('id_from must be from organism A and id_to must be from organism B.The valid zebrafish')
            print(f'Gene ID options are: {NCBI_ID}, {ZFIN_ID}, {ENS_ID}, and {SYMBOL}. The valid human')
            print(f'Gene ID type is: {HUMAN_ID}.')
            raise InvalidGeneTypeError


def _check_column_name_matches_id_choice(data: pd.DataFrame, id_from: str, column_name: str):
    data_cols = data.columns
    if column_name:
        if column_name not in data_cols:
            print('The column name you specified does not exist in your dataset.')
            raise InvalidGeneTypeError
    else:    
        if id_from not in data_cols:
            print('The Gene ID type you chose does not match a column name in the given')
            print('dataset. Please either change the column name in the dataset or provide')
            print('the matching column name in the column_name parameter.')
            raise InvalidGeneTypeError

if __name__ == '__main__':
    ids = ['ZDB-GENE-030131-2931',
    'ZDB-GENE-081031-55',
    'ZDB-GENE-081031-61',
    'ZDB-GENE-110420-1',
    'ZDB-GENE-200107-1',
    'ZDB-GENE-200107-2',
    'ZDB-GENE-071004-33',
    'ZDB-GENE-081105-47']
    out = convert_ids(ids, ZFIN_ID, NCBI_ID)
    print(out.dropna())