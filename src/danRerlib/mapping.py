'''
This script contains all the tools necessary to map zebrafish
gene ids from the varying identification types. The mapping 
data has been downloaded directly from NCBI and ZFIN. Check
out the database README for more information. 
'''

import pandas as pd
import numpy as np
import os.path 
from pathlib import Path

file_dir = Path(os.path.dirname(os.path.realpath(__file__)))
version_number = 1
database_dir = file_dir / Path('database/')
raw_data_dir = file_dir / Path('database/raw_data/')
master_mapping_file_path = database_dir / Path('master_gene_mapping_file_V' + str(version_number) + '.txt')
master_ortho_file_path = database_dir / Path('master_ortho_mapping_file_V' + str(version_number) + '.txt')


class DatabaseNotFoundError(Exception):
    "Raised when the the database directory is not found"
    pass

class InvalidGeneTypeError(Exception):
    "Raised when the the database directory is not found"
    pass


def build_gene_mapping():
    '''
    purpose:    the purpose of this function is to build the master gene mapping file 
                named 'master_gene_mapping_file_V<version_number>.txt
                - the files used to build the mapping file have been downloaded from
                  ZFIN and should be in the database/raw_data sub directory
    dependencies:   zfin_to_ncbi_V<version_number>.txt
                    zfin_to_ensembl_V<version_number>.txt
    '''

    try:
        # read in all the data files we have
        zfin_to_ncbi = raw_data_dir / Path('zfin_to_ncbi_V' + str(version_number) + '.txt')
        zfin_to_ens = raw_data_dir / Path('zfin_to_ensembl_V' + str(version_number) + '.txt')

        # check if database directory exists
        if not database_dir.is_dir():
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

        out_file_name = database_dir / Path('master_gene_mapping_file_V' + str(version_number) + '.txt')
        merged_df.to_csv(out_file_name, sep = '\t', index=False)
        
    except FileNotFoundError:
        print('To build the gene mapping file, we require these two files:')
        print(' - zfin_to_ncbi_V<version_number>.txt')
        print(' - zfin_to_ensembl_V<version_number>.txt')
        print('One (or both) of these files was not found in the raw_data folder')
        print('of the database directory: ', raw_data_dir)
        print('Therefore, the build cannot be completed.')
    except DatabaseNotFoundError:
        print('The database should be accessible via the current working directory.')
        print('It cannot be found. Therefore, the build cannot be completed.')

def build_ortho_mapping():
    '''
    purpose:    the purpose of this function is to build the master orthology mapping file 
                named 'master_ortho_mapping_file_V<version_number>.txt
                - the files used to build the orthology file have been downloaded from
                  ZFIN and should be in the database/raw_data sub directory
    dependencies:   zfin_to_ncbi_V<version_number>.txt
                    zfin_to_ensembl_V<version_number>.txt
    '''

    try:

        # read in the gene data and human ortholog information
        # data downloaded from ZFIN!
        zfish_human_orthology_path = raw_data_dir / Path('zfish_human_orthology_V' + str(version_number) + '.txt')
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
        out_file_name = database_dir / Path('master_ortho_mapping_file_V'+str(version_number)+'.txt')
        zfin_orthos.to_csv(out_file_name, sep = '\t', index=False)
    
    except FileNotFoundError:
        print('To build the orthology mapping file, we require the files:')
        print(' - zfish_human_orthology_V<version_number>.txt')
        print('This file was not found in the raw_data folder')
        print('of the database directory: ', raw_data_dir)
        print('Therefore, the build cannot be completed.')
    except DatabaseNotFoundError:
        print('The database should be accessible via the current working directory.')
        print('It cannot be found. Therefore, the build cannot be completed.')


def convert_ids_keep_mapping(gene_list, id_from, id_to):

    # read in the master mapping file
    master_mapping = pd.read_csv(file_dir / master_mapping_file_path, sep = '\t')
    
    # handle the case that NCBI Gene ID is not a string when it should be
    gene_list = gene_list.astype(str) if (id_from == 'NCBI Gene ID') and (gene_list.dtype == int) else gene_list

    # grab the columns we want
    df_desired_columns = master_mapping[[id_from, id_to]]
    filtered_df = df_desired_columns[df_desired_columns[id_from].isin(gene_list)]

    # get the list of mapped genes 
    return filtered_df.drop_duplicates()



def convert_ids(gene_list: any, id_from: str, id_to: str) -> pd.Series:
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

    try: 
        # some error handling
        # -------------------
        _check_valid_gene_id_type([id_from, id_to])
        gene_list = _make_sure_is_pandas_series(gene_list)
        # just in case the NCBI Gene IDs are strings
        gene_list = gene_list.astype(str) if (id_from == 'NCBI Gene ID') and (gene_list.dtype == int) else gene_list
        gene_list = gene_list.drop_duplicates()

        # read in the master mapping file
        master_mapping = pd.read_csv(file_dir / master_mapping_file_path, sep = '\t')
        master_mapping['NCBI Gene ID'] = master_mapping['NCBI Gene ID'].astype(str)
        df_desired_columns = master_mapping[[id_from, id_to]]

        # get the rows where our ids come from
        filtered_df = df_desired_columns[df_desired_columns[id_from].isin(gene_list)]

        # get the list of mapped genes 
        mapped_genes = filtered_df[id_to]

        return mapped_genes.drop_duplicates()
    
    except InvalidGeneTypeError:
        pass



def add_mapped_column(data: pd.DataFrame, id_from: str, id_to: str, keep_old_ids = True, position = 0):
    
    if type(data) != pd.DataFrame:
        raise TypeError
    
    # get the gene list from the given data
    gene_list = data[id_from]

    # convert the ids in the gene list
    converted_ids = convert_ids_keep_mapping(gene_list, id_from, id_to)

    # # add the converted genes to the dataset
    data = pd.merge(data, converted_ids, on=id_from, how='outer')

    if keep_old_ids == False:
        data = data.drop(id_from, axis = 1)

    return data.drop_duplicates()


def get_ortho_ids(gene_list, id_from, id_to):

    gene_list = _make_sure_is_pandas_series(gene_list)
    gene_list = gene_list.drop_duplicates()
    
    # read in the master mapping file
    master_mapping = pd.read_csv(master_ortho_file_path, sep = '\t', dtype=str)

    # grab the columns we want
    df_desired_columns = master_mapping[[id_from, id_to]]
    filtered_df = df_desired_columns[df_desired_columns[id_from].isin(gene_list)]

    # get the list of mapped genes 
    mapped_genes = filtered_df[id_to]
    
    return mapped_genes.drop_duplicates()


def get_ortho_ids_keep_mapping(gene_list, id_from, id_to): 
    # read in the master mapping file
    master_mapping = pd.read_csv(master_ortho_file_path, sep = '\t', dtype=str)
    if type(gene_list) == pd.DataFrame:
        gene_list = gene_list[id_from]
    # grab the columns we want
    df_desired_columns = master_mapping[[id_from, id_to]]
    filtered_df = df_desired_columns[df_desired_columns[id_from].isin(gene_list)]

    # get the list of mapped genes 
    return filtered_df.drop_duplicates()


def add_ortholog_column(data: pd.DataFrame, id_from: str, id_to: str, keep_old_ids = True):
    # TODO:  right now this deletes any row that does not have an ortholog?
    if type(data) != pd.DataFrame:
        raise TypeError
    
    # get the gene list from the given data
    gene_list = data[id_from]

    # convert the ids in the gene list
    converted_ids = get_ortho_ids_keep_mapping(gene_list, id_from, id_to)

    # # add the converted genes to the dataset
    data = pd.merge(data, converted_ids, on=id_from, how = 'right')

    if keep_old_ids == False:
        data = data.drop(id_from, axis = 1)

    return data.drop_duplicates()

#-------

def _make_sure_is_pandas_series(gene_list):
    
    if type(gene_list) != pd.Series:
        if type(gene_list) == pd.DataFrame:
            gene_list = gene_list.squeeze()
        else:
            gene_list = pd.Series(gene_list)
    return gene_list

def _check_valid_gene_id_type(gene_id_types):
    gene_id_options = ['NCBI Gene ID', 'ZFIN ID', 'Ensembl ID', 'Symbol', 'Human NCBI Gene ID']

    if type(gene_id_types) == str:
        gene_id_types = [gene_id_types]

    invalid_gene_id_types = [item for item in gene_id_types if item not in gene_id_options]
    if invalid_gene_id_types:
        print('One or more of the Gene ID types you gave is invalid. The valid Gene ID options')
        print('are: NCBI Gene ID, ZFIN ID, Ensembl ID, and Symbol. The ID(s) you gave that are')
        print('invalid are:')
        print('----------')
        for item in invalid_gene_id_types:
            print(item)
        print('----------')
        print('Reminder: Gene ID types are case and spelling sensitive.')
        raise InvalidGeneTypeError

def main():
    # gene_id_type = ['bka', 'NCBI Gene ID']
    # _check_valid_gene_id_type(gene_id_type)
    list_of_gene_ids = [ 
        100000252, 100000750, 100001198, 100001260, 100002225, 100002263, 
        100002756, 100003223, 100007521, 100149273, 100149794, 100170795,
        100321746, 100329897, 100330617,
    ]
    list_of_ids = np.array(list_of_gene_ids)
    out_ids = convert_ids(list_of_ids, 'NCBI Gene ID', 'ZFIN ID')
    print(out_ids)



if __name__ == '__main__':
    main()