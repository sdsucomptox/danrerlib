"""
Gene Mapping Module
===================

This module provides functions for gene ID mapping and orthology checks between different species, 
with a focus on Zebrafish and human gene IDs.

Functions:
    - ``convert_ids``: Convert a list of Gene IDs between different types.
    - ``add_mapped_column``: Add a new column with mapped Gene IDs to a pandas DataFrame.
    - ``convert_to_human``: Convert Zebrafish gene IDs to their corresponding human orthologs.
    - ``convert_to_zebrafish``: Convert human gene IDs to their corresponding Zebrafish orthologs.
    - ``add_mapped_ortholog_column``: Add a new column with mapped ortholog Gene IDs to a pandas DataFrame.

Constants:
    - ``NCBI_ID``: Identifier for NCBI Gene ID.
    - ``ZFIN_ID``: Identifier for ZFIN ID.
    - ``ENS_ID``: Identifier for Ensembl ID.
    - ``SYMBOL``: Identifier for gene Symbol.
    - ``HUMAN_ID``: Identifier for Human NCBI Gene ID.

Database Rebuild Functions:
    - ``build_gene_mapping``: Build a master gene mapping file based on Zebrafish gene ID conversions.
    - ``build_ortho_mapping``: Build a master orthology mapping file between Zebrafish and human genes.

Notes:
    - This module is designed for use in genomics research and gene ID mapping projects.
    - It provides functions to convert gene IDs between different types, build mapping files, and perform orthology checks.
    - If you would like to re-build mapping files, ensure that your dataset and mapping files are accurately named and organized as per the function requirements. It is recommended to use the current package database build. 

Example:
    To convert a list of Zebrafish gene IDs to human orthologs:
    ```
    gene_list = ['ENSDARG00000012345', 'ENSDARG00000067890']
    converted_genes = convert_to_human(gene_list, ENS_ID, HUMAN_ID)
    ```

For detailed information on each function and their usage, please refer to the documentation. For more examples of full functionality, please refer to tutorials.
"""

from danrerlib import utils
from danrerlib.settings import *
import os.path 
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Union, Optional, List


class DatabaseNotFoundError(Exception):
    "Raised when the the database directory is not found"
    pass

class InvalidGeneTypeError(Exception):
    "Raised when the Gene ID Type is invalid"
    pass

# GENE MAPPING FUNCTIONS
# ----------------------

def convert_ids(gene_list: Union[list, pd.Series, pd.DataFrame, np.array],
                id_from: str,
                id_to: str,
                keep_mapping: bool = False,
                out_format: Optional[str] = None
                ) -> Union[pd.Series, pd.DataFrame]:
    """
    Convert a list of zebrafish Gene IDs.

    Parameters:
        - ``gene_list (array-like)``: A list of Gene IDs with supported formats like list, pd.Series, pd.DataFrame, or np.array.
        - ``id_from (str)``: The current Gene ID type.
        - ``id_to (str)``: The Gene ID type to convert to.
        - ``keep_mapping (bool, optional)``: Whether to keep a mapping of the original IDs to the converted IDs. Default is False.
        - ``out_format (str, optional)``: The desired output format for the result. Default is None.

    Returns:
        - ``mapped_genes (pd.Series or pd.DataFrame)``: A pandas Series or DataFrame containing the converted Gene IDs.

    Other:
        - Gene ID Type Options: NCBI Gene ID, ZFIN ID, Symbol, Ensembl ID

    Notes:
        - If the order is important, it is recommended to have keep_mapping = True. Otherwise, 
          the mapping order is not guaranteed.
    """
    # last updated: July 7, 2023
    # last checked: July 7, 2023.]
    # tutorial? yes

    try: 
        # some error handling
        # -------------------
        id_from = utils.normalize_gene_id_type(id_from)
        id_to = utils.normalize_gene_id_type(id_to)
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

def add_mapped_column(data: Union[pd.DataFrame, list],
                      id_from: str,
                      id_to: str,
                      column_name_with_ids: str = None,
                      keep_old_ids: bool = True,
                      drop_na: bool = False
                      ) -> pd.DataFrame:
    """
    Add a new column to a pandas DataFrame with mapped zebrafish Gene IDs.

    Parameters:
        - ``data (pd.DataFrame, list)``: A pandas DataFrame containing a column that has Gene IDs of some type.
        - ``id_from (str)``: The current Gene ID type. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, or Symbol.
        - ``id_to (str)``: The Gene ID type to convert to. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, or Symbol.
        - ``column_name_with_ids (str, optional)``: The name of the column containing the Gene IDs if it doesn't match id_from.
        - ``keep_old_ids (bool, optional)``: Whether to keep the old Gene ID column. Default is True.
        - ``drop_na (bool, optional)``: Whether to drop rows with NA values in the resulting mapped column. Default is False.

    Returns:
        - ``data (pd.DataFrame)``: A pandas DataFrame containing the added mapped column.

    Notes:
        - This function adds a new column to the input DataFrame containing the mapped Gene IDs.
        - The new column will have the name 'Gene ID (Mapped)' unless specified otherwise.
    """

    # last updated: July 13, 2023
    # last checked: July 13, 2023
    # tutorial? yes

    try:
        # some error handling
        # -------------------    
        id_from = utils.normalize_gene_id_type(id_from)
        id_to = utils.normalize_gene_id_type(id_to)
        _check_valid_zebrafish_gene_id_type([id_from, id_to])
        if type(data) != pd.DataFrame:
            if type(data) == list:
                data = pd.DataFrame(data, columns=[id_from])
            # raise TypeError
        _check_column_name_matches_id_choice(data, id_from, column_name_with_ids)

        if column_name_with_ids:
            # rename data column to the id in options
            data = data.rename(columns={column_name_with_ids:id_from})

        gene_list = data[id_from]
        if id_from == NCBI_ID:
            gene_list = gene_list.astype(int)
            data[NCBI_ID] = gene_list.astype(str) if (id_from == NCBI_ID) and (gene_list.dtype == int) else gene_list
            data[NCBI_ID] = data[NCBI_ID].to_numpy()
            # data[NCBI_ID] = gene_list.replace('nan', np.nan) if (id_from == NCBI_ID) and (gene_list.dtype == int) else gene_list
            gene_list = data[id_from]
        
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

def convert_to_human(gene_list: List[str],
                     zfish_gene_type: str,
                     keep_mapping: bool = False,
                     keep_missing_orthos: bool = False
                     ) -> List[str]:
    """
    Convert a list of zebrafish gene IDs to their human orthologs.

    Parameters:
        - ``gene_list (list)``: A list of Zebrafish gene IDs to be converted to human orthologs.
        - ``zfish_gene_type (str)``: The current gene ID type for the Zebrafish genes. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, or Symbol.
        - ``keep_mapping (bool, optional)``: Whether to retain the mapping information. Default is False.
        - ``keep_missing_orthos (bool, optional)``: Whether to keep gene IDs with missing orthologs. Default is False.

    Returns:
        - ``human_ids (list)``: A list of human gene IDs corresponding to the Zebrafish gene IDs.

    Notes:
        - This function converts zebrafish gene IDs to their human orthologs using orthology mapping.
        - The ``zfish_gene_type`` parameter specifies the type of Zebrafish gene IDs. 
        - To retain the mapping information, set ``keep_mapping`` to True.
        - To keep Zebrafish gene IDs with missing orthologs, set `keep_missing_orthos` to True.
    """
    id_to = HUMAN_ID
    human_ids = get_ortho_ids(gene_list, zfish_gene_type, id_to, 
                              keep_mapping, keep_missing_orthos)
    return human_ids

def convert_to_zebrafish(gene_list: List[str],
                         zfish_gene_type: str,
                         keep_mapping: bool = False,
                         keep_missing_orthos: bool = False
                         ) -> List[str]:
    """
    Convert a list of human gene IDs to their zebrafish orthologs.

    Parameters:
        - ``gene_list (list)``: A list of human gene IDs to be converted to zebrafish orthologs.
        - ``zfish_gene_type (str)``: The target gene ID type for the zebrafish orthologs. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, or Symbol.
        - ``keep_mapping (bool, optional)``: Whether to retain the mapping information. Default is False.
        - ``keep_missing_orthos (bool, optional)``: Whether to keep gene IDs with missing orthologs. Default is False.

    Returns:
        - ``zebrafish_ids (list)``: A list of zebrafish gene IDs corresponding to the human gene IDs.

    Notes:
        - This function converts human gene IDs to their zebrafish orthologs using orthology mapping.
        - The ``zfish_gene_type`` parameter specifies the type of Zebrafish gene IDs.
        - To retain the mapping information, set ``keep_mapping`` to True.
        - To keep human gene IDs with missing Zebrafish orthologs, set ``keep_missing_orthos`` to True.
    """
    id_from = HUMAN_ID
    zebrafish_ids = get_ortho_ids(gene_list, id_from, zfish_gene_type, 
                                  keep_mapping, keep_missing_orthos)
    return zebrafish_ids

def get_ortho_ids(gene_list: List[str],
                  id_from: str,
                  id_to: str,
                  keep_mapping: bool = False,
                  keep_missing_orthos: bool = False
                  ) -> List[str]:
    """
    Retrieve orthologous gene IDs for a given list of genes.

    Parameters:
        - ``gene_list (list)``: A list of gene IDs.
        - ``id_from (str)``: The current gene ID type. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, Symbol, or Human NCBI Gene ID.
        - ``id_to (str)``: The target gene ID type to convert to. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, Symbol, or Human NCBI Gene ID.
        - ``keep_mapping (bool, optional)``: Whether to retain the mapping information. Default is False.
        - ``keep_missing_orthos (bool, optional)``: Whether to keep gene IDs with missing orthologs. Default is False.

    Returns:
        - ``mapped_genes (list)``: A list of orthologous gene IDs.

    Notes:
        - This function retrieves orthologous gene IDs for the provided gene list.
        - The mapping information can be retained by setting ``keep_mapping`` to True.
        - Gene IDs with missing orthologs can be retained by setting ``keep_missing_orthos`` to True.
    """

    # last updated: July 13, 2023
    # last checked: July 13, 2023
    # tutorial? yes

    try:
        # some error handling
        # -------------------
        id_from = utils.normalize_gene_id_type(id_from)
        id_to = utils.normalize_gene_id_type(id_to)    
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

def add_mapped_ortholog_column(data: pd.DataFrame,
                               id_from: str,
                               id_to: str,
                               column_name_with_ids: str = None,
                               keep_old_ids: bool = True,
                               drop_na: bool = False
                               ) -> pd.DataFrame:
    """
    Add a new column to a pandas DataFrame with mapped ortholog Gene IDs.

    Parameters:
        - ``data (pd.DataFrame)``: A pandas DataFrame containing a column with Gene IDs of the specified 'id_from' type.
        - ``id_from (str)``: The current Gene ID type to convert from. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, Symbol, or Human NCBI Gene ID.
        - ``id_to (str)``: The Gene ID type to convert to. Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, Symbol, or Human NCBI Gene ID.
        - ``column_name_with_ids (optional, str) ``: The name of the column containing the Gene IDs if it doesn't match 'id_from'.
        - ``keep_old_ids (optional, bool)``: Whether to keep the original Gene ID column. Default is True.
        - ``drop_na (optional, bool)``: Whether to drop rows with NA values in the resulting mapped column. Default is False.

    Returns:
        - data (pd.DataFrame): A pandas DataFrame containing the added mapped ortholog Gene IDs.

    Notes:
        - This function adds a new column to the input DataFrame containing the mapped ortholog Gene IDs.
        - The new column will have the name 'Ortholog Gene ID (Mapped)' unless specified otherwise.
    """

    # last updated: July 13, 2023
    # last checked: July 13, 2023.
    # tutorial? yes
    
    try:
        # some error handling
        # -------------------          
        if type(data) != pd.DataFrame:
            if type(data) == pd.Series:
                data = pd.DataFrame(data)
            else:
                raise TypeError
        id_from = utils.normalize_gene_id_type(id_from)
        id_to = utils.normalize_gene_id_type(id_to)
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

        return data.drop_duplicates().reset_index(drop=True)
    
    except InvalidGeneTypeError:
        pass

# DATABASE BUILDING FUNCTIONS
# ----------------------

def build_gene_mapping():
    """
    Build the master gene mapping file.

    Raises:
        - ``FileNotFoundError``: If one or both of the required data files are not found in the ``raw_data`` folder of the database directory.
        - ``DatabaseNotFoundError``: If the database directory is not found in the current working directory.

    Notes:
        - This function reads data files from the ``database/raw_data`` subdirectory to create the mapping file. The required data files includes ``zfin_to_ncbi_V<VERSION_NUM>.txt`` and ``zfin_to_ensembl_V<VERSION_NUM>.txt``.
        - The resulting merged DataFrame is saved as ``master_gene_mapping_file_V<VERSION_NUM>.txt`` in the database directory.
        - Make sure to have the required data files in the ``raw_data`` subdirectory before running this function.
    """

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
    """
    Build the master orthology mapping file.

    Raises:
        ``FileNotFoundError``: If the required data file is not found in the ``raw_data`` folder of the database directory.
        ``DatabaseNotFoundError``: If the database directory is not found in the current working directory.

    Notes:
        - This function reads data files from the ``database/raw_data`` subdirectory to create the orthology mapping file. The required data file includes ``zfish_human_orthology_V<VERSION_NUM>.txt``.
        - This function processes the 'zfish_human_orthology' dataset, removing unnecessary columns.
        - Duplicate entries are removed, and the DataFrame is saved as ``master_ortho_mapping_file_V<VERSION_NUM>.txt`` in the database directory.
        - Ensure that you have the required data file in the ``raw_data`` folder of the database directory before running this function.
    """

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

def _make_sure_is_pandas_series(gene_list: Union[pd.Series, pd.DataFrame, list], id_from: str
                                ) -> pd.Series:
    """
    Ensure that the input is a pandas Series.

    Parameters:
        - ``gene_list (pd.Series, pd.DataFrame, list)``: Input data, which can be a Series, DataFrame, or list.
        - ``id_from (str)``: Identifier for the data.

    Returns:
        - gene_list (pd.Series): A pandas Series based on the input data.

    Notes:
        - This function checks the type of the input data and converts it into a pandas Series if necessary.
        - If the input is a DataFrame, it will be squeezed into a Series.
        - If the input is a list, it will be converted into a Series with the specified name.
    """
    if type(gene_list) != pd.Series:
        if type(gene_list) == pd.DataFrame:
            gene_list = gene_list.squeeze()
        else:
            gene_list = pd.Series(gene_list, name = id_from)
    return gene_list

def _check_valid_zebrafish_gene_id_type(gene_id_types: Union[str, List[str]]) -> None:
    """
    Check the validity of Zebrafish gene ID types.

    This function checks if the provided Zebrafish gene ID types are valid options.

    Parameters:
        gene_id_types (str or list): A string or a list of Zebrafish gene ID types to be validated.

    Raises:
        InvalidGeneTypeError: If one or more of the provided gene ID types are invalid.

    Notes:
        - Valid Zebrafish gene ID types include: NCBI Gene ID, ZFIN ID, Ensembl ID, or Symbol.
        - The input can be a single gene ID type as a string or multiple types in a list.
        - Gene ID types are case and spelling sensitive.
    """
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
    
def _check_valid_gene_id_type_for_orthology(gene_id_from: str, 
                                            gene_id_to: str
                                            ) -> None:
    """
    Check the validity of gene ID types for orthology mapping.

    This function checks if the provided gene ID types are valid options for orthology mapping.
    It ensures that the gene ID types are not from the same organism.

    Parameters:
        gene_id_from (str): The source gene ID type.
        gene_id_to (str): The target gene ID type.

    Raises:
        InvalidGeneTypeError: If one or both of the provided gene ID types are invalid or from the same organism.

    Notes:
        - Valid Zebrafish gene ID types include: NCBI_ID, ZFIN_ID, ENS_ID, and SYMBOL.
        - Valid human gene ID type is: HUMAN_ID.
        - Gene ID types are case and spelling sensitive.
        - To perform orthology mapping, id_from must be from organism A, and id_to must be from organism B.
    """

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


def _check_column_name_matches_id_choice(data: pd.DataFrame, 
                                         id_from: str, 
                                         column_name: str
                                         ) -> None:
    """
    Check if the specified column name matches the chosen gene ID type.

    This function ensures that the specified column name exists in the DataFrame if provided,
    or that the chosen gene ID type matches a column name in the DataFrame.

    Parameters:
        data (pd.DataFrame): The DataFrame containing gene ID information.
        id_from (str): The chosen gene ID type.
        column_name (str): The name of the column containing gene IDs (if specified).

    Raises:
        InvalidGeneTypeError: If the specified column name does not exist in the dataset or
                             if the chosen gene ID type does not match any column name.

    Notes:
        - If the 'column_name' parameter is provided, this function checks if it exists in the DataFrame.
        - If 'column_name' is not provided, it checks if 'id_from' matches any column name in the DataFrame.
        - Ensure that the dataset and column names are accurate to perform gene ID conversion.
    """
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