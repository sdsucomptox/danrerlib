"""
KEGG Module
===========

This module provides functions for retrieving gene information associated with KEGG pathways and diseases
for various organisms, including human, zebrafish, and mapped zebrafish.

Functions:
    - ``get_genes_in_pathway``: Retrieve genes associated with a specific KEGG pathway.
    - ``get_genes_in_disease``: Retrieve genes associated with a disease for a specified organism.

Constants:
    - ``NCBI_ID``: Identifier for NCBI Gene ID.
    - ``ZFIN_ID``: Identifier for ZFIN ID.
    - ``ENS_ID``: Identifier for Ensembl ID.
    - ``SYMBOL``: Identifier for gene Symbol.
    - ``HUMAN_ID``: Identifier for Human NCBI Gene ID.

Database Rebuild Functions:
    - ``build_kegg_database``: Build or update the KEGG pathway and disease database.

Notes:
    - This module is designed for accessing gene information from KEGG pathways and diseases.
    - It provides functions to retrieve gene data associated with specific pathways and diseases.
    - Data is obtained either through direct API calls or by reading pre-processed files.
    - Organism options include 'dre' (Danio rerio), 'hsa' (Human), or 'dreM' (Mapped Danio rerio from Human).

Example:
    To retrieve genes associated with a specific KEGG pathway:
    ```
    pathway_genes = get_genes_in_pathway('hsa04010', 'hsa')
    ```

For detailed information on each function and their usage, please refer to the documentation. For more examples of full functionality, please refer to tutorials.
"""

from danrerlib.settings import *
from danrerlib import mapping, utils
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Union, Optional, List
import time
from urllib.error import HTTPError

zebrafish_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_dre_V' + str(VERSION_NUM) + '.txt')
human_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_hsa_V' + str(VERSION_NUM) + '.txt')
mapped_zebrafish_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_dreM_V' + str(VERSION_NUM) + '.txt')

zebrafish_pathways_dir = KEGG_DATA_DIR / Path('dre/')
human_pathways_dir = KEGG_DATA_DIR / Path('hsa/')
dre_mapped_dir = KEGG_DATA_DIR / Path('dreM/')

human_disease_path = KEGG_DATA_DIR / Path('disease_ids_V'+ str(VERSION_NUM) + '.txt')
human_disease_genes_path = KEGG_DATA_DIR / Path('disease_ids_and_genes_V'+ str(VERSION_NUM) + '.txt')
dre_disease_genes_path = KEGG_DATA_DIR / Path('disease_ids_dreM_and_genes_V'+ str(VERSION_NUM) + '.txt')
empty_disease_ids_path = KEGG_DATA_DIR / Path('empty_disease_ids_V'+ str(VERSION_NUM) + '.txt')
valid_disease_ids_path = KEGG_DATA_DIR / Path('disease_ids_valid_V'+ str(VERSION_NUM) + '.txt')

def get_genes_in_pathway(pathway_id: str, org=None, do_check = True):
    """
    Retrieve genes associated with a specific pathway from KEGG.

    Parameters:
        - ``pathway_id (str)``: The KEGG pathway ID for the pathway of interest.
        - ``org (str, optional)``: The organism for which to retrieve pathway information. Options include

                - 'dre': Danio rerio (Zebrafish) pathways from KEGG.
                - 'hsa': Human pathways from KEGG.
                - 'dreM': Mapped Danio rerio pathways from human.
                - (Default is None, which will use the provided 'org' or 'hsa' if 'org' is None.)

    Returns:
        - ``df (pd.DataFrame)``: A pandas DataFrame containing genes associated with the specified pathway.

    Notes:
        - This function retrieves gene information associated with a specific KEGG pathway.
        - The ``pathway_id`` parameter should be a valid KEGG pathway identifier.
        - The ``org`` parameter specifies the organism for which to retrieve pathway information.
        - Organism options include 'dre' (Danio rerio), 'hsa' (Human), or 'dreM' (Mapped Danio rerio from Human).
    """
    try:
        if do_check:
            if org:
                org = utils.normalize_organism_name(org)
            org, pathway_id = _check_for_organism(pathway_id, org)
            result = _check_if_pathway_id_exists(pathway_id, org, error_message=True)
            if result:
                file_name = KEGG_DATA_DIR / Path(org) / Path(pathway_id + '.txt')
                df = pd.read_csv(file_name, sep='\t')
                return df
            else:
                raise ValueError
        else:                 
            file_name = KEGG_DATA_DIR / Path(org) / Path(pathway_id + '.txt')
            df = pd.read_csv(file_name, sep='\t')
            return df
    except ValueError:
        pass

def get_genes_in_disease(disease_id: str, 
                         org: str, 
                         do_check = True,
                         API: bool = False
                         ) -> pd.DataFrame:
    """
    Retrieve genes associated with a disease for a specified organism.

    Parameters:
        - ``disease_id (str)``: The disease ID for which genes are to be retrieved.
        - ``org (str)``: The organism for which genes should be retrieved. Options include

                   - 'hsa': Returns human genes associated with the disease.
                   - 'dre' or 'dreM': Returns human genes mapped to Zebrafish genes (same result, as Zebrafish disease not characterized on KEGG).

        - ``API (bool, optional)``: Whether to use the KEGG REST API for retrieval. Default is False.

    Returns:
       - ``genes (pd.DataFrame)``: A pandas DataFrame containing genes associated with the disease for the specified organism.

    Notes:
        - If 'API' is set to True, the function fetches gene information using the KEGG REST API.
        - If 'API' is False, the function reads gene information from pre-processed files.
    """

    if API:
        url = 'https://rest.kegg.jp/link/hsa/' + disease_id
        column_names = ['trash', HUMAN_ID]
        gene_df = pd.read_csv(url, names=column_names, sep='\t')
        genes = gene_df[HUMAN_ID].str[4:]
        if org == 'hsa':
            if type(genes) != pd.DataFrame:
                genes = genes.to_frame()
                genes[HUMAN_ID] = genes[HUMAN_ID].values.astype(np.int64)
            return genes
        elif org == 'dre' or org == 'dreM':
            genes = mapping.convert_to_zebrafish(genes, NCBI_ID)
            if type(genes) != pd.DataFrame:
                genes = genes.to_frame()
                genes[NCBI_ID] = genes[NCBI_ID].values.astype(np.int64)
            return genes
    else:
        if do_check:
            if org:
                org = utils.normalize_organism_name(org)
            result = _check_if_disease_id_exists(disease_id, org, error_message=True)
            if not result:
                raise ValueError
            
        file_dict_by_org = {
            'hsa': human_disease_genes_path,
            'dre': dre_disease_genes_path,
            'dreM': dre_disease_genes_path
        }

        id_type_by_org = {
            'hsa': HUMAN_ID,
            'dre': NCBI_ID,
            'dreM': NCBI_ID
        }

        genes_and_disease_df = pd.read_csv(file_dict_by_org[org], sep = '\t')
        filtered_df = genes_and_disease_df[genes_and_disease_df['Disease ID'] == disease_id]

        genes = filtered_df[id_type_by_org[org]]
        if type(genes) != pd.DataFrame:
            genes = genes.to_frame()
            genes[id_type_by_org[org]] = genes[id_type_by_org[org]].values.astype(np.int64)

        return genes.reset_index(drop=True)

def _get_total_number_kegg(org, db):
    """
    Retrieve the total number of KEGG pathways or KEGG diseases available for a specified organism.

    Parameters:
    - `org (str)`: The organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish, 'hsa' for human).

    Returns:
    - `int`: The total number of KEGG pathways/diseases for the specified organism.

    Raises:
    - KeyError: If the provided organism code is not found in the known organism codes.

    """

    if db == 'pathway':
        org_path_dict = {
            'dre' : zebrafish_pathways_path,
            'dreM' :mapped_zebrafish_pathways_path,
            'hsa' : human_pathways_path
        }
    elif db == 'disease':
        org_path_dict = {
            'dre' : human_disease_genes_path,
            'dreM' :human_disease_genes_path,
            'hsa' : human_disease_genes_path
        } 


    path = org_path_dict[org]

    if path is None:
        raise KeyError(f"Invalid organism code: {org}")
    
    df = pd.read_csv(path, sep = '\t')
    total_num_pathways = len(df)
    return total_num_pathways
    
def _check_for_organism(pathway_id: str, 
                        org: Optional[str]
                        ) -> [str, str]:

    """
    Check for and validate the organism information in a KEGG pathway ID.

    Parameters:
        - pathway_id (str): The KEGG pathway ID, which may or may not include organism information.
        - org (str, optional): The specified organism code to validate against the pathway ID.

    Returns:
        tuple: A tuple containing the validated organism code and the modified pathway ID.

    Raises:
        ValueError: If the organism code in the pathway ID does not match the specified 'org' parameter.

    Notes:
        - This internal function is used to validate and extract organism information from a KEGG pathway ID.
        - If the pathway ID includes organism information, it validates it against the 'org' parameter.
        - If 'org' is None, it extracts and returns the organism code from the pathway ID.
        - Raises a ValueError if there is a mismatch between the organism code in the pathway ID and the 'org' parameter.
    """
    pathway_id = str(pathway_id)
    if pathway_id[0].isdigit():
        if org == None:
            print('ERROR - no organism specified or included in Pathway ID')
            raise ValueError
        pathway_id = org + pathway_id
    if org == None:
        org = pathway_id[:-5]
        return org, pathway_id
    else:
        org_from_id = pathway_id[:-5]
        if org != org_from_id:
            print('The organism you specified does not match the organism prefix')
            print('for the pathway ID.')
            raise ValueError
        return org, pathway_id

def _check_if_pathway_id_exists(pathway_id: str, 
                                org: str,
                                error_message = False
                                ) -> bool:
    """
    Check if a KEGG pathway ID exists for the specified organism.

    Parameters:
        pathway_id (str): The KEGG pathway ID to be validated.
        org (str): The specified organism code (e.g., 'dre', 'dreM', or 'hsa').

    Returns:
        bool: True if the pathway ID exists for the specified organism, False otherwise.

    Notes:
        - This internal function checks whether a given KEGG pathway ID exists for the specified organism.
        - It reads relevant pathway data files based on the organism code and checks if the pathway ID is present.
        - Returns True if the pathway ID exists, False otherwise.
        - Provides an error message if the pathway ID or organism code is invalid.
        - Organism code options include: 'dre' (Danio rerio), 'dreM' (Mapped Danio rerio from Human), or 'hsa' (Human).
        - Pathway IDs must be entered as strings if the organism prefix is omitted.
    """
    result = False
    if org == 'dreM':
        if not pathway_id.startswith('dreM'):
            pathway_id_nums = pathway_id[3:]
            pathway_id = 'dreM'+pathway_id_nums
        df = pd.read_csv(mapped_zebrafish_pathways_path, sep='\t')
        result = True if pathway_id in df['Pathway ID'].values else False
    elif org == 'hsa':
        df = pd.read_csv(human_pathways_path, sep='\t')
        result = True if pathway_id in df['Pathway ID'].values else False
    elif org == 'dre':
        df = pd.read_csv(zebrafish_pathways_path, sep='\t')
        result = True if pathway_id in df['Pathway ID'].values else False
    if error_message and not result:
        print('The Pathway ID you gave does not exist for the organism you specified.')
        print('A Pathway ID is identified by the combination of a 3-4 letter prefix')
        print('code and a 5 digit number.\n')
        print('Prefix options include: dre, dreM, or hsa\n')
        print('reminder - python does not support integers that start with a 0;')
        print('           therefore, Pathway IDs must be entered as strings if')
        print('           the organism prefix is omitted.')
    return result

def _check_if_disease_id_exists(disease_id: str, 
                                org: str,
                                error_message = False) -> bool:
    """
    Check if a KEGG disease ID exists for the specified organism.

    Parameters:
        disease_id (str): The KEGG disease ID to be validated.
        org (str): The specified organism code (e.g., 'dre', 'dreM', or 'hsa').

    Returns:
        bool: True if the disease ID exists for the specified organism, False otherwise.

    Notes:
        - This internal function checks whether a given KEGG disease ID exists for the specified organism.
        - It reads relevant disease data files based on the organism code and checks if the disease ID is present.
        - Returns True if the disease ID exists, False otherwise.
        - Provides an error message if the disease ID or organism code is invalid.
        - Organism code options include: 'dre' (Danio rerio), 'dreM' (Mapped Danio rerio from Human), or 'hsa' (Human).
        - Disease IDs must be entered as strings if the organism prefix is omitted.
    """
    result = False
    disease_id = disease_id.strip()
    orgs = ['dreM', 'dre', 'hsa']
    if org in ['dreM', 'dre']:
        df = pd.read_csv(dre_disease_genes_path, sep='\t')
        result = True if disease_id in df['Disease ID'].values else False
    elif org == 'hsa':
        df = pd.read_csv(human_disease_genes_path, sep='\t')
        result = True if disease_id in df['Disease ID'].values else False       
    if error_message and not result:
        print('The Disease ID you gave does not exist for the organism you specified.')
        print('A Disease ID is identified by the prefix H and a 5 digit number.\n')
        print('reminder - python does not support integers that start with a 0;')
        print('           therefore, Dathway IDs must be entered as strings if')
        print('           the prefix is omitted.')
    return result

# DATABASE BUILDING FUNCTIONS
# ---------------------------

def build_kegg_database():
    """
    Build or update the KEGG pathway and disease database.

    Notes:
        - This function is intended for database creation or update and should be run only during version updates.
        - The process may take some time, so exercise caution when running it.
        - The function performs the following steps:
        
            1. Downloads zebrafish pathway IDs and stores them in a specified file.
            2. Downloads human pathway IDs and stores them in a specified file.
            3. Creates mapped zebrafish pathway IDs from human pathways and stores them in a specified file.
            4. Downloads genes associated with true KEGG pathways for zebrafish.
            5. Downloads genes associated with true KEGG pathways for humans.
            6. Builds mapped zebrafish KEGG pathways from human pathways and stores them in a specified directory.
            7. Downloads KEGG disease IDs and stores them in a specified file.

        - Running this function should be done carefully, as it involves downloading and processing data from KEGG.
    """

    # pathway ids
    _download_zebrafish_pathway_ids(zebrafish_pathways_path)
    _download_human_pathway_ids(human_pathways_path)
    _create_mapped_zebrafish_pathway_ids(mapped_zebrafish_pathways_path, human_pathways_path)

    # get genes in true kegg pathways
    _download_zebrafish_pathway_genes(zebrafish_pathways_path, zebrafish_pathways_dir)
    _download_human_pathway_genes(human_pathways_path, human_pathways_dir)

    # get mapped human kegg pathways to zebrafish
    _build_dre_mapped_pathway(human_pathways_path, human_pathways_dir, dre_mapped_dir)

    # get kegg disease ids
    _download_disease_ids(human_disease_path)
    _download_human_disease_genes(human_disease_path)

    # get mapped human kegg pathways to zebrafish
    _build_dre_mapped_disease()

def _download_zebrafish_pathway_ids(zebrafish_pathways_path):
    """
    Download and store the list of zebrafish pathways and pathway names.

    Parameters:
        zebrafish_pathways_path (str): The file path to save the zebrafish pathway data.

    Notes:
        - This function is intended for database building and should be run only during version updates.
    """
    url = 'https://rest.kegg.jp/list/pathway/dre'
    names=['Pathway ID', 'Pathway Description']
    zebrafish_pathways = pd.read_csv(url, names=names, sep='\t')
    zebrafish_pathways['Pathway Description'] = zebrafish_pathways['Pathway Description'].str.replace(' - Danio rerio (zebrafish)', '', regex=False)
    zebrafish_pathways.to_csv(zebrafish_pathways_path, sep='\t', index=False)

def _download_human_pathway_ids(human_pathways_path):
    """
    Download and store the list of human pathways and pathway names.

    Parameters:
        human_pathways_path (str): The file path to save the human pathway data.

    Notes:
        - This function is intended for database building and should be run only during version updates.
    """
    url = 'https://rest.kegg.jp/list/pathway/hsa'
    names=['Pathway ID', 'Pathway Description']
    human_pathways = pd.read_csv(url, names=names, sep='\t')
    human_pathways['Pathway Description'] = human_pathways['Pathway Description'].str.replace(' - Homo sapiens (human)', '', regex=False)
    human_pathways.to_csv(human_pathways_path, sep='\t', index=False)

def _create_mapped_zebrafish_pathway_ids(mapped_zebrafish_pathways_path, human_pathways_path):
    """
    Create mapped zebrafish pathway IDs from human pathways and store them.

    Parameters:
        mapped_zebrafish_pathways_path (str): The file path to save the mapped zebrafish pathway data.
        human_pathways_path (str): The file path containing human pathway data.

    Notes:
        - This function is intended for database building and should be run only during version updates.
    """  
    df = pd.read_csv(human_pathways_path, sep='\t')
    df['Pathway ID'] = df['Pathway ID'].str.replace('hsa', 'dreM')
    df.to_csv(mapped_zebrafish_pathways_path, sep='\t', index=False)

def _download_disease_ids(human_disease_path):
    """
    Download and store the list of KEGG disease IDs and descriptions for humans.

    Parameters:
        human_disease_path (str): The file path to save the human disease data.

    Notes:
        - This function is intended for database building and should be run only during version updates.
    """
    url = 'https://rest.kegg.jp/list/disease'
    names=['Disease ID', 'Disease Description']
    human_diseases = pd.read_csv(url, names=names, sep='\t')
    human_disease_path = KEGG_DATA_DIR / Path('disease_ids_V' + str(VERSION_NUM) + '.txt')
    human_diseases.to_csv(human_disease_path, sep='\t', index=False)

def _download_human_pathway_genes(human_pathways_path, human_pathways_dir):
    """
    Download and store genes associated with true KEGG pathways for humans.

    Parameters:
        human_pathways_path (str): The file path containing human pathway data.
        human_pathways_dir (str): The directory to save gene data associated with human pathways.

    Notes:
        - This function is intended for database building and should be run only during version updates.
    """
    human_pathways = pd.read_csv(human_pathways_path, sep='\t')
    for pathway_id in human_pathways['Pathway ID']:
        url = 'https://rest.kegg.jp/link/hsa/'+pathway_id
        column_names = ['trash', HUMAN_ID]
        genes_df = pd.read_csv(url, names = column_names, sep='\t')
        genes = genes_df[HUMAN_ID].str[4:]
        file_path = human_pathways_dir / Path(str(pathway_id) + '.txt')
        genes.to_csv(file_path, sep='\t', index=False)

def _download_zebrafish_pathway_genes(zebrafish_pathways_path, zebrafish_pathways_dir):
    """
    Download and store genes associated with true KEGG pathways for zebrafish.

    Parameters:
        zebrafish_pathways_path (str): The file path containing zebrafish pathway data.
        zebrafish_pathways_dir (str): The directory to save gene data associated with zebrafish pathways.

    Notes:
        - This function is intended for database building and should be run only during version updates.
    """
    zebrafish_pathways = pd.read_csv(zebrafish_pathways_path, sep='\t')
    for pathway_id in zebrafish_pathways['Pathway ID']:
        url = 'https://rest.kegg.jp/link/dre/'+pathway_id
        column_names = ['trash', NCBI_ID]
        genes_df = pd.read_csv(url, names = column_names, sep='\t')
        genes = genes_df[NCBI_ID].str[4:]
        file_path = zebrafish_pathways_dir / Path(str(pathway_id) + '.txt')
        genes.to_csv(file_path, sep='\t', index=False)

def _download_human_disease_genes(disease_ids_path):
    """
    Download and store genes associated with human diseases.

    Parameters:
        disease_ids_path (str): The file path to the list of disease IDs.

    Returns:
        None

    Notes:
        - This function reads a list of human disease IDs from the specified file.
        - It iterates through each disease ID and retrieves genes associated with the disease.
        - Genes are stored in a DataFrame with columns 'Human NCBI Gene ID' and 'Disease ID'.
        - The resulting DataFrame is saved to a file named 'human_disease_genes.txt'.
        - This function is intended for use in building or updating a database of human disease genes.
    """
    # there is an API issue for the number of requests made..... right now you have to run this function
    # multiple times to actually get all the data. it needs to be modified 
    try:
        disease_data = pd.read_csv(disease_ids_path, sep='\t')
        m, n = disease_data.shape

        # check which ids we already have 
        result_df = pd.read_csv(human_disease_genes_path, sep = '\t')
        completed_disease_ids = result_df['Disease ID'].unique()

        print(f"num of diseases: {m}, \t num completed: {len(completed_disease_ids)}")

        # Create an empty DataFrame to store the results
        # result_df = pd.DataFrame(columns=[HUMAN_ID, 'Disease ID'])
        empty_disease_ids = pd.read_csv(empty_disease_ids_path, sep = '\t')
        empty_disease_ids = empty_disease_ids['Disease ID'].tolist()

        # Iterate through each Disease ID and retrieve genes
        for index, row in disease_data.iterrows():
            disease_id = row['Disease ID']

            if (disease_id not in completed_disease_ids 
                and disease_id not in empty_disease_ids):
                print(f"{index} / {m}: {disease_id}")
                genes = _get_genes_for_disease(disease_id)
                if genes.empty:
                    empty_disease_ids.append(disease_id)
                else:
                    genes_df = pd.DataFrame({HUMAN_ID: genes, 'Disease ID': disease_id})
                    result_df = pd.concat([result_df, genes_df], ignore_index=True)

        empty_disease_ids_df = pd.DataFrame(empty_disease_ids, columns=['Disease ID'])
        empty_disease_ids_df.to_csv(empty_disease_ids_path, sep = '\t')
        result_df.to_csv(human_disease_genes_path, index=False, sep = '\t')

    except HTTPError as e:
        if e.code == 403:
            # Handle error here
            # Save the result to a file
            result_df.to_csv(human_disease_genes_path, index=False, sep = '\t')
            empty_disease_ids_df = pd.DataFrame(empty_disease_ids, columns=['Disease ID'])
            empty_disease_ids_df.to_csv(empty_disease_ids_path, index=False, sep = '\t')

def _get_genes_for_disease(disease_id):
    url = 'https://rest.kegg.jp/link/hsa/' + disease_id
    column_names = ['trash', HUMAN_ID]
    gene_df = pd.read_csv(url, names=column_names, sep='\t')
    genes = gene_df[HUMAN_ID].str[4:]
    return genes

def _build_dre_mapped_pathway(human_pathways_path, human_pathways_dir, dre_mapped_dir):
    """
    Build mapped zebrafish KEGG pathways from human pathways and store them.

    Parameters:
        human_pathways_path (str): The file path containing human pathway data.
        human_pathways_dir (str): The directory containing gene data associated with human pathways.
        dre_mapped_dir (str): The directory to save mapped zebrafish pathway data.

    Notes:
        - This function is intended for database building and should be run only during version updates.
    """
    human_pathways_df = pd.read_csv(human_pathways_path, sep='\t')
    for pathway_id in human_pathways_df['Pathway ID']:
        in_file_name = human_pathways_dir / Path(pathway_id+'.txt')
        humman_genes = pd.read_csv(in_file_name, sep='\t')
        zebrafish_genes = mapping.convert_to_zebrafish(humman_genes, NCBI_ID, keep_missing_orthos=False)
        stripped_pathway_id = pathway_id[3:]
        out_file_name = dre_mapped_dir / Path(stripped_pathway_id+'.txt')
        zebrafish_genes.to_csv(out_file_name, sep='\t', index=False)

def _build_dre_mapped_disease(human_disease_genes_path: str, 
                              dre_disease_genes_path: str
                              ) -> None:
    """
    Build mapped Zebrafish KEGG disease genes from human genes and store them.

    Parameters:
        human_disease_genes_path (str): Path to the file containing human disease-associated genes.
        dre_disease_genes_path (str): Path to the file where the mapped Zebrafish disease genes will be stored.

    Notes:
        - This function takes a list of human genes associated with diseases, converts them to Zebrafish orthologs,
          and stores the resulting genes in the specified file.
    """

    # disease ids and associated genes
    human_disease_genes_df = pd.read_csv(human_disease_genes_path, sep='\t')

    zebrafish_genes = mapping.add_mapped_ortholog_column(human_disease_genes_df, HUMAN_ID, NCBI_ID, keep_old_ids=False, drop_na=True)
    zebrafish_genes.to_csv(dre_disease_genes_path, index=False, sep = '\t')