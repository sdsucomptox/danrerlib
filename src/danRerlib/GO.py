"""
GO Module
===========

The Gene Ontology (GO) module provides functions for retrieving gene information associated with 
Gene Ontology Biological Process (BP), Molecularfunctionsns (MF) and Cellular Components (CC)
for various organisms, including human, zebrafish, and mapped zebrafish.

Functions:
    - ``get_genes_in_GO_concept``: Retrieve genes associated with a specific GO concept ID.

Constants:
    - ``NCBI_ID``: Identifier for NCBI Gene ID.
    - ``ZFIN_ID``: Identifier for ZFIN ID.
    - ``ENS_ID``: Identifier for Ensembl ID.
    - ``SYMBOL``: Identifier for gene Symbol.
    - ``HUMAN_ID``: Identifier for Human NCBI Gene ID.

Database Rebuild Functions:
    - ``build_gene_ontology_database``: Build or update the gene ontology database.

Notes:
    - This module is designed for accessing gene information from Gene Ontology.
    - It provides functions to retrieve gene data associated with biological processes, cellular components, and molecular functions.
    - Data is obtained either through reading pre-processed files from the database build.
    - Organism options include 'dre' (Danio rerio), 'hsa' (Human), or 'dreM' (Mapped Danio rerio from Human).

Example:
    To retrieve genes associated with a specific GO concept ID:
    ```
    go_genes = get_genes_in_GO_concept('GO:0007582', 'dre')
    ```

For detailed information on each function and their usage, please refer to the documentation. For more examples of full functionality, please refer to tutorials.
"""

from danrerlib.settings import *
from danrerlib import mapping, utils
from typing import Optional
import pandas as pd

GO_IDS_PATH = GO_DATA_DIR / Path('GO_ids_V' + str(VERSION_NUM) + '.txt')
GO_PATH_dre = GO_DATA_DIR / Path('GO_dre_V' + str(VERSION_NUM) + '.txt')
GO_PATH_dreM = GO_DATA_DIR / Path('GO_dreM_V' + str(VERSION_NUM) + '.txt')
GO_PATH_hsa = GO_DATA_DIR / Path('GO_hsa_V' + str(VERSION_NUM) + '.txt')

GO_BASIC_URL = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
GO_ZFIN_URL = 'https://current.geneontology.org/annotations/zfin.gaf.gz'
GO_NCBI_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'

def get_genes_in_GO_concept(concept_id: str, 
                            organism: str, 
                            gene_id_type: Optional[str] = None,
                            out_data_type: Optional[type] = pd.DataFrame,
                            do_check = True
                            ) -> pd.Series:
    """
    Retrieve gene IDs associated with a Gene Ontology (GO) concept for a specified organism.

    Parameters:
        - ``concept_id (str)``: The Gene Ontology (GO) concept ID for which gene IDs are to be retrieved.
        - ``organism (str)``: The organism for which gene IDs should be retrieved. Options include

                       - 'hsa': Human.
                       - 'dre': Zebrafish.
                       - 'dreM': Mapped Zebrafish from Human (same as 'dre' as Zebrafish GO IDs are not characterized).

        - ``gene_id_type (str, optional)``: The desired gene ID type for the returned gene IDs. Default is None.

    Returns:
       - ``gene_ids_series (pd.Series)``: A pandas Series containing gene IDs associated with the specified GO concept.

    Raises:
        ValueError: If the GO concept ID does not exist for the given organism or if the organism is invalid.

    Notes:
        - This function retrieves gene IDs associated with a Gene Ontology (GO) concept for a specified organism.
        - The 'concept_id' parameter should be a valid GO concept identifier.
        - The 'organism' parameter specifies the organism for which gene IDs should be retrieved.
        - Options for 'organism' include 'hsa' (Human), 'dre' (Zebrafish), or 'dreM' (Mapped Zebrafish from Human).

    Example:
        To retrieve Zebrafish gene IDs associated with a GO concept:
        ```
        concept_id = 'GO:0001234'
        organism = 'dre'
        gene_ids = get_genes_in_GO_concept(concept_id, organism)
        ```
    """
    try:
        if do_check:
            # check given organism format and raise exception if invalid
            organism = utils.normalize_organism_name(organism)
            utils.check_valid_organism(organism)

            # check given ID format and raise exception if invalid
            if gene_id_type:
                gene_id_type = utils.normalize_gene_id_type(gene_id_type)
                if organism == 'hsa':
                    utils.check_valid_human_gene_id_type(gene_id_type)
                else:
                    utils.check_valid_zebrafish_gene_id_type(gene_id_type)

            # check if ID is in required format:
            concept_id = _check_id_format(concept_id)
            # check to see if the GO id exists for given organism
            if not id_exists_given_organism(concept_id, organism):
                raise ValueError('GO ID does not exist for given organism.')
            
        # the default gene id types for the GO database
        if organism == 'hsa':
            path = GO_PATH_hsa
            gene_id_type_from_db = HUMAN_ID
        elif organism == 'dre':
            path = GO_PATH_dre
            gene_id_type_from_db = ZFIN_ID
        elif organism == 'dreM':
            path = GO_PATH_dreM
            gene_id_type_from_db = ZFIN_ID
        else:
            raise ValueError('Invalid organism.')
        
        df = pd.read_csv(path, sep = '\t')
        filtered_df = df[df['GO ID'] == concept_id]

        gene_ids_series = filtered_df[gene_id_type_from_db]
        if gene_id_type and (organism != 'hsa'):
            # map to desired gene id type
            if gene_id_type != gene_id_type_from_db:
                gene_ids_series = mapping.convert_ids(gene_ids_series, 
                                gene_id_type_from_db, gene_id_type)
        gene_ids_series = gene_ids_series.reset_index(drop=True)

        if out_data_type == pd.Series:
            return gene_ids_series
        elif out_data_type == pd.DataFrame:
            return pd.DataFrame(gene_ids_series)
        elif out_data_type == list:
            return gene_ids_series.tolist()
        else:
            return gene_ids_series
        
    except utils.InvalidGeneTypeError as e:
        pass

def id_exists_given_organism(concept_id: str, organism: str):
    """
    Check if a Gene Ontology (GO) ID exists for the specified organism.

    Parameters:
        - ``concept_id (str)``: The GO ID to be validated.
        - ``organism (str)``: The specified organism code (e.g., 'hsa', 'dre', or 'dreM').

    Returns:
        - bool: True if the GO ID exists for the specified organism, False otherwise.

    Raises:
        - ValueError: If the organism code is invalid.

    Notes:
        - This function checks whether a given GO ID exists for the specified organism.
        - It reads a GO ID data file and looks for a match in the 'GO ID' column.
        - The 'organism' parameter specifies the organism code for which to perform the check.
        - Returns True if the GO ID exists for the organism, False otherwise.
        - Provides an error message if the organism code is invalid.
    """

    # check if ID exists for given organism:
    df = pd.read_csv(GO_IDS_PATH, sep = '\t')
    filtered_df = df[df['GO ID'] == concept_id]
    if organism == 'hsa':
        target_column = 'exists_hsa'
    elif organism == 'dre':
        target_column = 'exists_dre'
    elif organism == 'dreM':
        target_column = 'exists_hsa'
    else:
        raise ValueError('Invalid Organism')
    
    return filtered_df[target_column].iloc[0]

def build_gene_ontology_database():
    """
    Build or update the Gene Ontology (GO) database.

    Notes:
        - This function is intended for database creation or update and should be run only during version updates.
        - The process may take some time, so exercise caution when running it.
        - The function performs the following steps:
        
            1. Downloads all current GO IDs.
            2. Builds the human GO file and with IDs and associated genes.
            3. Builds the zebrafish GO file and with IDs and associated genes.
            5. Builds the mapped zebrafish GO file and with IDs and associated genes.
            6. Add a column to the database file to identify if the GO ID exists for human and/or zebrafish. 

        - Running this function should be done carefully, as it involves downloading and processing data fromGO.
    """

    _build_GO_IDs(GO_BASIC_URL, GO_IDS_PATH)
    _build_human_GO(GO_PATH_hsa, GO_NCBI_URL)
    _build_zebrafish_GO(GO_PATH_dre, GO_ZFIN_URL)
    _build_dre_mapped_GO(GO_PATH_hsa, GO_PATH_dreM)
    _add_to_GO_IDs(GO_IDS_PATH, GO_PATH_dre, GO_PATH_hsa)

def _build_GO_IDs(GO_BASIC_URL, GO_IDS_PATH):
    '''
    This function generates a list of all GO IDs, the Name of GO IDs, and the 
    ontology it belongs to. This is downloaded directly from the GO website using 
    the GO API.
    '''
    url = GO_BASIC_URL
    output_file_for_processed = GO_IDS_PATH
    content = utils.download_url(url)
    if content:
        go_terms = _process_GOBasic_content(content)
        _write_GOBasic_to_file(go_terms, output_file_for_processed)
    
def _add_to_GO_IDs():
    '''
    the purpose of this function is to add on two columns to identify if
    the GO concept exists for zebrafish and humans.
    '''
    df = pd.read_csv(GO_IDS_PATH, sep='\t')
    df_zfish = pd.read_csv(GO_PATH_dre, sep='\t')
    df_human = pd.read_csv(GO_PATH_hsa, sep='\t')

    df['exists_hsa'] = df['GO ID'].isin(df_human['GO ID'])
    df['exists_dre'] = df['GO ID'].isin(df_zfish['GO ID'])
    utils.save_data(df, GO_IDS_PATH)

def _build_zebrafish_GO():
    column_names = ['Database Designation',
                    'Marker ID',
                    'Gene Symbol',
                    'Qualifiers',
                    'GO Term ID',
                    'Reference ID',
                    'GO Evidence Code',
                    'Inferred From',
                    'Ontology',
                    'Marker Name',
                    'Marker Synonyms',
                    'Marker Type',
                    'Taxon',
                    'Modification Date',
                    'Assigned By',
                    'Annotation Extension', 
                    'Gene Product Form ID']
    df = pd.read_csv(GO_ZFIN_URL, sep='\t', comment='!', names=column_names, low_memory=False)
    desired_columns = ['Marker ID', 'GO Term ID', 'Ontology']
    df = df.loc[:, desired_columns].drop_duplicates()
    new_column_names = {'Marker ID': 'ZFIN ID',
                        'GO Term ID': 'GO ID'}
    df.rename(columns=new_column_names, inplace=True)
    utils.save_data(df, GO_PATH_dre)

def _build_human_GO():
    df = pd.read_csv(GO_NCBI_URL, sep='\t', compression='gzip')
    df = df[df['#tax_id']==9606]
    desired_columns = ['GeneID', 'GO_ID', 'Category']
    df = df.loc[:, desired_columns].drop_duplicates()
    new_column_names = {'GeneID': 'Human NCBI Gene ID',
                        'GO_ID': 'GO ID',
                        'Category' : 'Orthology'}
    df.rename(columns=new_column_names, inplace=True)
    df['Orthology'] = df['Orthology'].apply(_get_ontology)
    utils.save_data(df, GO_PATH_hsa)

def _build_dre_mapped_GO():
    df = pd.read_csv(GO_PATH_hsa, sep = '\t')
    ortho_df = mapping.add_mapped_ortholog_column(df, 'Human NCBI Gene ID', 'ZFIN ID', keep_old_ids=False, drop_na=True)
    utils.save_data(ortho_df, GO_PATH_dreM)

def _get_ontology(namespace):
    if (namespace == 'biological_process'
        or namespace == 'Process') :
        return 'P'
    elif (namespace == 'molecular_function'
          or namespace == 'Function'):
        return 'F'
    elif (namespace == 'cellular_component'
          or namespace == 'Component'):
        return 'C'
    else:
        return None
    
def _process_GOBasic_content(content: list) -> list:
    go_terms = []
    current_term = None

    for line in content:
        line = line.strip().decode()  # Decode bytes to string

        if line.startswith('[Term]'):
            if current_term is not None:
                # Check if the current term has all the required fields
                if all(field is not None for field in current_term.values()):
                    go_terms.append(current_term)

            current_term = {
                'id': None,
                'name': None,
                'ontology': None
            }
        elif line.startswith('id: '):
            current_term['id'] = line.split(': ')[1]
        elif line.startswith('name: '):
            current_term['name'] = line.split(': ')[1]
        elif line.startswith('namespace: '):
            namespace = line.split(': ')[1]
            current_term['ontology'] = _get_ontology(namespace)

    # Check and add the last term after the loop ends
    if current_term is not None and all(field is not None for field in current_term.values()):
        go_terms.append(current_term)

    return go_terms

def _write_GOBasic_to_file(go_terms: list, output_file: str) -> None:
    with open(output_file, 'w') as f_out:
        f_out.write('GO ID\tGO Name\tOntology\n')  # Header line

        for term in go_terms:
            f_out.write(f'{term["id"]}\t{term["name"]}\t{term["ontology"]}\n')

def _check_id_format(id):
    if _is_numeric(id):
        length_of_numeric = len(str(id))
        if length_of_numeric < 7:
            zeros_needed = 7-length_of_numeric
            id = 'GO:' + zeros_needed*'0' + str(id)
        elif length_of_numeric == 7:
            id = 'GO:' + str(id)
        else:
            raise ValueError('Unknown Gene Ontology ID')
    else:
        if len(id) != 10:
            raise ValueError('GO ID should be length 10')
        if id[0:3] != 'GO:':
            raise ValueError('The prefix should be \'GO:\'')
    return id

def _is_numeric(value):
    if type(value) == str:
        return value.isnumeric()
    else:
        return isinstance(value, (int, float))