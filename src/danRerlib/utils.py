import pandas as pd
import urllib.request
from danrerlib.settings import *
from typing import Union, Optional, List

class InvalidGeneTypeError(Exception):
    "Raised when the Gene ID Type is invalid"
    pass

class InvalidOrganismError(Exception):
    "Raised when the organism is invalid"
    pass

def check_valid_organism(org: str) -> None:
    """
    Check the validity of the chosen organism.

    This function checks if the provided organism is a valid option.

    Parameters:
        org (str): A string representing an organism to be validated.

    Raises:
        InvalidOrganismError: If one or more of the provided gene ID types are invalid.

    Notes:
        - Valid organisms include: hsa, dre, dreM.
    """
    org = ['dre', 'hsa', 'dreM']

    if type(org) == str:
        org = [org]

    invalid_orgs = [item for item in org if item not in org]
    if invalid_orgs:
        print('The organism choice you gave is invalid. The valid organisms for this library are:')
        print('library are: dre, hsa, and dreM. The organism you gave that is invalid is:')
        print('----------')
        for item in invalid_orgs:
            print(item)
        print('----------')
        print('Reminder: The organism is case and spelling sensitive.')
        raise InvalidOrganismError   
    
def normalize_organism_name(organism_name: str) -> str:
    """
    Normalize an organism name to a specified format.

    Parameters:
       - organism_name (str): The organism name to be normalized.

    Returns:
        - organism_name (str): The normalized organism name.
    """
    # Define mappings from common variations to the desired format
    organism_mappings = {
        'human': 'hsa',
        'homo sapiens': 'hsa',
        'hsa' : 'hsa',
        'zebrafish': 'dre',
        'zfish': 'dre',
        'dre' : 'dre',
        'danio rerio' : 'dre',
        'mapped zebrafish': 'dreM',
        'mapped': 'dreM',
    }

    # Strip the input organism name just in case
    stripped_organism_name = organism_name.strip()
    # Convert the input organism name to lowercase for case-insensitive matching
    lowercase_organism_name = stripped_organism_name.lower()

    # Check if the lowercase organism name exists in the mappings
    if lowercase_organism_name in organism_mappings:
        return organism_mappings[lowercase_organism_name]
    else:
        # If no mapping is found, return the original input
        return organism_name
    
def check_valid_zebrafish_gene_id_type(gene_id_types: Union[str, List[str]]) -> None:
    """
    Check the validity of Zebrafish gene ID types.

    This function checks if the provided Zebrafish gene ID types are valid options.

    Parameters:
        - ``gene_id_types (str or list)``: A string or a list of Zebrafish gene ID types to be validated.

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

def check_valid_human_gene_id_type(gene_id_types: Union[str, List[str]]) -> None:
    """
    Check the validity of Zebrafish gene ID types.

    This function checks if the provided Zebrafish gene ID types are valid options.

    Parameters:
        - ``gene_id_types (str or list)``: A string or a list of Zebrafish gene ID types to be validated.

    Raises:
        InvalidGeneTypeError: If one or more of the provided gene ID types are invalid.

    Notes:
        - Valid Zebrafish gene ID types include: NCBI Gene ID, ZFIN ID, Ensembl ID, or Symbol.
        - The input can be a single gene ID type as a string or multiple types in a list.
        - Gene ID types are case and spelling sensitive.
    """
    gene_id_options = [HUMAN_ID]

    if type(gene_id_types) == str:
        gene_id_types = [gene_id_types]

    invalid_gene_id_types = [item for item in gene_id_types if item not in gene_id_options]
    if invalid_gene_id_types:
        print('The human Gene ID type you gave is invalid. The only valid Gene ID option')
        print(f'is: {HUMAN_ID}. The ID(s) you gave that are invalid are:')
        print('----------')
        for item in invalid_gene_id_types:
            print(item)
        print('----------')
        print('Reminder: Gene ID types are case and spelling sensitive.')
        raise InvalidGeneTypeError
     
def normalize_gene_id_type(gene_id_type: str) -> str:
    """
    Normalize a gene ID type to a specified format.

    Parameters:
       - ``gene_id_type (str)``: The gene ID type to be normalized.

    Returns:
        - ``gene_id_type (str)``: The normalized gene ID type.
    """
    # Define mappings from common variations to the desired format
    id_type_mappings = {
        'symbol': SYMBOL,
        'sym': SYMBOL,
        'genesymbol': SYMBOL,
        'gene symbol': SYMBOL,
        'ncbi gene id': NCBI_ID,
        'ncbi_gene_id': NCBI_ID,
        'ncbi': NCBI_ID,
        'ncbi id': NCBI_ID,
        'ncbi_id': NCBI_ID,
        'ncbiid': NCBI_ID,
        'ncbigeneid': NCBI_ID,
        'entrez':NCBI_ID,
        'entrez id': NCBI_ID,
        'entrez gene id': NCBI_ID, 
        'entrezgeneid': NCBI_ID,
        'zebrafish ncbi': NCBI_ID,
        'zfin': ZFIN_ID,
        'zfin id': ZFIN_ID,
        'zfinid': ZFIN_ID,
        'zfingeneid': ZFIN_ID,
        'zfin_id': ZFIN_ID,
        'ensembl_id': ENS_ID,
        'ensembl id': ENS_ID,
        'ensembl gene id': ENS_ID,
        'ensemblid': ENS_ID,
        'ensemblgeneid': ENS_ID,
        'ens': ENS_ID,
        'ensid': ENS_ID,
        'ens gene id': ENS_ID,
        'ensgeneid': ENS_ID,
        'ensembl': ENS_ID,
        'human id': HUMAN_ID,
        'human ncbi gene id': HUMAN_ID, 
        'human ncbi': HUMAN_ID,
    }

    # Strip the gene id string just in case
    stripped_gene_id_type = gene_id_type.strip()
    # Convert the input gene ID type to lowercase for case-insensitive matching
    lowercase_gene_id_type = stripped_gene_id_type.lower() 

    # Check if the lowercase gene ID type exists in the mappings
    if lowercase_gene_id_type in id_type_mappings:
        return id_type_mappings[lowercase_gene_id_type]
    else:
        # If no mapping is found, return the original input
        return gene_id_type

def pretty_print_series(series: pd.Series):
    for element in series:
        print(element)

def save_data(data: pd.Series or pd.DataFrame, location):
    data.to_csv(location, sep = '\t', index=False)

def download_url(url_to_download: str) -> list:
    response = None

    try:
        request = urllib.request.Request(url_to_download)
        response = urllib.request.urlopen(request)

        return response.readlines()

    except urllib.error.HTTPError as e:
        print('Failed to download contents of URL')
        print('Status code: {}'.format(e.code))
        print()

    finally:
        if response != None:
            response.close()

    return []

def download_url_to_file(url_to_download: str, output_file_path: str) -> None:
    response = None
    file_to_save = None

    try:
        request = urllib.request.Request(url_to_download)
        response = urllib.request.urlopen(request)

        file_to_save = open(output_file_path, 'wb')
        file_to_save.write(response.read())

    except urllib.error.HTTPError as e:
        print('Failed to download contents of URL')
        print('Status code: {}'.format(e.code))
        print()

    finally:
        if file_to_save != None:
            file_to_save.close()
        
        if response != None:
            response.close()

def load_data(data_file_name):
    file_dict = {
        'raw human gene info': RAW_DATA_DIR / Path('Homo_sapiens.gene_info'),
        'raw orthology': RAW_DATA_DIR/ Path(f'zfish_human_orthology_V{VERSION_NUM}.txt'),
        'master orthology': DATABASE_DIR / Path(f'master_ortho_mapping_file_V{VERSION_NUM}.txt'),
        'master mapping': DATABASE_DIR / Path(f'master_gene_mapping_file_V{VERSION_NUM}.txt'),
        'GO dre': DATABASE_DIR / Path(f'GO/GO_dre_V{VERSION_NUM}.txt'),
        'GO dreM': DATABASE_DIR / Path(f'GO/GO_dreM_V{VERSION_NUM}.txt'),
        'GO hsa': DATABASE_DIR / Path(f'GO/GO_hsa_V{VERSION_NUM}.txt')
    }
    return pd.read_csv(file_dict[data_file_name], sep = '\t')