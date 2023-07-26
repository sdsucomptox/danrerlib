from settings import *
import mapping, utils

GO_IDS_PATH = GO_DATA_DIR / Path('GO_ids_V' + str(VERSION_NUM) + '.txt')
GO_PATH_dre = GO_DATA_DIR / Path('GO_dre_V' + str(VERSION_NUM) + '.txt')
GO_PATH_hsa = GO_DATA_DIR / Path('GO_hsa_V' + str(VERSION_NUM) + '.txt')

GO_BASIC_URL = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
GO_ZFIN_URL = 'https://current.geneontology.org/annotations/zfin.gaf.gz'
GO_NCBI_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'

def build_GO_IDs():
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
    
    
# add on the human and zebrafish columns
def add_to_GO_IDs():
    df = pd.read_csv(GO_IDS_PATH, sep='\t')
    df_zfish = pd.read_csv(GO_PATH_dre, sep='\t')
    df_human = pd.read_csv(GO_PATH_hsa, sep='\t')

    df['exists_hsa'] = df['GO ID'].isin(df_human['GO ID'])
    df['exists_dre'] = df['GO ID'].isin(df_zfish['GO ID'])
    utils.save_data(df, GO_IDS_PATH)

def build_zebrafish_GO():
    url = GO_ZFIN_URL
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
    df = pd.read_csv(url, sep='\t', comment='!', names=column_names, low_memory=False)
    desired_columns = ['Marker ID', 'GO Term ID', 'Ontology']
    df = df.loc[:, desired_columns].drop_duplicates()
    new_column_names = {'Marker ID': 'ZFIN ID',
                        'GO Term ID': 'GO ID'}
    df.rename(columns=new_column_names, inplace=True)
    utils.save_data(df, GO_PATH_dre)

def build_human_GO():
    url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'
    df = pd.read_csv(url, sep='\t', compression='gzip')
    df = df[df['#tax_id']==9606]
    desired_columns = ['GeneID', 'GO_ID', 'Category']
    df = df.loc[:, desired_columns].drop_duplicates()
    new_column_names = {'GeneID': 'Human NCBI Gene ID',
                        'GO_ID': 'GO ID',
                        'Category' : 'Orthology'}
    df.rename(columns=new_column_names, inplace=True)
    df['Orthology'] = df['Orthology'].apply(_get_ontology)
    utils.save_data(df, GO_PATH_hsa)

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

add_to_GO_IDs()