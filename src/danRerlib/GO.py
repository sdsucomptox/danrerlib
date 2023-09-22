from danRerLib.settings import *
import danRerLib.mapping as mapping
import danRerLib.utils as utils

GO_IDS_PATH = GO_DATA_DIR / Path('GO_ids_V' + str(VERSION_NUM) + '.txt')
GO_PATH_dre = GO_DATA_DIR / Path('GO_dre_V' + str(VERSION_NUM) + '.txt')
GO_PATH_dreM = GO_DATA_DIR / Path('GO_dreM_V' + str(VERSION_NUM) + '.txt')
GO_PATH_hsa = GO_DATA_DIR / Path('GO_hsa_V' + str(VERSION_NUM) + '.txt')

GO_BASIC_URL = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
GO_ZFIN_URL = 'https://current.geneontology.org/annotations/zfin.gaf.gz'
GO_NCBI_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'

def get_genes_in_GO_concept(concept_id, organism, gene_id_type = None):
    '''
    
    gene_id_type: the desired gene id type to be returned
    '''
    # check if ID is in required format:
    concept_id = check_id_format(concept_id)
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
            
    return gene_ids_series

def id_exists_given_organism(concept_id, organism):

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
    
def check_id_format(id):
    if is_numeric(id):
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

def is_numeric(value):
    if type(value) == str:
        return value.isnumeric()
    else:
        return isinstance(value, (int, float))

def build_all():

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
    the GO concept exists for zebrafish and himans.
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