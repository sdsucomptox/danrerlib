from settings import *

import mapping

zebrafish_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_dre_V' + str(VERSION_NUM) + '.txt')
human_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_hsa_V' + str(VERSION_NUM) + '.txt')
mapped_zebrafish_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_dreM_V' + str(VERSION_NUM) + '.txt')

zebrafish_pathways_dir = KEGG_DATA_DIR / Path('dre/')
human_pathways_dir = KEGG_DATA_DIR / Path('hsa/')
dre_mapped_dir = KEGG_DATA_DIR / Path('dreM/')

human_disease_path = KEGG_DATA_DIR / Path('disease_ids_V'+ str(VERSION_NUM) + '.txt')

def get_genes_in_pathway(pathway_id, org = None):
    '''
    organism options:
        dre : tru danio rerio pathways from KEGG
        hsa : true human pathways from KEGG
        dreM : mapped danio rerio pathways from human
    '''
    try:
        org, pathway_id = _check_for_organism(pathway_id, org)
        _check_if_pathway_id_exists(pathway_id, org)

        file_name = KEGG_DATA_DIR / Path(org) / Path(pathway_id + '.txt')
        df = pd.read_csv(file_name, sep='\t')
        return df
    except ValueError:
        pass

def get_genes_in_disease(disease_id, org = 'dreM'):
    '''
    organism options:
        hsa - will return human genes in the disease
        dre - will return human genes mapped to zebrafish in the disease
        dreM - will return human genes mapped to zebrafish in the disease
    note: dre and dreM produce the same result as human disease is not 
          characterized for zebrafish on KEGG
    There are so many disease pathways I opted to just use the API
    instead of saving all the genes in each pathway to a file
    '''
    API = True
    if API:
        url = 'https://rest.kegg.jp/link/hsa/' + disease_id
        column_names = ['trash', HUMAN_ID]
        gene_df = pd.read_csv(url, names=column_names, sep='\t')
        genes = gene_df[HUMAN_ID].str[4:]
        if org == 'hsa':
            return genes
        elif org == 'dre' or org == 'dreM':
            genes = mapping.convert_to_zebrafish(genes, NCBI_ID)
            return genes
    
def _check_for_organism(pathway_id, org):
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

def _check_if_pathway_id_exists(pathway_id, org):
    result = False
    if org == 'dreM':
        df = pd.read_csv(mapped_zebrafish_pathways_path, sep='\t')
        result = True if pathway_id in df['Pathway ID'].values else False
    elif org == 'hsa':
        df = pd.read_csv(human_pathways_path, sep='\t')
        result = True if pathway_id in df['Pathway ID'].values else False
    elif org == 'dre':
        df = pd.read_csv(zebrafish_pathways_path, sep='\t')
        result = True if pathway_id in df['Pathway ID'].values else False
    if not result:
        print('The Pathway ID you gave does not exist for the organism you specified.')
        print('A Pathway ID is identified by the combination of a 2-4 letter prefix')
        print('code and a 5 digit number.\n')
        print('Prefix options include: dre, dreM, or hsa\n')
        print('reminder - python does not support integers that start with a 0;')
        print('           therefore, Pathway IDs must be entered as strings if')
        print('           the organism prefix is omitted.')
    return result

# DATABASE BUILDING FUNCTIONS
# ---------------------------

def build_kegg_database():
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    this takes a while, be careful:)
    '''

    # pathway ids
    _download_zebrafish_pathway_ids(zebrafish_pathways_path)
    _download_human_pathway_ids(human_pathways_path)
    _create_mapped_zebrafish_pathway_ids(mapped_zebrafish_pathways_path, human_pathways_path)

    # get genes in true kegg pathways
    _download_zebrafish_pathway_genes(zebrafish_pathways_path, zebrafish_pathways_dir)
    _download_human_pathway_genes(human_pathways_path, human_pathways_dir)

    # get mapped human kegg pathways to zebrafish
    _build_dre_mapped(human_pathways_path, human_pathways_dir, dre_mapped_dir)

    # get kegg disease ids
    _download_disease_ids(human_disease_path)


def _download_zebrafish_pathway_ids(zebrafish_pathways_path):
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    download the list of zebrafish pathways and pathway names
    '''
    url = 'https://rest.kegg.jp/list/pathway/dre'
    names=['Pathway ID', 'Pathway Description']
    zebrafish_pathways = pd.read_csv(url, names=names, sep='\t')
    zebrafish_pathways['Pathway Description'] = zebrafish_pathways['Pathway Description'].str.replace(' - Danio rerio (zebrafish)', '', regex=False)
    zebrafish_pathways.to_csv(zebrafish_pathways_path, sep='\t', index=False)

def _download_human_pathway_ids(human_pathways_path):
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    '''
    url = 'https://rest.kegg.jp/list/pathway/hsa'
    names=['Pathway ID', 'Pathway Description']
    human_pathways = pd.read_csv(url, names=names, sep='\t')
    human_pathways['Pathway Description'] = human_pathways['Pathway Description'].str.replace(' - Homo sapiens (human)', '', regex=False)
    human_pathways.to_csv(human_pathways_path, sep='\t', index=False)

def _create_mapped_zebrafish_pathway_ids(mapped_zebrafish_pathways_path, human_pathways_path):
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    '''    
    df = pd.read_csv(human_pathways_path, sep='\t')
    df['Pathway ID'] = df['Pathway ID'].str.replace('hsa', 'dreM')
    df.to_csv(mapped_zebrafish_pathways_path, sep='\t', index=False)

def _download_disease_ids(human_disease_path):
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    '''
    url = 'https://rest.kegg.jp/list/disease'
    names=['Disease ID', 'Disease Description']
    human_diseases = pd.read_csv(url, names=names, sep='\t')
    human_disease_path = KEGG_DATA_DIR / Path('disease_ids_V' + str(VERSION_NUM) + '.txt')
    human_diseases.to_csv(human_disease_path, sep='\t', index=False)

def _download_human_pathway_genes(human_pathways_path, human_pathways_dir):
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    '''
    human_pathways = pd.read_csv(human_pathways_path, sep='\t')
    for pathway_id in human_pathways['Pathway ID']:
        url = 'https://rest.kegg.jp/link/hsa/'+pathway_id
        column_names = ['trash', HUMAN_ID]
        genes_df = pd.read_csv(url, names = column_names, sep='\t')
        genes = genes_df[HUMAN_ID].str[4:]
        file_path = human_pathways_dir / Path(str(pathway_id) + '.txt')
        genes.to_csv(file_path, sep='\t', index=False)

def _download_zebrafish_pathway_genes(zebrafish_pathways_path, zebrafish_pathways_dir):
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    '''
    zebrafish_pathways = pd.read_csv(zebrafish_pathways_path, sep='\t')
    for pathway_id in zebrafish_pathways['Pathway ID']:
        url = 'https://rest.kegg.jp/link/dre/'+pathway_id
        column_names = ['trash', NCBI_ID]
        genes_df = pd.read_csv(url, names = column_names, sep='\t')
        genes = genes_df[NCBI_ID].str[4:]
        file_path = zebrafish_pathways_dir / Path(str(pathway_id) + '.txt')
        genes.to_csv(file_path, sep='\t', index=False)

def _build_dre_mapped(human_pathways_path, human_pathways_dir, dre_mapped_dir):
    '''
    A DATABASE BUILDING FUNCTION
    ONLY RUN ON VERSION UPDATE!!
    '''
    human_pathways_df = pd.read_csv(human_pathways_path, sep='\t')
    for pathway_id in human_pathways_df['Pathway ID']:
        in_file_name = human_pathways_dir / Path(pathway_id+'.txt')
        humman_genes = pd.read_csv(in_file_name, sep='\t')
        zebrafish_genes = mapping.convert_to_zebrafish(humman_genes, NCBI_ID, keep_missing_orthos=False)
        stripped_pathway_id = pathway_id[3:]
        out_file_name = dre_mapped_dir / Path(stripped_pathway_id+'.txt')
        zebrafish_genes.to_csv(out_file_name, sep='\t', index=False)

# def is_dre_pathway(pathway_id):
#     id_str = str(pathway_id)
#     if not id_str.startswith('0'):
#         pathway_id = '0'+ id_str
#     pathway_id = 'dre'+str(pathway_id)
#     df = pd.read_csv('database/kegg_ids_dre.txt', sep='\t')
#     if pathway_id in df['id']:
#         return True
#     else:
#         return False

# def _get_only_id(df: pd.DataFrame) -> pd.Series:
#     id_column = df['Pathway ID']
#     ids = id_column.str[3:]
#     return ids

# def _download_url(url_to_download: str, file_path: str, column_names = None) -> None:
#     response = None
#     file_to_save = None

#     try:
#         request = urllib.request.Request(url_to_download)
#         response = urllib.request.urlopen(request)

#         file_to_save = open(file_path, 'wb')

#         if column_names != None:
#             header_string = '\t'.join(column_names) + '\n'
#             file_to_save.write(header_string.encode())

#         # Because response.read() returns a bytes object and because we
#         # opened the file with the 'wb' option, we can write those bytes
#         # directly to the file without first decoding them to a
#         # string.
#         file_to_save.write(response.read())

#     except urllib.error.HTTPError as e:
#         print('Failed to download contents of URL')
#         print('Status code: {}'.format(e.code))
#         print()

#     finally:
#         if file_to_save != None:
#             file_to_save.close()
        
#         if response != None:
#             response.close()



def testing():
    id = 'H00001'
    print(get_genes_in_disease(id, 'dre'))

if __name__ == '__main__':
    testing()