import pandas as pd
import urllib.error
import urllib.request
from settings import *
import mapping

def build_kegg_database():
    '''
    this takes a while, be careful:)
    '''
    zebrafish_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_dre_V' + str(VERSION_NUM) + '.txt')
    human_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_hsa_V' + str(VERSION_NUM) + '.txt')
    zebrafish_pathways_dir = KEGG_DATA_DIR / Path('dre/')
    human_pathways_dir = KEGG_DATA_DIR / Path('hsa/')
    human_disease_path = KEGG_DATA_DIR / Path('disease_ids_V'+ str(VERSION_NUM) + '.txt')
    dre_mapped_dir = KEGG_DATA_DIR / Path('dre_mapped/')

    # kegg pathway
    # _download_zebrafish_pathway_ids(zebrafish_pathways_path)
    # _download_human_pathway_ids(human_pathways_path)
    # _download_zebrafish_pathway_genes(zebrafish_pathways_path, zebrafish_pathways_dir)
    # _download_human_pathway_genes(human_pathways_path, human_pathways_dir)

    # map human kegg pathways to zebrafish
    _build_dre_mapped(human_pathways_path, human_pathways_dir, dre_mapped_dir)

    # # kegg disease
    # _download_disease_ids(human_disease_path)


def _download_zebrafish_pathway_ids(zebrafish_pathways_path):
    '''
    download the list of zebrafish pathways and pathway names
    '''
    url = 'https://rest.kegg.jp/list/pathway/dre'
    names=['Pathway ID', 'Pathway Description']
    zebrafish_pathways = pd.read_csv(url, names=names, sep='\t')
    zebrafish_pathways['Pathway Description'] = zebrafish_pathways['Pathway Description'].str.replace(' - Danio rerio (zebrafish)', '', regex=False)
    zebrafish_pathways.to_csv(zebrafish_pathways_path, sep='\t', index=False)

def _download_human_pathway_ids(human_pathways_path):
    url = 'https://rest.kegg.jp/list/pathway/hsa'
    names=['Pathway ID', 'Pathway Description']
    human_pathways = pd.read_csv(url, names=names, sep='\t')
    human_pathways_path = KEGG_DATA_DIR / Path('pathway_ids_hsa_V' + str(VERSION_NUM) + '.txt')
    human_pathways['Pathway Description'] = human_pathways['Pathway Description'].str.replace(' - Homo sapiens (human)', '', regex=False)
    human_pathways.to_csv(human_pathways_path, sep='\t', index=False)

def _download_disease_ids(human_disease_path):
    url = 'https://rest.kegg.jp/list/disease'
    names=['Disease ID', 'Disease Description']
    human_diseases = pd.read_csv(url, names=names, sep='\t')
    human_disease_path = KEGG_DATA_DIR / Path('disease_ids_V' + str(VERSION_NUM) + '.txt')
    human_diseases.to_csv(human_disease_path, sep='\t', index=False)

def _download_human_pathway_genes(human_pathways_path, human_pathways_dir):
    human_pathways = pd.read_csv(human_pathways_path, sep='\t')
    for pathway_id in human_pathways['Pathway ID']:
        url = 'https://rest.kegg.jp/link/hsa/'+pathway_id
        column_names = ['trash', HUMAN_ID]
        genes_df = pd.read_csv(url, names = column_names, sep='\t')
        genes = genes_df[HUMAN_ID].str[4:]
        file_path = human_pathways_dir / Path(str(pathway_id) + '.txt')
        genes.to_csv(file_path, sep='\t', index=False)

def _download_zebrafish_pathway_genes(zebrafish_pathways_path, zebrafish_pathways_dir):
    zebrafish_pathways = pd.read_csv(zebrafish_pathways_path, sep='\t')
    for pathway_id in zebrafish_pathways['Pathway ID']:
        url = 'https://rest.kegg.jp/link/dre/'+pathway_id
        column_names = ['trash', NCBI_ID]
        genes_df = pd.read_csv(url, names = column_names, sep='\t')
        genes = genes_df[NCBI_ID].str[4:]
        file_path = zebrafish_pathways_dir / Path(str(pathway_id) + '.txt')
        genes.to_csv(file_path, sep='\t', index=False)

def _build_dre_mapped(human_pathways_path, human_pathways_dir, dre_mapped_dir):
    human_pathways_df = pd.read_csv(human_pathways_path, sep='\t')
    for pathway_id in human_pathways_df['Pathway ID']:
        in_file_name = human_pathways_dir / Path(pathway_id+'.txt')
        humman_genes = pd.read_csv(in_file_name, sep='\t')
        zebrafish_genes = mapping.convert_to_zebrafish(humman_genes, NCBI_ID, keep_missing_orthos=False)
        stripped_pathway_id = pathway_id[3:]
        out_file_name = dre_mapped_dir / Path(stripped_pathway_id+'.txt')
        zebrafish_genes.to_csv(out_file_name, sep='\t', index=False)

def get_genes_in_disease(disease_id, API = True):
    '''
    There are so many disease pathways I opted to just use the API
    instead of saving all the genes in each pathway to a file
    '''
    if API:
        url = 'https://rest.kegg.jp/link/hsa/' + disease_id
        column_names = ['trash', HUMAN_ID]
        gene_df = pd.read_csv(url, names=column_names, sep='\t')
        gene_df[HUMAN_ID] = gene_df[HUMAN_ID].str[4:]
        return gene_df[HUMAN_ID]

def get_genes_in_pathway(pathway_id, org):
    file_name = KEGG_DATA_DIR / Path(org + pathway_id + '.txt')
    return pd.read_csv(file_name, sep='\t')

def is_dre_pathway(pathway_id):
    id_str = str(pathway_id)
    if not id_str.startswith('0'):
        pathway_id = '0'+ id_str
    pathway_id = 'dre'+str(pathway_id)
    df = pd.read_csv('database/kegg_ids_dre.txt', sep='\t')
    if pathway_id in df['id']:
        return True
    else:
        return False
    
def download_pathway(pathway_id):
    pathway_id = str(pathway_id)
    if not pathway_id.startswith('0'):
        pathway_id = '0'+ pathway_id

    if is_dre_pathway(pathway_id):
        human_only = False
        pathway_id = 'dre'+str(pathway_id)
        gene_list = download_zebrafish_pathway(pathway_id)
    else:
        human_only = True
        pathway_id = 'hsa'+str(pathway_id)
        gene_list = download_human_pathway(pathway_id)

    return gene_list, human_only

def _get_only_id(df: pd.DataFrame) -> pd.Series:
    id_column = df['Pathway ID']
    ids = id_column.str[3:]
    return ids

def _download_url(url_to_download: str, file_path: str, column_names = None) -> None:
    response = None
    file_to_save = None

    try:
        request = urllib.request.Request(url_to_download)
        response = urllib.request.urlopen(request)

        file_to_save = open(file_path, 'wb')

        if column_names != None:
            header_string = '\t'.join(column_names) + '\n'
            file_to_save.write(header_string.encode())

        # Because response.read() returns a bytes object and because we
        # opened the file with the 'wb' option, we can write those bytes
        # directly to the file without first decoding them to a
        # string.
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

def testing():
    # this is for testing purposed
    url = 'https://rest.kegg.jp/list/pathway/dre'
    file_path = KEGG_DATA_DIR/ Path('trial.txt')
    header = 'trash\thello\n'
    _download_url(url, file_path, header)
    return None


if __name__ == '__main__':
    build_kegg_database()