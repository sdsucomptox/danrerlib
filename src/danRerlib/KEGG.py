import pandas as pd

def download_human_pathway(pathway_id):
    '''
    input:  pathway_id
            - KEGG human pathway id, staring with hsa
    output: gene_list
            - list of genes in NCBI format title HumanGeneID
    '''
    # read in data from the KEGG API URL
    full_url = 'https://rest.kegg.jp/link/hsa/'+pathway_id
    gene_list = pd.read_csv(full_url, sep = '\t', names = ["Pathway ID", "Human NCBI Gene ID"])
    gene_list.drop(columns = ['Pathway ID'], inplace=True)
    gene_list['Human NCBI Gene ID'] = gene_list['Human NCBI Gene ID'].str.replace('hsa:', '')
#     gene_list['HumanGeneID'] = gene_list['HumanGeneID'].astype(int)
    return gene_list

def download_zebrafish_pathway(pathway_id):
    '''
    input:  pathway_id
            - KEGG zebrafish pathway id, staring with dre
    output: gene_list
            - list of genes in NCBI format title GeneID
    '''
    # read in data from the KEGG API URL
    full_url = 'https://rest.kegg.jp/link/dre/'+pathway_id
    gene_list = pd.read_csv(full_url, sep = '\t', names = ["Pathway ID", "NCBI Gene ID"])
    gene_list.drop(columns = ['Pathway ID'], inplace=True)
    gene_list['NCBI Gene ID'] = gene_list['NCBI Gene ID'].str.replace('dre:', '')
    # gene_list['NCBI Gene ID'] = gene_list['NCBI Gene ID'].astype(int)
    return gene_list

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