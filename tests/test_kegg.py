from danrerlib import KEGG
from danrerlib.settings import *
from pandas.testing import assert_frame_equal
from itertools import permutations

def test_get_genes_from_pathway():

    # human
    pathway_id = 'hsa04010'
    true_gene_ids = _format_from_kegg(pathway_id, 'hsa')
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id, 'hsa')
    assert_frame_equal(true_gene_ids, generated_gene_ids) 
    
    pathway_id = '04010'
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id, 'hsa')
    assert_frame_equal(true_gene_ids, generated_gene_ids) 

    pathway_id = 'hsa04010'
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id)
    assert_frame_equal(true_gene_ids, generated_gene_ids) 

    # zfish
    pathway_id = 'dre04910'
    true_gene_ids = _format_from_kegg(pathway_id, 'dre')
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id, 'dre')
    assert_frame_equal(true_gene_ids, generated_gene_ids) 
    
    pathway_id = '04910'
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id, 'dre')
    assert_frame_equal(true_gene_ids, generated_gene_ids) 

    pathway_id = 'dre04910'
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id)
    assert_frame_equal(true_gene_ids, generated_gene_ids) 

    # zfish mapped
    pathway_id = 'dreM00232'
    true_gene_ids = pd.read_csv('tests/data/out_data/kegg/dreM00232.txt', sep = '\t')
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id, 'dreM')
    assert_frame_equal(true_gene_ids, generated_gene_ids) 
    
    pathway_id = '00232'
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id, 'dreM')
    assert_frame_equal(true_gene_ids, generated_gene_ids) 

    pathway_id = 'dreM00232'
    generated_gene_ids = KEGG.get_genes_in_pathway(pathway_id)
    assert_frame_equal(true_gene_ids, generated_gene_ids) 

def _format_from_kegg(pathway_id, org):
    file_name = 'tests/data/out_data/kegg/'+pathway_id+'.txt'
    from_kegg_directly = pd.read_csv(file_name, sep = '\t', names = ['trash', 'ids'])
    
    strip_str = org+':'
    from_kegg_directly['ids'] = from_kegg_directly['ids'].str.strip(strip_str)
    from_kegg_directly['ids'] = from_kegg_directly['ids'].astype(np.int64)


    from_kegg_directly = from_kegg_directly.drop(columns = ['trash'])
    if org == 'hsa':
        from_kegg_directly = from_kegg_directly.rename(columns={'ids':'Human NCBI Gene ID'})
    else:
        from_kegg_directly = from_kegg_directly.rename(columns={'ids':'NCBI Gene ID'})
    return from_kegg_directly

def test_get_genes_from_disease():

    disease_id = 'H00001'
    file_name = 'tests/data/out_data/kegg/'+disease_id+'.txt'
    true_gene_ids = pd.read_csv(file_name, sep = '\t')

    generated_gene_ids = KEGG.get_genes_in_disease(disease_id, 'hsa')
    assert_frame_equal(true_gene_ids, generated_gene_ids) 

    file_name = 'tests/data/out_data/kegg/'+disease_id+'_dre.txt'
    true_gene_ids = pd.read_csv(file_name, sep = '\t').sort_values(by = NCBI_ID).reset_index(drop = True)

    generated_gene_ids = KEGG.get_genes_in_disease(disease_id, 'dre').sort_values(by = NCBI_ID).reset_index(drop = True)
    assert_frame_equal(true_gene_ids, generated_gene_ids) 



