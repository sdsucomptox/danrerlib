from danrerlib import enrichment
from danrerlib.settings import *
from pandas.testing import assert_frame_equal
from itertools import permutations
import json


def test_fishers():
    gene_universe = pd.read_csv('tests/data/in_data/example_diff_express_data.txt', sep = '\t')
    # background_gene_set = pd.DataFrame(gene_universe['NCBI Gene ID'])
    kegg_pathway = pd.read_csv('tests/data/in_data/enrich/dre04340.txt', sep = '\t')
    sig_genes_set = enrichment._get_sig_genes_set(gene_universe, 'both', 'NCBI Gene ID', 0.05, 0.0)
    out = enrichment.fishers(gene_universe, 
                             sig_genes_set,
                             kegg_pathway, 
                        'NCBI Gene ID', 
                        'KEGG Pathway', 
                        'dre04340',
                        'WNT Signaling Pathway', 
                        len(gene_universe))
    file_path1 = 'tests/data/out_data/enrich/fishers1.txt'
    with open(file_path1, 'r') as file:
        true_dictionary = json.load(file)
    assert out == true_dictionary

# Sample gene_universe DataFrame for testing
gene_universe_data = {
    'GeneID': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
    'PValue': [0.01, 0.05, 0.001, 0.1, 0.04],
    'logFC': [0.5, -0.8, 1.2, -0.3, -0.1]
}

gene_universe = pd.DataFrame(gene_universe_data)

def test_get_sig_genes_set_up():

    sig_genes_set = enrichment._get_sig_genes_set(gene_universe, 'up', 'GeneID', 0.05, 0.0)
    assert sig_genes_set == {'Gene1', 'Gene3'}

def test_get_sig_genes_set_down():
    sig_genes_set = enrichment._get_sig_genes_set(gene_universe, 'down', 'GeneID', 0.05, 0.0)
    assert sig_genes_set == {'Gene5'}

def test_get_sig_genes_set_both():
    sig_genes_set = enrichment._get_sig_genes_set(gene_universe, 'both', 'GeneID', 0.05, 0.0)
    assert sig_genes_set == {'Gene1', 'Gene3', 'Gene5'}

