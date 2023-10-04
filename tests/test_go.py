from danrerlib import GO
from danrerlib.settings import *
from pandas.testing import assert_frame_equal
from itertools import permutations

def test_get_genes_from_go_id():
    go_id = 'GO:0005173'
    org = 'dre'
    generated_genes = GO.get_genes_in_GO_concept(go_id, org, 'sym')
    generated_genes = generated_genes.to_frame().sort_values(by='Symbol', ascending=True).reset_index(drop=True)

    true_genes = pd.read_csv('tests/data/out_data/go/id_dre_0005173.txt')
    true_genes = true_genes.sort_values(by='Symbol', ascending=True).reset_index(drop=True)
    assert_frame_equal(true_genes, generated_genes) 

    org = 'human'
    generated_genes = GO.get_genes_in_GO_concept(go_id, org, 'human id')
    generated_genes = generated_genes.to_frame().sort_values(by='Human NCBI Gene ID', ascending=True).reset_index(drop=True)

    true_genes = pd.read_csv('tests/data/out_data/go/id_hsa_0005173.txt')
    true_genes = true_genes.sort_values(by='Human NCBI Gene ID', ascending=True).reset_index(drop=True)
    assert_frame_equal(true_genes, generated_genes) 
