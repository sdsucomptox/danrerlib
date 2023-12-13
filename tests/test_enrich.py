from danrerlib import enrichment
from danrerlib.settings import *
from pandas.testing import assert_frame_equal
from itertools import permutations
import json
import pytest


# Sample gene_universe DataFrame for testing
gene_universe_data = {
    'GeneID': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
    'PValue': [0.01, 0.05, 0.001, 0.1, 0.04],
    'log2FC': [0.5, -0.8, 1.2, -0.3, -0.1]
}

gene_universe_synthetic = pd.DataFrame(gene_universe_data)
gene_universe_sample = pd.read_csv('tests/data/in_data/example_diff_express_data.txt', sep = '\t')
gene_universe_full = pd.read_csv('tests/data/in_data/01_TPP.txt', sep = '\t')

def test_fishers():
    gene_set = pd.read_csv('tests/data/in_data/enrich/dre04340.txt', sep = '\t')
    sig_genes_df = enrichment._get_sig_genes_df(gene_universe_sample, 'NCBI Gene ID', 0.05, 0.0)
    gene_set = gene_set[gene_set[NCBI_ID].isin(gene_universe_sample[NCBI_ID])]
    out = enrichment._fishers(gene_universe_sample, 
                             sig_genes_df,
                             gene_set, 
                        'NCBI Gene ID', 
                        'KEGG Pathway', 
                        'dre04340',
                        'WNT Signaling Pathway', 
                        len(gene_universe_sample),
                        'two-sided')
    file_path1 = 'tests/data/out_data/enrich/fishers1.txt'
    with open(file_path1, 'r') as file:
        true_dictionary = json.load(file)
    assert out == true_dictionary

def test_logistic():
    gene_set = pd.read_csv('tests/data/in_data/enrich/dre04340.txt', sep = '\t')
    sig_genes_df = enrichment._get_sig_genes_df(gene_universe_sample, 'NCBI Gene ID', 0.05, 0.0)
    gene_set = gene_set[gene_set[NCBI_ID].isin(gene_universe_sample[NCBI_ID])]
    out = enrichment._logistic(gene_universe_sample, 
                             sig_genes_df,
                             gene_set, 
                        'NCBI Gene ID', 
                        'KEGG Pathway', 
                        'dre04340',
                        'WNT Signaling Pathway', 
                        len(gene_universe_sample),
                        'non-directional',
                        0.05,
                        0)
    file_path1 = 'tests/data/out_data/enrich/logistic1.txt'
    with open(file_path1, 'r') as file:
        true_dictionary = json.load(file)
    assert out == true_dictionary

def test_basic_kegg_enrichment():
    out = enrichment.enrich_logistic(gene_universe_full, 
                  gene_id_type='NCBI Gene ID', 
                  database='KEGG Pathway')
    true = pd.read_csv('tests/data/out_data/enrich/KEGG_enrich.txt', sep = '\t')
    assert_frame_equal(out,true)

def test_basic_kegg_enrichment_smaller():
    out = enrichment.enrich_logistic(gene_universe_sample, 
                  gene_id_type='NCBI Gene ID', 
                  database='KEGG Pathway',
                  directional_test=False)
    true = pd.read_csv('tests/data/out_data/enrich/KEGG_enrich_sub.txt', sep = '\t')
    assert_frame_equal(out,true)

def test_go_enrichment_with_concept_ids():
    concept_ids_path = 'tests/data/in_data/enrich/list_of_go_ids.txt'
    concept_ids = pd.read_csv(concept_ids_path, sep = '\t')
    out = enrichment.enrich_logistic(gene_universe_full, 
                    gene_id_type='NCBI Gene ID', 
                    database='GO MF', 
                    directional_test=False, 
                    concept_ids = concept_ids)
    true = pd.read_csv('tests/data/out_data/enrich/GOMF_four_ids.txt', sep = '\t')
    assert_frame_equal(out,true)

def test_get_sig_genes_set_up():

    sig_genes_set = enrichment._get_sig_genes_set(gene_universe_synthetic, 'up', 'GeneID', 0.05, 0.0)
    assert sig_genes_set == {'Gene1', 'Gene3'}

def test_get_sig_genes_set_down():
    sig_genes_set = enrichment._get_sig_genes_set(gene_universe_synthetic, 'down', 'GeneID', 0.05, 0.0)
    assert sig_genes_set == {'Gene5'}

def test_get_sig_genes_set_both():
    sig_genes_set = enrichment._get_sig_genes_set(gene_universe_synthetic, 'both', 'GeneID', 0.05, 0.0)
    assert sig_genes_set == {'Gene1', 'Gene3', 'Gene5'}

kegg_options = ['KEGG Pathway', 'KEGG Disease']
go_options = ['GO', 'GO BP', 'GO CC', 'GO MF']

def test_master_gene_id_no_change():
    gene_universe = pd.read_csv('tests/data/in_data/example_diff_express_data.txt', sep = '\t')
    gene_id_type_out, gene_universe_out = enrichment._map_to_master_geneid(gene_universe, NCBI_ID, 
                                                                       'dre', 'KEGG Pathway', 
                                                                       kegg_options, go_options)
    assert_frame_equal(gene_universe, gene_universe_out) 
    assert NCBI_ID == gene_id_type_out


sample_gene_universe = pd.DataFrame({
    'GeneID': [1, 2, 3, 4, 5],
    'PValue': [0.01, 0.02, 0.03, 0.04, 0.05],
    'Log2FC': [1.0, 0.8, -0.5, 0.2, -0.9]
})

sample_gene_id_type = NCBI_ID
sample_databases = ['KEGG Pathway']

def test_normalize_and_validate_database():
    # Test for a single valid database
    result = enrichment._normalize_and_validate_database('KEGG Pathway')
    assert result == 'KEGG Pathway'

    # Test for a single valid database
    result = enrichment._normalize_and_validate_database('Kegg pathway')
    assert result == 'KEGG Pathway'

    # Test for an invalid database
    with pytest.raises(ValueError):
        enrichment._normalize_and_validate_database('Invalid Database')

    # Add more test cases as needed

def test_expand_database_option():
    # Test for a single valid database
    result = enrichment._expand_database_option('KEGG Pathway')
    assert result == ['KEGG Pathway']

    # Test for 'GO' option
    result = enrichment._expand_database_option('GO')
    assert result == ['GO BP', 'GO CC', 'GO MF']


def test_process_database_options_string():
    input_db = "KEGG Pathway"
    result = enrichment._process_database_options(input_db)
    assert result == ["KEGG Pathway"]

def test_process_database_options_list():
    input_db = ["KEGG Pathway", "GO BP"]
    result = enrichment._process_database_options(input_db)
    assert result == ["KEGG Pathway", "GO BP"]

def test_process_database_options_list_expansion():
    input_db = "GO"
    result = enrichment._process_database_options(input_db)
    assert result == ["GO BP", "GO CC", "GO MF"]

def test_process_database_options_mixed():
    input_db = ["KEGG Pathway", "GO"]
    result = enrichment._process_database_options(input_db)
    assert result == ["KEGG Pathway", "GO BP", "GO CC", "GO MF"]

def test_process_database_options_invalid_type():
    input_db = 123  # Some invalid type
    with pytest.raises(TypeError):
       enrichment._process_database_options(input_db)


def test_get_enrichment_method_valid():
    method_name = 'logistic'
    result = enrichment._get_enrichment_method(method_name)
    assert result == enrichment._logistic

def test_get_enrichment_method_valid_case_insensitive():
    method_name = 'LoGiStIc'
    result = enrichment._get_enrichment_method(method_name)
    assert result == enrichment._logistic

def test_get_enrichment_method_valid_fishers():
    method_name = 'fishers'
    result = enrichment._get_enrichment_method(method_name)
    assert result == enrichment._fishers

def test_get_enrichment_method_invalid():
    with pytest.raises(ValueError, match=r"Invalid method: invalid_method. Supported methods are logistic, fishers."):
        enrichment._get_enrichment_method('invalid_method')