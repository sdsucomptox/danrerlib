import danrerlib as db
import pandas as pd
from scipy.stats import fisher_exact
import numpy as np


NCBI_ID = 'NCBI Gene ID'

def enrich(gene_universe,
           gene_set,
           concept_type: str, 
           concept_id: str, 
           concept_name: str, 
           sig_cutoff = 0.05,
           log2FC_cutoff = 1.0,
           background_gene_set = None
           ):
    
    if not background_gene_set:
        background_gene_set = pd.DataFrame(gene_universe[NCBI_ID])
        gene_universe = pd.merge(gene_universe, background_gene_set, on=NCBI_ID, how='right', suffixes=('_universe', '_background'))
    
    total_number_of_genes_in_universe = len(gene_universe)

    # Filter genes in the gene set that are in our gene universe
    gene_set = gene_set[gene_set[NCBI_ID].isin(gene_universe[NCBI_ID])]
    gene_set_set = set(gene_set[NCBI_ID])

    if log2FC_cutoff:
        sig_genes_set = set(gene_universe[NCBI_ID][(gene_universe['PValue'] < sig_cutoff) 
                                                     & (np.abs(gene_universe['logFC']) > log2FC_cutoff)])
    else:
        sig_genes_set = set(gene_universe[NCBI_ID][gene_universe['PValue'].lt(sig_cutoff)])

    # Number of genes that are both in the gene set and significantly expressed
    a = len(sig_genes_set.intersection(gene_set_set))

    # Number of genes in the gene set but not significantly expressed
    b = len(gene_set_set.difference(sig_genes_set))

    # Number of genes that are significantly expressed but not in the gene set
    c = len(sig_genes_set.difference(gene_set_set))

    # Number of genes neither in the gene set nor significantly expressed
    d = total_number_of_genes_in_universe - (a + b + c)

    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

    # Determine enrichment direction
    direction = 'enriched' if odds_ratio > 1 else 'depleted'

    # organize important stats
    # ------------------------
    data = {
        'Concept Type': concept_type,
        'Concept ID': concept_id,
        'Concept Name': concept_name,
        '# Genes in Concept in Universe': a+b,
        '# Sig Genes Belong to Concept': a,
        'Proportion of Sig Genes in Set': a/(a+b),
        'Odds Ratio': odds_ratio,
        'P-value': p_value,
        'Direction': direction
    }

    return data

if __name__ == "__main__":
    gene_universe = pd.read_csv('data/test_data/example_diff_express_data.txt', sep = '\t')
    gene_set = pd.read_csv('data/test_data/dre04910.txt', sep = '\t')
    out = enrich(gene_universe, gene_set, 'KEGG', 'blah', 'blah')
    print(out)
    pass