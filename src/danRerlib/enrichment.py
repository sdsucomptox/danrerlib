# Set up python environment
import numpy as np
import pandas as pd

import mapping_zfin
import KEGG
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import scipy.stats
from scipy.stats import chi2
from scipy.stats import fisher_exact

def logistic(gene_universe, gene_set, concept_type, concept_id, sig_cutoff = 0.05):

    # only include genes in the gene set that are in 
    # our gene universe (full gene de)
    gene_set = gene_set[gene_set['NCBI Gene ID'].isin(gene_universe['NCBI Gene ID'])]

    # create essential dataframe
    # --------------------------

    master_df = gene_universe[['NCBI Gene ID', 'PValue']]
    # determine which genes are in the gene set
    master_df['In Gene Set'] = master_df['NCBI Gene ID'].isin(gene_set['NCBI Gene ID']).astype(int)
    # determine which genes have significant differential expression
    master_df['Sig'] = np.where(master_df['PValue'] < sig_cutoff, 1, 0)
    # determine which genes are in gene set and also significant
    master_df['Sig and in Gene Set'] = np.where((master_df['Sig'] == 1) & (master_df['In Gene Set'] == 1), 1, 0)

    # add important for test to dataframe
    # -----------------------------------

    # Y is defined as 1 for genes in gene set, and 0 
    # for all other genes
    master_df['Y'] = master_df['In Gene Set']
    master_df['x'] = -np.log10(master_df['PValue'])

    # perform logistic regression
    # ---------------------------   

    log_reg = smf.logit("Y ~ x", data=master_df).fit(disp=0)
    # Access the summary results
    summary = log_reg.summary()

    # get important stats
    # -------------------

    # Get the beta coefficient (slope)
    beta = log_reg.params['x']

    # Get the p-value
    p_value = log_reg.pvalues['x']

    # Calculate FDR-adjusted p-values
    p_values_adjusted = smm.multipletests(p_value, method='fdr_bh')[1]

    # Get the odds ratio
    odds_ratio = np.exp(beta)

    if odds_ratio > 1:
        direction = 'enriched'
    else:
        direction = 'depleted'

    # organize important stats
    # ------------------------

    data = {
        'Concept Type': concept_type,
        'Concept ID': concept_id,
        '# Gene in Concept in Universe': len(gene_set),
        '# Sig Genes Belong to Concept': master_df['Sig and in Gene Set'].sum(),
        'Proportion of Genes': master_df['Sig and in Gene Set'].sum()/len(gene_set),
        'Coeff': beta,
        'P-value': p_value,
        'FDR': p_values_adjusted,
        'Odds Ratio': odds_ratio,
        'Direction': direction
    }

    df = pd.DataFrame(data)

    return df

def fishers(gene_universe, gene_set, concept_type, concept_id, sig_cutoff = 0.05):

    # only include genes in the gene set that are in 
    # our gene universe (full gene de)
    gene_set = gene_set[gene_set['NCBI Gene ID'].isin(gene_universe['NCBI Gene ID'])]

    # create essential dataframe
    # --------------------------

    master_df = gene_universe[['NCBI Gene ID', 'PValue']]
    # determine which genes are in the gene set
    master_df['In Gene Set'] = master_df['NCBI Gene ID'].isin(gene_set['NCBI Gene ID']).astype(int)
    # determine which genes have significant differential expression
    master_df['Sig'] = np.where(master_df['PValue'] < sig_cutoff, 1, 0)
    # determine which genes are in gene set and also significant
    master_df['Sig and in Gene Set'] = np.where((master_df['Sig'] == 1) & (master_df['In Gene Set'] == 1), 1, 0)

    # get key parameters
    # ------------------

    # N is the total number of genes tested
    N = len(master_df)

    # n is the number of significantly expressed genes
    n = master_df['Sig'].sum()

    # m is the number of genes in the gene set that are
    # in the gene universe
    m = len(gene_set)

    # k is the number of differentially expressed genes
    # in the gene set
    k = master_df['Sig and in Gene Set'].sum()

    # Run Statistical Test
    # --------------------

    # Create a DataFrame with the data
    data = {'DE': [k, n-k], 'non-DE': [m-k, N+k-n-m]}
    df = pd.DataFrame(data, index=['Inside Gene Set', 'Outside Gene Set'])

    # Create the contingency table
    contingency_table = df.values

    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')

    if odds_ratio > 1:
        direction = 'enriched'
    else:
        direction = 'depleted'

    # organize important stats
    # ------------------------

    data = {
        'Concept Type': concept_type,
        'Concept ID': concept_id,
        '# Gene in Concept in Universe': len(gene_set),
        '# Sig Genes Belong to Concept': master_df['Sig and in Gene Set'].sum(),
        'Proportion of Genes': master_df['Sig and in Gene Set'].sum()/len(gene_set),
        'P-value': p_value,
        'Odds Ratio': odds_ratio,
        'Direction': direction
    }

    df = pd.DataFrame(data, index = [0])

    return df

def perform_enrichment(concept_type: str, concept_ids: list, gene_universe_path: str, method = 'logistic'):
    # GENE UNIVERSE

    gene_universe = pd.read_csv(gene_universe_path, sep='\t')
    gene_universe['NCBI Gene ID'] = gene_universe['NCBI Gene ID'].astype(str)

    # GENE SET
    resuling_df_list = []    
    for concept_id in concept_ids:
        gene_set, human_only = KEGG.download_pathway(concept_id)
        if human_only:
            # get the zebrafish orthologs
            gene_set = mapping_zfin.add_ortholog_column(gene_set, 'Human NCBI Gene ID', 'ZFIN ID')
            gene_set = mapping_zfin.add_mapped_column(gene_set, 'ZFIN ID', 'NCBI Gene ID')
        else:
            gene_set = mapping_zfin.add_mapped_column(gene_set, 'ZFIN ID', 'NCBI Gene ID')

        if method == 'logistic':
            out = logistic(gene_universe, gene_set, concept_type, concept_id)
        elif method == 'fishers':
            out = fishers(gene_universe, gene_set, concept_type, concept_id)

        
        # Append the DataFrame to the list
        resuling_df_list.append(out)

    # Concatenate the list of DataFrames into a single DataFrame
    result = pd.concat(resuling_df_list, ignore_index=True)
    return result