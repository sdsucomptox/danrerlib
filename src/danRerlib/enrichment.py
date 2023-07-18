# Set up python environment
from settings import *

import mapping
import KEGG

import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import scipy.stats
from scipy.stats import chi2
from scipy.stats import fisher_exact

def logistic(gene_universe: pd.DataFrame, gene_set: pd.DataFrame, 
             concept_type: str, concept_id: str, 
             sig_cutoff = 0.05) -> pd.DataFrame:
    '''
    This function performs gene enrichment using the logistic regression method. A 
    p-value cutoff is required for this method and the default is taken to be 0.05
    '''
    # TODO:
    # - include a log2FC cutoff as well?
    
    gene_set = gene_set[gene_set[NCBI_ID].isin(gene_universe[NCBI_ID])]

    # create essential dataframe
    # --------------------------

    master_df = gene_universe[[NCBI_ID, 'PValue']]
    # determine which genes are in the gene set
    master_df['In Gene Set'] = master_df[NCBI_ID].isin(gene_set[NCBI_ID]).astype(int)
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
        enriched = 'enriched'
    else:
        enriched = 'depleted'

    # organize important stats
    # ------------------------

    data = {
        'Concept Type': concept_type,
        'Concept ID': concept_id,
        '# Genes in Concept in Universe': len(gene_set),
        '# Sig Genes Belong to Concept': master_df['Sig and in Gene Set'].sum(),
        'Proportion of Genes': master_df['Sig and in Gene Set'].sum()/len(gene_set),
        'Coeff': beta,
        'P-value': p_value,
        'FDR': p_values_adjusted,
        'Odds Ratio': odds_ratio,
        'Enriched': enriched
    }

    df = pd.DataFrame(data)

    return df

def fishers(gene_universe, gene_set, concept_type, concept_id, sig_cutoff = 0.05):
    '''
    gene_universe - all genes 
    gene_set - genes in the concept of interest
    concept_type - the database (KEGG, GO, etc)
    concept_id - the id (e.g. the KEGG pathway id)
    '''

    # only include genes in the gene set that are in 
    # our gene universe (full gene de)
    gene_set = gene_set[gene_set[NCBI_ID].isin(gene_universe[NCBI_ID])]

    # create essential dataframe
    # --------------------------

    master_df = gene_universe[[NCBI_ID, 'PValue']]
    # determine which genes are in the gene set
    master_df['In Gene Set'] = master_df[NCBI_ID].isin(gene_set[NCBI_ID]).astype(int)
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

def enrich_KEGG(gene_universe: str, gene_id_type = NCBI_ID,
                database =  'pathway', concept_ids = None, 
                org = 'dre', method = 'logistic',
                sig_gene_cutoff_pvalue = 0.05,
                sig_conceptID_cutoff_pvalue = None,
                sig_conceptID_cutoff_FDR = None,
                order_by_p_value = True):
    
    # TODO
    # - make it so that the first column is assumed to be the Gene ID
    #  and the second column is the p-value
    # - support different zebrafish Gene ID inputs
    
    # only include genes in the gene set that are in 
    # our gene universe (full gene de)
    if type(gene_universe) != pd.DataFrame:
        gene_universe = pd.DataFrame(gene_universe)
    
    # quality control TODO:
    _check_gene_universe(gene_universe)

    if gene_id_type != NCBI_ID:
        raise ValueError('NCBI Gene ID is the only currently supported type.')

    # quality control concept ids:
    _check_valid_org(org)

    # identify concept
    concept_dict = {
        'pathway': KEGG.get_genes_in_pathway,
        'disease': KEGG.get_genes_in_disease
    }

    get_genes_function = concept_dict[database]
    concept_type = 'KEGG ' + database

    # identify enrichment method
    methods_dict = {
        'logistic': logistic,
        'fishers': fishers
    }

    # Get the function based on the method name
    enrich_method_function = methods_dict[method]

    # Check if the method is valid
    if method not in methods_dict:
        raise ValueError(f"Invalid method: {method}")
    
    org_pathway_list_dict = {
        'dre' : KEGG.zebrafish_pathways_path,
        'dreM': KEGG.mapped_zebrafish_pathways_path,
        'hsa' : KEGG.human_pathways_path
    }
    
    if concept_ids is None:
        if database == 'pathway':
            path = org_pathway_list_dict[org]
        elif database == 'disease':
            raise ValueError('Testing the full list of diseases is not yet supported.')
        
        all_ids = pd.read_csv(path, sep='\t')
        concept_ids = all_ids['Pathway ID']
    else: 
        concept_ids = _check_concept_ids(concept_ids, org)

    resulting_df_list = []    
    for concept_id in concept_ids:
        gene_set =  get_genes_function(concept_id, org)
        if type(gene_set) != pd.DataFrame:
            gene_set = pd.DataFrame(gene_set)

        out = enrich_method_function(gene_universe, gene_set, 
                                     concept_type, concept_id,
                                     sig_cutoff = sig_gene_cutoff_pvalue)
        
        # Append the DataFrame to the list
        resulting_df_list.append(out)

    # Concatenate the list of DataFrames into a single DataFrame
    result = pd.concat(resulting_df_list, ignore_index=True)
    if sig_conceptID_cutoff_pvalue:
        result = result[result["P-value"] <= sig_conceptID_cutoff_pvalue]
    if method == 'logistic' and sig_conceptID_cutoff_FDR:
        result = result[result["FDR"] <= sig_conceptID_cutoff_FDR]
    if order_by_p_value:
        result = result.sort_values(by='P-value', ascending=True)  
    result = result.reset_index(drop=True)
    return result

def _check_gene_universe(gene_universe):
    if (NCBI_ID not in gene_universe.columns 
        and 'PValue' not in gene_universe.columns):
        raise ValueError('Required columns do not exist.')
    
def _check_concept_ids(concept_ids, org):
    if type(concept_ids) == pd.DataFrame:
        if concept_ids.shape[1] != 1:
            raise ValueError('Concept IDs given should be a 1 dimensional list.')
        else:
            column_name = concept_ids.columns[0]
            concept_ids = concept_ids[column_name].to_list()

    modified_list = []
    warn = False
    for id in concept_ids:
        new_id = str(id)
        if new_id[0].isdigit():
            id_len = len(new_id)
            if id_len < 5:
                new_id = org + '0' * (5-id_len) + new_id
            elif id_len == 5:
                pass
            else:
                raise ValueError('The ID given is greater than 5 numeric values.')
        elif new_id[0].isalpha():
            if len(new_id) >= 5 and new_id[-5:].isdigit():
                org_from_id = new_id[:-5]
                if org_from_id != org:
                    new_id = org + new_id[-5:]
                    warn = True
        else:
            raise ValueError
        modified_list.append(new_id)
    if warn:
        print('WARNING: The organism you specified does not match the prefix')
        print(f'         in the given concept IDs. The organism chosen, {org}, is')
        print('         used, not the organism in the concept IDs.')
    concept_ids = modified_list
    return concept_ids

def _check_valid_org(org):
    valid_orgs = ['dre', 'hsa', 'dreM']
    if org not in valid_orgs:
        raise ValueError('Invalid organism chosen. Options are: [\'dre\', \'hsa\', \'dreM\']')


def testing():

    # cwd = Path().absolute() 
    # test_data_path = cwd / Path('tutorials/data/test_data/example_diff_express_data.txt')
    # gene_univese_full_data = pd.read_csv(test_data_path, sep='\t')

    concept_ids = ['dreM00010', 'dreM00020']
    concept_ids = ['10', '20']
    org = 'dreM'
    df = _check_concept_ids(concept_ids, org)
    print(df)

    return None

if __name__ == '__main__':
    testing()