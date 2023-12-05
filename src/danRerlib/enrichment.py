"""
Enrichment Module
=================

The Enrichment module offers a collection of functions to perform gene enrichment analyses. These analyses allow you to identify overrepresented gene sets or concepts in a given set of genes compared to a background or universe of genes. The module supports various enrichment databases, methods, and organisms, making it a versatile tool for uncovering biological insights.

Functions:
    - ``enrich_GO``: Perform Gene Ontology (GO) enrichment analysis.
    - ``enrich_KEGG``: Conduct enrichment analysis using the Kyoto Encyclopedia of Genes and Genomes (KEGG) database.
    - ``logistic``: Perform gene enrichment using logistic regression.
    - ``fishers``: Perform gene enrichment using Fisher's exact test.

Constants:
    - ``NCBI_ID``: Identifier for NCBI Gene ID.
    - ``ZFIN_ID``: Identifier for ZFIN ID.
    - ``ENS_ID``: Identifier for Ensembl ID.
    - ``SYMBOL``: Identifier for gene Symbol.
    - ``HUMAN_ID``: Identifier for Human NCBI Gene ID.

Notes:
    - The Enrichment module is designed for conducting gene enrichment analyses using various databases and methods.
    - It provides functions to analyze gene sets in the context of Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) concepts.
    - Users can choose from different gene ID types and organisms for analysis, enhancing flexibility.
    - The module includes statistical methods such as logistic regression and Fisher's exact test for enrichment analysis.

Example:
    To perform GO enrichment analysis for a set of zebrafish genes:
    
    ``results = enrich_GO(gene_universe, gene_id_type=ZFIN_ID, database='BP', concept_ids=['GO:0007582'], org='dre')``

    This example retrieves genes associated with the specified GO concept ID ('GO:0007582') for zebrafish genes (org='dre').

For comprehensive details on each function and usage examples, please consult the documentation. You can also find tutorials demonstrating the full functionality of the Enrichment module.
"""

# Set up python environment
from danrerlib.settings import *
from danrerlib import mapping, KEGG, GO, utils

import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import scipy.stats
from scipy.stats import chi2, fisher_exact

def logistic(gene_universe: pd.DataFrame, gene_set: pd.DataFrame, gene_id_type: str,
             concept_type: str, concept_id: str, concept_name: str,
             sig_cutoff = 0.05) -> pd.DataFrame:
    """
    Perform gene enrichment analysis using the logistic regression method.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame representing the universe of genes.
        - ``gene_set (pd.DataFrame)``: A DataFrame containing the genes of interest.
        - ``gene_id_type (str)``: The type of gene identifier used in the DataFrames.
        - ``concept_type (str)``: The type of concept (e.g., GO term) being analyzed.
        - ``concept_id (str)``: The ID of the concept being analyzed.
        - ``concept_name (str)``: The name or description of the concept being analyzed.
        - ``sig_cutoff (float, optional)``: The significance cutoff for gene inclusion. Default is 0.05.

    Returns:
        - ``df (pd.DataFrame)``: A DataFrame containing enrichment analysis results, including concept details,
        the number of genes in the concept in the universe, the number of significant genes belonging
        to the concept, the proportion of genes in the concept, the coefficient, p-value, false discovery
        rate (FDR)-adjusted p-value, odds ratio, and enrichment status ('enriched' or 'depleted').

    Notes:
        - This function performs gene enrichment analysis using logistic regression.
        - It calculates enrichment statistics for a specified concept (e.g., GO term) by comparing a gene set
          of interest to a larger gene universe.
        - The 'gene_id_type' parameter specifies the type of gene identifiers used in the DataFrames.
        - The 'sig_cutoff' parameter sets the significance cutoff for gene inclusion based on p-values.
        - Enrichment results include the coefficient (slope) of the logistic regression, p-value,
          FDR-adjusted p-value, odds ratio, and enrichment status.
    """
    # TODO:
    # - include a log2FC cutoff as well?
    gene_set = gene_set[gene_set[gene_id_type].isin(gene_universe[gene_id_type])]

    # create essential dataframe
    # --------------------------

    master_df = gene_universe[[gene_id_type, 'PValue']]
    # determine which genes are in the gene set
    master_df['In Gene Set'] = master_df[gene_id_type].isin(gene_set[gene_id_type]).astype(int)
    # determine which genes have significant differential expression
    master_df['Sig'] = np.where(master_df['PValue'] < sig_cutoff, 1, 0)
    # determine which genes are in gene set and also significant
    master_df['Sig and in Gene Set'] = np.where((master_df['Sig'] == 1) & (master_df['In Gene Set'] == 1), 1, 0)
    if len(gene_set) == 0:
        proportion_of_genes = 0
    else:
        proportion_of_genes = master_df['Sig and in Gene Set'].sum()/len(gene_set)

    if proportion_of_genes != 0 and len(gene_set != 0):

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

        np.seterr(over='ignore')
        # Get the odds ratio
        odds_ratio = np.exp(beta)

        if odds_ratio > 1:
            enriched = 'enriched'
        else:
            enriched = 'depleted'
    
    else:
        # CURRENT DEFAULT PARAMS FOR INVALID ENRICHMENT
        # my thought is these should not even be included in the normal case, but 
        # for my case I wanted a place holder...
        # place holders to avoid errors when finding logistic regression if the
        # values don't exist
        beta = 2
        p_value = 1
        p_values_adjusted = 1
        odds_ratio = 1
        enriched = 'depleted'


    # organize important stats
    # ------------------------

    data = {
        'Concept Type': concept_type,
        'Concept ID': concept_id,
        'Concept Name': concept_name,
        '# Genes in Concept in Universe': len(gene_set),
        '# Sig Genes Belong to Concept': master_df['Sig and in Gene Set'].sum(),
        'Proportion of Genes': proportion_of_genes,
        'Coeff': beta,
        'P-value': p_value,
        'FDR': p_values_adjusted,
        'Odds Ratio': odds_ratio,
        'Enriched': enriched
    }
    df = pd.DataFrame(data, index = [0])
    return df

def fishers_old(gene_universe: pd.DataFrame,
            gene_set: pd.DataFrame, 
            concept_type: str, 
            concept_id: str, 
            concept_name: str, 
            background_gene_list: pd.DataFrame,
            sig_cutoff=0.05):
    """
    Perform gene enrichment analysis using Fisher's exact test.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame representing the universe of genes.
        - ``gene_set (pd.DataFrame)``: A DataFrame containing the genes of interest.
        - ``concept_type (str)``: The type of concept (e.g., KEGG pathway, GO term) being analyzed.
        - ``concept_id (str)``: The ID of the concept being analyzed (e.g., KEGG pathway ID).
        - ``concept_name (str)``: The name or description of the concept being analyzed.
        - ``sig_cutoff (float, optional)``: The significance cutoff for gene inclusion. Default is 0.05.

    Returns:
        - ``df (pd.DataFrame)``: A DataFrame containing enrichment analysis results, including concept details,
        the number of genes in the concept in the universe, the number of significant genes belonging
        to the concept, the proportion of genes in the concept, p-value, odds ratio, and enrichment direction.

    Notes:
        - This function performs gene enrichment analysis using Fisher's exact test.
        - It calculates enrichment statistics for a specified concept (e.g., KEGG pathway) by comparing a gene set
          of interest to a larger gene universe.
        - The 'sig_cutoff' parameter sets the significance cutoff for gene inclusion based on p-values.
        - Enrichment results include p-value, odds ratio, and enrichment direction ('enriched' or 'depleted').
    """

    # Input validation
    if gene_universe.empty or gene_set.empty:
        raise ValueError("Input DataFrames cannot be empty")

    # Filter genes in the gene set that are in our gene universe
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
        'Concept Name': concept_name,
        '# Gene in Concept in Universe': len(gene_set),
        '# Sig Genes Belong to Concept': master_df['Sig and in Gene Set'].sum(),
        'Proportion of Genes': master_df['Sig and in Gene Set'].sum()/len(gene_set),
        'P-value': p_value,
        'Odds Ratio': odds_ratio,
        'Direction': direction
    }

    df = pd.DataFrame([data])

    return df

def produce_gene_universe(gene_id_type, 
                          gene_universe, 
                          background_gene_set,
                          sig_cutoff,
                          log2FC_cutoff):

    if not background_gene_set:
        background_gene_set = pd.DataFrame(gene_universe[gene_id_type])
        gene_universe = pd.merge(gene_universe, background_gene_set, on=gene_id_type, how='right', suffixes=('_universe', '_background'))
    
    total_number_of_genes_in_universe = len(gene_universe)

    if log2FC_cutoff:
        sig_genes_set = set(gene_universe[gene_id_type][(gene_universe['PValue'] < sig_cutoff) 
                                                     & (np.abs(gene_universe['logFC']) > log2FC_cutoff)])
    else:
        sig_genes_set = set(gene_universe[gene_id_type][gene_universe['PValue'].lt(sig_cutoff)])
    
    return total_number_of_genes_in_universe,  sig_genes_set

def fishers(gene_universe: pd.DataFrame,
            sig_genes_set: pd.DataFrame,
           gene_set: pd.DataFrame,
           gene_id_type: str, 
           concept_type: str, 
           concept_id: str, 
           concept_name: str,
           total_number_of_genes_in_universe: int):
    

    # Filter genes in the gene set that are in our gene universe
    gene_set = gene_set[gene_set[gene_id_type].isin(gene_universe[gene_id_type])]
    gene_set_set = set(gene_set[gene_id_type])

    # Number of genes that are both in the gene set and significantly expressed
    a = len(sig_genes_set.intersection(gene_set_set))

    # Number of genes in the gene set but not significantly expressed
    b = len(gene_set_set.difference(sig_genes_set))

    # Number of genes that are significantly expressed but not in the gene set
    c = len(sig_genes_set.difference(gene_set_set))

    # Number of genes neither in the gene set nor significantly expressed
    d = total_number_of_genes_in_universe - (a + b + c)

    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='two-sided')

    # Determine enrichment direction
    direction = 'enriched' if odds_ratio > 1 else 'depleted'

    # Calculate the proportion of significant genes in the set (avoid division by zero)
    proportion_sig_genes_in_set = a / (a + b) if (a + b) > 0 else 0

    # organize important stats
    # ------------------------
    data = {
        'Concept Type': concept_type,
        'Concept ID': concept_id,
        'Concept Name': concept_name,
        '# Genes in Concept in Universe': a+b,
        '# Sig Genes Belong to Concept': a,
        'Proportion of Sig Genes in Set': proportion_sig_genes_in_set,
        'Odds Ratio': odds_ratio,
        'P-value': p_value,
        'Direction': direction
    }

    return data

def enrich(gene_universe: pd.DataFrame, 
            gene_id_type: str,
            database: str, 
            org: 'dre',
            method = 'logistic',
            direction = 'both',
            sig_gene_cutoff_pvalue = 0.05,
            log2FC_cutoff_value = 0,
            concept_ids = None, 
            background_gene_set = None,
            sig_conceptID_cutoff_pvalue = 0.05,
            order_by_p_value = True,
            include_all = False):
    """
    Perform gene enrichment analysis.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame containing gene information, including gene IDs, p-values, and log2FC.
        - ``gene_id_type (str)``: The type of gene ID in the gene universe. The recommended gene id type is NCBI Gene ID (NCBI_ID). Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, Symbol, or for human: Human NCBI Gene ID.
        - ``database (str)``: The functional annotation database. Options include:

                - 'KEGG Pathway': 
                - 'KEGG Disease': 
                - 'GO Biological Processes': 
                - 'GO Cellular Component': 
                - 'GO Molecular Function': 
                - 'all': all databases shown above

        - ``org (str, optional)``: The organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish, 'hsa' for human). Default is 'dre'.
        - ``method (str, optional):`` The enrichment analysis method ('logistic' or 'fishers'). Default is 'logistic'.
        - ``sig_gene_cutoff_pvalue (float, optional)``: The significance cutoff for gene inclusion based on p-values. Default is 0.05.
        - ``concept_ids (list, optional)``: A list of concept IDs (e.g., pathway IDs or disease IDs) to analyze. Default is None.
        - ``sig_conceptID_cutoff_pvalue (float, optional)``: The significance cutoff for concept IDs based on p-values. Default is None.
        - ``sig_conceptID_cutoff_FDR (float, optional)``: The significance cutoff for concept IDs based on FDR (only for 'logistic' method). Default is None.
        - ``order_by_p_value (bool, optional)``: Whether to order the results by p-value. Default is True.

    Returns:
        - ``result (pd.DataFrame)``: A DataFrame containing enrichment analysis results, including concept details, the number of genes in the concept in the universe,
        the number of significant genes belonging to the concept, the proportion of genes in the concept, p-value, odds ratio, and enrichment direction.
    """    

    # QUALITY CONTROL :
    # Provided gene universe must be a dataframe
    if type(gene_universe) != pd.DataFrame:
        gene_universe = pd.DataFrame(gene_universe)
    gene_id_type = utils.normalize_gene_id_type(gene_id_type)

    org = utils.normalize_organism_name(org)
    utils.check_valid_organism(org)

    if org == 'dre' or org == 'dreM':
        utils.check_valid_zebrafish_gene_id_type(gene_id_type)
        _check_gene_universe(gene_universe, gene_id_type, org)
    elif org == 'hsa':
        utils.check_valid_human_gene_id_type(gene_id_type)
        _check_gene_universe(gene_universe, gene_id_type, org)

    # -------------------------
    # DEAL WITH DATABASE CHOICE
    # -------------------------
    # identify concept function to use
    concept_dict = {
        'KEGG Pathway': KEGG.get_genes_in_pathway,
        'KEGG Disease': KEGG.get_genes_in_disease,
        'GO': GO.get_genes_in_GO_concept,
        'GO BP': GO.get_genes_in_GO_concept,
        'GO CC': GO.get_genes_in_GO_concept,
        'GO MF': GO.get_genes_in_GO_concept,
    }

    get_genes_function = concept_dict[database]
    concept_type = database
    all_ids, id_column_name, name_column_name = _get_pathway_ids_and_names(database, org)

    # -------------------------
    # DEAL WITH METHODS CHOICE
    # -------------------------
    # identify enrichment method (function)
    methods_dict = {
        'logistic': logistic,
        'fishers': fishers
    }

    # Check if the method is valid
    if method not in methods_dict:
        raise ValueError(f"Invalid method: {method}")
    
    # Get the function based on the method name
    enrich_method_function = methods_dict[method]

    # DEAL WITH CONCEPT_IDS CHOICE
    if concept_ids is None:
        concept_ids = all_ids[id_column_name]
    else: 
        concept_ids = _check_concept_ids_KEGG(concept_ids, org, database)
        print(concept_ids)

    # -------------------------
    # SIGNIFICANT GENE UNIVERSE
    # -------------------------  

    # the background gene set is all genes in the universe since the gene universe is 
    # all genes detected (not just the significant genes)
    # maybe I should actually do all the genes in the genome .....

    # DEAL WITH THIS EVENTUALLY - I THINK SOMETHING IS MESSED UP
    # if not background_gene_set:
    #     background_gene_set = pd.DataFrame(gene_universe[gene_id_type])
    #     gene_universe = pd.merge(gene_universe, background_gene_set, on=gene_id_type, how='right', suffixes=('_universe', '_background'))

    if background_gene_set:
        gene_universe = _preprocess_gene_universe(gene_universe, background_gene_set, gene_id_type)
        # another option would be to use all genes in the genome. 
    total_number_of_genes_in_universe = len(gene_universe) 
    
    sig_genes_set = _get_sig_genes_set(gene_universe, direction, gene_id_type, 
                                       sig_gene_cutoff_pvalue, log2FC_cutoff_value)

    # -------------------------
    # LAUNCH ENRICHMENT
    # -------------------------
    resulting_dictionary_list = []    
    for concept_id in concept_ids:
        gene_set =  get_genes_function(concept_id, org, do_check = False)
        concept_name = all_ids.loc[all_ids[id_column_name] == concept_id, name_column_name].values[0]
        out = enrich_method_function(gene_universe, 
                                     sig_genes_set,
                                     gene_set, 
                                     gene_id_type,
                                     concept_type, 
                                     concept_id, 
                                     concept_name,
                                     total_number_of_genes_in_universe,
                                     direction)
        
        # Append the dictionary to the list
        resulting_dictionary_list.append(out)
    
    result = pd.DataFrame(resulting_dictionary_list)
    if sig_conceptID_cutoff_pvalue:
        result = result[result["P-value"] <= sig_conceptID_cutoff_pvalue]
    if order_by_p_value:
        result = result.sort_values(by='P-value', ascending=True)
    return result.reset_index(drop=True)

def _preprocess_gene_universe(gene_universe, background_gene_set, gene_id_type, pvalue_col='PValue', logFC_col='logFC'):
    """
    Preprocess the gene universe by adding missing genes from a background set with 'NaN' values for p-value and log2FC.

    Parameters:
    - gene_universe: DataFrame containing information about the subset of genes.
    - background_gene_set: DataFrame containing information about the background genes.
    - gene_id_type: Name of the column containing gene IDs in gene_universe.
    - pvalue_col: Name of the column containing p-values (default: 'PValue').
    - logFC_col: Name of the column containing log2 fold changes (default: 'logFC').

    Returns:
    - DataFrame representing the preprocessed gene universe.
    """

    # Right merge gene_universe with background_gene_set, filling missing values with NaN
    preprocessed_gene_universe = pd.merge(background_gene_set, gene_universe, on=gene_id_type, how='right', suffixes=('_background', '_universe'))

    # Drop extra columns from the background_gene_set (if any)
    preprocessed_gene_universe = preprocessed_gene_universe.drop(columns=[gene_id_type + '_background'])

    return preprocessed_gene_universe

def _get_sig_genes_set(gene_universe, method, gene_id_type, pval_cutoff, log2FC_cutoff):
    # significant genes are defined by the pvalue cutoff choice and the log2FC choice. 

    # Filter genes based on the specified method
    if method == 'up':
        sig_genes_set = set(gene_universe[gene_id_type][(gene_universe['PValue'] < pval_cutoff)
                                                         & (gene_universe['logFC'] > log2FC_cutoff)])
    elif method == 'down':
        sig_genes_set = set(gene_universe[gene_id_type][(gene_universe['PValue'] < pval_cutoff)
                                                         & (gene_universe['logFC'] < -log2FC_cutoff)])
    elif method == 'both':
        sig_genes_set = set(gene_universe[gene_id_type][(gene_universe['PValue'] < pval_cutoff)
                                                        & (np.abs(gene_universe['logFC']) > log2FC_cutoff)])
    else:
        raise ValueError("Invalid method. Supported methods are 'up', 'down', or 'both'.")
    
    return sig_genes_set

def _get_pathway_ids_and_names(database, org):
    if database == 'KEGG Pathway':
        org_pathway_list_dict = {
            'dre' : KEGG.zebrafish_pathways_path,
            'dreM': KEGG.mapped_zebrafish_pathways_path,
            'hsa' : KEGG.human_pathways_path
        }
        path = org_pathway_list_dict[org]
        all_ids = pd.read_csv(path, sep='\t')
        id_column_name = 'Pathway ID'
        name_column_name = 'Pathway Description'
    elif database == 'KEGG Disease':
        path = KEGG.valid_disease_ids_path
        all_ids = pd.read_csv(path, sep='\t')
        id_column_name = 'Disease ID'
        name_column_name = 'Disease Description'
    elif database == 'GO':
        path = GO.GO_IDS_PATH
        all_ids = pd.read_csv(path, sep='\t')
        id_column_name = 'GO ID'
        name_column_name = 'GO Name'
    elif database == 'GO BP':
        path = GO.GO_IDS_PATH
        all_ids = pd.read_csv(path, sep='\t')
        org_col_name = 'exists_dre' if org == 'dre' else 'exists_hsa'
        all_ids = all_ids[(all_ids['Ontology'] == 'P') 
                          & (all_ids[org_col_name] == True)]
        id_column_name = 'GO ID'
        name_column_name = 'GO Name'
    elif database == 'GO CC':
        path = GO.GO_IDS_PATH
        all_ids = pd.read_csv(path, sep='\t')
        org_col_name = 'exists_dre' if org == 'dre' else 'exists_hsa'
        all_ids = all_ids[(all_ids['Ontology'] == 'C') 
                          & (all_ids[org_col_name] == True)]
        id_column_name = 'GO ID'
        name_column_name = 'GO Name'
    elif database == 'GO MF':
        path = GO.GO_IDS_PATH
        all_ids = pd.read_csv(path, sep='\t')
        org_col_name = 'exists_dre' if org == 'dre' else 'exists_hsa'
        all_ids = all_ids[(all_ids['Ontology'] == 'F') 
                          & (all_ids[org_col_name] == True)]
        id_column_name = 'GO ID'
        name_column_name = 'GO Name'
    
    pathway_list = all_ids[[id_column_name, name_column_name]].copy()
    return pathway_list, id_column_name, name_column_name

def enrich_KEGG(gene_universe: str, 
                gene_id_type = NCBI_ID,
                database =  'pathway', 
                concept_ids = None, 
                org = 'dre', 
                method = 'logistic',
                sig_gene_cutoff_pvalue = 0.05,
                sig_conceptID_cutoff_pvalue = None,
                sig_conceptID_cutoff_FDR = None,
                order_by_p_value = True):
    """
    Perform gene enrichment analysis using KEGG pathway or disease databases.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame containing gene information, including gene IDs, p-values, and log2FC.
        - ``gene_id_type (str, optional)``: The type of gene ID in the gene universe. Default is NCBI Gene ID (NCBI_ID).
        - ``database (str, optional)``: The KEGG database to use ('pathway' or 'disease'). Default is 'pathway'.
        - ``concept_ids (list, optional)``: A list of concept IDs (e.g., pathway IDs or disease IDs) to analyze. Default is None.
        - ``org (str, optional)``: The organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish, 'hsa' for human). Default is 'dre'.
        - ``method (str, optional):`` The enrichment analysis method ('logistic' or 'fishers'). Default is 'logistic'.
        - ``sig_gene_cutoff_pvalue (float, optional)``: The significance cutoff for gene inclusion based on p-values. Default is 0.05.
        - ``sig_conceptID_cutoff_pvalue (float, optional)``: The significance cutoff for concept IDs based on p-values. Default is None.
        - ``sig_conceptID_cutoff_FDR (float, optional)``: The significance cutoff for concept IDs based on FDR (only for 'logistic' method). Default is None.
        - ``order_by_p_value (bool, optional)``: Whether to order the results by p-value. Default is True.

    Returns:
        - ``result (pd.DataFrame)``: A DataFrame containing enrichment analysis results, including concept details, the number of genes in the concept in the universe,
        the number of significant genes belonging to the concept, the proportion of genes in the concept, p-value, odds ratio, and enrichment direction.

    Notes:
        - This function performs gene enrichment analysis using the KEGG pathway or disease database.
        - It supports two enrichment analysis methods: 'logistic' and 'fishers'.
        - The 'sig_gene_cutoff_pvalue' parameter sets the significance cutoff for gene inclusion based on gene p-values.
        - The 'sig_conceptID_cutoff_pvalue' parameter can be used to filter concept IDs based on their p-values.
        - The 'sig_conceptID_cutoff_FDR' parameter is applicable when using the 'logistic' method and filters concept IDs based on FDR.
        - The 'order_by_p_value' parameter determines whether to order the results by p-value.
    """
    
    # TODO
    # - make it so that the first column is assumed to be the Gene ID
    #  and the second column is the p-value
    # - support different zebrafish Gene ID inputs
    
    # QUALITY CONTROL :
    # Provided gene universe must be a dataframe
    if type(gene_universe) != pd.DataFrame:
        gene_universe = pd.DataFrame(gene_universe)
    gene_id_type = utils.normalize_gene_id_type(gene_id_type)

    org = utils.normalize_organism_name(org)
    utils.check_valid_organism(org)

    if org == 'dre' or org == 'dreM':
        utils.check_valid_zebrafish_gene_id_type(gene_id_type)
        _check_gene_universe(gene_universe, gene_id_type)
    elif org == 'hsa':
        utils.check_valid_human_gene_id_type(gene_id_type)
        _check_gene_universe(gene_universe, gene_id_type)

    # if gene_id_type != NCBI_ID:
    #     # I would map here
    #     gene_id_to = HUMAN_ID if org == 'hsa' else NCBI_ID
    #     gene_universe = mapping.add_mapped_column(gene_universe, gene_id_type, gene_id_to,
    #                                               keep_old_ids=False, drop_na=True)
    #     gene_id_type = gene_id_to

    # -------------------------
    # DEAL WITH DATABASE CHOICE
    # -------------------------
    # identify concept function to use
    concept_dict = {
        'pathway': KEGG.get_genes_in_pathway,
        'disease': KEGG.get_genes_in_disease
    }

    get_genes_function = concept_dict[database]
    concept_type = 'KEGG ' + database
    
    org_pathway_list_dict = {
        'dre' : KEGG.zebrafish_pathways_path,
        'dreM': KEGG.mapped_zebrafish_pathways_path,
        'hsa' : KEGG.human_pathways_path
    }

    if database == 'pathway':
        path = org_pathway_list_dict[org]
        all_ids = pd.read_csv(path, sep='\t')
        id_column_name = 'Pathway ID'
        name_column_name = 'Pathway Description'
    elif database == 'disease':
        path = KEGG.human_disease_path
        all_ids = pd.read_csv(path, sep='\t')
        id_column_name = 'Disease ID'
        name_column_name = 'Disease Description'
        if concept_ids is None:
            raise ValueError('Testing the full list of diseases is not yet supported.')

    # -------------------------
    # DEAL WITH METHODS CHOICE
    # -------------------------
    # identify enrichment method (function)
    methods_dict = {
        'logistic': logistic,
        'fishers': fishers
    }

    # Get the function based on the method name
    enrich_method_function = methods_dict[method]

    # Check if the method is valid
    if method not in methods_dict:
        raise ValueError(f"Invalid method: {method}")
    
    # DEAL WITH CONCEPT_IDS CHOICE
    if concept_ids is None:
        concept_ids = all_ids[id_column_name]
    else: 
        concept_ids = _check_concept_ids_KEGG(concept_ids, org, database)
   
    # -------------------------
    # LAUNCH ENRICHMENT
    # -------------------------
    resulting_df_list = []    
    for concept_id in concept_ids:
        gene_set =  get_genes_function(concept_id, org)
        if type(gene_set) != pd.DataFrame:
            gene_set = pd.DataFrame(gene_set)
        concept_name = all_ids.loc[all_ids[id_column_name] == concept_id, name_column_name].values[0]
        out = enrich_method_function(gene_universe, gene_set, gene_id_type,
                                     concept_type, concept_id, concept_name,
                                     sig_cutoff = sig_gene_cutoff_pvalue)
        
        # Append the DataFrame to the list
        resulting_df_list.append(out)

    # -------------------------
    # ORGANIZE OUTPUT
    # -------------------------
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

def enrich_GO(gene_universe: str, 
              gene_id_type = ZFIN_ID, 
              database = None, 
              concept_ids = None, 
              org = 'dre', 
              method = 'logistic',
              sig_gene_cutoff_pvalue = 0.05,
              sig_conceptID_cutoff_pvalue = None,
              sig_conceptID_cutoff_FDR = None,
              order_by_p_value = True,
              gene_id_col_name = ZFIN_ID,
              pval_col_name = 'PValue',
              log2fc_col_name = 'log2FC'):
    """
    Perform gene enrichment analysis using Gene Ontology (GO) databases.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame containing gene information, including gene IDs, p-values, and log2FC.
        - ``gene_id_type (str, optional)``: The type of gene ID in the gene universe. Default is ZFIN ID (ZFIN_ID).
        - ``database (str, optional)``: The GO database to use ('BP' for Biological Processes, 'CC' for Cellular Component, 'MF' for Molecular Function). Default is None.
        - ``concept_ids (list, optional)``: A list of GO concept IDs to analyze. Default is None.
        - ``org (str, optional)``: The organism code ('dre' for zebrafish, 'hsa' for human). Default is 'dre'.
        - ``method (str, optional)``: The enrichment analysis method ('logistic' or 'fishers'). Default is 'logistic'.
        - ``sig_gene_cutoff_pvalue (float, optional)``: The significance cutoff for gene inclusion based on p-values. Default is 0.05.
        - ``sig_conceptID_cutoff_pvalue (float, optional)``: The significance cutoff for concept IDs based on p-values. Default is None.
        - ``sig_conceptID_cutoff_FDR (float, optional)``: The significance cutoff for concept IDs based on FDR (only for 'logistic' method). Default is None.
        - ``order_by_p_value (bool, optional)``: Whether to order the results by p-value. Default is True.

    Returns:
        - ``result (pd.DataFrame)``: A DataFrame containing enrichment analysis results, including concept details, the number of genes in the concept in the universe,
        the number of significant genes belonging to the concept, the proportion of genes in the concept, p-value, odds ratio, and enrichment direction.

    Notes:
        - This function performs gene enrichment analysis using the Gene Ontology (GO) databases.
        - It supports two enrichment analysis methods: 'logistic' and 'fishers'.
        - The 'sig_gene_cutoff_pvalue' parameter sets the significance cutoff for gene inclusion based on gene p-values.
        - The 'sig_conceptID_cutoff_pvalue' parameter can be used to filter concept IDs based on their p-values.
        - The 'sig_conceptID_cutoff_FDR' parameter is applicable when using the 'logistic' method and filters concept IDs based on FDR.
        - The 'order_by_p_value' parameter determines whether to order the results by p-value.
    """ 
    # TODO
    # - make it so that the first column is assumed to be the Gene ID
    #  and the second column is the p-value
    # - support different zebrafish Gene ID inputs
    
    # only include genes in the gene set that are in 
    # our gene universe (full gene de)
    if type(gene_universe) != pd.DataFrame:
        gene_universe = pd.DataFrame(gene_universe)
    
    # quality control:
    gene_id_type = utils
    gene_universe = _check_gene_universe(gene_universe, gene_id_type)
    _check_valid_org(org)

    if gene_id_type != ZFIN_ID:
        gene_id_to = HUMAN_ID if org == 'hsa' else ZFIN_ID
        gene_universe = mapping.add_mapped_column(gene_universe, gene_id_type, gene_id_to,
                                                  keep_old_ids=False, drop_na=True)
        gene_id_type = gene_id_to

    # -------------------------
    # DEAL WITH DATABASE CHOICE
    # -------------------------

    all_ids = pd.read_csv(GO.GO_IDS_PATH, sep='\t')
    id_column_name = 'GO ID'
    name_column_name = 'GO Name'
    concept_type_dict = {
        'P': 'GO Biological Processes',
        'C': 'GO Cellular Component',
        'F': 'GO Molecular Function'
    }
    org_col_name = 'exists_dre' if org == 'dre' else 'exists_hsa'
    if database == 'BP':
        all_ids = all_ids[(all_ids['Ontology'] == 'P') 
                          & (all_ids[org_col_name] == True)]
        concept_type = concept_type_dict['P']
    elif database == 'CC':
        all_ids = all_ids[(all_ids['Ontology'] == 'C') 
                          & (all_ids[org_col_name] == True)]
        concept_type = concept_type_dict['C']
    elif database == 'MF':
        all_ids = all_ids[(all_ids['Ontology'] == 'F') 
                          & (all_ids[org_col_name] == True)]
        concept_type = concept_type_dict['F'] 
    elif database == None:
        all_ids = all_ids[all_ids[org_col_name] == True]
        concept_type = 'varies'
    else:
        raise ValueError('Invalid Database')   

    # identify concept function to use
    get_genes_function = GO.get_genes_in_GO_concept

    # # this needs to be adapted for the case when no database is chosen
    # concept_type = 'GO ' + database

    # -------------------------
    # DEAL WITH METHODS CHOICE
    # -------------------------
    # identify enrichment method (function)
    methods_dict = {
        'logistic': logistic,
        'fishers': fishers
    }

    # Get the function based on the method name
    enrich_method_function = methods_dict[method]

    # Check if the method is valid
    if method not in methods_dict:
        raise ValueError(f"Invalid method: {method}")

    # DEAL WITH CONCEPT_IDS CHOICE
    if concept_ids is None:
        concept_ids = all_ids[id_column_name]
    else: 
        concept_ids = _check_concept_ids_GO(concept_ids, org)

    # -------------------------
    # LAUNCH ENRICHMENT
    # -------------------------
    resulting_df_list = []    
    for concept_id in concept_ids:
        gene_set =  get_genes_function(concept_id, org)
        if type(gene_set) != pd.DataFrame:
            gene_set = pd.DataFrame(gene_set)
        if concept_type == 'varies':
            ontology = all_ids.loc[all_ids['GO ID'] == concept_id, 'Ontology'].values[0]
            concept_type = concept_type_dict[ontology]
        concept_name = all_ids.loc[all_ids[id_column_name] == concept_id, name_column_name].values[0]
        out = enrich_method_function(gene_universe, gene_set, gene_id_type,
                                     concept_type, concept_id, concept_name,
                                     sig_cutoff = sig_gene_cutoff_pvalue)
        
        # Append the DataFrame to the list
        resulting_df_list.append(out)

    # -------------------------
    # ORGANIZE OUTPUT
    # -------------------------
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

def _check_gene_universe(gene_universe, specified_gene_id_type, org):

    if org == 'dre' or org == 'dreM':
        required_columns = ['NCBI Gene ID', 'PValue', 'log2FC']
    elif org == 'hsa':
        required_columns = ['Human NCBI Gene ID', 'PValue', 'log2FC']

    existing_columns = gene_universe.columns

    gene_id_type_from_data = existing_columns[0]
    if specified_gene_id_type != gene_id_type_from_data:
        gene_universe = gene_universe.rename(columns={gene_id_type_from_data: specified_gene_id_type})
    
    if specified_gene_id_type != required_columns[0]:
        gene_universe = mapping.add_mapped_column(gene_universe, specified_gene_id_type, required_columns[0], keep_old_ids=False)

    missing_columns = [col for col in required_columns if col not in existing_columns]
    
    # If any required columns are missing, rename existing columns to required names
    # THIS ASSUMES THE FIRST COLUMN IS GENE ID, SECOND COLUMN PVAL, THIRD COLUMN LOG2FC
    if missing_columns:
        # Create a dictionary to map existing column names to the required names
        column_mapping = {col: required_columns[i] for i, col in enumerate(existing_columns)}

        # Rename the columns based on the mapping
        gene_universe = gene_universe.rename(columns=column_mapping)
    
    return gene_universe
    
def _check_concept_ids_KEGG(concept_ids, org, database) -> list:
    if type(concept_ids) == pd.DataFrame:
        if concept_ids.shape[1] != 1:
            raise ValueError('Concept IDs given should be a 1 dimensional list.')
        else:
            column_name = concept_ids.columns[0]
            concept_ids = concept_ids[column_name].to_list()
    if database == 'pathway':
        modified_list = []
        warn = False
        for id in concept_ids:
            new_id = str(id).strip()
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
    elif database == 'disease':
        modified_list = []
        for id in concept_ids:
            new_id = str(id).strip()
            if new_id[0].isdigit():
                id_length = len(new_id)
                if id_length > 5:
                    raise ValueError('Invalid ID')
                elif id_length == 5:
                    new_id = 'H' + new_id
                elif id_length < 5:
                    new_id = 'H' + '0'* 5-id_length + new_id
                else:
                    raise('Invalid ID')
                if id[0] != 'H':
                    raise ValueError('Invalid ID')
            modified_list.append(new_id)
            concept_ids = modified_list
    return concept_ids

def _check_concept_ids_GO(concept_ids, org):
    if type(concept_ids) == pd.DataFrame:
        if concept_ids.shape[1] != 1:
            raise ValueError('Concept IDs given should be a 1 dimensional list.')
        else:
            column_name = concept_ids.columns[0]
            concept_ids = concept_ids[column_name].to_list()
    df = pd.read_csv(GO.GO_IDS_PATH, sep = '\t')
    modified_list = []
    bad_ids = []
    warning = False
    for id in concept_ids:
        new_id = _check_GO_id_format(id)
        if not _id_exists_given_organism(id, org, df):
            bad_ids.append(id)
            warning = True
            # raise ValueError(f'GO ID {id} does not exist for given organism.')
        else:
            modified_list.append(new_id)
    if warning:
        print(f'WARNING: {len(bad_ids)} given GO IDs do not exist for given organism.')
        print()
        print('Omitted ids include:')
        for idx, id in enumerate(bad_ids, start = 1):
            print(id, end=" ")
            if idx % 4 == 0:
                print()  # Print a newline after every fourth ID
    concept_ids = modified_list
    return concept_ids

def _id_exists_given_organism(concept_id, organism, df):

    # check if ID exists for given organism:
    filtered_df = df[df['GO ID'] == concept_id]
    if organism == 'hsa':
        target_column = 'exists_hsa'
    elif organism == 'dre':
        target_column = 'exists_dre'
    elif organism == 'dreM':
        target_column = 'exists_hsa'
    else:
        raise ValueError('Invalid Organism')
    
    return filtered_df[target_column].iloc[0] 
def _check_GO_id_format(id):
    if _is_numeric(id):
        length_of_numeric = len(str(id))
        if length_of_numeric < 7:
            zeros_needed = 7-length_of_numeric
            id = 'GO:' + zeros_needed*'0' + str(id)
        elif length_of_numeric == 7:
            id = 'GO:' + str(id)
        else:
            raise ValueError('Unknown Gene Ontology ID')
    else:
        if len(id) != 10:
            raise ValueError('GO ID should be length 10')
        if id[0:3] != 'GO:':
            raise ValueError('The prefix should be \'GO:\'')
    return id
def _is_numeric(value):
    if type(value) == str:
        return value.isnumeric()
    else:
        return isinstance(value, (int, float))
# def _check_valid_org(org):
#     valid_orgs = ['dre', 'hsa', 'dreM']
#     if org not in valid_orgs:
#         raise ValueError('Invalid organism chosen. Options are: [\'dre\', \'hsa\', \'dreM\']')

# def test_KEGG_enrich(option1  = False, option2 = False, option3 = False):
#     # option 1: all pathways
#     test_data_dir = FILE_DIR / Path('../../tutorials/data/test_data/')
#     file_path = test_data_dir / Path('TPP.txt')
#     tpp_df = pd.read_csv(file_path, sep='\t')
#     if option1:
#         id_type = 'NCBI Gene ID'
#         out = enrich_KEGG(tpp_df, gene_id_type = id_type, org = 'dreM')
#         print(out.head(3))
#     if option2:
#         file_path = test_data_dir / Path('kegg_pathways.txt')
#         pathway_ids = pd.read_csv(file_path, sep='\t')
#         zfish_gene_type = 'NCBI Gene ID' # the gene type in the tpp_df
#         out = enrich_KEGG(tpp_df, 
#                              gene_id_type = zfish_gene_type,
#                              org = 'dreM',
#                              database = 'pathway',
#                              concept_ids = pathway_ids)
#         print(out.head(3))
#     # kegg disease
#     if option3:
#         file_path = test_data_dir / Path('kegg_disease.txt')
#         pathway_ids = pd.read_csv(file_path, sep='\t')
#         zfish_gene_type = 'NCBI Gene ID' # the gene type in the tpp_df
#         out = enrich_KEGG(tpp_df, 
#                         gene_id_type = zfish_gene_type,
#                         org = 'dreM',
#                         database = 'disease',
#                         concept_ids = pathway_ids)
#         print(out)

# def investigation():
#     # kegg_id = '04911'
#     # org = 'dreM'
#     # genes = KEGG.get_genes_in_pathway(kegg_id, org)
#     # print(genes)
#     # print(type(genes[NCBI_ID].values[0]))
#     concept_id = 'H00019'
#     genes = KEGG.get_genes_in_disease(concept_id, 'dreM')
#     # print(genes)
#     # print(type(genes[NCBI_ID].values[0]))


#     test_data_dir = FILE_DIR / Path('../../tutorials/data/test_data/')
#     file_path = test_data_dir / Path('TPP.txt')
#     tpp_df = pd.read_csv(file_path, sep='\t')
#     gene_universe = tpp_df


#     gene = genes[NCBI_ID].values[0]

#     is_in_universe = gene in gene_universe[NCBI_ID].values
#     print(is_in_universe)
    
# def test_GO_enrich(option1, option2, option3, option4, option5):
#     test_data_dir = FILE_DIR / Path('../../tutorials/data/test_data/')
#     file_path = test_data_dir / Path('TPP.txt')
#     tpp_df = pd.read_csv(file_path, sep='\t')
#     zfish_gene_type = 'NCBI Gene ID' # the gene type in the tpp_df
#     # note: these take forever
#     # option 1 all BP
#     if option1:
#         out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
#                                 org = 'dre', database='BP')
#     # option 2 all MF
#     if option2:
#         out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
#                                 org = 'dre', database='BP')
#     # option 3 all CC
#     if option3:
#         out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
#                                 org = 'dre', database='BP')
#     # option 4 all GO!
#     if option4:
#         out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
#                                 org = 'dre', database='BP')
#     # option 5 given list
#     if option4:
#         file_path = test_data_dir / Path('go_ids.txt')
#         pathway_ids = pd.read_csv(file_path, sep='\t')
#         out = enrich_GO(tpp_df, 
#                                 gene_id_type = zfish_gene_type,
#                                 org = 'dre',
#                                 concept_ids = pathway_ids)
#         print(out)


# def testing():

#     # cwd = Path().absolute() 
#     # test_data_path = cwd / Path('tutorials/data/test_data/example_diff_express_data.txt')
#     # gene_univese_full_data = pd.read_csv(test_data_path, sep='\t')

#     # concept_ids = ['dreM00010', 'dreM00020']
#     # concept_ids = ['10', '20']
#     # org = 'dreM'
#     # df = _check_concept_ids(concept_ids, org)
#     # print(df)
#     # test_KEGG_enrich(option3  = True)
#     # investigation()

#     return None

# if __name__ == '__main__':
#     testing()