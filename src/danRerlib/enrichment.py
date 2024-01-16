"""
Enrichment Module
=================

The Enrichment module offers a collection of functions to perform gene enrichment analyses. These analyses allow you to identify overrepresented gene sets or concepts in a given set of genes compared to a background or universe of genes. The module supports various enrichment databases, methods, and organisms, making it a versatile tool for uncovering biological insights.

Functions:
    - ``enrich_logistic``: Perform functional enrichment analysis using logistic regression for any of the supported annotation databases:
                
                - 'KEGG Pathway': Kyoto Encyclopedia of Genes and Genomes Pathways
                - 'KEGG Disease': Kyoto Encyclopedia of Genes and Genomes Diseases
                - 'GO BP': Gene Ontology Biological Processes
                - 'GO CC': Gene Ontology Cellular Components
                - 'GO MF': Gene Ontology Molecular Function
                - 'GO': Gene Ontology Biological Processes, Cellular Components, and Molecular Function
    
    - ``enrich_fishers``: Perform functional enrichment analysis using Fisher's exact test for any of the supported annotation databases:

                - 'KEGG Pathway': Kyoto Encyclopedia of Genes and Genomes Pathways
                - 'KEGG Disease': Kyoto Encyclopedia of Genes and Genomes Diseases
                - 'GO BP': Gene Ontology Biological Processes
                - 'GO CC': Gene Ontology Cellular Components
                - 'GO MF': Gene Ontology Molecular Function
                - 'GO': Gene Ontology Biological Processes, Cellular Components, and Molecular Function
    
Private Functions (used by enrich_logistic and enrich_fishers):

    - ``_enrich``: Launch enrichment and perform quality control.
    - ``_enrich_variable_org``: Launch enrichment and perform quality for a variable organism.
    - ``_logistic``: Perform statistical test using logistic regression.
    - ``_fishers``: Perform statistical test using Fisher's exact test.

Constants:
    - ``NCBI_ID``: Identifier for NCBI Gene ID.
    - ``ZFIN_ID``: Identifier for ZFIN ID.
    - ``ENS_ID``: Identifier for Ensembl ID.
    - ``SYMBOL``: Identifier for gene Symbol.
    - ``HUMAN_ID``: Identifier for Human NCBI Gene ID.

Notes:
    - The Enrichment module is designed for conducting functional enrichment analyses using various databases and methods.
    - It provides functions to analyze gene sets in the context of Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) concepts.
    - Users can choose from different gene ID types and organisms for analysis, enhancing flexibility.
    - The module includes statistical methods such as logistic regression and Fisher's exact test for enrichment analysis.

Example:
    To perform gene enrichment analysis for a set of zebrafish genes:
    
    ``results = enrich_logistic(gene_universe, gene_id_type=ZFIN_ID, database='KEGG Pathway', org='dre')``

    This example performs functional enrichment analysis using logistic regression on the specified gene set for zebrafish genes (org='dre').

For comprehensive details on each function and usage examples, please consult the documentation. You can also find tutorials demonstrating the full functionality of the Enrichment module.
"""

# Set up python environment
from danrerlib.settings import *
from danrerlib import mapping, KEGG, GO, utils

import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
from statsmodels.tools.sm_exceptions import ConvergenceWarning, PerfectSeparationWarning

from scipy.stats import chi2, fisher_exact
import warnings

def enrich_fishers(gene_universe: pd.DataFrame, 
            database: list[str], 
            gene_id_type: str,
            org = 'dre',
            direction = 'non-directional',
            sig_gene_cutoff_pvalue = 0.05,
            log2FC_cutoff_value = 0,
            concept_ids = None, 
            background_gene_set = None,
            sig_conceptID_cutoff_pvalue = 0.05,
            order_by_p_value = True,
            min_num_genes_in_concept = 10,
            include_all = False,
            orthology_base = 'dre') -> pd.DataFrame:
    """
    Perform functional enrichment analysis using the logistic method, a cut-off free method.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame containing gene information, including gene IDs, p-values, and log2FC.
        - ``database (str or list[str])``: A list of functional annotation databases to test. Options include:

                - 'KEGG Pathway': Kyoto Encyclopedia of Genes and Genomes Pathways
                - 'KEGG Disease': Kyoto Encyclopedia of Genes and Genomes Diseases
                - 'GO BP': Gene Ontology Biological Processes
                - 'GO CC': Gene Ontology Cellular Components
                - 'GO MF': Gene Ontology Molecular Function
                - 'GO': Gene Ontology Biological Processes, Cellular Components, and Molecular Function
                - 'all': all databases shown above
                - a list of any combination of the databases shown above. eg. databases = [KEGG Pathway, KEGG Disease]

        - ``gene_id_type (str)``: The type of gene ID in the gene universe. The recommended gene id type is NCBI Gene ID (NCBI_ID). Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, Symbol, or for human: Human NCBI Gene ID.
        - ``org (str)``: The organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish, 'hsa' for human, 'variable' for varied).
        - ``directional_test (str, optional)``: 'up' to test for up-regulation, 'down' to test for down-regulation, 'non-directional' for enrichment/depletion. Default is 'non-directional'
        - ``sig_gene_cutoff_pvalue (float, optional)``: The significance cutoff for gene inclusion based on p-values. Default is 0.05.
        - ``log2FC_cutoff_value (float, optional)``: The log2 fold change cutoff value for gene inclusion. Default is 0.
        - ``concept_ids (list, optional)``: A list of concept IDs (e.g., pathway IDs or disease IDs) to analyze. Default is None.
        - ``background_gene_set (pd.DataFrame, optional)``: A DataFrame representing a background gene set. Default is None.
        - ``sig_conceptID_cutoff_pvalue (float, optional)``: The significance cutoff for concept IDs based on p-values. Default is 0.05.
        - ``order_by_p_value (bool, optional)``: Whether to order the results by p-value. Default is True.
        - ``min_num_genes_in_concept (int, optional)``: The minimum number of genes in a concept for it to be considered. Default is 10.
        - ``include_all (bool, optional)``: Include all results without filtering based on significance. Default is False.
        - `orthology_base (str, optional)`: The base organism code for orthology comparisons. Default is 'dre'.

    Returns:
        - ``result (pd.DataFrame)``: A DataFrame containing enrichment analysis results, including concept details, the number of genes in the concept in the universe, the number of significant genes belonging to the concept, the proportion of genes in the concept, p-value, odds ratio, and enrichment direction.

    Note:
        - If you are providing a list of concept ids, they must come from one database only.
        - Fisher's Exact test is cutoff dependent. 
    """    

    direction_options = ['non-directional', 'up', 'down']
    if direction not in direction_options:
        raise ValueError('Invalid Direction Choice.')

    if org != 'variable':
        out = _enrich(gene_universe, database, gene_id_type, org, 'fishers', direction, 
                    sig_gene_cutoff_pvalue, log2FC_cutoff_value, concept_ids, 
                    background_gene_set, sig_conceptID_cutoff_pvalue, order_by_p_value, 
                    min_num_genes_in_concept, include_all)
    else:
        out = _enrich_variable_org(gene_universe, database, gene_id_type, org, 
                                            'logistic', direction, 
                    sig_gene_cutoff_pvalue, log2FC_cutoff_value, concept_ids, 
                    background_gene_set, sig_conceptID_cutoff_pvalue, order_by_p_value, 
                    min_num_genes_in_concept, include_all, orthology_base)
    return out

def enrich_logistic(gene_universe: pd.DataFrame, 
                    database: list[str], 
                    gene_id_type: str,
                    org: str = 'dre',
                    directional_test: bool = True,
                    sig_gene_cutoff_pvalue: float = 0.05,
                    log2FC_cutoff_value: float = 0,
                    concept_ids: list[str] = None, 
                    background_gene_set: pd.DataFrame = None,
                    sig_conceptID_cutoff_pvalue: float = 0.05,
                    order_by_p_value: bool = True,
                    min_num_genes_in_concept: int = 10,
                    include_all: bool = False,
                    orthology_base: str = 'dre') -> pd.DataFrame:
    """
    Perform functional enrichment analysis using the logistic method, a cut-off free method.

    Parameters:
      - `gene_universe (pd.DataFrame)`: A DataFrame containing gene information, including gene IDs, p-values, and log2FC.
      - `database (str or list[str])`: A list of functional annotation databases to test.
      - `gene_id_type (str)`: The type of gene ID in the gene universe.
      - `org (str, optional)`: The organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish, 'hsa' for human, 'variable' for varied).
      - `directional_test (bool, optional)`: True for directional test (up/down regulation), False for non-directional (enrichment/depletion). Default is True.
      - `sig_gene_cutoff_pvalue (float, optional)`: The significance cutoff for gene inclusion based on p-values. Default is 0.05.
      - `log2FC_cutoff_value (float, optional)`: The log2 fold change cutoff value for gene inclusion. Default is 0.
      - `concept_ids (list[str], optional)`: A list of concept IDs (e.g., pathway IDs or disease IDs) to analyze. Default is None.
      - `background_gene_set (pd.DataFrame, optional)`: A DataFrame representing a background gene set. Default is None.
      - `sig_conceptID_cutoff_pvalue (float, optional)`: The significance cutoff for concept IDs based on p-values. Default is 0.05.
      - `order_by_p_value (bool, optional)`: Whether to order the results by p-value. Default is True.
      - `min_num_genes_in_concept (int, optional)`: The minimum number of genes in a concept for it to be considered. Default is 10.
      - `include_all (bool, optional)`: Include all results without filtering based on significance. Default is False.
      - `orthology_base (str, optional)`: The base organism code for orthology comparisons. Default is 'dre'.

    Returns:
      - `result (pd.DataFrame)`: A DataFrame containing enrichment analysis results.

    Note:
      - If you are providing a list of concept IDs, they must come from one database only.
      - The logistic regression method does not depend on the p-value cutoff for significance.
    """

    if not isinstance(directional_test, bool):
        raise ValueError("directional_test must be a boolean.")
    
    if directional_test:
        direction = 'directional'
    else:
        direction = 'non-directional'

    if org != 'variable':
        out = _enrich(gene_universe, database, gene_id_type, org, 'logistic', direction, 
                    sig_gene_cutoff_pvalue, log2FC_cutoff_value, concept_ids, 
                    background_gene_set, sig_conceptID_cutoff_pvalue, order_by_p_value, 
                    min_num_genes_in_concept, include_all)
    else:
        out = _enrich_variable_org(gene_universe, database, gene_id_type, org, 
                                            'logistic', direction, 
                    sig_gene_cutoff_pvalue, log2FC_cutoff_value, concept_ids, 
                    background_gene_set, sig_conceptID_cutoff_pvalue, order_by_p_value, 
                    min_num_genes_in_concept, include_all, orthology_base)
    return out

def _enrich_variable_org(gene_universe: pd.DataFrame, 
                                  database: str or list[str], 
                                  gene_id_type: str,
                                  org: str = 'dre',
                                  method: str = 'logistic',
                                  direction: str = 'non-directional',
                                  sig_gene_cutoff_pvalue: float = 0.05,
                                  log2FC_cutoff_value: float = 0,
                                  concept_ids: list[str] = None, 
                                  background_gene_set: pd.DataFrame = None,
                                  sig_conceptID_cutoff_pvalue: float = 0.05,
                                  order_by_p_value: bool = True,
                                  min_num_genes_in_concept: int = 10,
                                  include_all: bool = False,
                                  orthology_base: str = 'dre') -> pd.DataFrame:
    """
    Perform gene enrichment analysis for a variable organism.

    Parameters:
      - `gene_universe (pd.DataFrame)`: A DataFrame containing gene information, including gene IDs, p-values, and log2FC.
      - `database (Union[str, List[str]])`: A list of functional annotation databases to test.
      - `gene_id_type (str)`: The type of gene ID in the gene universe.
      - `org (str, optional)`: The organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish, 'hsa' for human).
      - `method (str, optional)`: The enrichment analysis method ('logistic' or 'fishers'). Default is 'logistic'.
      - `direction (str, optional)`: The direction of statistical test for enrichment (Fishers options: 'up', 'down', or 'non-directional', Logistic options: 'directional', 'non-directional'). Default is 'non-directional'.
      - `sig_gene_cutoff_pvalue (float, optional)`: The significance cutoff for gene inclusion based on p-values. Default is 0.05.
      - `log2FC_cutoff_value (float, optional)`: The log2 fold change cutoff value for gene inclusion. Default is 0.
      - `concept_ids (List[str], optional)`: A list of concept IDs (e.g., pathway IDs or disease IDs) to analyze. Default is None.
      - `background_gene_set (pd.DataFrame, optional)`: A DataFrame representing a background gene set. Default is None.
      - `sig_conceptID_cutoff_pvalue (float, optional)`: The significance cutoff for concept IDs based on p-values. Default is 0.05.
      - `order_by_p_value (bool, optional)`: Whether to order the results by p-value. Default is True.
      - `min_num_genes_in_concept (int, optional)`: The minimum number of genes in a concept for it to be considered. Default is 10.
      - `include_all (bool, optional)`: Include all results without filtering based on significance. Default is False.
      - `orthology_base (str, optional)`: The base organism code for orthology comparisons. Default is 'dre'.

    Returns:
      - `result (pd.DataFrame)`: A DataFrame containing enrichment analysis results.

    Note:
      - If providing a list of concept IDs, they must come from one database only.
    """
    
    # quality control database choice
    database = _process_database_options(database)
    if concept_ids is None:
        if orthology_base == 'dre':
            # get all concept ids for dre given db
            concept_id_dict = {}
            for db in database:
                dre_concept_ids, dreM_concept_ids = _organize_concept_ids(db, 'dre', 'dreM')
                concept_id_dict[('dre', db)] = dre_concept_ids
                concept_id_dict[('dreM', db)] = dreM_concept_ids
        elif orthology_base == 'dreM' or orthology_base == 'hsa':
            # get all concept ids for dre given db
            concept_id_dict = {}
            for db in database:
                dreM_concept_ids, dre_concept_ids = _organize_concept_ids(db, 'dreM', 'dre')
                concept_id_dict[('dreM', db)] = dreM_concept_ids
                concept_id_dict[('dre', db)] = dre_concept_ids 
    else:
        if len(database) > 1:
            raise ValueError('Only one database allowed for provided concept ids.')
        else:
            db = database[0]
        if orthology_base == 'dre':
            # get all concept ids for dre given db
            concept_id_dict = {}
            dre_concept_ids, dreM_concept_ids = _organize_concept_ids(db, 'dre', 'dreM', concept_ids)
            concept_id_dict[('dre', db)] = dre_concept_ids
            concept_id_dict[('dreM', db)] = dreM_concept_ids
        elif orthology_base == 'dreM' or orthology_base == 'hsa':
            dreM_concept_ids, dre_concept_ids = _organize_concept_ids(db, 'dreM', 'dre', concept_ids)
            concept_id_dict[('dreM', db)] = dreM_concept_ids
            concept_id_dict[('dre', db)] = dre_concept_ids            
    
    out_dataframe_list = []
    for key, value in concept_id_dict.items():   
        org_inner = key[0]
        db = key[1]
        concept_ids = value
        if len(concept_ids) == 0:
            continue
        out_inner = _enrich(gene_universe, db, gene_id_type, org_inner, method, direction, 
                sig_gene_cutoff_pvalue, log2FC_cutoff_value, concept_ids, 
                background_gene_set, sig_conceptID_cutoff_pvalue, order_by_p_value, 
                min_num_genes_in_concept, include_all)
        out_dataframe_list.append(out_inner)
    out = pd.concat(out_dataframe_list, ignore_index=True)
    if order_by_p_value:
        out = out.sort_values(by='P-value', ascending=True)
    out = out.reset_index(drop=True)
    return out

def _enrich(gene_universe: pd.DataFrame, 
            database: list[str], 
            gene_id_type: str,
            org = 'dre',
            method = 'logistic',
            direction = 'non-directional',
            sig_gene_cutoff_pvalue = 0.05,
            log2FC_cutoff_value = 0,
            concept_ids = None, 
            background_gene_set = None,
            sig_conceptID_cutoff_pvalue = 0.05,
            order_by_p_value = True,
            min_num_genes_in_concept = 10,
            include_all = False) -> pd.DataFrame:
    """
    Perform gene enrichment analysis.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame containing gene information, including gene IDs, p-values, and log2FC.
        - ``database (str or list[str])``: A list of functional annotation databases to test. Options include:

                - `KEGG Pathway`
                - `KEGG Disease`
                - `GO BP`: Gene Ontology Biological Processes
                - `GO CC`: Gene Ontology Cellular Components
                - `GO MF`: Gene Ontology Molecular Function
                - `GO`: Gene Ontology Biological Processes, Cellular Components, and Molecular Function
                - `all`: all databases shown above
                - a list of any combination of the databases shown above. eg. databases = [KEGG Pathway, KEGG Disease]

        - ``gene_id_type (str)``: The type of gene ID in the gene universe. The recommended gene id type is NCBI Gene ID (NCBI_ID). Must be one of: NCBI Gene ID, ZFIN ID, Ensembl ID, Symbol, or for human: Human NCBI Gene ID.
        - ``org (str)``: The organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish, 'hsa' for human).
        - ``method (str, optional)``: The enrichment analysis method ('logistic' or 'fishers'). Default is 'logistic'.
        - ``direction (str, optional)``: The direction of statistical test for enrichment (Fishers options: 'up', 'down', or 'non-directional', Logistic options: 'directional', 'non-directional'). .
        - ``sig_gene_cutoff_pvalue (float, optional)``: The significance cutoff for gene inclusion based on p-values. Default is 0.05.
        - ``log2FC_cutoff_value (float, optional)``: The log2 fold change cutoff value for gene inclusion. Default is 0.
        - ``concept_ids (list, optional)``: A list of concept IDs (e.g., pathway IDs or disease IDs) to analyze. Default is None.
        - ``background_gene_set (pd.DataFrame, optional)``: A DataFrame representing a background gene set. Default is None.
        - ``sig_conceptID_cutoff_pvalue (float, optional)``: The significance cutoff for concept IDs based on p-values. Default is 0.05.
        - ``order_by_p_value (bool, optional)``: Whether to order the results by p-value. Default is True.
        - ``min_num_genes_in_concept (int, optional)``: The minimum number of genes in a concept for it to be considered. Default is 10.
        - ``include_all (bool, optional)``: Include all results without filtering based on significance. Default is False.

    Returns:
        - ``result (pd.DataFrame)``: A DataFrame containing enrichment analysis results, including concept details, the number of genes in the concept in the universe, the number of significant genes belonging to the concept, the proportion of genes in the concept, p-value, odds ratio, and enrichment direction.

    Note:
        - If you are providing a list of concept ids, they must come from one database only.
    """      

    # QUALITY CONTROL :
    if type(gene_universe) != pd.DataFrame:
        gene_universe = pd.DataFrame(gene_universe)
    gene_id_type = utils.normalize_gene_id_type(gene_id_type)

    org = utils.normalize_organism_name(org)
    utils.check_valid_organism(org)

    if org == 'dre' or org == 'dreM':
        utils.check_valid_zebrafish_gene_id_type(gene_id_type)
        gene_universe = _check_gene_universe(gene_universe, gene_id_type, org)
    elif org == 'hsa':
        utils.check_valid_human_gene_id_type(gene_id_type)
        gene_universe = _check_gene_universe(gene_universe, gene_id_type, org)
    
    # quality control database choice
    database = _process_database_options(database)

    kegg_options = ['KEGG Pathway', 'KEGG Disease']
    go_options = ['GO', 'GO BP', 'GO CC', 'GO MF']

    fishers_direction_options = ['up', 'down', 'non-directional']
    logistic_direction_options = ['directional', 'non-directional']
    if method == 'fishers' and direction not in fishers_direction_options:
        raise ValueError('Invalid direction for Fisher\'s Method')
    elif method == 'logistic' and direction not in logistic_direction_options:
        raise ValueError('Invalid direction for Logistic Method')

    # -------------------------
    # DEAL WITH METHODS CHOICE
    # -------------------------
    enrich_method_function = _get_enrichment_method(method)
    # -------------------------
    # DEAL WITH DATABASES CHOICE
    # -------------------------
    if concept_ids is None:
        original_concept_ids = concept_ids
    else:
        original_concept_ids = concept_ids.copy()  # Create a copy of the original concept_ids
        if len(database) > 1: 
            raise ValueError('Provided Concept IDs must come from one database only.')
    resulting_dataframe_list = []  # List to store results for each database
    
    for db in database:
        # identify concept function to use
        concept_dict = {
            'KEGG Pathway': KEGG.get_genes_in_pathway,
            'KEGG Disease': KEGG.get_genes_in_disease,
            'GO': GO.get_genes_in_GO_concept,
            'GO BP': GO.get_genes_in_GO_concept,
            'GO CC': GO.get_genes_in_GO_concept,
            'GO MF': GO.get_genes_in_GO_concept,
        }

        get_genes_function = concept_dict[db]
        concept_type = db
        all_ids, id_column_name, name_column_name = _get_pathway_ids_and_names(db, org)
        gene_id_type, gene_universe = _map_to_master_geneid(gene_universe, gene_id_type, org, 
                                                            db, kegg_options, go_options)

        # DEAL WITH CONCEPT_IDS CHOICE
        if original_concept_ids is None:
            current_concept_ids = all_ids[id_column_name]
        else:
            if db in kegg_options:
                current_concept_ids = _check_concept_ids_KEGG(original_concept_ids, org, db, all_ids, id_column_name)
            elif db in go_options:
                current_concept_ids = _check_concept_ids_GO(original_concept_ids, org, db, all_ids, id_column_name)

        # -------------------------
        # SIGNIFICANT GENE UNIVERSE
        # -------------------------  

        if background_gene_set:
            gene_universe = _preprocess_gene_universe(gene_universe, background_gene_set, gene_id_type)
            # another option would be to use all genes in the genome. 
        total_number_of_genes_in_universe = len(gene_universe) 
        
        
        # to test for overrepresentation, you would want only the over expressed genes
        # to test for either enrichment or depletion you would want both
        # to test downregulated genes you would do down 
        if method == 'fishers':
            if direction == 'up':
                test_direction = 'greater'
            elif direction == 'down':
                test_direction = 'less'
            elif direction == 'non-directional':
                test_direction = 'two-sided'
        if method == 'logistic':
            if direction == 'directional' or direction == 'non-directional':
                test_direction = direction
        
        if method == 'fishers':
            sig_genes_set = _get_sig_genes_set(gene_universe, direction, gene_id_type, sig_gene_cutoff_pvalue, log2FC_cutoff_value)
        else:
            sig_genes_set = None
        # -------------------------
        # LAUNCH ENRICHMENT
        # -------------------------
        resulting_dictionary_list = []    
        for concept_id in current_concept_ids:
            gene_set =  get_genes_function(concept_id, org, do_check = False)
            gene_set = gene_set[gene_set[gene_id_type].isin(gene_universe[gene_id_type])]
            num_genes = len(gene_set)
            if num_genes >= min_num_genes_in_concept:
                concept_name = all_ids.loc[all_ids[id_column_name] == concept_id, name_column_name].values[0]
                out = enrich_method_function(gene_universe, 
                                            sig_genes_set,
                                            gene_set, 
                                            gene_id_type,
                                            concept_type, 
                                            concept_id, 
                                            concept_name,
                                            total_number_of_genes_in_universe,
                                            test_direction,
                                            sig_gene_cutoff_pvalue,
                                            log2FC_cutoff_value)
                if out:
                    # Append the dictionary to the list
                    resulting_dictionary_list.append(out)
    
        result = pd.DataFrame(resulting_dictionary_list)
        if sig_conceptID_cutoff_pvalue and not include_all:
            result = result[result["P-value"] <= sig_conceptID_cutoff_pvalue]
        if order_by_p_value:
            result = result.sort_values(by='P-value', ascending=True)
        if method == 'fishers':
            if direction == 'up':
                result = result[result['Direction'] == 'upregulated']
            elif direction == 'down':
                result = result[result['Direction'] == 'downregulated']
        # Append the result DataFrame to the list
        resulting_dataframe_list.append(result)

    # Concatenate results for all databases into a single DataFrame
    final_result = pd.concat(resulting_dataframe_list, ignore_index=True)
    if order_by_p_value:
        final_result = final_result.sort_values(by='P-value', ascending=True)
    
    return final_result.reset_index(drop=True)

def _logistic(gene_universe_in: pd.DataFrame,
             sig_genes_set: None,
             gene_set: pd.DataFrame, 
             gene_id_type: str, 
             concept_type: str, 
             concept_id: str, 
             concept_name: str,
             total_number_of_genes_in_universe: int,
             test_direction,
             pval_cutoff, 
             log2FC_cutoff):
    """
    Perform functional enrichment analysis using logistic regression.

    Parameters:
        - ``gene_universe_in (pd.DataFrame)``: A DataFrame representing the universe of genes.
        - ``sig_genes_set (None)``: Placeholder.
        - ``gene_set (pd.DataFrame)``: A DataFrame containing the genes of interest.
        - ``gene_id_type (str)``: The type of gene identifier used in the DataFrames.
        - ``concept_type (str)``: The type of concept (e.g., pathway) being analyzed.
        - ``concept_id (str)``: The ID of the concept being analyzed.
        - ``concept_name (str)``: The name or description of the concept being analyzed.
        - ``total_number_of_genes_in_universe (int)``: The total number of genes in the universe.
        - ``test_direction (str)``: The directionality of the test ('greater', 'less', or 'two-sided').
        
    Returns:
        - ``data (dict)``: A dictionary containing enrichment analysis results, including concept details,
          the number of genes in the concept in the universe, the number of significant genes belonging
          to the concept, the proportion of genes in the concept, odds ratio, p-value, and enrichment direction.

    Notes:
        - This function performs functional enrichment analysis using logistic regression.
        - It calculates enrichment statistics for a specified concept (e.g., pathway) by comparing a gene set
          of interest to a larger gene universe.
        - The 'gene_id_type' parameter specifies the type of gene identifiers used in the DataFrames.
        - The 'test_direction' parameter determines the directionality of the test ('greater', 'less', or 'two-sided').
        - Enrichment results include the odds ratio, p-value, and enrichment direction.
    """    
    gene_universe = gene_universe_in.copy()
    
    # determine which genes are in the gene set
    gene_universe['InGeneSet'] = gene_universe[gene_id_type].isin(gene_set[gene_id_type]).astype(int)

    # test direction
    if test_direction == 'directional':
        gene_universe['direction'] = np.where(gene_universe['log2FC'] > 0, 1, -1)
        gene_universe['NegLogPValue'] = -np.log10(gene_universe['PValue'])
        gene_universe['NegLogPValue'] *= gene_universe['direction']
    elif test_direction == 'non-directional':
        gene_universe['NegLogPValue'] = -np.log10(gene_universe['PValue'])
    else:
        raise ValueError('Invalid Test Direction.')
    
    try:
        # perform logistic regression
        # ---------------------------   
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ConvergenceWarning)
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            warnings.filterwarnings ("ignore", category = PerfectSeparationWarning)
            formula = "InGeneSet ~ NegLogPValue"

            log_reg = smf.logit(formula, data=gene_universe).fit(disp=0)
            beta = log_reg.params['NegLogPValue']
            p_value = log_reg.pvalues['NegLogPValue']
            odds_ratio = np.exp(beta)

        # Determine enrichment direction based on the test direction
        if test_direction == 'directional':
            direction = 'upregulated' if beta > 0 else 'downregulated'
        else:
            # For two-sided or invalid test directions, use the sign of the coefficient
            direction = 'enriched' if beta > 0 else 'depleted' if beta < 0 else 'neutral'

        if direction == 'upregulated':
            sig_genes_set = _get_sig_genes_set(gene_universe, 'up', gene_id_type, pval_cutoff, log2FC_cutoff)
        elif direction == 'downregulated':
            sig_genes_set = _get_sig_genes_set(gene_universe, 'down', gene_id_type, pval_cutoff, log2FC_cutoff)
        else:
            sig_genes_set = _get_sig_genes_set(gene_universe,'non-directional', gene_id_type, pval_cutoff, log2FC_cutoff)

        gene_set_set = set(gene_set[gene_id_type])
        # Number of genes that are both in the gene set and significantly expressed
        a = len(sig_genes_set.intersection(gene_set_set))

        # Number of genes in the gene set but not significantly expressed
        b = len(gene_set_set.difference(sig_genes_set))

        # Calculate the proportion of significant genes in the set (avoid division by zero)
        proportion_sig_genes_in_set = a / (a + b) if (a + b) > 0 else 0

        # organize important stats
        # ------------------------
        if odds_ratio == np.inf:
            data = None
        else:
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
    
    except np.linalg.LinAlgError as e:
        data = None
        return data

def _fishers(gene_universe: pd.DataFrame,
            sig_genes_set: pd.DataFrame,
           gene_set: pd.DataFrame,
           gene_id_type: str, 
           concept_type: str, 
           concept_id: str, 
           concept_name: str,
           total_number_of_genes_in_universe: int,
           test_direction,
           pval_cutoff = None,
           log2fc_cutoff = None):
    """
    Perform functional enrichment analysis using Fisher's exact test.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame representing the universe of genes.
        - ``sig_genes_set (set)``: A set containing the significantly expressed genes.
        - ``gene_set (pd.DataFrame)``: A DataFrame containing the genes of interest.
        - ``gene_id_type (str)``: The type of gene identifier used in the DataFrames.
        - ``concept_type (str)``: The type of concept (e.g., pathway) being analyzed.
        - ``concept_id (str)``: The ID of the concept being analyzed.
        - ``concept_name (str)``: The name or description of the concept being analyzed.
        - ``total_number_of_genes_in_universe (int)``: The total number of genes in the universe.
        - ``test_direction (str)``: The directionality of the test. 
          Options: 'two-sided', 'greater', 'less'.
        
    Returns:
        - ``data (dict)``: A dictionary containing enrichment analysis results, including concept details,
          the number of genes in the concept in the universe, the number of significant genes belonging
          to the concept, the proportion of genes in the concept, odds ratio, p-value, and enrichment direction.

    Notes:
        - This function performs functional enrichment analysis using Fisher's exact test.
        - It calculates enrichment statistics for a specified concept (e.g., pathway) by comparing a gene set
          of interest to a larger gene universe.
        - The 'gene_id_type' parameter specifies the type of gene identifiers used in the DataFrames.
        - The 'test_direction' parameter determines the directionality of the test ('two-sided', 'greater', 'less').
        - Enrichment results include the odds ratio, p-value, and enrichment direction.
    """

    # Filter genes in the gene set that are in our gene universe
    # gene_set = gene_set[gene_set[gene_id_type].isin(gene_universe[gene_id_type])]
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

    alt = 'two-sided'
    if test_direction != 'two-sided':
        alt = 'greater'

    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative=alt)

    # Determine enrichment direction
    if test_direction == 'two-sided':
        direction = 'enriched' if odds_ratio > 1 else 'depleted'
    elif test_direction == 'greater':
        direction = 'upregulated' if odds_ratio > 1 else 'neutral'
    elif test_direction == 'less':
        direction = 'downregulated' if odds_ratio > 1 else 'neutral'


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

def _organize_concept_ids(db: str, 
                          base_org: str, 
                          additional_org: str, 
                          given_ids: pd.DataFrame = None):
    """
    Organize concept IDs for gene enrichment analysis.

    Parameters:
      - `db (str)`: The database type ('KEGG Pathway', 'KEGG Disease', 'GO BP', 'GO CC', 'GO MF').
      - `base_org (str)`: The base organism code ('dre' for zebrafish, 'dreM' for mapped zebrafish).
      - `additional_org (str)`: The additional organism code for comparison ('dre' for zebrafish, 'dreM' for mapped zebrafish).
      - `given_ids (Optional[pd.DataFrame])`: A DataFrame containing concept IDs or None. If provided, the function checks the existence of IDs for the specified database and organisms.

    Returns:
      - `base_list (List[str])`: A list of concept IDs for the base organism.
      - `additional_list (List[str])`: A list of concept IDs for the additional organism.

    Note:
      - If `given_ids` is None, the function retrieves concept IDs from the specified databases and organisms, removes common numeric portions, and formats the IDs accordingly.
      - If `given_ids` is provided, it checks the existence of each ID for the specified database and organisms and categorizes them into base and additional lists.
    """

    if given_ids is None:
        all_ids_base, id_col_name, id_desc_col_name = _get_pathway_ids_and_names(db, base_org)
        concept_ids_base = all_ids_base[id_col_name]
        all_ids_additional, id_col_name, id_desc_col_name = _get_pathway_ids_and_names(db, additional_org)
        concept_ids_additional = all_ids_additional[id_col_name]

        # Extract numeric portions and convert to integers for both series
        numeric_additional = concept_ids_additional.str.extract('(\d+)')
        numeric_base = concept_ids_base.str.extract('(\d+)')

        # Extract numeric portions as sets
        numeric_additional_set = set(numeric_additional.squeeze())
        numeric_base_set = set(numeric_base.squeeze())

        # Identify common numeric portions using sets
        common_numeric_set = numeric_base_set.intersection(numeric_base_set)

        # Remove common elements from numeric_additional_set
        filtered_numeric_additional_set = list(numeric_additional_set - common_numeric_set)
        if db == 'KEGG Pathway':
            prefixed_list = [additional_org + entry for entry in filtered_numeric_additional_set]
            base_list = concept_ids_base.squeeze().to_list()
            additional_list = prefixed_list
        elif db == 'KEGG Disease':
            prefixed_list = ['H' + entry for entry in filtered_numeric_additional_set]
            base_list = concept_ids_base.squeeze().to_list()
            additional_list = prefixed_list            
        else:
            prefixed_list = ['GO:' + entry for entry in filtered_numeric_additional_set]
            base_list = concept_ids_base.squeeze().to_list()
            additional_list = prefixed_list

    else:
        if type(given_ids) == pd.DataFrame:
            if given_ids.shape[1] != 1:
                raise ValueError('Concept IDs given should be a 1 dimensional list.')
            else:
                column_name = given_ids.columns[0]
                given_ids = given_ids[column_name].to_list()   
        func_dict = {
            'KEGG Pathway': KEGG._check_if_pathway_id_exists,
            'KEGG Disease': KEGG._check_if_disease_id_exists,
            'GO BP': GO.id_exists_given_organism,
            'GO CC': GO.id_exists_given_organism,
            'GO MF': GO.id_exists_given_organism
        }
        if db == 'KEGG Disease':
            base_list = []
            additional_list = []
            bad_ids = []
            fun = func_dict[db]
            for id in given_ids:
                if fun(id, 'dre'):
                    base_list.append(id)
                else:
                    bad_ids.append(id)
            if len(bad_ids) >0:
                print(f'WARNING: {len(bad_ids)} given KEGG Disease ids are invalid or do not contain')
                print(f'         any genes reported on the KEGG Database website.')
                print()
                print('Omitted ids include:')
                for idx, id in enumerate(bad_ids, start = 1):
                    print(id, end=" ")
                    if idx % 4 == 0:
                        print()  # Print a newline after every fourth ID
                print()
        else:
            base_list = []
            additional_list = []
            fun = func_dict[db]
            for id in given_ids:
                if fun(id, 'dre'):
                    base_list.append(id)
                elif fun(id, 'dreM'):
                    additional_list.append(id)
                else:
                    continue

    return base_list, additional_list

def _map_to_master_geneid(gene_universe, gene_id_type, org, database, kegg_options, go_options):
    """
    Map gene identifiers in the gene universe DataFrame to a common master gene ID.

    Parameters:
        - ``gene_universe (pd.DataFrame)``: A DataFrame containing gene information, including gene IDs.
        - ``gene_id_type (str)``: The type of gene ID in the gene universe.
        - ``org (str)``: The organism code ('dre' for zebrafish, 'hsa' for human).
        - ``database (str)``: The functional annotation database.
        - ``kegg_options (List[str])``: Options for KEGG databases.
        - ``go_options (List[str])``: Options for GO databases.

    Returns:
        - ``gene_id_type (str)``: The updated gene ID type after mapping to the master gene ID.
        - ``gene_universe (pd.DataFrame)``: The DataFrame with mapped gene identifiers.
    """
    if org == 'dre' or org == 'dreM':
        if database in kegg_options:
            master_gene_id = NCBI_ID
        elif database in go_options:
            master_gene_id = ZFIN_ID
    elif org == 'hsa':
        master_gene_id = HUMAN_ID

    if gene_id_type != master_gene_id:
        gene_universe = mapping.add_mapped_column(gene_universe, gene_id_type, master_gene_id, gene_id_type, keep_old_ids=False, drop_na=True)
        gene_id_type = master_gene_id
    
    return gene_id_type, gene_universe

def _get_enrichment_method(method: str):
    """
    Get the enrichment method function based on the provided method name.

    Parameters:
       - `method`: The name of the enrichment method.

    Returns:
        - Enrichment method function.

    Raises:
        - ValueError: If an invalid method is provided.
    """
    methods_dict = {
        'logistic': _logistic,
        'fishers': _fishers
    }

    # Normalize method name to lowercase for case-insensitive matching
    normalized_method = method.lower()

    # Check if the method is valid
    if normalized_method not in methods_dict:
        raise ValueError(f"Invalid method: {method}. Supported methods are {', '.join(methods_dict.keys())}.")

    # Get the function based on the method name
    return methods_dict[normalized_method]

def _process_database_options(databases):
    """
    Process the input databases to ensure they are in the correct format.

    Parameters:
       - `databases`: The input databases, which can be a string or a list.

    Returns:
        - List of expanded databases.
    """
    if not (isinstance(databases, str) or isinstance(databases, list)):
        raise TypeError("Input databases must be either a string or a list.")
    
    if isinstance(databases, str):
        databases = _normalize_and_validate_database(databases)
        databases = _expand_database_option(databases)
    elif isinstance(databases, list):
        expanded_databases = []
        for db_option in databases:
            normalized_db_option = _normalize_and_validate_database(db_option)
            expanded_databases.extend(_expand_database_option(normalized_db_option))
        databases = expanded_databases
    return databases

def _normalize_and_validate_database(database: str) -> str:
    """
    Normalize and validate a database option to a specified format.

    Parameters:
       - database (str): The database option to be normalized and validated.

    Returns:
        - normalized_database (str): The normalized and validated database option.
    """

    database_mappings = {
        'kegg pathway': 'KEGG Pathway',
        'kegg pathways': 'KEGG Pathway',
        'kegg': 'KEGG Pathway',
        'kegg disease': 'KEGG Disease',
        'kegg diseases': 'KEGG Disease',
        'go biological processes': 'GO BP',
        'go bp': 'GO BP',
        'biological processes': 'GO BP',
        'go cellular component': 'GO CC',
        'go cc': 'GO CC',
        'cellular component': 'GO CC',
        'go molecular function': 'GO MF',
        'go mf': 'GO MF',
        'molecular function': 'GO MF',
        'go': 'GO',
        'all': 'all',
    }

    # Strip and convert the input database option to lowercase for case-insensitive matching
    lowercase_database = database.strip().lower()

    # Check if the lowercase database option exists in the mappings
    if lowercase_database in database_mappings:
        return database_mappings[lowercase_database]
    else:
        # If no mapping is found, raise a ValueError
        raise ValueError(f"Invalid database option: {database}. Please choose from: {', '.join(database_mappings.values())}")
    
def _expand_database_option(database_option: str) -> list:
    """
    Expand a database option to a list of corresponding databases.

    Parameters:
        - `database_option` (str): The selected database option.

    Returns:
        - `databases` (list): A list of corresponding databases.
    
    Note: this function is only ran after the database options have been checked. 
    """
    database_mappings = {
        'kegg pathway': ['KEGG Pathway'],
        'kegg disease': ['KEGG Disease'],
        'go bp': ['GO BP'],
        'go cc': ['GO CC'],
        'go mf': ['GO MF'],
        'go': ['GO BP', 'GO CC', 'GO MF'],
        'all': ['KEGG Pathway', 'KEGG Disease', 'GO BP', 'GO CC', 'GO MF'],
    }

    # Convert the input database option to lowercase for case-insensitive matching
    lowercase_database_option = database_option.lower()

    # Check if the lowercase database option exists in the mappings
    if lowercase_database_option in database_mappings:
        return database_mappings[lowercase_database_option]
    else:
        # If no mapping is found, return a list with the original input
        return [database_option]
    
def _preprocess_gene_universe(gene_universe, background_gene_set, gene_id_type, pvalue_col='PValue', log2FC_col='log2FC'):
    """
    Preprocess the gene universe by adding missing genes from a background set with 'NaN' values for p-value and log2FC.

    Parameters:
    - gene_universe: DataFrame containing information about the subset of genes.
    - background_gene_set: DataFrame containing information about the background genes.
    - gene_id_type: Name of the column containing gene IDs in gene_universe.
    - pvalue_col: Name of the column containing p-values (default: 'PValue').
    - log2FC_col: Name of the column containing log2 fold changes (default: 'log2FC').

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
                                                         & (gene_universe['log2FC'] > log2FC_cutoff)])
    elif method == 'down':
        sig_genes_set = set(gene_universe[gene_id_type][(gene_universe['PValue'] < pval_cutoff)
                                                         & (gene_universe['log2FC'] < -log2FC_cutoff)])
    elif method == 'non-directional':
        sig_genes_set = set(gene_universe[gene_id_type][(gene_universe['PValue'] < pval_cutoff)
                                                        & (np.abs(gene_universe['log2FC']) > log2FC_cutoff)])
    else:
        raise ValueError("Invalid method. Supported methods are 'up', 'down', or 'non-directional'.")
    
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
        org_col_name = 'exists_dre' if org == 'dre' else 'exists_hsa'
        all_ids = all_ids[(all_ids[org_col_name] == True)]
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
    
    if specified_gene_id_type in ['NCBI Gene ID', 'Human NCBI Gene ID']:
        if type(gene_universe[specified_gene_id_type].values[0]) == str:
            gene_universe[specified_gene_id_type] = gene_universe[specified_gene_id_type].astype('int64')

    return gene_universe
    
def _check_concept_ids_KEGG(concept_ids, org, database, all_ids, id_column_name) -> list:
    if type(concept_ids) == pd.DataFrame:
        if concept_ids.shape[1] != 1:
            raise ValueError('Concept IDs given should be a 1 dimensional list.')
        else:
            column_name = concept_ids.columns[0]
            concept_ids = concept_ids[column_name].to_list()
    if database == 'KEGG Pathway':
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
    elif database == 'KEGG Disease':
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
    if any(concept_id not in all_ids[id_column_name].to_list() for concept_id in concept_ids):
        raise ValueError('Invalid ID')

    return concept_ids

def _check_concept_ids_GO(concept_ids, org, db, all_ids, id_column_name):
    if type(concept_ids) == pd.DataFrame:
        if concept_ids.shape[1] != 1:
            raise ValueError('Concept IDs given should be a 1 dimensional list.')
        else:
            column_name = concept_ids.columns[0]
            concept_ids = concept_ids[column_name].to_list()
    df = pd.read_csv(GO.GO_IDS_PATH, sep = '\t')
    modified_list = []
    bad_ids = []
    bad_ids_diff_db = []
    warning = False
    warning_ont = False
    for id in concept_ids:
        new_id = _check_GO_id_format(id)
        if not _id_exists_given_organism(id, org, df):
            bad_ids.append(id)
            warning = True
            # raise ValueError(f'GO ID {id} does not exist for given organism.')
        elif not _go_id_exists_given_db(id, db, df):
            bad_ids_diff_db.append(id)
            warning_ont = True
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
        print()
    if warning_ont:
        print(f'WARNING: {len(bad_ids_diff_db)} given GO IDs do not exist for given database.')
        print()
        print('Omitted ids include:')
        for idx, id in enumerate(bad_ids_diff_db, start = 1):
            print(id, end=" ")
            if idx % 4 == 0:
                print()  # Print a newline after every fourth ID
        print()
    concept_ids = modified_list
    if any(concept_id not in all_ids[id_column_name].to_list() for concept_id in concept_ids):
        invalid =  [concept_id for concept_id in concept_ids if concept_id not in all_ids[id_column_name].to_list()]
        raise ValueError(f'Invalid ID(s): {invalid}')
    return concept_ids

def _combine_results(dre_df, dreM_df, truth_base = 'dre'):

    dreM_df = dreM_df.copy()
    dre_df = dre_df.copy()

    # match concept ids
    dreM_df['Concept ID'] = dreM_df['Concept ID'].apply(lambda x: x.replace('dreM', 'dre'))

    if truth_base == 'dre':
        concept_ids_in_dre = set(dre_df['Concept ID'])

        # Filter out rows from results_dreM where Concept ID is in results_dre
        filtered_results_dreM = dreM_df[~dreM_df['Concept ID'].isin(concept_ids_in_dre)]

        final_merged_df = pd.concat([dre_df, filtered_results_dreM], ignore_index=True)
        
    elif truth_base == 'dreM':

        concept_ids_in_dreM = set(dreM_df['Concept ID'])
        filtered_results_dre = dre_df[~dre_df['Concept ID'].isin(concept_ids_in_dreM)]
        final_merged_df = pd.concat([dreM_df, filtered_results_dre], ignore_index=True)

    return final_merged_df

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

def _go_id_exists_given_db(concept_id, db, df):
    # check if ID exists for given db:
    if db == 'GO BP':
        target_column = 'P'
    elif db == 'GO MF':
        target_column = 'F'
    elif db == 'GO CC':
        target_column = 'C'
    else:
        raise ValueError('db')
    filtered_df = df[(df['GO ID'] == concept_id) & (df['Ontology'] == target_column)]
    return not filtered_df.empty

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
    