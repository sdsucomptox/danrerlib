# Set up python environment
from settings import *

import mapping
import KEGG
import GO

import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import scipy.stats
from scipy.stats import chi2
from scipy.stats import fisher_exact

def logistic(gene_universe: pd.DataFrame, gene_set: pd.DataFrame, gene_id_type: str,
             concept_type: str, concept_id: str, concept_name: str,
             sig_cutoff = 0.05) -> pd.DataFrame:
    '''
    This function performs gene enrichment using the logistic regression method. A 
    p-value cutoff is required for this method and the default is taken to be 0.05
    '''
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

def fishers(gene_universe, gene_set, concept_type, concept_id, concept_name, sig_cutoff = 0.05):
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

def enrich_KEGG(gene_universe: str, gene_id_type = NCBI_ID,
                database =  'pathway', concept_ids = None, 
                org = 'dre', method = 'logistic',
                sig_gene_cutoff_pvalue = 0.05,
                sig_conceptID_cutoff_pvalue = None,
                sig_conceptID_cutoff_FDR = None,
                order_by_p_value = True):
    '''
    gene_universe:  the given dataset containing genes, pvalues, and log2FC (in that order)
    '''
    
    # TODO
    # - make it so that the first column is assumed to be the Gene ID
    #  and the second column is the p-value
    # - support different zebrafish Gene ID inputs
    
    # only include genes in the gene set that are in 
    # our gene universe (full gene de)
    if type(gene_universe) != pd.DataFrame:
        gene_universe = pd.DataFrame(gene_universe)
    
    # quality control TODO:
    _check_gene_universe(gene_universe, gene_id_type)
    _check_valid_org(org)

    if gene_id_type != NCBI_ID:
        # I would map here
        gene_id_to = HUMAN_ID if org == 'hsa' else NCBI_ID
        gene_universe = mapping.add_mapped_column(gene_universe, gene_id_type, gene_id_to,
                                                  keep_old_ids=False, drop_na=True)
        gene_id_type = gene_id_to

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

def enrich_GO(gene_universe: str, gene_id_type = ZFIN_ID, 
              database = None, concept_ids = None, 
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
    _check_gene_universe(gene_universe, gene_id_type)
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

def _check_gene_universe(gene_universe, gene_id_type):
    if (gene_id_type not in gene_universe.columns 
        and 'PValue' not in gene_universe.columns):
        raise ValueError('Required columns do not exist.')  
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
def _check_valid_org(org):
    valid_orgs = ['dre', 'hsa', 'dreM']
    if org not in valid_orgs:
        raise ValueError('Invalid organism chosen. Options are: [\'dre\', \'hsa\', \'dreM\']')

def test_KEGG_enrich(option1  = False, option2 = False, option3 = False):
    # option 1: all pathways
    test_data_dir = FILE_DIR / Path('../../tutorials/data/test_data/')
    file_path = test_data_dir / Path('TPP.txt')
    tpp_df = pd.read_csv(file_path, sep='\t')
    if option1:
        id_type = 'NCBI Gene ID'
        out = enrich_KEGG(tpp_df, gene_id_type = id_type, org = 'dreM')
        print(out.head(3))
    if option2:
        file_path = test_data_dir / Path('kegg_pathways.txt')
        pathway_ids = pd.read_csv(file_path, sep='\t')
        zfish_gene_type = 'NCBI Gene ID' # the gene type in the tpp_df
        out = enrich_KEGG(tpp_df, 
                             gene_id_type = zfish_gene_type,
                             org = 'dreM',
                             database = 'pathway',
                             concept_ids = pathway_ids)
        print(out.head(3))
    # kegg disease
    if option3:
        file_path = test_data_dir / Path('kegg_disease.txt')
        pathway_ids = pd.read_csv(file_path, sep='\t')
        zfish_gene_type = 'NCBI Gene ID' # the gene type in the tpp_df
        out = enrich_KEGG(tpp_df, 
                        gene_id_type = zfish_gene_type,
                        org = 'dreM',
                        database = 'disease',
                        concept_ids = pathway_ids)
        print(out)

def investigation():
    # kegg_id = '04911'
    # org = 'dreM'
    # genes = KEGG.get_genes_in_pathway(kegg_id, org)
    # print(genes)
    # print(type(genes[NCBI_ID].values[0]))
    concept_id = 'H00019'
    genes = KEGG.get_genes_in_disease(concept_id, 'dreM')
    # print(genes)
    # print(type(genes[NCBI_ID].values[0]))


    test_data_dir = FILE_DIR / Path('../../tutorials/data/test_data/')
    file_path = test_data_dir / Path('TPP.txt')
    tpp_df = pd.read_csv(file_path, sep='\t')
    gene_universe = tpp_df


    gene = genes[NCBI_ID].values[0]

    is_in_universe = gene in gene_universe[NCBI_ID].values
    print(is_in_universe)
    
def test_GO_enrich(option1, option2, option3, option4, option5):
    test_data_dir = FILE_DIR / Path('../../tutorials/data/test_data/')
    file_path = test_data_dir / Path('TPP.txt')
    tpp_df = pd.read_csv(file_path, sep='\t')
    zfish_gene_type = 'NCBI Gene ID' # the gene type in the tpp_df
    # note: these take forever
    # option 1 all BP
    if option1:
        out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
                                org = 'dre', database='BP')
    # option 2 all MF
    if option2:
        out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
                                org = 'dre', database='BP')
    # option 3 all CC
    if option3:
        out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
                                org = 'dre', database='BP')
    # option 4 all GO!
    if option4:
        out = enrich_GO(tpp_df, gene_id_type = zfish_gene_type,
                                org = 'dre', database='BP')
    # option 5 given list
    if option4:
        file_path = test_data_dir / Path('go_ids.txt')
        pathway_ids = pd.read_csv(file_path, sep='\t')
        out = enrich_GO(tpp_df, 
                                gene_id_type = zfish_gene_type,
                                org = 'dre',
                                concept_ids = pathway_ids)
        print(out)


def testing():

    # cwd = Path().absolute() 
    # test_data_path = cwd / Path('tutorials/data/test_data/example_diff_express_data.txt')
    # gene_univese_full_data = pd.read_csv(test_data_path, sep='\t')

    # concept_ids = ['dreM00010', 'dreM00020']
    # concept_ids = ['10', '20']
    # org = 'dreM'
    # df = _check_concept_ids(concept_ids, org)
    # print(df)
    # test_KEGG_enrich(option3  = True)
    # investigation()

    return None

if __name__ == '__main__':
    testing()