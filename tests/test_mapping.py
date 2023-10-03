from danrerlib import mapping
from danrerlib.settings import *
from pandas.testing import assert_frame_equal
from itertools import permutations

# MAPPING COMBINATIONS ARE:
# ZFIN_ID and NCBI_ID  done
# ZFIN_ID and SYMBOL 
# ZFIN_ID and ENS_ID 
# NCBI_ID and SYMBOL 
# NCBI_ID and ENS_ID 
# SYMBOL and ENS_ID done

# ------------------------
# test convert_ids
# ------------------------

def test_mapping_list_of_ids_ncbi_and_zfin():

    ncbi_ids = ['562552', '100150233', '405818', '564006', '796163']
    
    zfin_ids = ['ZDB-GENE-001222-1', 'ZDB-GENE-030131-4309',
        'ZDB-GENE-040426-2356', 'ZDB-GENE-081001-1', 'ZDB-GENE-131127-614']
    
    out_ids = mapping.convert_ids(ncbi_ids, NCBI_ID, ZFIN_ID, out_format = list)
    assert sorted(zfin_ids) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(zfin_ids, ZFIN_ID, NCBI_ID, out_format = list)
    assert sorted(ncbi_ids) ==  sorted(out_ids)

def test_mapping_list_of_ids_zfin_and_symbol():

    zfin_ids = ['ZDB-GENE-060503-764', 'ZDB-GENE-010817-2', 'ZDB-GENE-070705-89',
                'ZDB-GENE-030131-296', 'ZDB-GENE-040608-1']
    
    symbols = ['si:ch211-152p11.4', 'angpt2b', 'trhde.2', 'ttpal', 'harbi1']
    
    out_ids = mapping.convert_ids(zfin_ids, ZFIN_ID, SYMBOL, out_format = list)
    assert sorted(symbols) ==  sorted(out_ids)
    
    out_ids = mapping.convert_ids(symbols, SYMBOL, ZFIN_ID, out_format = list)
    assert sorted(zfin_ids) ==  sorted(out_ids)


def test_mapping_list_of_ids_zfin_and_ens():

    zfin_ids = ['ZDB-GENE-040120-4', 'ZDB-LINCRNAG-110914-87', 'ZDB-GENE-141216-345',
                'ZDB-GENE-050809-56', 'ZDB-GENE-091204-383']
    ens_ids = ['ENSDARG00000008785', 'ENSDARG00000092478', 'ENSDARG00000100251', 
               'ENSDARG00000102792', 'ENSDARG00000096078']
    
    out_ids = mapping.convert_ids(zfin_ids, ZFIN_ID, ENS_ID, out_format = list)
    assert ens_ids ==  out_ids

    out_ids = mapping.convert_ids(ens_ids, ENS_ID, ZFIN_ID, out_format = list)
    assert sorted(zfin_ids) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(ens_ids, 'ens', 'zfin', out_format = list)
    assert sorted(zfin_ids) ==  sorted(out_ids)


def test_mapping_list_of_ids_ncbi_and_symbol():

    ncbi_ids = ['402853', '65231', '565228', '436656', '324103']
    symbols = ['uts2b', 'mdkb', 'si:ch211-215a10.4', 'krt97', 'wu:fc18g07']
    
    out_ids = mapping.convert_ids(ncbi_ids, NCBI_ID, SYMBOL, out_format = list)
    assert sorted(symbols) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(ncbi_ids, 'ncbi id', 'sym', out_format = list)
    assert sorted(symbols) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(symbols, SYMBOL, NCBI_ID, out_format = list)
    assert sorted(ncbi_ids) ==  sorted(out_ids)

def test_mapping_list_of_ids_ncbi_and_ens():

    ncbi_ids = ['550419', '794572', '101885923', '678616', '100535123']
    ens_ids = ['ENSDARG00000039730', 'ENSDARG00000091990', 'ENSDARG00000087413', 
                'ENSDARG00000058372', 'ENSDARG00000103395']
    
    out_ids = mapping.convert_ids(ncbi_ids, NCBI_ID, ENS_ID, out_format = list)
    assert sorted(ens_ids) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(ens_ids, ENS_ID, NCBI_ID, out_format = list)
    assert sorted(ncbi_ids) ==  sorted(out_ids)

def test_mapping_list_of_ids_ens_and_symbol():

    ens_ids = ['ENSDARG00000029019', 'ENSDARG00000070913', 'ENSDARG00000079440',
                'ENSDARG00000090873', 'ENSDARG00000097710']
    symbols = ['epb41b', 'sox2', 'coro2ba', 'ccl34a.4', 'si:dkey-58f6.3']
    
    out_ids = mapping.convert_ids(ens_ids, ENS_ID, SYMBOL, out_format = list)
    assert sorted(symbols) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(symbols, SYMBOL, ENS_ID, out_format = list)
    assert sorted(ens_ids) ==  sorted(out_ids)


# ------------------------
# test add_mapped_column
# ------------------------

def test_mapping_list_of_ids_to_mapped_df():

    ids = [ZFIN_ID, ENS_ID, SYMBOL, NCBI_ID]

    # Generate pairs where order matters
    pairs = list(permutations(ids, 2))

    for pair in pairs:

        generated_df, true_df = generate_mapped_df(pair[0], pair[1])
        assert_frame_equal(generated_df, true_df) 

def test_mapping_df_of_ids_to_mapped_df():

    ids = [ZFIN_ID, ENS_ID, SYMBOL, NCBI_ID]

    # Generate pairs where order matters
    pairs = list(permutations(ids, 2))

    for pair in pairs:

        generated_df, true_df = generate_mapped_df(pair[0], pair[1])
        assert_frame_equal(generated_df, true_df, False) 

def test_non_normalized_id_type():
    # TODO
    pass

def generate_mapped_df(in_id, out_id, list = True):

    id_string_dict = {
        NCBI_ID: 'ncbi',
        ZFIN_ID: 'zfin',
        ENS_ID: 'ens',
        SYMBOL: 'symbol',
    }

    in_id_str = id_string_dict[in_id]
    out_id_str = id_string_dict[out_id]

    in_path = 'tests/data/in_data/mapping/'+in_id_str+'_genes.txt'
    true_data_path = 'tests/data/out_data/mapping/'+in_id_str+'_to_'+out_id_str+'.txt'

    in_df = pd.read_csv(in_path, sep = '\t', dtype = str)

    in_list = in_df[in_id].to_list()

    true_data_out_df = pd.read_csv(true_data_path, sep = '\t', dtype = str).sort_values(
        by = in_id, ascending=True).reset_index(drop = True).dropna()
    
    if list == True:
        in_list = in_df[in_id].to_list()
        out_ids = mapping.add_mapped_column(in_list, in_id, out_id)
    else:
        out_ids = mapping.add_mapped_column(in_df, in_id, out_id)
        
    out_ids_sorted = out_ids.sort_values(by = in_id, ascending=True).reset_index(drop = True).dropna()

    return out_ids_sorted, true_data_out_df

# ------------------------
# test convert_to_human and convert_toz_zebrafish
# ------------------------

def test_mapping_list_of_ids_to_ortho_df():

    ids = [ZFIN_ID, ENS_ID, SYMBOL, NCBI_ID]
    for option in ids:
        
        generated_df, true_df = generate_ortho_df(HUMAN_ID, option)
        assert_frame_equal(generated_df, true_df[option].to_frame()) 

        generated_df, true_df = generate_ortho_df(option, HUMAN_ID)
        assert_frame_equal(generated_df, true_df[HUMAN_ID].to_frame()) 

def test_mapping_df_of_ids_to_ortho_df():

    ids = [ZFIN_ID, ENS_ID, SYMBOL, NCBI_ID]
    for option in ids:
        
        generated_df, true_df = generate_ortho_df(HUMAN_ID, option, False)
        assert_frame_equal(generated_df, true_df[option].to_frame()) 

        generated_df, true_df = generate_ortho_df(option, HUMAN_ID, False)
        assert_frame_equal(generated_df, true_df[HUMAN_ID].to_frame()) 

# ------------------------
# add_mapped_ortholog_column
# ------------------------

def test_mapping_df_of_ids_to_ortho_df_keep_mapping():

    ids = [ZFIN_ID, ENS_ID, SYMBOL, NCBI_ID]
    for option in ids:
        
        generated_df, true_df = generate_ortho_df(HUMAN_ID, option, False, True)
        assert_frame_equal(generated_df, true_df) 

        generated_df, true_df = generate_ortho_df(HUMAN_ID, option, False, True)
        assert_frame_equal(generated_df, true_df) 

        generated_df, true_df = generate_ortho_df(option, HUMAN_ID, False, True)
        assert_frame_equal(generated_df, true_df) 

def generate_ortho_df(in_id, out_id, list = True, add_column = False):

    id_string_dict = {
        NCBI_ID: 'ncbi',
        ZFIN_ID: 'zfin',
        ENS_ID: 'ens',
        SYMBOL: 'symbol',
        HUMAN_ID: 'human'
    }

    in_id_str = id_string_dict[in_id]
    out_id_str = id_string_dict[out_id]

    in_path = 'tests/data/in_data/mapping/'+in_id_str+'_genes.txt'
    true_data_path = 'tests/data/out_data/mapping/'+in_id_str+'_to_'+out_id_str+'.txt'

    in_df = pd.read_csv(in_path, sep = '\t', dtype = str)

    in_list = in_df[in_id].to_list()
    true_data_out_df = pd.read_csv(true_data_path, sep = '\t', dtype = str).sort_values(
        by = out_id, ascending=True).reset_index(drop = True).dropna()

    if add_column:
        out_ids = mapping.add_mapped_ortholog_column(in_df, in_id, out_id)
    else:
        if out_id == HUMAN_ID:
            if list == True:
                in_list = in_df[in_id].to_list()
                out_ids = mapping.convert_to_human(in_list, in_id)
            else:
                out_ids = mapping.convert_to_human(in_list, in_id)
        else:
            if list == True:
                in_list = in_df[in_id].to_list()
                out_ids = mapping.convert_to_zebrafish(in_list, out_id)
            else:
                out_ids = mapping.convert_to_zebrafish(in_list, out_id)
    
    if type(out_ids) == pd.Series:
        out_ids = out_ids.to_frame()

    out_ids_sorted = out_ids.sort_values(by = out_id, ascending=True).reset_index(drop = True).dropna()

    return out_ids_sorted, true_data_out_df

