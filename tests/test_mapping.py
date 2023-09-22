from danRerLib import mapping

# MAPPING COMBINATIONS ARE:
# ZFIN_ID and NCBI_ID done
# ZFIN_ID and SYMBOL done
# ZFIN_ID and ENS_ID done
# NCBI_ID and SYMBOL done
# NCBI_ID and ENS_ID done
# SYMBOL and ENS_ID

def test_mapping_list_of_ids_ncbi_and_zfin():

    ncbi_ids = ['562552', '100150233', '405818', '564006', '796163']
    
    zfin_ids = ['ZDB-GENE-001222-1', 'ZDB-GENE-030131-4309',
        'ZDB-GENE-040426-2356', 'ZDB-GENE-081001-1', 'ZDB-GENE-131127-614']
    
    out_ids = mapping.convert_ids(ncbi_ids, mapping.NCBI_ID, mapping.ZFIN_ID, out_format = list)
    assert sorted(zfin_ids) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(zfin_ids, mapping.ZFIN_ID, mapping.NCBI_ID, out_format = list)
    assert sorted(ncbi_ids) ==  sorted(out_ids)

def test_mapping_list_of_ids_zfin_and_symbol():

    zfin_ids = ['ZDB-GENE-060503-764', 'ZDB-GENE-010817-2', 'ZDB-GENE-070705-89',
                'ZDB-GENE-030131-296', 'ZDB-GENE-040608-1']
    
    symbols = ['si:ch211-152p11.4', 'angpt2b', 'trhde.2', 'ttpal', 'harbi1']
    
    out_ids = mapping.convert_ids(zfin_ids, mapping.ZFIN_ID, mapping.SYMBOL, out_format = list)
    assert sorted(symbols) ==  sorted(out_ids)
    
    out_ids = mapping.convert_ids(symbols, mapping.SYMBOL, mapping.ZFIN_ID, out_format = list)
    assert sorted(zfin_ids) ==  sorted(out_ids)


def test_mapping_list_of_ids_zfin_and_ens():

    zfin_ids = ['ZDB-GENE-040120-4', 'ZDB-LINCRNAG-110914-87', 'ZDB-GENE-141216-345',
                'ZDB-GENE-050809-56', 'ZDB-GENE-091204-383']
    ens_ids = ['ENSDARG00000008785', 'ENSDARG00000092478', 'ENSDARG00000100251', 
               'ENSDARG00000102792', 'ENSDARG00000096078']
    
    out_ids = mapping.convert_ids(zfin_ids, mapping.ZFIN_ID, mapping.ENS_ID, out_format = list)
    assert ens_ids ==  out_ids

    out_ids = mapping.convert_ids(ens_ids, mapping.ENS_ID, mapping.ZFIN_ID, out_format = list)
    assert sorted(zfin_ids) ==  sorted(out_ids)


def test_mapping_list_of_ids_ncbi_and_symbol():

    ncbi_ids = ['402853', '65231', '565228', '436656', '324103']
    symbols = ['uts2b', 'mdkb', 'si:ch211-215a10.4', 'krt97', 'wu:fc18g07']
    
    out_ids = mapping.convert_ids(ncbi_ids, mapping.NCBI_ID, mapping.SYMBOL, out_format = list)
    assert sorted(symbols) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(symbols, mapping.SYMBOL, mapping.NCBI_ID, out_format = list)
    assert sorted(ncbi_ids) ==  sorted(out_ids)

def test_mapping_list_of_ids_ncbi_and_ens():

    ncbi_ids = ['550419', '794572', '101885923', '678616', '100535123']
    ens_ids = ['ENSDARG00000039730', 'ENSDARG00000091990', 'ENSDARG00000087413', 
                'ENSDARG00000058372', 'ENSDARG00000103395']
    
    out_ids = mapping.convert_ids(ncbi_ids, mapping.NCBI_ID, mapping.ENS_ID, out_format = list)
    assert sorted(ens_ids) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(ens_ids, mapping.ENS_ID, mapping.NCBI_ID, out_format = list)
    assert sorted(ncbi_ids) ==  sorted(out_ids)

def test_mapping_list_of_ids_ens_and_symbol():

    ens_ids = ['ENSDARG00000029019', 'ENSDARG00000070913', 'ENSDARG00000079440',
                'ENSDARG00000090873', 'ENSDARG00000097710']
    symbols = ['epb41b', 'sox2', 'coro2ba', 'ccl34a.4', 'si:dkey-58f6.3']
    
    out_ids = mapping.convert_ids(ens_ids, mapping.ENS_ID, mapping.SYMBOL, out_format = list)
    assert sorted(symbols) ==  sorted(out_ids)

    out_ids = mapping.convert_ids(symbols, mapping.SYMBOL, mapping.ENS_ID, out_format = list)
    assert sorted(ens_ids) ==  sorted(out_ids)
        
