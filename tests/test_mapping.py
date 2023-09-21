import unittest
from ..src.danRerLib import mapping

class TestConvertIDs(unittest.TestCase):
    def test_mapping_list_of_ids_zfin_to_ens(self):

        zfin_ids = ['ZDB-GENE-040426-1218',
            'ZDB-GENE-030616-582',
            'ZDB-GENE-050327-38',
            'ZDB-GENE-050302-29',
            'ZDB-MIRNAG-090929-162',
            'ZDB-GENE-041111-305',
            'ZDB-GENE-040426-1159',
            'ZDB-GENE-050227-21',
            'ZDB-GENE-060526-320',
            'ZDB-GENE-040718-307']
        
        ens_ids = ['ENSDARG00000054292',
            'ENSDARG00000038780',
            'ENSDARG00000002406',
            'ENSDARG00000089310',
            'ENSDARG00000100596',
            'ENSDARG00000030881',
            'ENSDARG00000063570',
            'ENSDARG00000099828',
            'ENSDARG00000094489',
            'ENSDARG00000080838']
        
        out_ids = mapping.convert_ids(zfin_ids, mapping.ZFIN_ID, mapping.ENS_ID, out_format = list)
        self.assertListEqual(ens_ids, out_ids)

if __name__ == '__main__':
    unittest.main()


        

        
