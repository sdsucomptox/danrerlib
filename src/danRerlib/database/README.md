# Database Description

In this document, you will find all information about the database and the relevant information. 

## Downloaded Mapping/Orthology Files

All downloaded files are treated as raw data as they are taken straight from the ZFIN or NCBI website. All data files listed in this section should be in the raw_data sub-directory.

| File Name | Original File Name | File Description | Database | Downloaded From | Notes | 
|--|--|--| --|--|--|
| `zfin_to_ncbi_V1.txt` | `gene_2023.06.27.txt` | mapping of ZFIN gene ids to NCBI gene ids | ZFIN | [ZFIN downloads](https://zfin.org/downloads) -> Sequence Data -> ZFIN Marker associations to NCBI Gene data | need to remove date at top of file | 
| `zfin_to_ensembl_V1.txt` | `ensembl_1_to_1_2023.06.27.txt` | mapping of ZFIN gene ids to ensembl gene ids | ZFIN | [ZFIN downloads](https://zfin.org/downloads) -> Sequence Data -> ZFIN Marker associations to Ensembl IDs | need to remove date at top of file | 
| `zfish_human_orthology_V1.txt` | `human_orthos_2023.06.27.txt` | mapping of ZFIN ids to human NCBI ids | ZFIN | [ZFIN downloads](https://zfin.org/downloads) -> Orthology Data -> Human and Zebrafish Orthology | need to remove date at top of file |
| `Homo_sapiens.gene_info_V1` | `Homo_sapiens.gene_info.gx` | full list of all human genes | NCBI | [NCBI API](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/) and downloaded homo sapiens file [Homo_sapiens.gene_info.gz](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz) | needs to be unzipped | 

## Generated Files

| File Name | Original File Name | File Description | buildig function | 
|--|--|--|--|
| `master_gene_mapping_file_V1.txt` | `master_gene_mapping_file_zfin.txt` | the generated gene mapping file from ZFIN to symbol to Ensembl ID to NCBI ID | 
| `master_ortho_mapping_file_V1.txt` | `master_ortho_mapping_file.txt` | the generated orthology mapping file from ZFIN to symbol to Ensembl ID to NCBI ID to Human NCBI Gene ID | 

## KEGG file

All data files for KEGG are located in the KEGG sub-directory. 

| File Name | Original File Name | File Description | Database | Downloaded From | Notes | 
|--|--|--| --|--|--|
| `kegg_ids_dre_V1.txt` | `kegg_ids_dre.txt` | list of all KEGG ids for zebrafish | KEGG | [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html) using link` https://rest.kegg.jp/list/pathway/dre` | using API | 
| `kegg_ids_hsa_V1.txt` | `kegg_ids_hsa.txt` | list of all KEGG ids for humans | KEGG | [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html) using link` https://rest.kegg.jp/list/pathway/hsa` | using API | 

## Important Variables 
| Variable | Description | Example Gene | 
| ZFIN ID	| ZFIN ID as listed on the ZFIN website | 
| Symbol | Gene symbol as defined using ZFIN naming convention | 
| NCBI Gene ID | integer Gene ID from NCBI | 
| Ensembl ID | ensembl, starts with ENSDAR | 
| Human NCBI Gene ID | integer Gene ID from NCBI | 3630 | 

## Additional Information 

### Downloading data from ZFIN

The ZFIN downloads page has many different datasets to download. Each of the files are in a `.tsv` format. There does not seem to be a way to download the data using an API or a link address so you have to download the files directly from the webpage. To do so, you fine the desired file on the webpage and click the `tsv` button. This should download the file directly to your downloads folder. It is a tab delimited file but it is saved as a `.txt` file. 

#### To update a dataset 

- download the desire `tsv` file
- move the file to the database directory
- remove the top two lines in the `.txt` file that have the date
- rename the file to `[base_file_name]_V[version_number].txt`

### Downloading data from KEGG

I have utilized the KEGG API to download all the necessary data. If you visit the link in the table above, you will see all the data on your webpage. I used the pandas package to retrieve the data and save it to a .txt file.

```python
import pandas as pd
link_address = # your link, ex: 'https://rest.kegg.jp/list/pathway/dre'
df = pd.read_csv(link_address, sep = '\t')
file_name =     # your desired file name, ex: kegg_ids_dre.txt 
                # make sure the full path is included if needed
df.to_csv(file_name, sep = '\t', index=False)
```

#### To update a dataset

- follow the code block above and insert the desired version number

### Downloading data from NCBI

I have utilized the NCBI API to download all the necessary data. You can either retrieve as it was done for KEGG, or just download it directly. I downloaded it directly because it is zipped and I needed to unzip it. Yes, I could have done this in python too..

#### To update a dataset

- visit the correct API location, ex: `https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/`
- click your desired file, ex:` Homo_sapiens.gene_info.gz`
- move to database directory and name the file accordingly including version number

## On going wants/ todos

- Bigger want it to create a sql database. 
- Develop a safeguard for database building if the downloaded data has different column names than expected. I think I should change the database building file to grab the columns that I want rather than drop the columns that I don't want. This in theory should stop the potential issue or random columns being added or name changes. 