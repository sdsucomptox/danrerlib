# Usage Principles

## Workflow

To import and begin using `danRerLib` in your Python project, import `danRerLib` as:

```
import danRerLib as drl
```

To import a specific module, for example the mapping module, you can import as:

```
from danRerLib import mapping
```

A typical workflow of `danRerLib` depends on the feature you wish to use (see [Introduction](introduction.md)). Please the {ref}`tutorial page <tutorials>` for some typical workflows.

## Zebrafish Gene ID Types

This toolkit supports conversions between a variety of gene id types for the zebrafish. To properly identify the correct gene id, please utilize the following reference. The currently supported zebrafish Gene ID types are:

| Gene ID Type | Referenced as | Description | More Information |
| --| --| --|  -- | 
| NCBI Gene ID | `'NCBI Gene ID' (str)` or `danRerLib.NCBI_ID` | Integer Gene IF managed by NCBI, also known as Entrez Gene ID | [NCBI Website](https://www.ncbi.nlm.nih.gov/) |
| ZFIN ID | `'ZFIN ID' (str)` or `danRerLib.ZFIN_ID` | Each ID starts with 'ZFIN-' and is managed by the Zebrafish Information Network | [ZFIN Website](https://zfin.org/) | 
| Ensembl ID | `'Ensembl ID' (str)` or `danRerLib.ENS_ID` | Each ID starts with 'ENSDAR' and is managed by ensembl | [Ensembl Website](https://useast.ensembl.org/index.html) | 
| Symbol| `'Symbol' (str)` or `danRerLib.SYMBOL` | A alphabetic gene id type where the nomenclature is managed by ZFIN and further supported by NCBI | [ZFIN Website](https://zfin.org/) | 

Mor information can also be found in {ref}`Tutorials <tutorials>` and the {ref}`API Reference <APIRef>`. 

## Human Gene ID Type

This toolkit supports orthology mapping between zebrafish and human. It is possible to map to and from any of the above zebrafish gene id types listed and the supported human gene id type. The currently supported human gene ID type:

| Gene ID Type | Referenced as | Description | More Information |
| --| --| --|  -- | 
| NCBI Gene ID | 'NCBI Gene ID' (_str_) or NCBI_ID | Integer Gene IF managed by NCBI, also known as Entrez Gene ID | [NCBI Website](https://www.ncbi.nlm.nih.gov/) |

More information can also be found in {ref}`Tutorials <tutorials>` and the {ref}`API Reference <APIRef>` for the mapping module. 

## Pathway Databases

This toolkit supports the retrieval and usage of a variety of pathway and annotation databases. Each database has been developed using unique modules within the toolkit. The currently supported databases are:

| Database | Module | Short Description | More Information | 
| -- | -- | -- | -- | 
| KEGG Pathway | `danRerLib.KEGG` | Collection of curated biochemical pathways | [KEGG Pathway Website](https://www.genome.jp/kegg/pathway.html) | 
| KEGG Disease | `danRerLib.KEGG` | Database cataloging human diseases and genes | [KEGG Disease Website](https://www.genome.jp/kegg/disease/)|
| Gene Ontology Biological Processes (GO BP) | `danRerLib.GO` | Standardized vocabulary for biological processes | [Gene Ontology Website](https://geneontology.org/) | 
| Gene Ontology Molecular Function (GO MF) | `danRerLib.GO` | Description of molecular activities of genes| [Gene Ontology Website](https://geneontology.org/) | 
| Gene Ontology Cellular Components (GO CC) | `danRerLib.GO` | Definition of subcellular structures and locations	| [Gene Ontology Website](https://geneontology.org/) | 

## Enrichment Types

This toolkit supports enrichment testing for the pathway and annotation databases listed above. Enrichment testing, also known as pathway enrichment analysis or functional enrichment analysis, is used to determine if a set of genes associated with a specific biological process, function, or pathway is statistically overrepresented or underrepresented compared to what would be expected by chance. A variety of statistical tests can be used to determine over or under representation. The currently supported methods are:

| Enrichment Analysis Method| Function | Short Description | More Information|
|-------|---------|--------|---------|
| Fisher's Exact Test | `danRerLib.enrichment.fishers()`  | Statistical test for pathway enrichment  | [More Info](#)|
| Logistic Regression  | `danRerLib.enrichment.logistic()`| Regression-based method for enrichment | [More Info](#)  |