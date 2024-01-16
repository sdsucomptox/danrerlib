# Danio rerio Library, danRerLib

## Description

danRerLib, short for Danio rerio library, is a comprehensive toolkit designed in Python specifically for zebrafish researchers, focusing on transcriptomics. 

## Documentation

The documentation for the library is currently being hosted on [San Diego State University's Computational Toxicology Laboratory website](https://sdsucomptox.github.io/danrerlib/index.html). 

### Installation

The package requires Python>3.9 and < 3.13; we recommend Python 3.9. Please see the documentation[installation guide](https://sdsucomptox.github.io/danrerlib/userguide/installation.html) for more information.

Via pip:
~~~~
pip install 
~~~~

Via conda: 
~~~~
conda install 
~~~~

### Examples

- Check out our in depth [tutorials](https://sdsucomptox.github.io/danrerlib/tutorials) on the documentation page for step-by-step instructions on how to use danRerLib to its full potential. 
- You can also look at our results directory to replicate a stand-alone workflow. 

## Project History

This project repository was started in July 2023 by [Ashley Valentina Schwartz](https://github.com/ashleyvsch) with the first version of `danRerLib` being developed in September 2023. The project is maintained by [Ashley Valentina Schwartz](https://github.com/ashleyvsch) and the [SDSU Computational Toxicology Laboratory](https://github.com/sdsucomptox). The release of danRerLib occurred in early 2024 and is an evolving repository. 

## Folder Descriptions

`results/`
- This folder contains a reproducible workflow for the package that replicates the published findings in the associated publication. 

`src/danrerlib/`
- This folder contains the source code and database for the package.

`tests/`
- This folder contains tests, developed with pytest, to test the functionality of the package to maintain a stable repository.

`sphinx/`
- This folder contains the workable documentation for the repository including tutorials and the API reference.
- This is build using [spinx documentation builder](https://www.sphinx-doc.org/en/master/).

`docs/`
- This folder contains the documentation that is present on the website.
- This is a duplicate of `sphinx/_build/html/`.

## File Map (partial)

```
danRerLib/
│   README.md
|   pyproject.toml
|   poetry.tock
|   Makefile
│   LICENSE
|   CONTRIBUTING.md
|   CONDUCT.md
|   CHANGELOG.md
└───src/
│   └───danRerLib
│       │   __init__.py
│       │   mapping.py
│       │   KEGG.py
│       │   GO.py
│       │   enrichment.py
│       │   enrichplots.py
│       │   utils.py
│       │   settings.py
|       └───database/
└───tests/
└───sphinx/
└───docs/

```