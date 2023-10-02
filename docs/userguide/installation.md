# Installation

`danRerLib` was built with some dependencies. Therefore, to install `danRerLib`, you will need a proper version of Python and the required dependencies. The current dependencies are:

- python >=3.9,<3.13
- numpy 1.26.0
- pandas 2.1.1

Please see respective package information if dependencies are not already installed.

## Quick Install

### Conda

If you use `conda`, you can install `danRerLib` as:

```
conda install danrerlib
```

### PyPi/pip

If you use `pip`, you can install `danRerLib` with:

```
pip install danrerlib
```

### Development Version

If you would like to make changes to the library, either for your own purposes or to contribute, you should install the development version of this package. Please see the [Contributing Guidelines](../contributing.md) in the Developer's Guide to follow the recommended steps.

## Installation Guide

It is recommended to use a python virtual environment or conda environment to manage all python packages and projects. If you are a new Python user, installing and managing Python packages can be complicated. We recommend using a package management software such as [Anaconda](https://www.anaconda.com/). We refer the user elsewhere to make further installation decisions about Python packages in general. 

### Beginning users

On all of Windows, macOS, and Linux:

- Install Anaconda (it installs all packages you need and all other tools mentioned below).
- Utilize the Anaconda App to manage environments and package downloads

### Advanced users (recommended method)

**Conda**

- Create a conda environment specifically for the task or project you’re working on.

```
# Best practice, use an environment rather than install in the base env
conda create -n my-env
conda activate my-env
# install dependencies
conda install numpy
conda install pandas
# install library
conda install danrerlib
```

**pip**

- create a python virtual environment for the task or project you are working on. 

```
# Create a virtual environment
python -m venv .
mkdir src
# Activate virtual environment and install
cd bin
source activate
pip install numpy
pip install pandas
pip install danrerlib
```