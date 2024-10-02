# iSAPIENS
An in silico tool to optimize protein interaction energies

## Installation

iSAPIENS has been successfully tested on Linux systems


### Environment
- Download and install conda environment for 64-bit linux or Mac (https://conda.io/miniconda.html) using Python 3.7
  - Install in command line: ./Miniconda3-latest-Linux-x86_64.sh (Linux)
  - Open a new terminal window or source your bash with: source ~/.bashrc
- Create iSAPIENS environment in command line: conda create --name iSAPIENS
- Activate iSAPIENS environment in command line: conda activate iSAPIENS
- Install required packages by running the following in command line:
  - conda install -c r r
  - conda install -c bioconda r-bio3d
  - conda install -c r r-doparallel
  - Download R's Interpol package (https://cran.r-project.org/src/contrib/Archive/Interpol/Interpol_1.3.1.tar.gz) and install manually:
    - Within R terminal execute: install.packages("Interpol_1.3.1.tar.gz", source = TRUE, repos = NULL )
  - 
