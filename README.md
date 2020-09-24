# Linking genes to GWAS loci in keratinocytes

## Welcome!
This repository contains part of the code to reproduce the main results of *"Chromatin looping links target genes with genetic risk loci for dermatological traits"*.
This can also be used with HiChIP produced in different cell types or with different GWAS results, such as loci not reaching genome wide significance or credible SNP sets.

The preprint of the manuscript is available at https://www.biorxiv.org/content/10.1101/2020.03.05.973271v2

Please note that this paper is still a work in progress and some parts of the analysis might change in the future.
Moreover, this code has been edited from the original to allow easier usage by other people but it produces the same results as presented in the paper.

## Requirements and installation

To use the scripts from this repository you will need to clone this into your system (linux only). You can do this by running this code.

'''
git clone https://github.com/ChenfuShi/keratinocyte_gene_link.git
'''

Now you need to set up the environment with the python requirements:
- bedtools
- python
- pybedtools
- pandas
- numpy
- scipy
- statistics
- jupyterlab

For a working enviromnet we suggest using conda and the environment.yml file provided.

'''
conda env create --file environment.yml --name keratinocyte_gene_link
conda activate keratinocyte_gene_link
'''

## Preprocessing of your data



## Running the scripts

The main script for linking genes to GWAS loci si located in ./scripts/Main_Linker.ipynb
You can set all your files and settings in this script but the default uses all the data from ./datasets/ which are all the files needed to reproduce the results from the paper.



to do
preprocessing of data, hichip, gwas etc.
conda env
documentation