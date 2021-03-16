# Linking genes to GWAS loci in keratinocytes

## Welcome!
This repository contains part of the code to reproduce the main results of *"Chromatin looping links target genes with genetic risk loci for dermatological traits"*.
This can also be used with HiChIP produced in different cell types or with different GWAS results, such as loci not reaching genome wide significance or credible SNP sets.

The preprint of the manuscript is available at https://www.biorxiv.org/content/10.1101/2020.03.05.973271v2

Please note that this paper is still a work in progress and some parts of the analysis might change in the future.
Moreover, this code has been edited from the original to allow easier usage by other people but it produces the same results as presented in the paper.

We also provide images for all loci studied in the paper in *./all_loci_images/*.

## Requirements and installation

To use the scripts from this repository you will need to clone this into your system (linux only). You can do this by running this code.

```
git clone https://github.com/ChenfuShi/keratinocyte_gene_link.git
```

Now you need to set up the environment with the python requirements:
- bedtools
- python=3.7.3
- pybedtools
- pandas
- numpy
- scipy
- statistics
- jupyterlab
- xlsxwriter

For a working enviromnet we suggest using conda and the environment.yml file provided.

```
conda env create --file environment.yml --name keratinocyte_gene_link
conda activate keratinocyte_gene_link
```

## Preprocessing of your data

If you want to use different data from the ones provided in this repository you need to process your data using the standard tools for HiChIP and Hi-C analysis. More detailed instructions and settings can be found in the methods section of the paper.

- HiChIP and Hi-C data is first preprocessed using Hi-C pro (https://github.com/nservant/HiC-Pro).

- To obtain interactions from HiChIP data use FitHiChIP (https://github.com/ay-lab/FitHiChIP).

- To obtain the peaks and the bedgraphs from HiChIP data use HiChIP-Peaks (https://github.com/ChenfuShi/HiChIP_peaks).

- To obtain TADs from Hi-C data you can use OnTAD (https://github.com/anlin00007/OnTAD) or you can disable this step in the script.

- RNA-seq data must be provided in a CSV file with each row being a ensemble gene id and each column the TPM values for each sample.

- Finally GWAS SNPs must be reported in the following format:
chr  start  end(start+1)  score  snpid_LOCI

If you have different formats feel free to modify the code or ask for help.

## Running the scripts

The main script for linking genes to GWAS loci si located in *./scripts/Main_Linker.ipynb*. you can open this by using jupyterlab.

You can set all your files and settings in the top of this script but the default uses all the data from *./datasets/* which are all the files needed to reproduce the results from the paper.

## citation

```
Chenfu Shi, Helen Ray-Jones, James Ding, Kate Duffus, Yao Fu, Vasanthi Priyadarshini Gaddi, Oliver Gough, Jenny Hankinson, Paul Martin, Amanda McGovern, Annie Yarwood, Patrick Gaffney, Steve Eyre, Magnus Rattray, Richard B. Warren, Gisela Orozco
Chromatin Looping Links Target Genes with Genetic Risk Loci for Dermatological Traits
Journal of Investigative Dermatology, 2021, ISSN 0022-202X,
https://doi.org/10.1016/j.jid.2021.01.015.
(https://www.sciencedirect.com/science/article/pii/S0022202X2100141X)
```

