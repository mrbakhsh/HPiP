<!-- badges: start -->
[![](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/github/last-commit/mrbakhsh/HPiP.svg)](https://github.com/mrbakhsh/HPiP/commits/main)
[![License: MIT (&gt;=
3)](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![DOI:10.1093/bioadv/vbac038](http://img.shields.io/badge/DOI-10./bioadv/vbac038-5680C1.svg)](https://doi.org/10.1093/bioadv/vbac038)
<!-- badges: end -->


# Host-Pathogen Interaction Prediction (HPiP)

HPiP (host-pathogen interaction prediction) is an R package for automated prediction of host-pathogen protein-protein interactions (HP-PPIs) using structural and physicochemical descriptors computed from amino acid-composition of host and pathogen proteins. Briefly, HPiP extracts gold-standard of experimentally verified HP-PPIs (i.e., positive interactions) from public repository, construct negative interactions via negative sampling, retrieve and convert protein sequences to numerical representation via various descriptors, applies multivariate feature selection based on correlation and recursive feature elimination (RFE)-embedded, and finally applies ensemble averaging to predict interactions. Taken together, we hope that the HPiP package not only contributes a useful predictor to accelerate the exploration of host-pathogen PPIs, but also provides some meaningful insights into host-pathogen relationships.
## Installation

You can install the `HPiP` from bioconductor using:

```r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
}
BiocManager::install("HPiP")
```

To view documentation for the version of this package installed in your system, start R and enter:

```r
browseVignettes("HPiP")
```

To install the development version in `R`, run:
  
```r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
devtools::install_github("mrbakhsh/HPiP")
```

## Adapting the pipeline for different datasets

The utility and performance of the proposed package were demonstrated using three different case studies, and data analysis codes are available from https://github.com/BabuLab-UofR/HPiP_pub, where guidelines and sample datasets are also offered for testing purposes.

## Contribute

Check the github page for [source code](https://github.com/BabuLab-UofR/HPiP)

## License
This project is licensed under the MIT License - see the LICENSE.md file for more details.

If using these scripts in your data analyses pipelines, please cite our paper:

Rahmatbakhsh,M. et al. (2022) HPiP: an R/Bioconductor package for predicting host–pathogen protein–protein interactions from protein sequences using ensemble machine learning approach. Bioinforma. Adv., 2, vbac038.
