
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

To install the development version in `R`, run:
  
```r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
devtools::install_github("mrbakhsh/HPiP")
```


