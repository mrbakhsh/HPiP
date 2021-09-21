# R package: HPiP

HPiP (host-pathogen interaction prediction) is an R package for automated prediction of host-pathogen protein-protein interactions (HP-PPIs) using structural and physicochemical descriptors computed from amino acid-composition of host and pathogen proteins. 
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

