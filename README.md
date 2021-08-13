## HPiP: Host-Pathogen Interaction Prediction

![HPiP ematic](https://github.com/mrbakhsh/HPiP/blob/main/vignettes/Figures/HPiP_PackagePipline.png)

The HPiP (host-pathogen interaction prediction) is an R package for automated prediction of HP-PPIs using structural and physicochemical descriptors computed from amino acid-composition of host and pathogen proteins. The proposed package can effectively address data shortages and data unavailability for HP-PPI network reconstructions. Moreover, establishing computational frameworks in that regard will reveal mechanistic insights into infectious diseases and suggest potential HP-PPI targets, thus narrowing down the range of possible candidates for subsequent wet-lab experimental validations.

### Installation from github:
```{r}
## install.packages("devtools")
install_github("mrbakhsh/HPiP")
library(HPiP)
```

