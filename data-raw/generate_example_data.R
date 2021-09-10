# R script used to generate the sample data for pred_ensembel function


#Loading the required packages
library(dplyr)
library(HPiP)




###load the unlabeled HP-PPIs
data('unlabel_data')


###Construct labeled data
#1.viral data
data('viral_se')
#extract features from viral_se
counts_v <- assays(viral_se)$counts
#extract row.names from viral_Se
rnames_v <- row.names(counts_v)
#Construct labeed data

#2.host data
data('host_se')
#extract features from host_se
counts_h <- assays(host_se)$counts
#extract row.names from viral_Se
rnames_h <- row.names(counts_h)


#3.loading gold-standard data
data('Gold_ReferenceSet')
gd <- Gold_ReferenceSet

x1_viral <- matrix(NA, nrow = nrow(gd), ncol = ncol(counts_v))
for (i in 1:nrow(gd))
  x1_viral[i, ] <- counts_v[which(gd$Pathogen_Protein[i] == rnames_v), ]

x1_host <- matrix(NA, nrow = nrow(gd), ncol = ncol(counts_h))
for (i in 1:nrow(gd))
  x1_host[i, ] <- counts_h[which(gd$Host_Protein[i] == rnames_h), ]

#construct gd data
x <- getHPI(x1_viral,x1_host, type = "combine")
x <- as.data.frame(x)
x <- cbind(gd$PPI, gd$class, x)
colnames(x)[1:2] <- c("PPI", "class")

sampledata <-
  dplyr::bind_rows(x, unlabel_data)

set.seed(100)
sampledata <-
  sample_n(sampledata, 500)
#only select the first column
example_data <-
  sampledata[, c(1:100)]

# save
use_data(example_data, overwrite = TRUE)
