# R script used to generate the data input for enrichment analysis
library(HPiP)


data('example_data')
features <- example_data[, -2]
gd <- example_data[, c(1,2)]
 gd <- na.omit(gd)
ppi <-pred_ensembel(features,gd,classifier = c("avNNet", "svmRadial", "rf"),
                    resampling.method = "cv",ncross = 2,verboseIter = FALSE,
                    plots = FALSE,filename = "plots.pdf")
#extract predicted interactions
pred_interaction <- ppi[["predicted_interactions"]]
predicted_PPIs <-pred_interaction[, c(2,3)]


use_data(predicted_PPIs, overwrite = TRUE)
