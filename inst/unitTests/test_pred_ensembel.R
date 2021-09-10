test_pred_ensembel <- function() {
  data('example_data')
  features <- example_data[, -2]
  gd <- example_data[, c(1,2)]
  gd <- na.omit(gd)
  ppi <-pred_ensembel(features,gd,
                      classifier = c("avNNet", "svmRadial", "ranger"),
                      resampling.method = "cv",ncross = 5,
                      verboseIter = FALSE,plots = FALSE,
                      filename = "plots.pdf")
  #extract predicted interactions
  pred_interaction <- ppi[["predicted_interactions"]]
  checkTrue(is.data.frame(pred_interaction) == TRUE)
}
