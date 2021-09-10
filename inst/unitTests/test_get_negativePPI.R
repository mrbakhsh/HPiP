test_get_negativePPI = function() {

  prot1 <- c("P0DTC4", "P0DTC5", "P0DTC9")
  prot2 <- c("Q9Y679", "Q9NW15", "Q9NXF8")
  TPset <- c("P0DTC4~P31948", "P0DTC8~Q13438")
  TN_PPI <- get_negativePPI(prot1, prot2, TPset)
  checkEquals(nrow(TN_PPI), 9L)

}
