test_getpositivePPI = function() {
  local = tempdir()
  df <- get_positivePPI(organism.taxID = 2697049,
  access.key = '81bb3b5a6bd9a8084a7be71f0963ab1e',
  filename = "PositiveInt.RData", path = local)
  checkTrue(is.data.frame(df) == TRUE)
  checkEquals(ncol(df), 6L)

}
