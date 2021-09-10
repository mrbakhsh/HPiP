test_calculateBE = function() {
  data('UP000464024_df')
  x_df <- calculateBE(UP000464024_df)
  x_df <- x_df[, -1]
  checkEqualsNumeric(sum(as.matrix(x_df)), 304, tolerance = 1e-2)
}
