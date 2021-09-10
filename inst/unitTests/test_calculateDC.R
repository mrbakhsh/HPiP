test_calculateDC = function() {
  data('UP000464024_df')
  x_df <- calculateDC(UP000464024_df)
  x_df <- x_df[, -1]
  checkEqualsNumeric(sum(as.matrix(x_df)), 17, tolerance = 1e-2)
}
