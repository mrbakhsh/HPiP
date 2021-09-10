test_getHPI = function() {
  x <- matrix(c(1, 2, 3, 1), nrow = 2, ncol = 2, byrow = TRUE)
  y <- matrix(c(0, 3, 2, 1), nrow = 2, ncol = 2, byrow = TRUE)
  RUnit::checkEquals(ncol(getHPI(x, y, 'combine')), 4L)
  RUnit::checkEquals(ncol(getHPI(x, y, 'kron.prod')), 4L)

}


