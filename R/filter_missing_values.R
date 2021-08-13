#' filter_missing_values
#' @title Drop the Missing Values Above a Certain Threshold
#' @param x A numeric matrix as input.
#' @param max_miss_rate Maximal missing rate allow for a feature;default is 20.
#' @return A dataframe with features with missingness rate of more than user-defined
#' threshold.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @description Given an input matrix, compute the missingness rate for each
#' features and keep only features with missing rate more than user-defined
#' percentage.
#' @export
#' @examples
#' x <- matrix(1:10, ncol = 2)
#' x[, 2] <- NA
#' filter_missing_values(x, 30)
filter_missing_values <-
    function(x, max_miss_rate = 20) {
        if (is.matrix(x)) x <- as.data.frame(x)
        if (missing(max_miss_rate)) stop("Must provide percentage value for features removal")
        if (length(colnames(x)) == 0) stop("Column (feature) names are missing...")

        # remove column with more than 20% missing values
        miss_prec <-
            colSums(is.na(x)) / nrow(x) * 100
        print(miss_prec[miss_prec > 20])
        col_miss <-
            names(miss_prec[miss_prec > 20])

        if (length(col_miss) == 0) {
            x
        } else {

            # Removing the columns with more than 20% missing value
            x[, c(col_miss)] <- NULL
        }
        return(x)
    }
