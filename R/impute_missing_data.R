#' impute_missing_data
#' @title Impute missing Values per Features (i.e., Columns)
#' @param x A numeric matrix as input.
#' @param method Imputation method for missing values (mean, median or zero).
#' @return Imputed matrix.
#' @author Matineh Rahmatbakhsh \email{matinerb.94@gmail.com}
#' @importFrom stats median
#' @description Given an input matrix, impute the missing values via three approaches
#' including mean, median or zero.
#' @export
#' @examples
#' x <- matrix(1:10, ncol = 2)
#' x[1:3, 2] <- NA
#' row.names(x) <- c("A", "B", "C", "D", "E")
#' colnames(x) <- c("col1", "col2")
#' impute_missing_data(x, method = "mean")
#' impute_missing_data(x, method = "median")
#' impute_missing_data(x, method = "zero")
impute_missing_data <-
    function(x, method = c("mean", "median", "zero")) {
        if (!is.matrix(x)) x <- as.matrix(x)


        if (length(rownames(x)) == 0) {
            stop("Row names are missing...")
        }
        if (length(colnames(x)) == 0) {
            stop("Column names are missing...")
        }


        impute.f <- function(vec, method = c("mean", "median", "zero")) {
            if (missing(method)) stop("Must provide at least one imputation approach")

            if (!method %in% c("mean", "median", "zero")) {
                stop("Imputation method must be one of mean or median...")
            }

            if (!is.numeric(vec)) vec <- as.numeric(vec)

            if (method == "mean") {
                vec[is.na(vec)] <- mean(vec, na.rm = TRUE)
            }
            if (method == "median") {
                vec[is.na(vec)] <- median(vec, na.rm = TRUE)
            }
            if (method == "zero") {
                vec[is.na(vec)] <- 0
            }
            return(vec)
        }
        names <- rownames(x)
        x.imputed <-
            apply(x, 2, impute.f, method = method)
        row.names(x.imputed) <- names
        return(x.imputed)
    }
