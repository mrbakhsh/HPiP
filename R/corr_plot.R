#' corr_plot
#' @title Plot Correlation Matrix between Input Features
#' @param cormat A correlation matrix.
#' @param method The visualization method of correlation matrix; defaults to number.
#' See \code{\link[corrplot]{corrplot}} for more details.
#' @param cex The size of x/y axis label.
#' @return A correlation plot.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}.
#' @importFrom corrplot corrplot
#' @description A graphical display of a correlation matrix.
#' @export
#' @examples
#' x <- readRDS(system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP"))
#' x <- calculateAAC(x)
#' x <- x[, -1]
#' x_cor <- cor(as.matrix(x))
#' corr_plot(x_cor, method = "square", cex = 0.9)
corr_plot <- function(cormat, method = "number", cex = 0.9) {
    if (!is.matrix(cormat)) cormat <- as.matrix(cormat)

    # plotting corr matrix
    plot.cor <- corrplot(cormat,
        order = "hclust",
        tl.col = "black",
        method = method,
        tl.cex = cex,
        cl.cex = 1.0
    ) # Make plot
    return(invisible(plot.cor))
}
