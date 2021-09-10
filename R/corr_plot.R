    #' corr_plot
    #' @title Plot Correlation Matrix between Input Features
    #' @param cormat A correlation matrix.
    #' @param method The visualization method of correlation matrix;
    #' defaults to number.
    #' See \code{\link[corrplot]{corrplot}} for more details.
    #' @param cex The size of x/y axis label.
    #' @return A correlation plot.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}.
    #' @importFrom corrplot corrplot
    #' @description A graphical display of a correlation matrix.
    #' @export
    #' @examples
    #' data('example_data')
    #' x <- na.omit(example_data)
    #' #perform feature selection
    #' s <- FSmethod(x, type = 'both',
    #' cor.cutoff = 0.7, resampling.method = "repeatedcv",
    #' iter = 5, repeats = 3, metric = "ROC", verbose = TRUE)
    #' corr_plot(s$cor.result$corProfile, method = 'square' , cex = 0.5)
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
