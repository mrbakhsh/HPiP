    #' var_imp
    #' @title Variable Importance Plot
    #' @param x A list of elements returned from RFE analysis.
    #' See \code{\link[caret]{rfe}} for more details
    #' @param cex.x The size of x axis label.
    #' @param cex.y The size of y axis label.
    #' @return Variable Importance Plot.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}.
    #' @importFrom caret varImp
    #' @importFrom ggplot2 ggplot
    #' @importFrom ggplot2 aes
    #' @importFrom ggplot2 geom_point
    #' @importFrom ggplot2 labs
    #' @importFrom ggplot2 theme_classic
    #' @importFrom ggplot2 theme
    #' @importFrom ggplot2 geom_text
    #' @importFrom ggplot2 ggtitle
    #' @importFrom ggplot2 element_text
    #' @importFrom stats reorder
    #' @importFrom ggplot2 unit
    #' @description A graphical display of variable importance of selected
    #' features.
    #' @export
    #' @examples
    #' data('example_data')
    #' x <- na.omit(example_data)
    #' #perform feature selection
    #' s <- FSmethod(x, type = 'both',
    #' cor.cutoff = 0.7, resampling.method = "repeatedcv",
    #' iter = 5, repeats = 3, metric = "ROC", verbose = TRUE)
    #' var_imp(s$rf.result$rfProfile, cex.x = 10, cex.y = 10)
    var_imp <- function(x, cex.x = 1, cex.y = 2) {
      feature <- NULL
      importance <- NULL

      # examine variable importance for the selected features
      varimp_data <-
        data.frame(
          feature = row.names(varImp(x))[seq_len(length(predictors(x)))],
          importance = varImp(x)[seq_len(length(predictors(x))), 1]
        )

      ggplot(
        data = varimp_data,
        aes(x = reorder(feature, -importance), y = importance, color = feature)
      ) +
        geom_point(size = 2) +
        labs(x = "Features", y = "Variable Importance") +
        theme_classic() +
        theme(legend.position = "none") +
        geom_text(aes(label = round(importance, 2)),
                  vjust = 1.6,
                  color = "white", size = 4
        ) +
        ggtitle("Variable importance for the selected features") +
        theme(plot.title = element_text(
          color = "red", size = 14,
          face = "bold.italic"
        )) +
        theme(
          text = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = cex.x, colour = "black"),
          axis.text.y = element_text(size = cex.y, colour = "black")
        ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
        theme(axis.ticks.length = unit(.2, "cm"))
    }
