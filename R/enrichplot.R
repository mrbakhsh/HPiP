#' enrichplot
#' @title Plot the Enrichment Reuslt
#' @param x A data.frame with the enrichment analysis results.
#' @param low Colours for low.
#' @param high Colours for high.
#' @param cex.size Text size.
#' @return An enrichment plot.
#' @seealso See \code{\link{enrichfindP}} for functional enrichment analysis.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom ggplot2 scale_colour_gradient2
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 element_rect
#' @description This function plots the enrichment result.
#' @export
#' @examples
#' enrich.result <- readRDS(system.file("examples/enirchment_result.rds", package = "HPiP"))
#' enrichplot(enrich.result,
#'     low = "blue",
#'     high = "red",
#'     cex.size = 12
#' )
enrichplot <-
    function(x,
             low = "blue",
             high = "red",
             cex.size = 15) {
        xaxis <- NULL
        term_name <- NULL
        intersection_size <- NULL
        p_value <- NULL

        x$xaxis <- 1

        print(ggplot(x, aes(xaxis, term_name, color = -log10(p_value))) +
            geom_point(aes(xaxis, term_name, size = intersection_size), stroke = 1, show.legend = TRUE) +
            scale_colour_gradient2(low = low, high = high) +
            theme_bw() +
            scale_x_continuous(expand = c(0, 0), limits = c(1, 1)) +
            theme(text = element_text(size = cex.size, color = "#000000")) +
            theme(axis.text.y = element_text(color = "#000000")) +
            theme(axis.ticks.x = element_blank()) +
            theme(axis.text.x = element_blank()) +
            labs(title = "", x = "", y = "term_name"))
    }
