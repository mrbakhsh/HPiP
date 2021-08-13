#' enrichfindP
#' @title Functional Enrichment Analysis for Pathogen Interactors in the High-Confidence Network.
#' @param ppi A data.frame containing pathogen proteins in the first column and host
#' proteins in the second column.
#' @param threshold Custom p-value threshold for significance.
#' @param sources A vector of data sources to use. See \code{\link[gprofiler2]{gost}} for more details.
#' @param p.corrction.method The algorithm used for multiple testing correction;defaults to 'bonferroni'.
#' See \code{\link[gprofiler2]{gost}} for more details.
#' @param org organism name;defaults to 'hsapiens'.
#' See \code{\link[gprofiler2]{gost}} for more details.
#' @return A data.frame with the enrichment analysis results.
#' @seealso See \code{\link{enrichplot}} for plotting enrichment analysis.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom purrr map_df
#' @importFrom stats aggregate
#' @importFrom stringr str_split
#' @description This function uses \code{\link[gprofiler2]{gost}} function in \code{gprofiler2} package
#' to perfrom functional enrichment analysis for pathogen interactors in the high-confidence network.
#' @export
#' @examples
#' ppi <- readRDS(
#'     system.file("examples/HC_predPPIs.rds", package = "HPiP")
#' )
#' enrich.df <- enrichfindP(ppi,
#'     threshold = 0.05,
#'     sources = c("GO", "KEGG"),
#'     p.corrction.method = "bonferroni",
#'     org = "hsapiens"
#' )
enrichfindP <-
    function(ppi, threshold = 0.05,
             sources = c("GO", "KEGG"),
             p.corrction.method = "bonferroni",
             org = "hsapiens") {
        p_value <- NULL
        ppi <- as.data.frame(ppi)

        agg.data <-
            aggregate(list(ppi[, 2]),
                by = list(ppi[, 1]),
                FUN = function(x) {
                      paste(unique(x),
                          collapse = ";"
                      )
                  }
            )
        colnames(agg.data)[2] <- "interactors"
        agg.data$interactors <-
            vapply(lapply(strsplit(agg.data$interactors, ";"), unique),
                paste, character(1L),
                collapse = ";"
            )


        indcpx <-
            str_split(agg.data$interactors, ";")
        names(indcpx) <- agg.data$Group.1

        annotCov <- lapply(indcpx, function(x) {
              gprofiler2::gost(x,
                  significant = TRUE,
                  exclude_iea = TRUE,
                  evcodes = TRUE,
                  user_threshold = threshold,
                  sources = sources,
                  correction_method = p.corrction.method,
                  organism = org
              )
          })

        df <-
            lapply(annotCov, function(x) as.data.frame(x[[1]]))
        ans <-
            map_df(df, ~ as.data.frame(.x), .id = "id")

        ans <-
            arrange(ans, p_value)

        return(ans)
    }
