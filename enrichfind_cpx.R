    #' enrichfind_cpx
    #' @title Functional Enrichment Analysis for Predicted Modules
    #' @param predcpx Predicted modules resulted from
    #' \code{\link[HPiP]{run_clustering}}.
    #' @param threshold Custom p-value threshold for significance.
    #' @param sources A vector of data sources to use.
    #' See \code{\link[gprofiler2]{gost}} for more details.
    #' @param p.corrction.method The algorithm used for multiple testing
    #' correction;defaults to 'bonferroni'.
    #' See \code{\link[gprofiler2]{gost}} for more details.
    #' @param org An organism name;defaults to 'hsapiens'.
    #' See \code{\link[gprofiler2]{gost}} for more details.
    #' @return A data.frame with the enrichment analysis results.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @description This function uses \code{\link[gprofiler2]{gost}} function
    #' in \code{gprofiler2} package to perform functional enrichment analysis
    #' for predicted modules.
    #' @export

    enrichfind_cpx <-
        function(predcpx,
                 threshold = 0.05,
                 sources = c("GO", "KEGG"),
                 p.corrction.method = "bonferroni",
                 org = "hsapiens") {

            p_value <- NULL
            predcpx <-
                as.data.frame(predcpx)

            agg.data <-
                aggregate(list(predcpx[, 1]),
                          by = list(predcpx[, 2]),
                          FUN = function(x) {
                              paste(unique(x),
                                    collapse = ";"
                              )
                          }
                )
            colnames(agg.data)[2] <- "members"
            agg.data$members <-
                vapply(lapply(strsplit(agg.data$members, ";"), unique),
                       paste, character(1L),
                       collapse = ";"
                )


            indcpx <-
                str_split(agg.data$members, ";")
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
