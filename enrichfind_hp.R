    #' enrichfind_hp
    #' @title Functional Enrichment Analysis of all Host Proteins
    #' @param ppi A data.frame containing pathogen proteins in the
    #' first column and host proteins in the second column.
    #' @param threshold Custom p-value threshold for significance.
    #' @param sources A vector of data sources to use.
    #' See \code{\link[gprofiler2]{gost}} for more details.
    #' @param p.corrction.method The algorithm used for multiple testing
    #' correction;defaults to 'bonferroni'.
    #' See \code{\link[gprofiler2]{gost}} for more details.
    #' @param org An organism name;defaults to 'hsapiens'.
    #' See \code{\link[gprofiler2]{gost}} for more details.
    #' @return A data.frame with the enrichment analysis results.
    #' @seealso See \code{\link{enrichplot}} for plotting enrichment analysis.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @description This function uses \code{\link[gprofiler2]{gost}} function
    #' in \code{gprofiler2} package to perfrom functional enrichment analysis
    #' for all predicted host proteins in the high-confidence network.
    #' @export
    #' @examples
    #' data('predicted_PPIs')
    #' #perform enrichment
    #' enrich.df <- enrichfind_hp(predicted_PPIs,
    #' threshold = 0.05,
    #' sources = c("GO", "KEGG"),
    #' p.corrction.method = "bonferroni",
    #' org = "hsapiens")
    
    
    enrichfind_hp <-
        function(ppi, threshold = 0.05,
                 sources = c("GO", "KEGG"),
                 p.corrction.method = "bonferroni",
                 org = "hsapiens") {
            p_value <- NULL
            ppi <- as.data.frame(ppi)
            
            prot <- unique(ppi[,2])
    
            annotCov <- 
                gprofiler2::gost(prot,
                                 significant = TRUE,
                                 exclude_iea = TRUE,
                                 evcodes = TRUE,
                                 user_threshold = threshold,
                                 sources = sources,
                                 correction_method = p.corrction.method,
                                 organism = org
                )
            
            df <-
                annotCov[["result"]]
            ans <-
                arrange(df, p_value)
            
            return(ans)
        }