#' get_positivePPI
#' @title Get Positive Reference Host-Pathogen Protein-Protein Interactions (HP-PPIs)
#' @param organism.taxID Taxonomy identifier for the pathogen. All identifiers are
#' available at \url{http://webservice.thebiogrid.org/organisms/?accessKey=81bb3b5a6bd9a8084a7be71f0963ab1e}.
#' @return A Data.frame containing true positive protein-protein interactions for the selected
#' pathogen.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @seealso See \code{\link{get_negativePPI}} for generating negative
#' protein-protein interaction.
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom readr cols
#' @importFrom readr col_double
#' @importFrom readr col_character
#' @description Get positive reference host-pathogen protein-protein interactions from BioGRID.
#' @export
#' @examples
#' positive_reflists <- get_positivePPI(organism.taxID = 2697049)
get_positivePPI <- function(organism.taxID = 2697049) {
    `Official Symbol Interactor A` <- NULL
    `Official Symbol Interactor B` <- NULL

    # construct query to retrive all interactiosn for the requested organism
    url <-
        "http://webservice.thebiogrid.org/interactions/?"
    query <- paste(url,
        "accessKey=81bb3b5a6bd9a8084a7be71f0963ab1e",
        "interSpeciesExcluded=false",
        "selfInteractionsExcluded=true",
        "includeHeader=true",
        paste0("taxId=", organism.taxID),
        sep = "&"
    )


    # number of results
    n <-
        httr::GET(paste(query, "format=count", sep = "&"))
    httr::stop_for_status(n)


    n <- as.numeric(
        httr::content(n, type = "text/csv", encoding = "UTF-8", col_types = "i")
    )

    # retrieve results in batches
    tpINT <- lapply(seq(1, n, 10000), function(m) {
        l <- paste(query, paste0("start=", m), sep = "&")
        request <-
            httr::GET(l)
        httr::stop_for_status(request)

        # parse the biogird results
        biog_result <- httr::content(request,
            type = "text/tab-separated-values",
            encoding = "UTF-8",
            col_types = "cccccccccccccccccccccccc"
        )
        q <- subset(biog_result, select = c(
            "Official Symbol Interactor A",
            "Official Symbol Interactor B",
            "Organism Interactor A",
            "Organism Interactor B",
            "Experimental System Type"
        ))

        return(q)
    })

    tpINT <- do.call("rbind", tpINT)
    dataTP <-
        tpINT
    dataTP <- # remove duplicate
        tpINT[!duplicated(apply(tpINT, 1, function(x) paste(sort(x), collapse = ""))), ]

    dataTP <-
        dataTP %>%
        mutate(PPI = paste(`Official Symbol Interactor A`,
            `Official Symbol Interactor B`,
            sep = "~"
        )) %>%
        dplyr::select(6, seq_len(5))


    return(dataTP)
}
