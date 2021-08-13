#' get_negativePPI
#' @title Construct Negative Reference Host-Pathogen Protein-Protein Interactions (HP-PPIs)
#' @param prot1 A character vector containing pathogen proteins.
#' @param prot2 A character vector containing host proteins.
#' @param TPset A character vector containing positive reference interactions.
#' @return A Data.frame containing true negative interactions.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @references Eid, F.-E., ElHefnawi, M., and Heath, L. S. (2016).
#' DeNovo: virus-host sequence-based protein–protein interaction prediction. Bioinformatics 32, 1144–1150.
#' @seealso See \code{\link{get_positivePPI}} for generating positive
#' protein-protein interaction.
#' @description Construct true negative protein-protein interactions
#' from the positive interactions. In the context of PPI prediction, a negative interaction
#' is a pair of proteins that unlikely to interact. Since there is no experimentally verified
#' non-interacting pair, the negative sampling can be used to construct the negative reference set.
#' The negative sampling can be constructed from a set of host proteins, a set of pathogen proteins, and
#' a list of positive reference interactions between members of host and pathogen proteins (Eid et al., 2016).
#' @export
#' @examples
#' prot1 <- c("P0DTC4", "P0DTC5", "P0DTC9")
#' prot2 <- c("Q9Y679", "Q9NW15", "Q9NXF8")
#' TPset <- c("P0DTC4~P31948", "P0DTC8~Q13438")
#' TN_PPI <- get_negativePPI(prot1, prot2, TPset)
#' head(TN_PPI)
get_negativePPI <-
    function(prot1, prot2, TPset) {
        if (!is.character(prot1)) {
            stop("prot1 must be a character vector")
        }

        if (!is.character(prot2)) {
            stop("prot2 must be a character vector")
        }

        if (!is.character(TPset)) {
            stop("TPset must be a character vector")
        }

        TN_PPI <-
            apply(expand.grid(prot1, prot2), 1, paste, collapse = "~")

        PPI <- setdiff(TN_PPI, TPset)
        TN_PPI <- as.data.frame(PPI)


        return(TN_PPI)
    }
