#' calculateTC_Sm
#' @title Calculate Tripeptide Composition (TC) Descriptor from Biochemical Similarity Classes
#' @param x A data.frame containing gene/protein names and their fasta sequences.
#' @return A length 216  named vector for the data input.
#' @seealso See \code{\link{calculateTC}} for Tripeptide Composition descriptor.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#' @importFrom magrittr set_rownames
#' @description This function calculates Tripeptide Composition (TC) descriptor from biochemical similarity classes.
#' @export
#' @references
#' Ahmed, I., Witbooi, P., and Christoffels, A. (2018). Prediction of human-Bacillus
#' anthracis protein–protein interactions using multi-layer neural network.
#' \emph{Bioinformatics} 34, 4159–4164.
#'
#' Cui, G., Fang, C., and Han, K. (2012). Prediction of protein-protein
#' interactions between viruses and human by an SVM model.
#' \emph{BMC bioinformatics}, 1–10.
#'
#' @examples
#' x <- readRDS(
#'     system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP")
#' )
#' x_df <- calculateTC_Sm(x)
#' head(x_df, n = 2L)
calculateTC_Sm <- function(x) {
    if (!is.data.frame(x)) x <- is.data.frame(x)

    . <- NULL
    V1 <- NULL

    # convert data frame to list
    fastalist <-
        as.list(unlist(x[, 2]))
    names(fastalist) <-
        unlist(x[, 1])

    # check if there is any unrecognized amino acid
    fastalist.check <-
        lapply(fastalist, function(x) protr::protcheck(x))
    fastalist.check <-
        do.call(rbind, fastalist.check) %>%
        as.data.frame(.) %>%
        rownames_to_column("ID") %>%
        filter(V1 == FALSE)

    if (nrow(fastalist.check) > 0) {
        stop("Fastalist has unrecognized amino acid type")
    }


    classes <- vector("list", 6)
    classes[[1]] <- c("I", "V", "L", "M")
    classes[[2]] <- c("F", "Y", "W")
    classes[[3]] <- c("H", "K", "R")
    classes[[4]] <- c("D", "E")
    classes[[5]] <- c("Q", "N", "T", "P")
    classes[[6]] <- c("A", "C", "G", "S")


    vspace <- expand.grid(seq_len(6), seq_len(6), seq_len(6))

    BSlict <- vector("list", 216)

    for (i in seq_len(216)) {
        tmp <- as.vector(outer(classes[[vspace[i, 1]]],
            classes[[vspace[i, 2]]], paste,
            sep = ""
        ))
        BSlict[[i]] <- as.vector(outer(tmp,
            classes[[vspace[i, 3]]], paste,
            sep = ""
        ))
    }
    names(BSlict) <- paste(1:length(BSlict))



    # compute TC_Sm
    p <- function(x, BSlict) {
        f <-
            x %>%
            lapply(., function(x) strsplit(x[[1]], "")) %>%
            lapply(., function(x) {
                paste(paste(x[[1]][-c(length(x[[1]]), length(x[[1]]) - 1)],
                    x[[1]][-c(1, length(x[[1]]))],
                    sep = ""
                ), x[[1]][-c(1, 2)],
                sep = ""
                )
            })

        inter <- unlist(lapply(BSlict, function(X) {
            lapply(f, function(Y) {
                length(which(Y %in% X))
            })
        }), recursive = FALSE)

        s <-
            do.call(cbind, inter) %>%
            set_rownames(names(x))
        s <- t(apply(s, 1, function(x) (x - min(x)) / (max(x) - min(x))))
        return(s)
    }

    TC_Sm <- lapply(fastalist, function(x) p(x[[1]], BSlict))
    TC_Sm_df <-
        do.call(rbind, TC_Sm)

    TC_Sm_df <-
        TC_Sm_df %>%
        as.data.frame(.) %>%
        set_rownames(names(TC_Sm)) %>%
        rownames_to_column("identifier")

    return(TC_Sm_df)
}
