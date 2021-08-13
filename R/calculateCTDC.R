#' calculateCTDC
#' @title Calculate CTD Descriptors - Composition (C)
#' @param x A data.frame containing gene/protein names and their fasta sequences.
#' @return A length 21 named vector for the data input.
#' @seealso See \code{\link{calculateCTDT}} and \code{\link{calculateCTDD}}
#' for Transition and Distribution descriptors.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom dplyr rename
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @description This function calculates Composition (C) descriptor for data input.
#' @export
#' @references
#' Dubchak, I., Muchnik, I., Holbrook, S. R., and Kim, S.-H. (1995).
#' Prediction of protein folding class using global description of amino acid sequence.
#' \emph{Proc. Natl. Acad. Sci.} 92, 8700–8704.
#' @examples
#' x <- readRDS(
#'     system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP")
#' )
#' x_df <- calculateCTDC(x)
#' head(x_df, n = 2L)
calculateCTDC <- function(x) {
    . <- NULL
    V1 <- NULL
    value <- NULL
    ID <- NULL
    group <- NULL


    if (!is.data.frame(x)) {
        stop("Input data must be data.frame")
    }

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

    # three-group classification of the 20 amino acids by each attribute
    # retrieved from "https://CRAN.R-project.org/package=protr"

    group1 <- list(
        "hydrophobicity" = c("R", "K", "E", "D", "Q", "N"),
        "normwaalsvolume" = c("G", "A", "S", "T", "P", "D", "C"),
        "polarity" = c("L", "I", "F", "W", "C", "M", "V", "Y"),
        "polarizability" = c("G", "A", "S", "D", "T"),
        "charge" = c("K", "R"),
        "secondarystruct" = c("E", "A", "L", "M", "Q", "K", "R", "H"),
        "solventaccess" = c("A", "L", "F", "C", "G", "I", "V", "W")
    )

    group2 <- list(
        "hydrophobicity" = c("G", "A", "S", "T", "P", "H", "Y"),
        "normwaalsvolume" = c("N", "V", "E", "Q", "I", "L"),
        "polarity" = c("P", "A", "T", "G", "S"),
        "polarizability" = c("C", "P", "N", "V", "E", "Q", "I", "L"),
        "charge" = c(
            "A", "N", "C", "Q", "G", "H", "I", "L",
            "M", "F", "P", "S", "T", "W", "Y", "V"
        ),
        "secondarystruct" = c("V", "I", "Y", "C", "W", "F", "T"),
        "solventaccess" = c("R", "K", "Q", "E", "N", "D")
    )

    group3 <- list(
        "hydrophobicity" = c("C", "L", "V", "I", "M", "F", "W"),
        "normwaalsvolume" = c("M", "H", "K", "F", "R", "Y", "W"),
        "polarity" = c("H", "Q", "R", "K", "N", "E", "D"),
        "polarizability" = c("K", "M", "H", "F", "R", "Y", "W"),
        "charge" = c("D", "E"),
        "secondarystruct" = c("G", "N", "P", "S", "D"),
        "solventaccess" = c("M", "S", "P", "T", "H", "Y")
    )


    # nchar of each fasta sequemce
    nchar_count <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], "")) %>%
        lapply(., function(x) length(x[[1]]))

    # convert the list to the df
    nchar_df <-
        bind_rows(nchar_count) %>%
        gather("identifier", "nchar", seq_len(ncol(.))) %>%
        mutate(nchar = nchar)

    fastaSplitted <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], ""))


    g1.intersect <- unlist(lapply(group1, function(X) {
        lapply(fastaSplitted, function(Y) {
            length(which(Y[[1]] %in% X))
        })
    }), recursive = FALSE)

    g1.intersect <-
        do.call(rbind, g1.intersect) %>%
        as.data.frame(.) %>%
        rename(value = V1) %>%
        rownames_to_column("ID") %>%
        separate(ID, c("group", "identifier"), sep = "\\.") %>%
        mutate(group = paste0("group1", ".", group)) %>%
        left_join(., nchar_df, by = "identifier") %>%
        mutate(value = value / nchar) %>%
        select(2, 1, 3)

    g2.intersect <- unlist(lapply(group2, function(X) {
        lapply(fastaSplitted, function(Y) {
            length(which(Y[[1]] %in% X))
        })
    }), recursive = FALSE)

    g2.intersect <-
        do.call(rbind, g2.intersect) %>%
        as.data.frame(.) %>%
        rename(value = V1) %>%
        rownames_to_column("ID") %>%
        separate(ID, c("group", "identifier"), sep = "\\.") %>%
        mutate(group = paste0("group2", ".", group)) %>%
        left_join(., nchar_df, by = "identifier") %>%
        mutate(value = value / nchar) %>%
        select(2, 1, 3)

    g3.intersect <- unlist(lapply(group3, function(X) {
        lapply(fastaSplitted, function(Y) {
            length(which(Y[[1]] %in% X))
        })
    }), recursive = FALSE)

    g3.intersect <-
        do.call(rbind, g3.intersect) %>%
        as.data.frame(.) %>%
        rename(value = V1) %>%
        rownames_to_column("ID") %>%
        separate(ID, c("group", "identifier"), sep = "\\.") %>%
        mutate(group = paste0("group3", ".", group)) %>%
        left_join(., nchar_df, by = "identifier") %>%
        mutate(value = value / nchar) %>%
        select(2, 1, 3)

    dfOut <-
        rbind(g1.intersect, g2.intersect, g3.intersect) %>%
        spread(group, value)
    dfOut[is.na(dfOut)] <- 0

    return(dfOut)
}
