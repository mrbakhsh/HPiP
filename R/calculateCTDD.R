#' calculateCTDD
#' @title Calculate CTD Descriptors - Distribution (D)
#' @param x A data.frame containing gene/protein names and their fasta sequences.
#' @return A length 105 named vector for the data input.
#' @seealso See \code{\link{calculateCTDC}} and \code{\link{calculateCTDT}}
#' for Composition and Transition descriptors.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom dplyr rename
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr funs
#' @importFrom dplyr row_number
#' @importFrom dplyr summarise_each
#' @importFrom dplyr left_join
#' @importFrom dplyr full_join
#' @importFrom dplyr select
#' @description This function calculates Distribution (D) descriptor for data input.
#' @export
#' @references
#' Dubchak, I., Muchnik, I., Holbrook, S. R., and Kim, S.-H. (1995).
#' Prediction of protein folding class using global description of amino acid sequence.
#' \emph{Proc. Natl. Acad. Sci.} 92, 8700–8704.
#' @examples
#' x <- readRDS(
#'     system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP")
#' )
#' x_df <- calculateCTDD(x)
#' head(x_df, n = 1L)
calculateCTDD <- function(x) {
    . <- NULL
    Freq <- NULL
    ID <- NULL
    V1 <- NULL
    V2 <- NULL
    V3 <- NULL
    V4 <- NULL
    V5 <- NULL
    identifier <- NULL
    category <- NULL
    d <- NULL
    key <- NULL
    value <- NULL

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
        tibble::rownames_to_column("ID") %>%
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

    fastaSplitted <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], ""))


    # nchar of each fasta sequence
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


    # intersect with group 1
    g1.intersect <- unlist(lapply(group1, function(X) {
        lapply(fastaSplitted, function(Y) {
            ifelse(Y[[1]] %in% X, "G1", NA)
        })
    }), recursive = FALSE)

    # intersect with group 2
    g2.intersect <- unlist(lapply(group2, function(X) {
        lapply(fastaSplitted, function(Y) {
            ifelse(Y[[1]] %in% X, "G2", NA)
        })
    }), recursive = FALSE)

    # intersect with group 3
    g3.intersect <- unlist(lapply(group3, function(X) {
        lapply(fastaSplitted, function(Y) {
            ifelse(Y[[1]] %in% X, "G3", NA)
        })
    }), recursive = FALSE)

    list.g1g2 <- Map(list, g1.intersect, g2.intersect)
    list.g1g2 <-
        lapply(list.g1g2, function(x) {
            ifelse(is.na(x[[1]]), x[[2]], x[[1]])
        })

    list.g1g2g3 <- Map(list, list.g1g2, g3.intersect)

    list.g1g2g3 <-
        lapply(list.g1g2g3, function(x) {
            ifelse(is.na(x[[1]]), x[[2]], x[[1]])
        })

    D <- vector("list", length(list.g1g2g3))
    names(D) <-
        names(list.g1g2g3)
    for (i in seq_along(list.g1g2g3)) D[[i]] <- matrix(ncol = 5, nrow = 3)

    for (i in seq_along(list.g1g2g3)) {
        for (j in 1:3) {
            inds <- which(list.g1g2g3[[i]] == paste0("G", j))
            quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
            quartiles[which(quartiles <= 0)] <- 1
            D[[i]][j, ] <- if (length(inds) > 0) {
                (inds[c(1, quartiles, length(inds))]) * 100
            } else {
                0
            }
        }
    }
    tbl2 <- do.call(rbind, D)

    names.tbl2 <-
        as.data.frame(matrix(unlist(names(D)), ncol = 1, byrow = TRUE))

    names.tbl3 <-
        names.tbl2 %>%
        mutate(g1 = "g1") %>%
        mutate(g2 = "g2") %>%
        mutate(g3 = "g3") %>%
        gather("key", "value", -V1) %>%
        mutate(names = paste(V1, value, sep = ".")) %>%
        separate(V1, c("group", "identifier"), sep = "\\.") %>%
        select(1, 2, 5)



    # sort based on the attribute
    d1 <-
        names(fastalist)
    Final_df <-
        names.tbl3[order(match(names.tbl3$identifier, d1)), ]


    d2 <- c(
        "hydrophobicity", "normwaalsvolume", "polarity", "polarizability",
        "charge", "secondarystruct", "solventaccess"
    )
    Final_df <-
        Final_df[order(match(Final_df$group, d2)), ]

    row.names(tbl2) <-
        Final_df$names

    dfOut <-
        tbl2 %>%
        as.data.frame(.) %>%
        rownames_to_column("ID") %>%
        separate(ID, c("group", "identifier"), "\\.", remove = FALSE) %>%
        rename(residue0 = V1) %>%
        rename(residue25 = V2) %>%
        rename(residue50 = V3) %>%
        rename(residue75 = V4) %>%
        rename(residue100 = V5) %>%
        select(-2) %>%
        gather("key", "value", 3:7) %>%
        left_join(., nchar_df, by = "identifier") %>%
        mutate(Freq = value / nchar) %>%
        separate(ID, c("c", "p", "d"), sep = "\\.") %>%
        mutate(category = paste(c, d, key, sep = ".")) %>%
        arrange(d) %>%
        select(4, 9, 8) %>%
        group_by(identifier) %>%
        mutate(row = row_number()) %>%
        spread(category, Freq) %>%
        select(-row) %>%
        group_by(identifier) %>%
        replace(is.na(.), 0) %>%
        summarise_each(funs(max(.)))


    return(dfOut)
}
