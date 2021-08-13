# Calculate maximum of the sum of squared length of SARs in a window of 6 residues
F2_calc <-
    function(x) {
        . <- NULL
        V1 <- NULL
        V2 <- NULL
        AAcid <- NULL

        f <-
            x %>%
            lapply(., function(x) strsplit(x[[1]], "")) %>%
            lapply(., function(x) {
                paste(paste(x[[1]][-c(length(x[[1]]), length(x[[1]]) - 1)],
                    x[[1]][-c(1, length(x[[1]]))], x[[1]][-c(seq_len(2), length(x[[1]]))],
                    x[[1]][-c(seq_len(3), length(x[[1]]))],
                    x[[1]][-c(seq_len(4), length(x[[1]]))],
                    x[[1]][-c(seq_len(5), length(x[[1]]))],
                    sep = ""
                ), x[[1]][-c(1, 2, 3, 4, 5, 6)],
                sep = ""
                )
            })


        x <- unlist(f)

        x <- paste(x, collapse = " ")


        x <- unlist(strsplit(x, split = ""))


        rl <- rle(x)
        lst <- lapply(split(
            seq_along(x),
            rep(seq_along(rl$values), rl$lengths)
        ), range)

        names(lst) <- rl$values
        rl <-
            do.call(rbind, lst) %>%
            as.data.frame(.) %>%
            tibble::rownames_to_column("AAcid") %>%
            mutate(AAcid = str_remove_all(AAcid, ".\\d+")) %>%
            mutate(length = ((V2 - V1) + 1)^2) %>%
            group_by(AAcid) %>%
            summarise(Freq = max(length)) %>%
            filter(AAcid != "X.")
        return(rl)
    }



#' calculateF2
#' @title Calculate Maximum of the Sum of Squared Length of Single Amino Acid Repeats (SARs)
#' @param x A data.frame containing gene/protein names and their fasta sequences.
#' @return A length 20 named vector for the data input.
#' @seealso See \code{\link{calculateF1}}
#' for Sum of Squared length of Single Amino Acid Repeats (SARs) descriptor.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr bind_rows
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom stringr str_remove_all
#' @description This function calculates maximum of the sum of Single Amino Acid Repeats (SARs)
#' in a window of 6 residue.
#' @export calculateF2
#' @references
#' Alguwaizani, S., Park, B., Zhou, X., Huang, D.-S., and Han, K. (2018).
#' Predicting interactions between virus and host proteins using repeat patterns and composition of amino acids.
#' \emph{J. Healthc. Eng.} 2018.
#' @examples
#' x <- readRDS(
#'     system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP")
#' )
#' x_df <- calculateF2(x)
#' head(x_df, n = 2L)
calculateF2 <- function(x) {
    if (!is.data.frame(x)) {
        stop("Input data must be data.frame")
    }

    . <- NULL
    V1 <- NULL
    AAcid <- NULL
    Freq <- NULL

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


    p <- # apply a function
        lapply(fastalist, function(x) F2_calc(x[[1]]))
    F2_calc_df <-
        bind_rows(p, .id = "identifier") %>%
        tidyr::spread(AAcid, Freq) %>%
        replace(is.na(.), 0)


    return(F2_calc_df)
}
