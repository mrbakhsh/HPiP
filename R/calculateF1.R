# calculate sum of squared length of single amino acid repeats (SARs) in the entire protein sequence
F1_calc <-
    function(x) {
        . <- NULL
        V1 <- NULL
        V2 <- NULL
        AAcid <- NULL

        string <- strsplit(x, "")[[1]]
        rl <- rle(string)
        lst <- lapply(split(
            seq_along(string),
            rep(seq_along(rl$values), rl$lengths)
        ), range)
        names(lst) <- rl$values
        rl <-
            do.call(rbind, lst) %>%
            as.data.frame(.) %>%
            rownames_to_column("AAcid") %>%
            mutate(AAcid = str_remove_all(AAcid, ".\\d+")) %>%
            mutate(length = ((V2 - V1) + 1)^2) %>%
            group_by(AAcid) %>%
            summarise(Freq = sum(length))
        return(rl)
    }

#' calculateF1
#' @title Calculate Sum of Squared length of Single Amino Acid Repeats (SARs)
#' @param x A data.frame containing gene/protein names and their fasta sequences.
#' @return A length 20 named vector for the data input.
#' @seealso See \code{\link{calculateF2}}
#' for maximum of the sum of squared length of SARs in a window of 6 residues descriptor.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr bind_rows
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr separate
#' @importFrom stringr str_remove_all
#' @description This function calculates sum of quared length of Single Amino Acid Repeats (SARs)
#' in the entire protein sequence.
#' @export calculateF1
#' @references
#' Alguwaizani, S., Park, B., Zhou, X., Huang, D.-S., and Han, K. (2018).
#' Predicting interactions between virus and host proteins using repeat patterns and composition of amino acids.
#' \emph{J. Healthc. Eng.} 2018.
#' @examples
#' x <- readRDS(
#'     system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP")
#' )
#' x_df <- calculateF1(x)
#' head(x_df, n = 2L)
calculateF1 <- function(x) {
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
        rownames_to_column("ID") %>%
        filter(V1 == FALSE)

    if (nrow(fastalist.check) > 0) {
        stop("Fastalist has unrecognized amino acid type")
    }

    p <- # apply a function
        lapply(fastalist, function(x) F1_calc(x[[1]]))
    F1_calc_df <-
        bind_rows(p, .id = "identifier") %>%
        spread(AAcid, Freq) %>%
        replace(is.na(.), 0)

    return(F1_calc_df)
}
