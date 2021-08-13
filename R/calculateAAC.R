#' calculateAAC
#' @title Calculate Amino Acid Composition (AAC) Descriptor
#' @param x A data.frame containing gene/protein names and their fasta sequences.
#' @return A length 20 named vector for the data input.
#' @seealso See \code{\link{calculateDC}} and \code{\link{calculateTC}}
#' for Dipeptide Composition and Tripeptide Composition descriptors.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom dplyr bind_rows
#' @importFrom tidyr separate
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @description This function calculates Amino Acid Composition (AAC) descriptor for the data input.
#' @export
#' @references
#' Dey, L., Chakraborty, S., and Mukhopadhyay, A. (2020). Machine learning techniques
#' for sequence-based prediction of viral–host interactions between SARS-CoV-2 and human proteins.
#' \emph{Biomed. J.} 43, 438–450.
#' @examples
#' x <- readRDS(
#'     system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP")
#' )
#' x_df <- calculateAAC(x)
#' head(x_df, n = 2L)
calculateAAC <- function(x) {
    . <- NULL
    V1 <- NULL
    AAC <- NULL
    AACfreq <- NULL
    count_AAC <- NULL



    if (!is.data.frame(x)) x <- is.data.frame(x)

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


    # 20 Amino Acid Abbrevation Dictionary from
    # https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties

    AADict <- c(
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    )

    # nchar of each fasta sequence
    nchar_count <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], "")) %>%
        lapply(., function(x) length(x[[1]]))

    # convert the list to the df
    nchar_df <-
        bind_rows(nchar_count) %>%
        gather("identifier", "nchar", seq_len(ncol(.)))

    # calculate amino acid composition descriptor
    AAC_calc <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], "")) %>%
        lapply(., function(x) {
            summary(factor(x[[1]], levels = AADict), maxsum = 20)
        })




    AAC_calc_df <-
        bind_rows(AAC_calc, .id = "identifier") %>%
        gather("AAC", "count_AAC", 2:ncol(.)) %>%
        left_join(., nchar_df, by = "identifier") %>%
        mutate(AACfreq = count_AAC / nchar) %>%
        select(1, 2, 5) %>%
        spread(AAC, AACfreq) %>%
        as.data.frame(.)

    return(AAC_calc_df)
}
