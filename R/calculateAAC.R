    .checkFASTA <- function(x) {
      . <- NULL
      V1 <- NULL
      fastalist.check <-
        lapply(x, function(x) protr::protcheck(x))
      fastalist.check <-
        do.call(rbind, fastalist.check) %>%
        as.data.frame(.) %>%
        rownames_to_column("ID") %>%
        filter(V1 == FALSE)
    }

    #' calculateAAC
    #' @title Calculate Amino Acid Composition (AAC) Descriptor
    #' @param x A data.frame containing gene/protein names and their
    #' fasta sequences.
    #' @return A length 20 named vector for the data input.
    #' @seealso See \code{\link{calculateDC}} and
    #' \code{\link{calculateTC}} for Dipeptide Composition and
    #' Tripeptide Composition descriptors.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom dplyr filter
    #' @importFrom dplyr left_join
    #' @importFrom dplyr bind_rows
    #' @importFrom tidyr gather
    #' @importFrom tibble rownames_to_column
    #' @description This function calculates Amino Acid Composition
    #' (AAC) descriptor for the data input.
    #' @export
    #' @references
    #' Dey, L., Chakraborty, S., and Mukhopadhyay, A. (2020).
    #' Machine learning techniques for sequence-based prediction of
    #' viral–host interactions between SARS-CoV-2 and human proteins.
    #' \emph{Biomed. J.} 43, 438–450.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateAAC(UP000464024_df)
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

      f <- .checkFASTA(fastalist)
      if (nrow(f) > 0) {
        stop("Fastalist has unrecognized amino acid type")
      }


      # 20 Amino Acid Abbreviation
      AADict <- c(
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
      )


      # calculate amino acid composition descriptor
      AAC_calc <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], "")) %>%
        lapply(., function(x) {
          summary(factor(x[[1]], levels = AADict),
            maxsum = 20
          ) / lengths(x)
        })


      AAC_calc_df <-
        bind_rows(AAC_calc, .id = "identifier") %>%
        gather("AAC", "AACfreq", 2:ncol(.)) %>%
        spread(AAC, AACfreq) %>%
        as.data.frame(.)

      return(AAC_calc_df)
    }
