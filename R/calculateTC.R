    #' calculateTC
    #' @title Calculate Tripeptide Composition (TC) Descriptor
    #' @param x A data.frame containing gene/protein names and their fasta
    #' sequences.
    #' @return A length 8,000  named vector for the data input.
    #' @seealso See \code{\link{calculateAAC}},\code{\link{calculateDC}} and
    #' \code{\link{calculateTC_Sm}} for
    #' Amino Acid Composition, Dipeptide Composition and Tripeptide
    #' Composition (TC) Descriptor from Biochemical Similarity Classes.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @description This function calculates Tripeptide Composition (TC)
    #' descriptor for data input.
    #' @export
    #' @references
    #' Liao, B., Jiang, J.-B., Zeng, Q.-G., and Zhu, W. (2011).
    #' Predicting apoptosis protein subcellular location with
    #' PseAAC by incorporating tripeptide composition.
    #' \emph{Protein Pept. Lett.} 18, 1086–1092
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateTC(UP000464024_df)
    #' head(x_df, n = 2L)
    calculateTC <- function(x) {
      . <- NULL
      V1 <- NULL
      AAC <- NULL
      TCfreq <- NULL
      count_TC <- NULL
      TC <- NULL

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


      vect1 <- c(
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
      )
      vect2 <- c(
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
      )
      vect3 <- c(
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
      )



      TClist <-
        apply(expand.grid(vect1, vect2, vect3), 1, paste, collapse = "")



      # calculate dipeptide composition descriptor
      TC_calc <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], "")) %>%
        lapply(., function(x) {
          summary(factor(paste(paste(x[[1]][-c(
            length(x[[1]]),
            length(x[[1]]) - 1
          )],
          x[[1]][-c(1, length(x[[1]]))],
          sep = ""
          ), x[[1]][-c(1, 2)],
          sep = ""
          ), levels = TClist), maxsum = 8000) / (lengths(x) - 2)
        })


      TC_calc_df <-
        bind_rows(TC_calc, .id = "identifier") %>%
        gather("TC", "TCfreq", 2:ncol(.)) %>%
        spread(., TC, TCfreq) %>%
        as.data.frame(.) %>%
        replace(is.na(.), 0)


      return(TC_calc_df)
    }
