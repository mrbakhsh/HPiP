    #' calculateDC
    #' @title Calculate Dipeptide Composition (DC) Descriptor
    #' @param x A data.frame containing gene/protein names and their fasta
    #' sequences.
    #' @return A length 400 named vector for the data input.
    #' @seealso See \code{\link{calculateAAC}} and \code{\link{calculateTC}}
    #' for Amino Acid Composition and Tripeptide Composition descriptors.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @description This function calculates Dipeptide Composition (DC)
    #' descriptor for data input.
    #' @export
    #' @references
    #' Bhasin, M., and Raghava, G. P. S. (2004).
    #' Classification of nuclear receptors based on amino acid composition
    #' and dipeptide composition. \emph{J. Biol. Chem.} 279, 23262–23266.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateDC(UP000464024_df)
    #' head(x_df, n = 2L)
    calculateDC <- function(x) {
      . <- NULL
      V1 <- NULL
      DCfreq <- NULL
      count_DC <- NULL
      DC <- NULL

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


      DClist <-
        apply(expand.grid(vect1, vect2), 1, paste, collapse = "")

      # calculate dipeptide composition descriptor
      DC_calc <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], "")) %>%
        lapply(., function(x) {
          summary(factor(paste(x[[1]][-length(x[[1]])], x[[1]][-1],
            sep = ""
          ), levels = DClist), maxsum = 400) / (lengths(x) - 1)
        })

      DC_calc_df <-
        bind_rows(DC_calc, .id = "identifier") %>%
        gather("DC", "DCfreq", 2:ncol(.)) %>%
        spread(., DC, DCfreq) %>%
        as.data.frame(.)

      return(DC_calc_df)
    }
