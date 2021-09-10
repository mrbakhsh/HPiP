    #' calculateBE
    #' @title Tranform a Seqeuence into Binary Encoding (BE)
    #' @param x A data.frame containing gene/protein names and their
    #' fasta sequences.
    #' @return A length 400 named vector for the data input.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom tibble add_column
    #' @description This function transform each residue in a peptide
    #' into 20 coding values.
    #' @export
    #' @references
    #' Al-Barakati, H. J., Saigo, H., and Newman, R. H. (2019).
    #' RF-GlutarySite: a random forest based predictor for glutarylation sites.
    #' \emph{Mol. Omi.} 15, 189–204.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateBE(UP000464024_df)
    #' head(x_df, n = 2L)
    calculateBE <- function(x) {
      . <- NULL
      V1 <- NULL


      if (!is.data.frame(x)) {
        stop("Input data must be data.frame")
      }

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

      dict <- list(
        "A" = 1, "C" = 2, "D" = 3, "E" = 4, "F" = 5, "G" = 6, "H" = 7, "I" = 8,
        "K" = 9, "L" = 10, "M" = 11, "N" = 12, "P" = 13, "Q" = 14, "R" = 15,
        "S" = 16, "T" = 17, "V" = 18, "W" = 19, "Y" = 20
      )


      bin_m <- function(x) {
        fp <- matrix(0L, nrow = length(x), ncol = 400)
        vect <- c(
          "A", "C", "D", "E", "F", "G", "H", "I",
          "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
        )
        n <- rep(vect, 20)
        colnames(fp) <- n

        charSeq <- unlist(strsplit(x[[1]], split = ""))
        pos <- as.numeric(dict[unique(charSeq)])
        lenSeqs <- length(pos)
        rng <- (0:(lenSeqs[1] - 1)) * 20
        pos1 <- rng + pos
        for (i in 1) fp[i, pos1] <- 1L
        return(fp)
      }


      BM_list <- lapply(fastalist, function(x) bin_m(x[[1]]))


      # generate df from list
      tbl <- do.call(rbind, BM_list)
      rownames(tbl) <- names(BM_list)

      tbl <-
        tbl %>%
        as.data.frame(.) %>%
        add_column(.name_repair = c("minimal")) %>%
        rownames_to_column("identifier")



      return(tbl)
    }
