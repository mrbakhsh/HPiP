    #' calculateKSAAP
    #' @title Calculate k-spaced Amino Acid Pairs (KSAAP) Descriptor
    #' @param x A data.frame containing gene/protein names and their fasta
    #' sequences.
    #' @param spc A number of spaces separating two adjacent residues by a
    #' distance of spc, which can be any number up to two less than
    #' the length of the peptide; default to 3.
    #' @return A length 400 named vector for the data input.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @description This function calculates k-spaced Amino Acid Pairs
    #' (KSAAP) Descriptor for data input.
    #' This function is adapted from the \code{\link[ftrCOOL]{CkSAApair}}
    #'  function in the ftrCOOL package.
    #' @export
    #' @references
    #' Kao, H.-J., Nguyen, V.-N., Huang, K.-Y., Chang, W.-C., and Lee, T.-Y.
    #' (2020).SuccSite: incorporating amino acid composition and informative
    #' k-spaced amino acid pairs to identify protein succinylation sites.
    #' \emph{Genomics. Proteomics Bioinformatics} 18, 208–219.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateKSAAP(UP000464024_df)
    #' head(x_df, n = 2L)
    calculateKSAAP <- function(x, spc = 3) {
      if (!is.data.frame(x)) {
        stop("Input data must be data.frame")
      }

      . <- NULL
      V1 <- NULL

      AADict <-
        list(
          "A" = 1, "C" = 2, "D" = 3, "E" = 4, "F" = 5, "G" = 6,
          "H" = 7, "I" = 8, "K" = 9, "L" = 10, "M" = 11, "N" = 12,
          "P" = 13, "Q" = 14, "R" = 15, "S" = 16, "T" = 17, "V" = 18,
          "W" = 19, "Y" = 20
        )


      vect <- c(
        "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
      )


      DClist <-
        apply(expand.grid(vect, vect), 1, paste, collapse = " ")

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

      seqs <- unlist(x[, 2])
      label <- names(fastalist)
      nSeqs <- length(fastalist)
      len <- 1

      attributesM <- matrix(0, ncol = 400, nrow = nSeqs)
      row.names(attributesM) <- label

      attributeName <- vector()
      for (i in seq_len(len)) {
        attributeName <- c(attributeName, gsub(
          " ", strrep("s", spc[i]),
          DClist
        ))
      }
      colnames(attributesM) <- attributeName


      for (n in seq_len(nSeqs)) {
        seq <- seqs[n]
        seqChars <- unlist(strsplit(seq, split = ""))
        lenSeq <- length(seqChars)

        for (i in seq_len(len)) {
          temp1 <- seqChars[seq_len(lenSeq - spc[i] - 1)]
          temp2 <- seqChars[((spc[i] + 1) + 1):(lenSeq)]
          kmers <- paste(temp1, temp2, sep = "")
          tbkmers <- table(kmers)
          nmtbkmers <- names(tbkmers)
          for (j in seq_along(tbkmers)) {
            tmp <- unlist(strsplit(nmtbkmers[j], split = ""))
            index <- (as.numeric(AADict[tmp[1]]) - 1) * 20 +
              as.numeric(AADict[tmp[2]])
            index <- index + ((i - 1) * 400)
            attributesM[n, index] <- tbkmers[j]
          }
        }
      }

      # normalize
      seqLen <- vapply(seqs, nchar, FUN.VALUE = numeric(1))
      attributesM <- attributesM / seqLen

      attributesdf <-
        attributesM %>%
        as.data.frame(.) %>%
        rownames_to_column("identifier")

      return(attributesdf)
    }
