    .nchar <- function(x) {
      . <- NULL
      # nchar of each fasta sequence
      nchar_count <-
        x %>%
        lapply(., function(x) strsplit(x[[1]], "")) %>%
        lapply(., function(x) length(x[[1]]))

      # convert the list to the df
      nchar_df <-
        bind_rows(nchar_count) %>%
        gather("identifier", "nchar", seq_len(ncol(.))) %>%
        mutate(nchar = nchar)

      return(nchar_df)
    }



    .cld <- function(x, g, nchar_df) {
      . <- NULL
      V1 <- NULL
      value <- NULL
      ID <- NULL
      group <- NULL
      df <- do.call(rbind, x) %>%
        as.data.frame(.) %>%
        rename(value = V1) %>%
        rownames_to_column("ID") %>%
        separate(ID, c("group", "identifier"), sep = "\\.") %>%
        mutate(group = paste0(g, ".", group)) %>%
        left_join(., nchar_df, by = "identifier") %>%
        mutate(value = value / nchar) %>%
        select(2, 1, 3)
      return(df)
    }




    #' calculateCTDC
    #' @title Calculate CTD Descriptors - Composition (C)
    #' @param x A data.frame containing gene/protein names and
    #' their fasta sequences.
    #' @return A length 21 named vector for the data input.
    #' @seealso See \code{\link{calculateCTDT}} and \code{\link{calculateCTDD}}
    #' for Transition and Distribution descriptors.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom dplyr rename
    #' @description This function calculates Composition (C)
    #' descriptor for data input.
    #' @export
    #' @references
    #' Dubchak, I., Muchnik, I., Holbrook, S. R., and Kim, S.-H. (1995).
    #' Prediction of protein folding class using global description of
    #' amino acid sequence.\emph{Proc. Natl. Acad. Sci.} 92, 8700–8704.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateCTDC(UP000464024_df)
    #' head(x_df, n = 2L)
    calculateCTDC <- function(x) {
      . <- NULL
      V1 <- NULL
      value <- NULL
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
      f <- .checkFASTA(fastalist)
      if (nrow(f) > 0) {
        stop("Fastalist has unrecognized amino acid type")
      }

      # three-group classification of the 20 amino acids by each attribute

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


      # nchar of each fasta sequence
      nchar_df <- .nchar(fastalist)


      fastaSplitted <-
        fastalist %>% lapply(., function(x) strsplit(x[[1]], ""))


      g1 <- unlist(lapply(group1, function(X) {
        lapply(fastaSplitted, function(Y) {
          length(which(Y[[1]] %in% X))
        })
      }), recursive = FALSE)


      g2 <- unlist(lapply(group2, function(X) {
        lapply(fastaSplitted, function(Y) {
          length(which(Y[[1]] %in% X))
        })
      }), recursive = FALSE)

      g3 <- unlist(lapply(group3, function(X) {
        lapply(fastaSplitted, function(Y) {
          length(which(Y[[1]] %in% X))
        })
      }), recursive = FALSE)


      g1 <- .cld(g1, "G1", nchar_df)
      g2 <- .cld(g2, "G2", nchar_df)
      g3 <- .cld(g3, "G3", nchar_df)


      dfOut <-
        rbind(g1, g2, g3) %>% spread(group, value)
      dfOut[is.na(dfOut)] <- 0

      return(dfOut)
    }
