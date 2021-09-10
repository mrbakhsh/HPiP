    #' calculateCTDD
    #' @title Calculate CTD Descriptors - Distribution (D)
    #' @param x A data.frame containing gene/protein names and their fasta
    #' sequences.
    #' @return A length 105 named vector for the data input.
    #' @seealso See \code{\link{calculateCTDC}} and \code{\link{calculateCTDT}}
    #' for Composition and Transition descriptors.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom dplyr rename
    #' @importFrom dplyr arrange
    #' @importFrom dplyr full_join
    #' @importFrom dplyr row_number
    #' @importFrom dplyr across
    #' @importFrom dplyr everything
    #' @description This function calculates Distribution (D) descriptor
    #' for data input.
    #' @export
    #' @references
    #' Dubchak, I., Muchnik, I., Holbrook, S. R., and Kim, S.-H. (1995).
    #' Prediction of protein folding class using global description of amino
    #' acid sequence.\emph{Proc. Natl. Acad. Sci.} 92, 8700–8704.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateCTDD(UP000464024_df)
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
      f <- .checkFASTA(fastalist)
      if (nrow(f) > 0) {
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


      nchar_df <- .nchar(fastalist)

      fastaSplitted <-
        fastalist %>%
        lapply(., function(x) strsplit(x[[1]], ""))


      g1 <- .gintersect(group1, fastaSplitted, gn = "G1")
      g2 <- .gintersect(group2, fastaSplitted, gn = "G2")
      g3 <- .gintersect(group3, fastaSplitted, gn = "G3")
      list.g1g2g3 <- .mapL(g1, g2, g3)


      D <- vector("list", length(list.g1g2g3))
      names(D) <-
        names(list.g1g2g3)
      for (i in seq_along(list.g1g2g3)) D[[i]] <- matrix(ncol = 5, nrow = 3)

      for (i in seq_along(list.g1g2g3)) {
        for (j in seq_len(3)) {
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
        separate(ID, c("class", "identifier", "g"), "\\.", remove = FALSE) %>%
        rename(
          residue0 = V1, residue25 = V2, residue50 = V3, residue75 = V4,
          residue100 = V5
        ) %>%
        select(-c(2, 4)) %>%
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
        summarise(across(everything(), ~ max(.)))


      return(dfOut)
    }
