    .gintersect <- function(x, y, gn = "G1") {
      g1 <- unlist(lapply(x, function(X) {
        lapply(y, function(Y) {
          ifelse(Y[[1]] %in% X, gn, NA)
        })
      }), recursive = FALSE)
    }


    .mapL <- function(x, y, z) {
      list.g1g2 <- Map(list, x, y)

      list.g1g2 <-
        lapply(list.g1g2, function(x) {
          ifelse(is.na(x[[1]]), x[[2]], x[[1]])
        })

      list.g1g2g3 <-
        Map(list, list.g1g2, z)

      list.g1g2g3 <-
        lapply(list.g1g2g3, function(x) {
          ifelse(is.na(x[[1]]), x[[2]], x[[1]])
        })

      return(list.g1g2g3)
    }


    #' calculateCTDT
    #' @title Calculate CTD Descriptors - Transition (T)
    #' @param x A data.frame containing gene/protein names and
    #' their fasta sequences.
    #' @return A length 21 named vector for the data input.
    #' @seealso See \code{\link{calculateCTDC}} and \code{\link{calculateCTDD}}
    #' for Composition and Distribution descriptors.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom dplyr full_join
    #' @importFrom dplyr select
    #' @description This function calculates Transition (T) descriptor
    #' for data input.
    #' @export
    #' @references
    #' Dubchak, I., Muchnik, I., Holbrook, S. R., and Kim, S.-H. (1995).
    #' Prediction of protein folding class using global description of
    #' amino acid sequence.
    #' \emph{Proc. Natl. Acad. Sci.} 92, 8700–8704.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateCTDT(UP000464024_df)
    #' head(x_df, n = 2L)
    calculateCTDT <- function(x) {
      . <- NULL
      V1 <- NULL
      ID <- NULL
      group <- NULL
      Tr1221.prob <- NULL
      Tr1331.prop <- NULL
      Tr2332.prob <- NULL
      Tclass <- NULL
      count_CTDT <- NULL




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


      vect1 <- c("G1", "G2", "G3")
      vect2 <- c("G1", "G2", "G3")
      TClist <-
        apply(expand.grid(vect1, vect2), 1, paste, collapse = "")


      nchar_df <- .nchar(fastalist)


      fastaSplitted <-
        fastalist %>% lapply(., function(x) strsplit(x[[1]], ""))

      g1 <- .gintersect(group1, fastaSplitted, gn = "G1")
      g2 <- .gintersect(group2, fastaSplitted, gn = "G2")
      g3 <- .gintersect(group3, fastaSplitted, gn = "G3")
      list.g1g2g3 <- .mapL(g1, g2, g3)


      CTDT_calc <-
        list.g1g2g3 %>%
        lapply(., function(x) paste(x[-length(x)], x[-1], sep = "")) %>%
        lapply(., function(x) factor(x, levels = TClist)) %>%
        lapply(., function(x) table(x))

      G1G2.G2G1 <-
        CTDT_calc %>%
        lapply(., function(x) sum(x[c("G1G2", "G2G1")])) %>%
        do.call(rbind, .) %>%
        as.data.frame(.) %>%
        rename(G1G2.G2G1 = V1) %>%
        rownames_to_column("ID")

      G1G3.G3G1 <-
        CTDT_calc %>%
        lapply(., function(x) sum(x[c("G1G3", "G3G1")])) %>%
        do.call(rbind, .) %>%
        as.data.frame(.) %>%
        rename(G1G3.G3G1 = V1) %>%
        rownames_to_column("ID")


      G2G3.G3G2 <-
        CTDT_calc %>%
        lapply(., function(x) sum(x[c("G2G3", "G3G2")])) %>%
        do.call(rbind, .) %>%
        as.data.frame(.) %>%
        rename(G2G3.G3G2 = V1) %>%
        rownames_to_column("ID")

      dfOut <-
        full_join(G1G2.G2G1, G1G3.G3G1, by = "ID") %>%
        full_join(., G2G3.G3G2, by = "ID") %>%
        separate(ID, c("group", "identifier"), sep = "\\.") %>%
        left_join(., nchar_df, by = "identifier") %>%
        mutate(nchar = nchar - 1) %>%
        mutate(Tr1221.prob = G1G2.G2G1 / nchar) %>%
        mutate(Tr1331.prop = G1G3.G3G1 / nchar) %>%
        mutate(Tr2332.prob = G2G3.G3G2 / nchar) %>%
        select(2, 1, 7:9) %>%
        gather("Tclass", "count_CTDT", 3:ncol(.)) %>%
        mutate(Tclass = paste(group, Tclass, sep = ".")) %>%
        select(-2) %>%
        spread(Tclass, count_CTDT)

      return(dfOut)
    }
