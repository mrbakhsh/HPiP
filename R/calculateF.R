    .F1_calc <-
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


    .F2_calc <-
      function(x) {
        . <- NULL
        V1 <- NULL
        V2 <- NULL
        AAcid <- NULL

        f <-
          x %>%
          lapply(., function(x) strsplit(x[[1]], "")) %>%
          lapply(., function(x) {
            paste(paste(x[[1]][-c(length(x[[1]]), length(x[[1]]) - 1)],
              x[[1]][-c(1, length(x[[1]]))], x[[1]][-c(
                seq_len(2),
                length(x[[1]])
              )],
              x[[1]][-c(seq_len(3), length(x[[1]]))],
              x[[1]][-c(seq_len(4), length(x[[1]]))],
              x[[1]][-c(seq_len(5), length(x[[1]]))],
              sep = ""
            ), x[[1]][-c(1, 2, 3, 4, 5, 6)],
            sep = ""
            )
          })


        x <- unlist(f)

        x <- paste(x, collapse = " ")
        x <- unlist(strsplit(x, split = ""))


        rl <- rle(x)
        lst <- lapply(split(
          seq_along(x),
          rep(seq_along(rl$values), rl$lengths)
        ), range)

        names(lst) <- rl$values
        rl <-
          do.call(rbind, lst) %>%
          as.data.frame(.) %>%
          tibble::rownames_to_column("AAcid") %>%
          mutate(AAcid = str_remove_all(AAcid, ".\\d+")) %>%
          mutate(length = ((V2 - V1) + 1)^2) %>%
          group_by(AAcid) %>%
          summarise(Freq = max(length)) %>%
          filter(AAcid != "X.")
        return(rl)
      }


    #' calculateF
    #' @title Calculate F1 or F2 Descriptors
    #' @param x A data.frame containing gene/protein names and their
    #' fasta sequences.
    #' @param type The descriptor type:
    #' \code{F1} or \code{F2}.
    #' @return A length 20 named vector for the data input.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom dplyr mutate
    #' @importFrom magrittr %>%
    #' @importFrom dplyr group_by
    #' @importFrom dplyr summarise
    #' @importFrom tibble rownames_to_column
    #' @importFrom tidyr separate
    #' @importFrom stringr str_remove_all
    #' @importFrom tidyr spread
    #' @description This function calculates F1 or F2 descriptors:
    #' \itemize{
    #' \item \code{F1} - sum of squared length of Single Amino Acid Repeats
    #' (SARs) in the entire protein sequence.
    #' \item \code{F2} - maximum of the sum of Single Amino Acid Repeats (SARs)
    #' in a window of 6 residues.
    #' }
    #' @export calculateF
    #' @references
    #' Alguwaizani, S., Park, B., Zhou, X., Huang, D.-S., and Han, K. (2018).
    #' Predicting interactions between virus and host proteins using repeat
    #' patterns and composition of amino acids.
    #' \emph{J. Healthc. Eng.} 2018.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateF(UP000464024_df, type = "F1")
    #' head(x_df, n = 2L)
    calculateF <- function(x, type = c("F1", "F2")) {
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
      f <- .checkFASTA(fastalist)
      if (nrow(f) > 0) {
        stop("Fastalist has unrecognized amino acid type")
      }

      if (type == "F1") {
        p <-
          lapply(fastalist, function(x) .F1_calc(x[[1]]))
      }

      if (type == "F2") {
        p <-
          lapply(fastalist, function(x) .F2_calc(x[[1]]))
      }

      F_calc_df <-
        bind_rows(p, .id = "identifier") %>%
        spread(AAcid, Freq) %>%
        replace(is.na(.), 0)

      return(F_calc_df)
    }
