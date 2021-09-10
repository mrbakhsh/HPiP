    #' calculateCTriad
    #' @title Calculate Conjoint Triad Descriptor
    #' @param x A data.frame containing gene/protein names and their
    #' fasta sequences.
    #' @return A length 343  named vector for the data input.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @description This function calculates Conjoint Triad descriptor
    #' for data input.
    #' @export
    #' @references
    #' Shen, J., Zhang, J., Luo, X., Zhu, W., Yu, K., Chen, K., et al.
    #' (2007). Predicting protein–protein interactions based only on
    #' sequences information. \emph{Proc. Natl. Acad. Sci.} 104, 4337–4341.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateCTriad(UP000464024_df)
    #' head(x_df, n = 2L)
    calculateCTriad <- function(x) {
      if (!is.data.frame(x)) x <- is.data.frame(x)

      . <- NULL
      V1 <- NULL

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


      classes <- vector("list", 7)
      classes[[1]] <- c("A", "G", "V")
      classes[[2]] <- c("I", "L", "F", "P")
      classes[[3]] <- c("Y", "M", "T", "S")
      classes[[4]] <- c("H", "N", "Q", "W")
      classes[[5]] <- c("R", "K")
      classes[[6]] <- c("D", "E")
      classes[[7]] <- c("C")


      vspace <- expand.grid(seq_len(7), seq_len(7), seq_len(7))

      CTDict <- vector("list", 343)

      for (i in seq_len(343)) {
        tmp <- as.vector(outer(classes[[vspace[i, 1]]],
                               classes[[vspace[i, 2]]], paste,
                               sep = ""
        ))
        CTDict[[i]] <- as.vector(outer(tmp,
                                       classes[[vspace[i, 3]]], paste,
                                       sep = ""
        ))
      }
      names(CTDict) <- paste(seq_along(CTDict))



      # compute TC_Sm
      p <- function(x, CTDict) {
        f <-
          x %>%
          lapply(., function(x) strsplit(x[[1]], "")) %>%
          lapply(., function(x) {
            paste(paste(x[[1]][-c(length(x[[1]]), length(x[[1]]) - 1)],
                        x[[1]][-c(1, length(x[[1]]))],
                        sep = ""
            ), x[[1]][-c(1, 2)],
            sep = ""
            )
          })

        inter <- unlist(lapply(CTDict, function(X) {
          lapply(f, function(Y) {
            length(which(Y %in% X))
          })
        }), recursive = FALSE)

        s <-
          do.call(cbind, inter) %>%
          set_rownames(names(x))
        s <- t(apply(s, 1, function(x) (x - min(x)) / max(x)))
        return(s)
      }

      CT <- lapply(fastalist, function(x) p(x[[1]], CTDict))
      CT_df <-
        do.call(rbind, CT)

      CT_df <-
        CT_df %>%
        as.data.frame(.) %>%
        set_rownames(names(CT)) %>%
        rownames_to_column("identifier")

      return(CT_df)
    }
