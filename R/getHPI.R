    .getKronecker <- function(A, B) {
      m <- matrix(0, nrow = nrow(A) * nrow(B), ncol = ncol(A) * ncol(B))

      for (i in seq_len(nrow(A))) {
        cur_rows <- seq_len(nrow(B)) + (i - 1) * nrow(B)
        for (j in seq_len(ncol(A))) {
          m[cur_rows, seq_len(ncol(B)) + (j - 1) * ncol(B)] <- A[i, j] * B
        }
      }
      m
    }

    #' getHPI
    #' @title Generating Host-Pathogen Protein-Protein Interaction (HP-PPI)
    #' Descriptors
    #' @param pathogenData The pathogen descriptor matrix.
    #' @param hostData The host descriptor matrix.
    #' @param type The interaction type, one or two of
    #' \code{"combine"} and \code{"kron.prod"}.
    #' @return A matrix containing the Host-Pathogen Protein-Protein Interaction
    #' (HP-PPI) descriptors.
    #' @author Matineh Rahmatbakhsh \email{matinerb.94@gmail.com}
    #' @description This function calculates Host-Pathogen Protein-Protein
    #' Interaction (HP-PPI) descriptors via two approaches
    #' \itemize{
    #' \item \code{combine} - combine the two descriptor matrix,
    #' result has \code{(p1 + p2)} columns
    #' \item \code{kron.prod} - if A has m x n matrix and B is q x p matrix,
    #' then the Kronecker product is the code{(pm × qn)} block matrix
    #' }
    #' @export
    #' @examples
    #' x <- matrix(c(1, 2, 3, 1), nrow = 2, ncol = 2, byrow = TRUE)
    #' y <- matrix(c(0, 3, 2, 1), nrow = 2, ncol = 2, byrow = TRUE)
    #' getHPI(x, y, "combine")
    #' getHPI(x, y, "kron.prod")
    getHPI <- function(pathogenData, hostData, type = c("combine", "kron.prod"))
      {
      if (!is.matrix(pathogenData)) pathogenData <- as.matrix(pathogenData)
      if (!is.matrix(hostData)) hostData <- as.matrix(hostData)

      pathogenRow <- nrow(pathogenData)
      hostRow <- nrow(hostData)

      if (pathogenRow != hostRow) stop("Matrix row count must match")
      if (missing(type)) stop("Must provide at least one interaction type")


      if (type == "combine") {
        result <- cbind(pathogenData, hostData)
      } else if (type == "kron.prod") {
        result <- .getKronecker(pathogenData, hostData)
      } else {
        stop("Interaction type must be  kron or combine")
      }


      return(result)
    }
