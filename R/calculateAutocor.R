    .Moran_Autocor <- function(x, AADict, nlag, n.props, Pr) {
      fastaSplitted <- strsplit(x, split = "")[[1]]

      if (nchar(x) <= nlag) {
        warning("length of the amino acid sequence is <= nlag;
                    NAs will be generated")
      }

      list <- fastaSplitted
      L <- vector("list", n.props)
      for (i in seq_len(n.props)) L[[i]] <- list

      for (i in seq_len(n.props)) {
        for (j in AADict) {
          L[[i]][which(L[[i]] == j)] <- Pr[i, j]
        }
      }
      L <- lapply(L, as.numeric)

      Moran <- vector("list", n.props)
      N <- length(list)
      Pbar <- vapply(L, mean, FUN.VALUE = numeric(1))

      for (i in seq_len(n.props)) {
        for (j in seq_len(nlag)) {
          Moran[[i]][j] <- ifelse(
            N - j > 0,
            (N / (N - j)) * ((sum((L[[i]][seq_len(N - j)] - Pbar[i]) *
              (L[[i]][(seq_len(N - j)) + j] - Pbar[i]))) /
              (sum((L[[i]] - Pbar[i])^2))),
            NA
          )
          cMoran <- Moran
          cMoran <- unlist(cMoran)
        }
      }
      return(cMoran)
    }


    .MB_Autocor <- function(x, AADict, nlag, n.props, Pr) {
      fastaSplitted <- strsplit(x, split = "")[[1]]
      if (nchar(x) <= nlag) {
        warning("length of the amino acid sequence is <= nlag;
                    NAs will be generated")
      }

      list <- fastaSplitted
      L <- vector("list", n.props)
      for (i in seq_len(n.props)) L[[i]] <- list

      for (i in seq_len(n.props)) {
        for (j in AADict) {
          L[[i]][which(L[[i]] == j)] <- Pr[i, j]
        }
      }
      L <- lapply(L, as.numeric)

      MB <- vector("list", n.props)
      N <- length(list)

      for (i in seq_len(n.props)) {
        for (j in seq_len(nlag)) {
          MB[[i]][j] <- ifelse(
            N - j > 0,
            sum(L[[i]][seq_len(N - j)] * L[[i]][(seq_len(N - j)) + j]) /
              (N - j), NA
          )
          cMB <- MB
          cMB <- unlist(cMB)
        }
      }
      return(cMB)
    }


    .Geary_Autocor <- function(x, AADict, nlag, n.props, Pr) {
      fastaSplitted <- strsplit(x, split = "")[[1]]
      if (nchar(x) <= nlag) {
        warning("length of the amino acid sequence is <= nlag;
                    NAs will be generated")
      }


      list <- fastaSplitted
      L <- vector("list", n.props)
      for (i in seq_len(n.props)) L[[i]] <- list

      for (i in seq_len(n.props)) {
        for (j in AADict) {
          L[[i]][which(L[[i]] == j)] <- Pr[i, j]
        }
      }
      L <- lapply(L, as.numeric)

      Geary <- vector("list", n.props)
      N <- length(list)
      Pbar <- vapply(L, mean, FUN.VALUE = numeric(1))


      for (i in seq_len(n.props)) {
        for (j in seq_len(nlag)) {
          Geary[[i]][j] <- ifelse(
            N - j > 0,
            ((N - 1) / (2 * (N - j))) *
              ((sum((L[[i]][seq_len(N - j)] - L[[i]][(seq_len(N - j)) + j])^2))
              / (sum((L[[i]] - Pbar[i])^2))),
            NA
          )
          cGeary <- Geary
          cGeary <- unlist(cGeary)
        }
      }
      return(cGeary)
    }



    #' calculateAutocor
    #' @title Calculate Autocorrelation Descriptors
    #' @param x A data.frame containing gene/protein names and their
    #' fasta sequences.
    #' @param target.props A character vector, specifying the
    #' accession number of the target properties.
    #' 8 properties are used by default, as listed below:
    #' \describe{
    #' \item{AccNo. CIDH920105}{Normalized average hydrophobicity scales
    #'  (Cid et al., 1992)}
    #' \item{AccNo. BHAR880101}{Average flexibility indices
    #' (Bhaskaran-Ponnuswamy, 1988)}
    #' \item{AccNo. CHAM820101}{Polarizability parameter
    #' (Charton-Charton, 1982)}
    #' \item{AccNo. CHAM820102}{Free energy of solution in water, kcal/mole
    #' (Charton-Charton, 1982)}
    #' \item{AccNo. CHOC760101}{Residue accessible surface area in
    #'  tripeptide (Chothia, 1976)}
    #' \item{AccNo. BIGC670101}{Residue volume (Bigelow, 1967)}
    #' \item{AccNo. CHAM810101}{Steric parameter (Charton, 1981)}
    #' \item{AccNo. DAYM780201}{Relative mutability (Dayhoff et al., 1978b)}}
    #' @param nlag Maximum value of the lag parameter. Default is \code{30}.
    #' @return A length \code{nlag} named vector for data input.
    #' @param type The autocorrelation type:
    #' \code{moran}, \code{geary}, or \code{moreaubroto}.
    #'
    #' @author Matineh Rahmatbakhsh <\email{matinerb.94@gmail.com}>,
    #'         Nan Xiao
    #' @importFrom stats sd
    #' @importFrom utils read.csv
    #' @importFrom dplyr select
    #' @description This function calculates autocorrelation descriptors:
    #' \itemize{
    #' \item \code{moran} - moran autocorrelation,
    #' (Dim: \code{length(target.props) * nlag}).
    #' \item \code{geary} - geary autocorrelation,
    #' (Dim: \code{length(target.props) * nlag}).
    #' \item \code{moreaubroto} - moreau-broto autocorrelation,
    #' (Dim: \code{length(target.props) * nlag}).
    #' }
    #' @export
    #' @references
    #' AAindex: Amino acid index database.
    #' \url{http://www.genome.ad.jp/dbget/aaindex.html}
    #'
    #' Feng, Z.P. and Zhang, C.T. (2000)
    #' Prediction of membrane protein types based on the hydrophobic
    #' index of amino acids.
    #' \emph{Journal of Protein Chemistry}, 19, 269-275.
    #'
    #' Horne, D.S. (1988)
    #' Prediction of protein helix content from
    #' an autocorrelation analysis of sequence hydrophobicities.
    #' \emph{Biopolymers}, 27, 451-477.
    #'
    #' Sokal, R.R. and Thomson, B.A. (2006)
    #' Population structure inferred by local spatial autocorrelation:
    #' an Usage from an Amerindian tribal population.
    #' \emph{American Journal of Physical Anthropology}, 129, 121-131.
    #' @examples
    #' data(UP000464024_df)
    #' x_df <- calculateAutocor(UP000464024_df,type = 'moran')
    #' head(x_df, n = 2L)




    calculateAutocor <-
      function(x,
               target.props =
                 c(
                   "CIDH920105", "BHAR880101", "CHAM820101",
                   "CHAM820102", "CHOC760101", "BIGC670101",
                   "CHAM810101", "DAYM780201"
                 ),
               nlag = 30L, type = c("moran", "geary", "moreaubroto")) {
        . <- NULL
        Pr <- NULL
        V1 <- NULL
        n.props <- NULL


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

        ## this has to change later
        AAidx <- read.csv(
          system.file("extdata/AAidx.csv", package = "HPiP"),
          header = TRUE
        )


        # Compute Pr for each type of property
        aaidx <-
          AAidx %>%
          magrittr::set_rownames(.$AccNo) %>%
          select(-1)

        n.props <-
          length(target.props)

        pmean <-
          rowMeans(aaidx[target.props, ])

        psd <-
          apply(aaidx[target.props, ], 1, sd) * sqrt((20 - 1) / 20)

        Pr <- data.frame(matrix(ncol = 20, nrow = n.props))
        for (i in seq_len(n.props)) {
          Pr[i, ] <- (aaidx[target.props[i], ] - pmean[i]) / psd[i]
        }


        AADict <- c(
          "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
          "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
        )
        names(Pr) <- AADict


        if (type == "moran") {

          # apply the function to multiple lists
          autocor_lists <-
            lapply(fastalist, function(x) {
              .Moran_Autocor(x[[1]],
                AADict = AADict,
                nlag = nlag,
                Pr = Pr, n.props = n.props
              )
            })
        }

        if (type == "geary") {

          # apply the function to multiple lists
          autocor_lists <-
            lapply(fastalist, function(x) {
              .Geary_Autocor(x[[1]],
                AADict = AADict,
                nlag = nlag,
                Pr = Pr, n.props = n.props
              )
            })
        }

        if (type == "moreaubroto") {

          # apply the function to multiple lists
          autocor_lists <-
            lapply(fastalist, function(x) {
              .MB_Autocor(x[[1]],
                AADict = AADict,
                nlag = nlag,
                Pr = Pr, n.props = n.props
              )
            })
        }

        # generate df from list
        tbl <- do.call(rbind, autocor_lists)

        # change the colnames
        colnames(tbl) <- as.vector(t(outer(target.props,
          paste(".lag", seq_len(nlag), sep = ""),
          paste,
          sep = ""
        )))

        tbl <-
          tbl %>%
          as.data.frame(.) %>%
          rownames_to_column("identifier")


        return(tbl)
      }
