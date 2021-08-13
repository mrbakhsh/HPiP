#' calculateMoran
#' @title Caclulate Moran  Autocorrelation Descriptor
#' @param x A data.frame containing gene/protein names and their fasta sequences.
#' @param target.props A character vector, specifying the
#' Accession Number of the target properties.
#' #'8 properties are used by default, as listed below:
#'              \describe{
#'              \item{AccNo. CIDH920105}{Normalized average hydrophobicity scales (Cid et al., 1992)}
#'              \item{AccNo. BHAR880101}{Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)}
#'              \item{AccNo. CHAM820101}{Polarizability parameter (Charton-Charton, 1982)}
#'              \item{AccNo. CHAM820102}{Free energy of solution in water, kcal/mole (Charton-Charton, 1982)}
#'              \item{AccNo. CHOC760101}{Residue accessible surface area in tripeptide (Chothia, 1976)}
#'              \item{AccNo. BIGC670101}{Residue volume (Bigelow, 1967)}
#'              \item{AccNo. CHAM810101}{Steric parameter (Charton, 1981)}
#'              \item{AccNo. DAYM780201}{Relative mutability (Dayhoff et al., 1978b)}}
#' @param nlag Maximum value of the lag parameter. Default is \code{30}.
#' @return A length \code{nlag} named vector for data input.
#' @seealso See \code{\link{calculateMoreauBroto}} and
#' \code{\link{calculateGeary}} for Normalized Moreau-Broto and Geary autocorrelation descriptors.
#'
#' @author Matineh Rahmatbakhsh <\email{matinerb.94@gmail.com}>,
#'         Nan Xiao
#' @importFrom stats sd
#' @importFrom utils read.csv
#' @description This function calculates Moran  autocorrelation Descriptor for data input.
#' (Dim: \code{length(target.props) * nlag})
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
#' x <- readRDS(
#'     system.file("protseq/UP000464024_FASTAList.rds", package = "HPiP")
#' )
#' x_df <- calculateMoran(x)
#' head(x_df, n = 2L)
calculateMoran <-
    function(x,
             target.props =
                 c(
                     "CIDH920105", "BHAR880101", "CHAM820101",
                     "CHAM820102", "CHOC760101", "BIGC670101",
                     "CHAM810101", "DAYM780201"
                 ),
             nlag = 30L) {
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
        fastalist.check <-
            lapply(fastalist, function(x) protr::protcheck(x))
        fastalist.check <-
            do.call(rbind, fastalist.check) %>%
            as.data.frame(.) %>%
            tibble::rownames_to_column("ID") %>%
            filter(V1 == FALSE)

        if (nrow(fastalist.check) > 0) {
            stop("Fastalist has unrecognized amino acid type")
        }


        ## this has to change later
        AAidx <- read.csv(
            system.file("sysdata/AAidx.csv", package = "HPiP"),
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
            apply(aaidx[target.props, ], 1, sd) * sqrt((20 - 1) / 20) # sd() uses (n-1)

        Pr <- data.frame(matrix(ncol = 20, nrow = n.props))
        for (i in seq_len(n.props)) Pr[i, ] <- (aaidx[target.props[i], ] - pmean[i]) / psd[i]


        AADict <- c(
            "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
            "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
        )
        names(Pr) <- AADict


        # Compute Moran Autocorrelation Descriptor
        Moran_Autocor <- function(x) {
            fastaSplitted <- strsplit(fastalist$P0DTC2, split = "")[[1]]

            if (nchar(x) <= nlag) {
                warning("length of the amino acid sequence is <= nlag; NAs will be generated")
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
                        (N / (N - j)) * ((sum((L[[i]][1:(N - j)] - Pbar[i]) * (L[[i]][(1:(N - j)) + j] - Pbar[i]))) / (sum((L[[i]] - Pbar[i])^2))),
                        NA
                    )
                    cMoran <- Moran
                    cMoran <- unlist(cMoran)
                }
            }
            return(cMoran)
        }


        # apply the function to multiple lists
        Moran_lists <-
            lapply(fastalist, function(x) Moran_Autocor(x[[1]]))

        # generate df from list
        tbl <- do.call(rbind, Moran_lists)

        # change the colnames
        colnames(tbl) <- as.vector(t(outer(target.props,
            paste(".lag", seq_len(nlag), sep = ""),
            paste,
            sep = ""
        )))

        tbl <-
            tbl %>%
            as.data.frame(.) %>%
            tibble::rownames_to_column("identifier")


        return(tbl)
    }
