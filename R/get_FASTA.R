#' get_FASTA
#' @title Retrieve FASTA Sequences from UniProt Databse
#' @param custombg A character vector of UniprotKB identifiers under study.
#' @param proteomeIDs A proteome identifier for the organism under study.The full list of
#' proteome IDs are available at  \url{https://www.uniprot.org/proteomes/}.
#' @return A Data.frames containing protein FASTA sequences from the UniProt database.
#' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom utils read.delim
#' @description This function retrieves FASTA sequences for the organism of interest
#' from the UniProt database and checks if the amino-acid composition of protein sequences
#' is in the 20 default types.
#' @export
#' @examples
#' # get fasta seqeunces for three proteins of SARS-Cov-2
#' x <- c("P0DTC4", "P0DTC5", "P0DTC9")
#' fasta_df <- get_FASTA(custombg = x, proteomeIDs = "UP000464024")
#' head(fasta_df)
get_FASTA <-
    function(custombg = NULL, proteomeIDs = "UP000464024") {
        Entry <- NULL
        . <- NULL
        V1 <- NULL


        `%nin%` <- Negate(`%in%`)


        ## retrive proteome information of desired organism
        url <- "http://www.uniprot.org/uniprot/?query=proteome:"
        idStr <- paste(proteomeIDs, collapse = "+or+")
        reviewed <- "%20reviewed:yes"
        format <- "&format=tab"
        query <- read.delim(paste0(url, idStr, reviewed, format))



        if (is.character(custombg) == TRUE) {
            data_filt <-
                query %>%
                filter(Entry %in% custombg) %>%
                dplyr::select(1)

            # extract fasta sequence for query proteins
            names <- data_filt$Entry
            len <- nrow(data_filt)
            fastalist.query <-
                vector(mode = "list", length = len)

            for (i in seq_len(len)) {
                fastalist.query[[i]] <-
                    protr::readFASTA(paste("https://www.uniprot.org/uniprot/",
                        data_filt[i, ], ".fasta",
                        sep = ""
                    ))[[1]]
            }

            names(fastalist.query) <-
                names
            # check if protein sequence's amino acid types re in the 20 default types
            s1.check <-
                lapply(fastalist.query, function(x) protr::protcheck(x))
            s1.check <-
                do.call(rbind, s1.check) %>%
                as.data.frame(.) %>%
                rownames_to_column("Uniprot") %>%
                filter(V1 == FALSE) %>%
                .$Uniprot
            fastalist.query.out <-
                fastalist.query[names(fastalist.query) %nin% s1.check]
            fastalist.query.out <-
                do.call(rbind, fastalist.query.out) %>%
                as.data.frame(.) %>%
                tibble::rownames_to_column("UniprotKBID") %>%
                rename(FASTASEQ = V1)
        }


        if (is.character(custombg) == FALSE) {
            df_unfilt <-
                query %>%
                select(1)

            # extract fasta sequence for query proteins
            names <- df_unfilt$Entry
            len <- nrow(df_unfilt)
            fastalist.query <-
                vector(mode = "list", length = len)

            for (i in seq_len(len)) {
                fastalist.query[[i]] <-
                    protr::readFASTA(paste("https://www.uniprot.org/uniprot/",
                        df_unfilt[i, ], ".fasta",
                        sep = ""
                    ))[[1]]
            }

            names(fastalist.query) <-
                names
            # check if protein sequence's amino acid types re in the 20 default types
            s1.check <-
                lapply(fastalist.query, function(x) protr::protcheck(x))
            s1.check <-
                do.call(rbind, s1.check) %>%
                as.data.frame(.) %>%
                rownames_to_column("Uniprot") %>%
                filter(V1 == FALSE) %>%
                .$Uniprot
            fastalist.query.out <-
                fastalist.query[names(fastalist.query) %nin% s1.check]
            fastalist.query.out <-
                do.call(rbind, fastalist.query.out) %>%
                as.data.frame(.) %>%
                rownames_to_column("UniprotKBID") %>%
                rename(FASTASEQ = V1)
        }


        return(fastalist.query.out)
    }
