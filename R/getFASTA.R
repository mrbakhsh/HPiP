    #' getFASTA
    #' @title Fetch FASTA Sequence from the UniProt Database
    #' @param uniprot.id A character vector of UniProt identifiers.
    #' @param filename A character string, indicating the output filename as an
    #' RData object to store the retrieved sequences.
    #' @param path A character string indicating the path to the project
    #' directory that contains the interaction data.
    #' If the directory is missing, it will be stored in the current directory.
    #' Default is FASTASeq.
    #' @return A list containing protein FASTA sequences.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom stats na.omit
    #' @importFrom utils read.delim
    #' @description This function retrieves protein sequences in FASTA format
    #' directly from the UniProt database via UniProt protein IDs. This function
    #' also checks if the amino-acid composition of protein sequences
    #' is in the 20 default types.
    #' @export
    #' @examples
    #' # get fasta sequences for three proteins of SARS-Cov-2
    #' local = tempdir()
    #' uniprot.id <- c("P0DTC4", "P0DTC5", "P0DTC9")
    #' fasta_df <- getFASTA(uniprot.id, filename = 'FASTA.RData', path = local)
    #' head(fasta_df)
    
    getFASTA <- function(uniprot.id, filename = "FASTA.RData",
                         path = "FASTASeq") {
      . <- NULL
      V1 <- NULL
      
      `%nin%` <- Negate(`%in%`)
      
      names <- unique(uniprot.id)
      if (!is.character(names)) names <- as.character(names)
      len <- length(names)
      fastalist.query <- vector("list", len)
      
      for (i in seq_len(len)) {
        fastalist.query[[i]] <-
          tryCatch(suppressWarnings(protr::readFASTA(paste(
            "https://www.uniprot.org/uniprot/", names[i], ".fasta",
            sep = ""
          )))[[1]], error = function(err) {
            cat("Error detected --> NA substituted\n")
            return(NA)
          })
      }
      
      names(fastalist.query) <- names
      
      # check if protein sequence's amino acid types re in the 20 default types
      s1.check <-
        lapply(fastalist.query, function(x) protr::protcheck(x))
      s1.check <-
        do.call(rbind, s1.check) %>%
        as.data.frame(.) %>%
        rownames_to_column("Uniprot") %>%
        filter(V1 == FALSE) %>%
        .$Uniprot
      fastalist.query <-
        fastalist.query[names(fastalist.query) %nin% s1.check]
      
      
      fastalist.query.out <-
        do.call(rbind, fastalist.query) %>%
        as.data.frame(.) %>%
        rownames_to_column("UniprotKBID") %>%
        rename(FASTASEQ = V1) %>%
        na.omit(.)
      
      if (!file.exists(path)) {
        dir.create(path, recursive = TRUE)
        localPath <- tools::file_path_as_absolute(path)
        fname <- file.path(localPath, filename)
        save(fastalist.query.out, file = fname)
      } else {
        fname <- file.path(path, filename)
        save(fastalist.query.out, file = fname)
      }
      
      return(fastalist.query)
    }
