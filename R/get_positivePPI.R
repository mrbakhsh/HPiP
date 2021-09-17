    #' get_positivePPI
    #' @title Fetch Positive Reference Host-Pathogen Protein-Protein
    #' Interactions (HP-PPIs) from the BioGRID Database
    #' @param organism.taxID Taxonomy identifier for the pathogen.
    #' @param access.key Access key for using BioGRID webpage. To retrieve
    #' interactions from the BioGRID database, the users are first required to
    #' register for access key at
    #' \url{https://webservice.thebiogrid.org/}.
    #' @param filename A character string, indicating the output filename as an
    #' RData object to store the retrieved interactions.
    #' @param path A character string indicating the path to the project
    #' directory that contains the interaction data. If the directory is
    #' missing, it will be stored in the current directory.
    #' Default is PositiveInt.
    #' @return A Data.frame containing true positive protein-protein
    #' interactions for the selected pathogen.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @seealso See \code{\link{get_negativePPI}} for generating negative
    #' protein-protein interaction.
    #' @importFrom readr cols
    #' @importFrom readr col_double
    #' @importFrom readr col_character
    #' @description This function retrieves positive reference host-pathogen
    #' protein-protein interactions directly from BioGRID database.
    #' @export
    #' @examples
    #' \donttest{
    #' local = tempdir()
    #' try(get_positivePPI(organism.taxID = 2697049,
    #' access.key = 'XXXX',
    #' filename = "PositiveInt.RData",
    #' path = local))
    #' }


    get_positivePPI <- function(organism.taxID = 2697049,
                                access.key = "81bb3b5a6bd9a8084a7be71f0963ab1e",
                                filename = "PositiveInt.RData",
                                path = "PositiveInt") {
      `Official Symbol Interactor A` <- NULL
      `Official Symbol Interactor B` <- NULL

      # construct query to retrive all interactiosn for the requested organism
      url <-
        "http://webservice.thebiogrid.org/interactions/?"
      query <- paste(url,
        paste0("accessKey=", access.key),
        "interSpeciesExcluded=false",
        "selfInteractionsExcluded=true",
        "includeHeader=true",
        paste0("taxId=", organism.taxID),
        sep = "&"
      )


      # number of results
      n <-
        httr::GET(paste(query, "format=count", sep = "&"))
      httr::stop_for_status(n)


      n <- as.numeric(
        httr::content(n, type = "text/csv", encoding = "UTF-8", col_types = "i")
      )

      # retrieve results in batches
      tpINT <- lapply(seq(1, n, 10000), function(m) {
        l <- paste(query, paste0("start=", m), sep = "&")
        request <-
          httr::GET(l)
        httr::stop_for_status(request)

        # parse the biogird results
        biog_result <- httr::content(request,
          type = "text/tab-separated-values",
          encoding = "UTF-8",
          col_types = "cccccccccccccccccccccccc"
        )
        q <- subset(biog_result, select = c(
          "Official Symbol Interactor A",
          "Official Symbol Interactor B",
          "Organism Interactor A",
          "Organism Interactor B",
          "Experimental System Type"
        ))

        return(q)
      })

      tpINT <- do.call("rbind", tpINT)
      dataTP <-
        tpINT
      dataTP <- # remove duplicate
        tpINT[!duplicated(apply(tpINT, 1, function(x) {
          paste(sort(x), collapse = "")
        })), ]

      dataTP <-
        dataTP %>%
        mutate(PPI = paste(`Official Symbol Interactor A`,
          `Official Symbol Interactor B`,
          sep = "~"
        )) %>%
        dplyr::select(6, seq_len(5))

      if (!file.exists(path)) {
        dir.create(path, recursive = TRUE)
        localPath <- tools::file_path_as_absolute(path)
        fname <- file.path(localPath, filename)
        save(dataTP, file = fname)
      } else {
        fname <- file.path(path, filename)
        save(dataTP, file = fname)
      }

      return(dataTP)
    }
