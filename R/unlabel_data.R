    #' HP-PPIs with Unknown Class Labels
    #' @description This dataset consists of interactions between
    #' SARS-CoV-2 and human proteins, achieved by AP-MS
    #' (affinity purification mass spectrometry).
    #' @details
    #' To construct this dataset, data (supplementary table 1) containing
    #' SARS-CoV-2-human PPIs was retrieved from (Gordon et al., 2020) and 700
    #' pairs were randomly selected from total pairs, followed by converting
    #' protein sequences of host or viral proteins to numerical features and
    #' finally concatenating the computed features in order to construct
    #' host-pathogen PPIs.
    #' @name unlabel_data
    #' @usage data(unlabel_data)
    #' @docType data
    #' @format A data.frame containing 700 SARS-CoV-2-Human protein-protein
    #' interactions (PPIs) with pre-computed numerical features using
    #' CTD (composition/transition/distribution) descriptors.
    #' @source \url{https://www.nature.com/articles/s41586-020-2286-9#Sec36}
    #' @references
    #' Gordon,D.E. et al. (2020) A SARS-CoV-2 protein interaction map reveals
    #' targets for drug repurposing. Nature, 583, 459–468.
    NULL
