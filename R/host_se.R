    #' Host SummarizedExperiment object
    #' @description SummarizedExperiment object of numerical features for host
    #' proteins.
    #' @details
    #' To construct this object, first protein sequences were converted to
    #' numerical features using (CTD) descriptors provided in the HPiP package,
    #' followed by converting each numerical features matrix to
    #' SummarizedExperiment object.Each object is then merged into one object
    #' using `cbind()`.
    #' @name host_se
    #' @usage data(host_se)
    #' @docType data
    NULL
