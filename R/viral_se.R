    #' Viral SummarizedExperiment object
    #' @description SummarizedExperiment object of numerical features for
    #' SARS-CoV-2 proteins.
    #' @details
    #' To construct this object, first protein sequences were converted to
    #' numerical features using (CTD) descriptors provided in the HPiP package,
    #' followed by converting each numerical features matrix to
    #' SummarizedExperiment object.
    #' Each object is then merged into one object using `cbind()`.
    #' @name viral_se
    #' @usage data(viral_se)
    #' @docType data
    NULL
