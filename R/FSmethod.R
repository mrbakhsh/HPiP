    .filter.corr <-
      function(features, class, cor.cutoff = 0.7) {
        . <- NULL
        x <- NULL
        if (!is.matrix(features)) features <- as.matrix(features)
        descrCorr <- cor(features)

        # find correlated features
        highCorr <-
          caret::findCorrelation(descrCorr, cor.cutoff, names = TRUE)

        # find the remining features
        clcol <- which(names(as.data.frame(features)) %in% highCorr)
        opt_var <- colnames(features[, -clcol])


        # remove correlated features
        reduced_Data <-
          features %>%
          as.data.frame(.) %>%
          select(-all_of(highCorr)) %>%
          rownames_to_column("PPI") %>%
          cbind(class, .)

        ret <- list()
        ret$corProfile <- descrCorr
        ret$corSelectedFeatures <- opt_var
        ret$cordf <- reduced_Data

        return(ret) # return RF profile
      }

    .rfeFS <- function(features, class, cor.cutoff = 0.7,
                       resampling.method = "cv", iter = 2, repeats = 2,
                       metric = "Accuracy", verbose = TRUE) {
      if (!is.data.frame(x)) x <- as.data.frame(x)

      . <- NULL

      caretFuncs <- rfFuncs
      if (metric == "ROC") {
        caretFuncs$summary <- twoClassSummary
      }

      sizeF <- ncol(features)
      #### recursive feature elimination-random forest
      rfProfile <- caret::rfe(
        x = features,
        y = as.factor(class),
        sizes = seq_len(sizeF), # all possible features to be tested
        metric = metric,
        rfeControl = rfeControl(
          functions = caretFuncs,
          method = resampling.method,
          repeats = repeats,
          number = iter,
          allowParallel = TRUE,
          verbose = verbose
        )
      )

      opt_var <- predictors(rfProfile)
      dataset_rfe <- features[, c(opt_var)]
      dataset_rfe <-
        dataset_rfe %>%
        as.data.frame(.) %>%
        tibble::rownames_to_column("PPI")
      dataset_rfe <- cbind(class, dataset_rfe)

      ret <- list()
      ret$rfProfile <- rfProfile
      ret$rfSelectedFeatures <- opt_var
      ret$rfdf <- dataset_rfe

      return(ret) # return RF profile
    }

    #' FSmethod
    #' @title Feature Selection via Matrix Correlation and
    #' Recursive Feature Elimination (RFE)
    #' @param x A data.frame containing protein-protein interactions,
    #' class labels and features.
    #' @param type The feature selection type, one or two of
    #' \code{filter.corr} and \code{rfeFS}.
    #' @param cor.cutoff Correlation coefficient cutoff used for filtering.
    #' See \code{filter.corr} for more details.
    #' @param resampling.method The resampling method for RFE :'boot',
    #' 'boot632', optimism_boot',boot_all', 'cv', 'repeatedcv', 'LOOCV',
    #' 'LGOCV';defaults to cv. See \code{rfeFS} and
    #' \code{\link[caret]{rfeControl}} for more details.
    #' @param iter Number of partitions for cross-validation;
    #' defaults to 2. See \code{rfeFS} and \code{\link[caret]{rfeControl}}
    #' for more details.
    #' @param repeats For repeated k-fold cross validation only;
    #'  defaults to 3.See \code{rfeFS} and \code{\link[caret]{rfeControl}}
    #'  for more details.
    #' @param metric  A string that specifies what summary metric will be used
    #' to select the optimal feature ; default to ROC.See \code{rfeFS} and
    #' \code{\link[caret]{rfe}} for more details.
    #' @param verbose Make the output verbose.See \code{rfeFS} and
    #' \code{\link[caret]{rfeControl}} for more details.

    #' @return
    #' If the type set to \code{filter.corr} , the output includes
    #' the following elements:
    #' \itemize{
    #' \item{corProfile} - A correlation matrix.
    #' \item{corSelectedFeatures} - Name of features that retained after
    #' the correlation analysis.
    #' \item{cordf} - A data.frame filtered.
    #' }
    #' If the type set to \code{rfeFS} , the output includes the following
    #' elements:
    #' \itemize{
    #' \item{rfProfile} - A list of elements. See \code{\link[caret]{rfe}} for
    #' more details.
    #' \item{rfSelectedFeatures} - Name of features that retained in the feature
    #' selection process.
    #' \item{rfdf} - A data.frame filtered.
    #' }
    #' If type set to \code{both} the output includes the following elements:
    #' \itemize{
    #' \item{rfdf} - The final data.frame that includes the selected features
    #' retained after both \code{filter.corr} and \code{rfeFS} analysis.
    #' }

    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}.
    #' @importFrom stats cor
    #' @importFrom caret twoClassSummary
    #' @importFrom dplyr select
    #' @importFrom dplyr all_of
    #' @importFrom tibble rownames_to_column
    #' @importFrom caret rfeControl
    #' @importFrom caret predictors
    #' @importFrom caret rfFuncs
    #' @description This function performs feature selections
    #' via two approaches
    #' \itemize{
    #' \item \code{filter.corr} - compute matrix correlation between features
    #' and filter using a threshold.
    #' \item \code{rfeFS} - perform recursive feature elimination (RFE) method
    #' wrapped with a Random Forest (RF) algorithm for feature importance
    #' evaluation.
    #' }
    #' @export
    #' @examples
    #' data('example_data')
    #' x <- na.omit(example_data)
    #' s <- FSmethod(x, type = 'both',
    #' cor.cutoff = 0.7, resampling.method = "repeatedcv",
    #' iter = 5, repeats = 3, metric = "ROC", verbose = TRUE)


    FSmethod <- function(x, type = c("cor", "rfe", "both"),
                         cor.cutoff = 0.7, resampling.method = "cv",
                         iter = 2, repeats = 3, metric = "Accuracy",
                         verbose = TRUE) {
      PPI <- NULL
      . <- NULL


      if (!is.data.frame(x)) x <- as.data.frame(x)

      if (all(colnames(x) != "class") == TRUE) {
        stop("class attribute is absent from the data.frame")
      }

      if (all(colnames(x) != "PPI") == TRUE) {
        stop("PPI attribute is absent from the data.frame")
      }

      if (missing(type)) stop("Must provide at least one feature selection
                          approach")

      if (!type %in% c("cor", "rfe", "both")) {
        stop("Feature selection type must be one of cor,rfe or both...")
      }
      features <-
        x %>%
        as.data.frame(.) %>%
        magrittr::set_rownames(.$PPI) %>%
        dplyr::select(-PPI, -class)

      class <-
        x$class

      if (type == "cor") {
        result <- .filter.corr(features, class, cor.cutoff)
      } else if (type == "rfe") {
        result <- .rfeFS(features, class,
          resampling.method = resampling.method,
          iter = iter, repeats = repeats,
          metric = metric, verbose = verbose
        )
      } else if (type == "both") {
        cor.result <-
          .filter.corr(features, class, cor.cutoff)
        df <-
          cor.result[["cordf"]]
        df.features <-
          df %>%
          as.data.frame(.) %>%
          magrittr::set_rownames(.$PPI) %>%
          dplyr::select(-PPI, -class)

        df.class <-
          df$class

        rf.result <-
          .rfeFS(df.features, df.class,
            resampling.method = resampling.method,
            iter = iter, repeats = repeats,
            metric = metric, verbose = verbose
          )

        result <- list()
        result$cor.result <- cor.result
        result$rf.result <- rf.result
      } else {
        stop('Feature selection type must be in "cor" and "rfe" or "both"')
      }

      return(result)
    }
