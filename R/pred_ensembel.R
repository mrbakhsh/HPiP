      .MajVot <- function(x) {
        tabulatedOutcomes <- table(x)
        sortedOutcomes <- sort(tabulatedOutcomes, decreasing = TRUE)
        mostCommonLabel <- names(sortedOutcomes)[1]
        mostCommonLabel
      }

      #' pred_ensembel
      #' @title Predict Interactions via Ensemble Learning Method
      #' @param features A data frame with host-pathogen protein-protein
      #' interactions (HP-PPIs) in the first column, and features to be passed to
      #' the classifier in the remaining columns.
      #' @param gold_standard A data frame with gold_standard HP-PPIs and class
      #' label indicating if such PPIs are positive or negative.
      #' @param classifier The type of classifier to use. See \code{caret}
      #' for the availbale classifiers.
      #' @param resampling.method The resampling method:'boot', 'boot632',
      #' 'optimism_boot', boot_all', 'cv', 'repeatedcv', 'LOOCV', 'LGOCV';
      #' defaults to cv. See \code{\link[caret]{trainControl}} for more details.
      #' @param ncross  Number of partitions for cross-validation;
      #' defaults to 5.See \code{\link[caret]{trainControl}} for more details.
      #' @param repeats for repeated k-fold cross validation only;
      #'  defaults to 3.See \code{\link[caret]{rfeControl}}
      #'  for more details.
      #' @param verboseIter Logical value, indicating whether to check the status
      #' of training process;defaults to FALSE.
      #' @param plots Logical value, indicating whether to plot the performance
      #' of ensemble learning algorithm as compared to individual classifiers;
      #' defaults to TRUE.If the argument set to TRUE, plots will be saved in
      #' the current working directory.
      #' These plots are :
      #' \itemize{
      #' \item{pr_plot} - Precision-recall plot of ensemble classifier vs
      #' selected individual classifiers.
      #' \item{roc_plot} - ROC plot of ensemble classifier vs selected individual
      #' classifiers.
      #' \item{points_plot} - Plot accuracy, F1-score ,positive predictive value
      #' (PPV),sensitivity (SE), and Matthews correlation coefficient (MCC) of
      #' ensemble classifier vs selected individual classifiers.
      #' }
      #' @param filename A character string, indicating the output filename as an
      #' pdf object.
      #' @return
      #' Ensemble_training_output
      #' \itemize{
      #' \item{prediction score} - Prediction scores for whole dataset from each
      #' individual classifier.
      #' \item{Best} - Selected hyper parameters.
      #' \item{Parameter range} - Tested hyper parameters.
      #' \item{prediction_score_test} - Scores probabilities for test data from
      #' each individual classifier.
      #' \item{class_label} -  Class probabilities for test data from each
      #' individual classifier.
      #' }
      #' classifier_performance
      #' \itemize{
      #' \item{cm} - A confusion matrix.
      #' \item{ACC} - Accuracy.
      #' \item{SE} - Sensitivity.
      #' \item{SP} - Specificity.
      #' \item{PPV} - Positive Predictive Value.
      #' \item{F1} - F1-score.
      #' \item{MCC} - Matthews correlation coefficient.
      #' \item{Roc_Object} - A list of elements.
      #' See \code{\link[pROC]{roc}} for more details.
      #' \item{PR_Object} - A list of elements.
      #' See \code{\link[PRROC]{pr.curve}} for more details.
      #' }
      #' predicted_interactions - The input data frame of pairwise interactions,
      #' including classifier scores averaged across all models.
      #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
      #' @importFrom caret trainControl
      #' @importFrom caret train
      #' @importFrom dplyr inner_join
      #' @importFrom pROC roc
      #' @importFrom stats predict
      #' @importFrom pROC plot.roc
      #' @importFrom graphics legend
      #' @importFrom PRROC pr.curve
      #' @importFrom pROC ggroc
      #' @importFrom ggplot2 geom_jitter
      #' @importFrom ggplot2 ylab
      #' @importFrom ggplot2 xlab
      #' @importFrom ggplot2 element_line
      #' @importFrom ggplot2 scale_y_continuous
      #' @importFrom ggplot2 geom_line
      #' @importFrom ggplot2 element_blank
      #' @importFrom ggplot2 annotation_custom
      #' @importFrom ggplot2 theme_bw
      #' @importFrom ggplot2 coord_equal
      #' @importFrom ggplot2 geom_abline
      #' @importFrom grDevices dev.off
      #' @importFrom grDevices pdf
      #' @description This function uses an ensemble of classifiers to
      #' predict interactions from the sequence-based dataset. This ensemble
      #' algorithm combines different results generated from individual
      #' classifiers within the ensemble via average to enhance prediction.
      #' @export
      #' @examples
      #' data('example_data')
      #' features <- example_data[, -2]
      #' gd <- example_data[, c(1,2)]
      #' gd <- na.omit(gd)
      #' ppi <-pred_ensembel(features,gd,
      #' classifier = c("avNNet", "svmRadial", "ranger"),
      #' resampling.method = "cv",ncross = 2,verboseIter = FALSE,plots = FALSE,
      #' filename = "plots.pdf")
      #' #extract predicted interactions
      #' pred_interaction <- ppi[["predicted_interactions"]]
      pred_ensembel <-
        function(features,
                 gold_standard,
                 classifier = c("avNNet", "svmRadial", "ranger"),
                 resampling.method = "cv",
                 ncross = 2,
                 repeats = 2,
                 verboseIter = TRUE,
                 plots = TRUE,
                 filename = "plots.pdf") {
          PPI <- NULL
          . <- NULL
          model <- NULL
          x <- NULL
          y <- NULL
          X1 <- NULL
          X2 <- NULL

          t.data <-
            features

          `%nin%` <- Negate(`%in%`)
          p.data <-
            features %>%
            filter(PPI %nin% gold_standard$PPI)


          if (!is.data.frame(features)) {
            stop("Input data should be data.frame")
          }
          if (all(colnames(features) != "PPI") == TRUE) {
            stop("PPI is missing from feature...")
          }
          if (length(colnames(features)) == 0) {
            stop("Column names are missing...")
          }

          if (all(colnames(gold_standard) != "class") == TRUE) {
            stop("class attribute is absent from the gold_standard data")
          }

          if (all(colnames(gold_standard) != "PPI") == TRUE) {
            stop("PPI is missing from gold_standard...")
          }
          if (!is.character(gold_standard$class)) {
            stop("Class column must include character vector")
          }


          # training data preparation

          x_train <-
            inner_join(t.data, gold_standard, by = "PPI") %>%
            select(-PPI)

          # change charactor to factor and reverse the factor level
          x_train$class <- # chanage charactor to factor
            as.factor(as.character(x_train$class))
          x_train$class <-
            factor(x_train$class, levels = rev(levels(x_train$class)))

          # check tghe factor level
          x <- unique(x_train$class)
          x <- data.frame(x, as.integer(x))
          x <-
            x %>%
            filter(x == "Positive") %>%
            .$as.integer.x.
          if (x == 2) {
            x_train$class <- factor(x_train$class,
                                    levels = rev(levels(x_train$class))
            )
          }

          # train control
          trControl <- trainControl(
            method = resampling.method,
            number = ncross, # number of resampling iterations
            repeats = repeats,
            verboseIter = verboseIter,
            classProbs = TRUE,
            savePredictions = "final",
            allowParallel = TRUE
          )

          Ensemble_output <-
            list()
          for (ml in classifier) {
            message(ml)
            models <-
              train(class ~ .,
                    data = x_train,
                    method = "ranger",
                    metric = "Accuracy",
                    trControl = trControl
              )

            clcol <- which(names(p.data) == "PPI")
            preds <- # predict on new data
              predict(models, newdata = p.data[, -clcol], type = "prob")
            Ensemble_output[[ml]][["prediction_score"]] <-
              preds[, "Positive"]

            # get the best parameter
            Ensemble_output[[ml]][["Best"]] <-
              models$bestTune
            # list of tested parameter
            Ensemble_output[[ml]][["Parameter range"]] <-
              models[["results"]]
            # prediction score on test set during cv
            Ensemble_output[[ml]][["prediction_score_test"]] <-
              models[["pred"]]$Positive

            # class label
            Ensemble_output[[ml]][["class_label"]] <-
              models[["pred"]]$obs
          }

          # get prediction interactions data
          Pred_interactions <-
            list()
          for (ml in classifier) {
            Pred_interactions[[ml]] <-
              Ensemble_output[[ml]][["prediction_score"]]
          }

          Pred_interactions <-
            do.call(cbind, Pred_interactions)
          ensemble <-
            rowMeans(Pred_interactions)
          datScore <-
            as.data.frame(cbind(p.data$PPI, ensemble))
          colnames(datScore)[1] <- "PPI"
          datScore$ensemble <- as.numeric(as.character(datScore$ensemble))
          datScore <-
            separate(datScore, PPI, c(
              "Pathogen_protein",
              "Host_protein"
            ), sep = "~", remove = FALSE)

          ## class label output
          label_output <-
            list()
          for (ml in classifier) {
            label_output[[ml]] <-
              Ensemble_output[[ml]]$class_label
          }
          ss <-
            do.call(cbind, label_output)

          ss <- apply(ss, 2, function(x) {
            as.factor(ifelse(x == 2, "Negative", "Positive"))
          })

          ensemble <- apply(ss, 1, function(x) .MajVot(x))
          ss <- cbind(ss, ensemble)
          ss <- as.data.frame(ss)


          # classifier prediction output
          classifier_output <-
            list()
          for (ml in classifier) {
            classifier_output[[ml]] <-
              Ensemble_output[[ml]]$prediction_score_test
          }

          classifier_output <-
            do.call(cbind, classifier_output)

          ensemble <-
            rowMeans(classifier_output)

          df_perf <-
            cbind(classifier_output, ensemble)
          df_perf <- as.data.frame(df_perf)

          output <- list()
          for (i in names(df_perf)) {
            for (j in names(ss)) {
              df_perf.response <- df_perf
              df_perf.response[, i] <-
                as.factor(ifelse(df_perf.response[, i] > 0.5,
                                 "Positive", "Negative"
                ))
              x <- unique(df_perf.response[, i])
              x <- data.frame(x, as.integer(x))
              x <-
                x %>%
                filter(x == "Positive") %>%
                .$as.integer.x.
              if (x == 2) {
                df_perf.response[, i] <-
                  factor(df_perf.response[, i],
                         levels = rev(levels(df_perf.response[, i]))
                  )
              }

              ss1 <- ss
              ss1[, j] <- as.factor(as.character(ss1[, j]))
              x <- unique(ss1[, j])
              x <- data.frame(x, as.integer(x))
              x <-
                x %>%
                filter(x == "Positive") %>%
                .$as.integer.x.
              if (x == 2) {
                ss1[, j] <- factor(ss1[, j],
                                   levels = rev(levels(ss1[, j]))
                )
              }
              cm <-
                table(df_perf.response[, i], ss1[, j])

              TP <- as.double(cm[1, 1])
              TN <- as.double(cm[2, 2])
              FP <- as.double(cm[1, 2])
              FN <- as.double(cm[2, 1])

              ACC <- (TP + TN) / (TP + TN + FP + FN)
              SE <- TP / (TP + FN)
              SP <- TN / (FP + TN)
              PPV <- TP / (TP + FP)
              F1 <- 2 * TP / (2 * TP + TP + FN)
              MCC <- (TP * TN - FP * FN) /
                sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

              output[[i]][["cm"]] <- cm
              output[[i]][["ACC"]] <- ACC
              output[[i]][["SE"]] <- SE
              output[[i]][["SP"]] <- SP
              output[[i]][["PPV"]] <- PPV
              output[[i]][["F1"]] <- F1
              output[[i]][["MCC"]] <- MCC



              output[[i]][["Roc_Object"]] <- roc(ss1[, j], df_perf[, i])

              output[[i]][["PR_Object"]] <- pr.curve(
                scores.class0 =
                  df_perf[, i][ss1[, j] == "Positive"],
                scores.class1 =
                  df_perf[, i][ss1[, j] == "Negative"],
                curve = TRUE
              )
            }
          }


          if (plots) {

            ################## plot pr.curve curve

            pr.auc.list <- list()
            for (i in names(output)) {
              pr.auc.list[[i]] <- output[[i]][["PR_Object"]][["auc.integral"]]
            }
            pr.auc.list <-
              do.call(rbind, pr.auc.list) %>%
              as.data.frame(.) %>%
              rownames_to_column("model")
            colnames(pr.auc.list)[2] <- "auc"
            pr.auc.list$auc <- round(pr.auc.list$au, 2)



            mytheme <- gridExtra::ttheme_default(
              core = list(fg_params = list(cex = 0.6)),
              colhead = list(fg_params = list(cex = 0.6)),
              rowhead = list(fg_params = list(cex = 0.6))
            )


            PR_Object <- list()
            for (i in names(output)) {
              PR_Object[[i]] <- output[[i]][["PR_Object"]]$curve
            }

            PR_m <-
              do.call(rbind, lapply(PR_Object, data.frame)) %>%
              # as.data.frame(.) %>%
              rownames_to_column("model") %>%
              mutate(model = str_remove_all(model, "\\.\\d+"))

            for (i in seq_along(PR_m)) {
              PR_plot <-
                ggplot(PR_m, aes(x = X1, y = X2, color = model)) +
                geom_line() +
                theme(plot.title = element_text(hjust = 0.5)) +
                theme_bw() +
                scale_y_continuous(limits = c(0, 1)) +
                ggtitle("pr_plot") +
                theme(plot.title = element_text(hjust = 0.5)) +
                theme(
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                ) +
                theme(axis.ticks.length = unit(.25, "cm")) +
                theme(legend.position = "top", legend.box = "horizontal") +
                theme(legend.title = element_blank()) +
                theme(text = element_text(size = 14, color = "black")) +
                theme(axis.text = element_text(size = 12, color = "black")) +
                xlab("Recall") +
                ylab("Percision") +
                annotation_custom(gridExtra::tableGrob(pr.auc.list,
                                                       rows = NULL,
                                                       theme = mytheme
                ),
                xmin = unit(0.9, "npc"),
                xmax = unit(0.9, "npc"), ymin = 0.22, ymax = 0.22
                )
            }


            ######## plot the result from confusion matrix
            cm_output <- list()
            for (i in names(output)) {
              cm_output[[i]][["ACC"]] <- output[[i]][["ACC"]]
              cm_output[[i]][["SE"]] <- output[[i]][["SE"]]
              cm_output[[i]][["PPV"]] <- output[[i]][["PPV"]]
              cm_output[[i]][["F1"]] <- output[[i]][["F1"]]
              cm_output[[i]][["MCC"]] <- output[[i]][["MCC"]]
            }

            cm_output <-
              do.call(rbind, cm_output) %>%
              as.data.frame(.) %>%
              rownames_to_column("model")

            point_plot <-
              cm_output %>%
              gather(x, y, 2:6) %>%
              ggplot(aes(x = x, y = as.numeric(y), color = model)) +
              geom_jitter(width = 0.2, size = 2) +
              theme_bw() +
              ggtitle("points_plot") +
              theme(legend.position = "top", legend.box = "horizontal") +
              theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
              labs(
                x = "Performance metrics", y = "value", title = "",
                colour = ""
              ) +
              theme(text = element_text(size = 14, color = "black")) +
              theme(axis.text = element_text(
                color = "black", # label color and size
                size = 12
              )) +
              theme(axis.ticks.length = unit(.25, "cm")) +
              theme(plot.title = element_text(hjust = 0.5))



            ################### plot ROC curve
            auc.list <- list()
            for (i in names(output)) {
              auc.list[[i]] <- output[[i]][["Roc_Object"]][["auc"]]
            }
            auc.list <-
              do.call(rbind, auc.list) %>%
              as.data.frame(.) %>%
              rownames_to_column("model")
            colnames(auc.list)[2] <- "auc"
            auc.list$auc <- round(auc.list$au, 2)



            Roc_Object <- list()
            for (i in names(output)) {
              Roc_Object[[i]] <- output[[i]][["Roc_Object"]]
            }

            ROC_plot <-
              ggroc(Roc_Object, legacy.axes = TRUE) +
              theme(plot.title = element_text(hjust = 0.5)) +
              geom_abline(
                slope = 1, intercept = 0,
                linetype = "dashed", alpha = 0.7, color = "black"
              ) +
              coord_equal() +
              theme_bw() +
              ggtitle("roc_plot") +
              theme(plot.title = element_text(hjust = 0.5)) +
              theme(
                axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
              ) +
              theme(axis.ticks.length = unit(.25, "cm")) +
              theme(legend.position = "top", legend.box = "horizontal") +
              theme(legend.title = element_blank()) +
              theme(text = element_text(size = 14, color = "black")) +
              theme(axis.text = element_text(size = 12, color = "black")) +
              xlab("False Positive Rate (1-Specificity)") +
              ylab("True Positive Rate (Sensitivity)") +
              annotation_custom(gridExtra::tableGrob(auc.list,
                                                     rows = NULL,
                                                     theme = mytheme
              ),
              xmin = unit(0.8, "npc"),
              xmax = unit(0.8, "npc"), ymin = 0.22, ymax = 0.22
              )


            # Make plots.
            plot_list <- list()
            plot_list$roc <- ROC_plot
            plot_list$pr <- PR_plot
            plot_list$point <- point_plot

            pdf(filename)
            for (i in seq_along(plot_list)) {
              print(plot_list[[i]])
            }
            dev.off()
          }

          final_output <- list()
          final_output[["Ensemble_training_output"]] <- Ensemble_output
          final_output[["classifier_performance"]] <- output
          final_output[["predicted_interactions"]] <- datScore

          return(final_output)
        }
