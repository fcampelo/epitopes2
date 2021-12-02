#'  fits a single random forest and assess the performance on a holdout split
#'  @importFrom dplyr %>%
#'  @importFrom rlang .data
fit_RF <- function(df,
                   sample.rebalancing,
                   holdout.split = NULL,
                   threshold,
                   ncpus,
                   min.node.size,
                   ntrees, ...){

  # isolate the holdout split (if present)
  if(!is.null(holdout.split)){
    message("Fit model excluding split: ", holdout.split)
    df.tr <- df[which(df$Info_split != holdout.split), ]
    df.ho <- df[which(df$Info_split == holdout.split), ]
  } else {
    df.tr <- df
  }

  # Determine how class imbalance is treated
  case.weights <- NULL
  if(sample.rebalancing == "by_tree"){
    case.weights <- ifelse(df.tr$Class == "1",
                           sum(df.tr$Class == "-1") / nrow(df.tr),
                           sum(df.tr$Class == "1") / nrow(df.tr))
  } else if (sample.rebalancing == "undersampling"){
    df.tr <- df.tr %>%
      group_by(Class) %>%
      sample_n(size = min(table(df.tr$Class)), replace = FALSE) %>%
      ungroup()
  }

  message("Fitting Random Forest...")
  myRF <- ranger::ranger(Class ~ .,
                         data =  dplyr::select(df.tr,
                                               dplyr::starts_with("feat_"),
                                               "Class"),
                         classification = TRUE,
                         probability    = TRUE,
                         case.weights   = case.weights,
                         num.threads    = ncpus,
                         min.node.size  = min.node.size,
                         oob.error      = FALSE,
                         ...)

  # Assess performance
  if (!is.null(holdout.split)){
    preds <- stats::predict(myRF,
                            data = dplyr::select(df.ho,
                                                 dplyr::starts_with("feat_"),
                                                 Class = .data$Class))

    perf <-  calc_performance(truth = df.ho$Class,
                              pred  = ifelse(preds$predictions[, 2] >= threshold,
                                             "1", "-1"),
                              prob  = preds$predictions[, 2],
                              ret.as.list = TRUE,
                              posValue = "1",
                              negValue = "-1",
                              ncpus = ncpus)
  } else {
    perf <- NULL
  }

  # Assemble outlist
  outlist <- list(RF.model = myRF,
                  perf     = perf)

  return(outlist)

}
