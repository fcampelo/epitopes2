#'  Fits a single random forest and assess the performance on a holdout split.
#'
#'  Not to be used separately, this function is called as part of calls to
#'  [fit_model()].
#'
#' @param df dataframe containing a `Class` attribute and one or more feature
#' attributes (columns starting with "feat_")
#' @param threshold probability threshold for attributing a prediction as
#' *positive*.
#' @param sample.rebalancing character: should the model try to compensate class
#' imbalances during training? See **Dealing with class imbalance** for details.
#' @param ncpus number of cores to use.
#' @param holdout.split name of split to be used as a holdout set.
#' @param ... other options to be passed down to [ranger::ranger()].
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export

fit_RF <- function(df,
                   sample.rebalancing,
                   holdout.split = NULL,
                   threshold,
                   ncpus,
                   ...){

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
  inbag        <- NULL
  if(sample.rebalancing == "by_tree"){
    message("Sample rebalancing by modified sampling probabilities")
    case.weights <- ifelse(df.tr$Class == "1",
                           sum(df.tr$Class == "-1") / nrow(df.tr),
                           sum(df.tr$Class == "1") / nrow(df.tr))

  } else if (sample.rebalancing == "undersampling"){
    message("Sample rebalancing by stratified undersampling")
    # Balance classes in such a way that all *peptides* of the majority
    # class retain some representation in the sub-sampled training set.
    minclass <- names(which.min(table(df.tr$Class)))
    nmin     <- sum(df.tr$Class == minclass)
    nmaj     <- sum(df.tr$Class != minclass)
    df.tmp   <- dplyr::filter(df.tr, .data$Class != minclass) %>%
      dplyr::group_by(.data$Info_PepID) %>%
      dplyr::sample_n(ceiling(dplyr::n() * nmin / nmaj)) %>%
      dplyr::ungroup()
    if (nrow(df.tmp) > nmin){
      df.tmp <- dplyr::sample_n(df.tmp, nmin)
    }
    df.tr <- dplyr::filter(df.tr, .data$Class == minclass) %>%
      dplyr::bind_rows(df.tmp)

    rm(df.tmp)
  }

  message("Fitting Random Forest...")
  myRF <- ranger::ranger(Class ~ .,
                         data =  dplyr::select(df.tr,
                                               dplyr::starts_with("feat_"),
                                               "Class"),
                         classification  = TRUE,
                         probability     = TRUE,
                         inbag           = inbag,
                         case.weights    = case.weights,
                         num.threads     = ncpus,
                         oob.error       = FALSE,
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
