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
  inbag        <- NULL
  if(sample.rebalancing == "by_tree"){
    message("Sample rebalancing by modified sampling probabilities")
    case.weights <- ifelse(df.tr$Class == "1",
                           sum(df.tr$Class == "-1") / nrow(df.tr),
                           sum(df.tr$Class == "1") / nrow(df.tr))

    # inbag  <- mypblapply(1:ntrees,
    #                      function(i, clvec, idvec){
    #                        minclass <- names(which.min(table(clvec)))
    #
    #                        minpool  <- data.frame(rown = which(clvec == minclass),
    #                                               idx = idvec[which(clvec == minclass)])
    #                        maxpool  <- data.frame(rown = which(clvec != minclass),
    #                                               idx = idvec[which(clvec != minclass)])
    #
    #                        minsmp <- dplyr::group_by(minpool, .data$idx) %>%
    #                          dplyr::sample_n(1)
    #                        maxsmp <- dplyr::group_by(maxpool, .data$idx) %>%
    #                          dplyr::sample_n(1) %>%
    #                          dplyr::ungroup() %>%
    #                          dplyr::sample_n(nrow(minsmp))
    #
    #                        return(c(minsmp$rown, maxsmp$rown))
    #                      },
    #                      clvec = df.tr$Class,
    #                      idvec = df.tr$Info_PepID,
    #                      ncpus = ncpus)

  } else if (sample.rebalancing == "undersampling"){
    message("Sample rebalancing by stratified undersampling")
    # Balance classes in such a way that all *peptides* of the majority
    # class retain some representation in the sub-sampledtraining set.
    minclass <- names(which.min(table(df.tr$Class)))
    nmin     <- sum(df.tr$Class == minclass)
    nmaj     <- sum(df.tr$Class != minclass)
    df.tr    <- dplyr::filter(df.tr, .data$Class != minclass) %>%
      dplyr::group_by(.data$Info_PepID) %>%
      dplyr::sample_frac(nmin/nmaj) %>%
      dplyr::ungroup() %>%
      dplyr::bind_rows(dplyr::filter(df.tr, .data$Class == minclass))
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
                         min.node.size   = min.node.size,
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
