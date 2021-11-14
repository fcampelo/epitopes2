# fits a single random forest and assess the performance on a holdout split
fit_RF_model <- function(df,
                         sample.rebalancing,
                         holdout.split,
                         threshold,
                         ncpus,
                         ...){


  # Get some Ranger parameters (if passed)
  ntrees <- ifelse("num.trees" %in% ...names(),
                   ...elt(which(...names() == "num.trees")),
                   formals(ranger::ranger)$num.trees)

  min.node.size <- ifelse("min.node.size" %in% ...names(),
                          ...elt(which(...names() == "min.node.size")),
                          50)

  # isolate the holdout split (if present)
  if(!is.null(holdout.split)){
    message("Model fitting with hold-out split = ", holdout.split)
    df.tr <- df[which(df$Info_split != holdout.split), ]
    df.ho <- df[which(df$Info_split == holdout.split), ]
  }

  # Determine sampled cases for each tree
  # (ignoring samples from holdout.split since they have weight = 0)
  if (sample.rebalancing){
    message("Sampling class-balanced observations for each tree...")
    inbag  <- mypblapply(1:ntrees,
                         function(i, clvec){
                           nsmp <- ceiling(min(table(clvec))* 2/3)
                           smp  <- c(sample(x       = which(clvec == 1),
                                            size    = nsmp,
                                            replace = FALSE),
                                     sample(x       = which(clvec == -1),
                                            size    = nsmp,
                                            replace = FALSE))
                           return(as.numeric(seq_along(clvec) %in% smp))},
                         clvec = df.tr$Class,
                         ncpus = ncpus)

  } else {
    message("Sampling observations for each tree...")
    inbag  <- mypblapply(1:ntrees,
                         function(i, clvec){
                           nsmp <- ceiling(length(clvec)* 2/3)
                           smp  <- sample.int(length(clvec),
                                              size = nsmp,
                                              replace = FALSE)
                           return(as.numeric(seq_along(clvec) %in% smp))},
                         clvec = df.tr$Class,
                         ncpus = ncpus)
  }

  message("Fitting Random Forest...")
  myRF <- ranger::ranger(Class ~ .,
                         data =  dplyr::select(df.tr,
                                               dplyr::starts_with("feat_"),
                                               "Class"),
                         classification = TRUE,
                         probability    = TRUE,
                         class.weights  = nrow(df.tr)/table(df.tr$Class),
                         inbag          = inbag,
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
