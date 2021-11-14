# fits a single random forest and assess the performance on a holdout split
fit_RF_model <- function(df,
                         sample.rebalancing,
                         holdout.split,
                         threshold,
                         ncpus,
                         ...){


  # Get number of trees (for inbag setup)
  ntrees <- ifelse("num.trees" %in% ...names(),
                   ...elt(which(...names() == "num.trees")),
                   formals(ranger::ranger)$num.trees)

  # remove the holdout split (if present)
  if(!is.null(holdout.split)){
    df.tr <- df[which(df$Info_split != holdout.split), ]
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
                         min.node.size  = 100,
                         oob.error      = FALSE,
                         ...)

  return(myRF)

}
