#' Calculate classification performance indices
#'
#' This function is used to calculate several distinct BINARY classification
#' performance indicators. `NA` values are not considered in the calculations.
#'
#' @param truth vector of reference values
#' @param pred  vector of predicted values
#' @param prob  optional, vector of predicted probabilities
#' @param posValue value that indicates the "positive" class.
#' @param negValue value that indicates the "negative" class.
#' @param ret.as.list optional, should the result be returned as a list (if
#'        TRUE the routine also returns the `tpr` and `fpr` vectors used to
#'        calculate AUC, as a data frame "roc")
#' @param ncpus number of cores to use.
#' @param roc.points maximum number of points to use for the calculation of the
#' ROC curve and AUC value.
#'
#' @return list containing several performance indicators
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data

calc_performance <- function(truth, pred, prob = NULL,
                             posValue = 1, negValue = -1,
                             ret.as.list = FALSE,
                             ncpus = 1,
                             roc.points = 200){

  assertthat::assert_that(length(posValue) == 1,
                          length(negValue) == 1,
                          length(truth) == length(pred),
                          all(stats::na.omit(truth) %in% c(posValue, negValue)),
                          all(stats::na.omit(pred)  %in% c(posValue, negValue)),
                          is.null(prob) | all(prob >= 0),
                          is.null(prob) | all(prob <= 1),
                          is.null(prob) | length(prob) == length(pred),
                          is.logical(ret.as.list),
                          assertthat::is.count(roc.points),
                          roc.points >= 10)

  idx     <- which(is.na(truth) | is.na(pred))
  nPos    <- sum(pred == posValue)
  nNeg    <- sum(pred == negValue)
  truePos <- sum(truth == posValue, na.rm = TRUE)
  trueNeg <- sum(truth == negValue, na.rm = TRUE)

  if(length(idx) > 0){
    truth <- truth[-idx]
    pred  <- pred[-idx]
    prob  <- prob[-idx]
  }

  TN <- as.numeric(sum(truth == negValue & pred == negValue))
  FN <- as.numeric(sum(truth == posValue & pred == negValue))
  FP <- as.numeric(sum(truth == negValue & pred == posValue))
  TP <- as.numeric(sum(truth == posValue & pred == posValue))

  mccNum  <- TP * TN - FP * FN
  mccDen  <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  sens    <- TP / ifelse(TP + FN == 0, 1, TP + FN)
  spec    <- TN / ifelse(TN + FP == 0, 1, TN + FP)
  ppv     <- TP / ifelse(TP + FP == 0, 1, TP + FP)
  npv     <- TN / ifelse(TN + FN == 0, 1, TN + FN)
  f1      <- 2 * TP / ifelse(2 * TP + FP + FN == 0, 1, 2 * TP + FP + FN)

  if (truePos == 0) sens <- ppv <- f1 <- NA
  if (trueNeg == 0) spec <- npv <- NA


  if(!is.null(prob)){
    tr <- sort(unique(prob), decreasing = TRUE)
    if (length(unique(prob)) > roc.points){
      ii <- round(seq(from = 1, to = length(tr), length.out = roc.points))
      tr <- tr[ii]
    }

    #message("Calculating ROC curve:")
    roc <- lapply(tr,
                  function(tri, prob, truth, posValue, negValue){
                    tpr  <- sum(prob >= tri & truth == posValue) / sum(truth == posValue)
                    fpr  <- sum(prob >= tri & truth == negValue) / sum(truth == negValue)
                    return(data.frame(tpr = tpr, fpr = fpr))
                  },
                  prob = prob, truth = truth,
                  posValue = posValue, negValue = negValue) %>%
      dplyr::bind_rows()

    auc <- roc %>%
      dplyr::mutate(tprl = dplyr::lag(.data$tpr),
                    fprl = dplyr::lag(.data$fpr)) %>%
      dplyr::summarise(auc = sum((.data$fpr - .data$fprl) * (.data$tpr + .data$tprl) / 2,
                                 na.rm = TRUE)) %>%
      as.numeric()

    # # Implementation validated against pROC::roc
    # df  <- data.frame(prob = prob, class = truth)
    # tr  <- sort(unique(prob), decreasing = TRUE)
    # tpr <- fpr <- numeric(length(tr))
    # auc <- 0
    # for (i in seq_along(tr)){
    #   df$pred <- ifelse(df$prob >= tr[i],
    #                     yes = posValue,
    #                     no  = negValue)
    #   tpr[i] <- sum(df$pred == posValue & df$class == posValue) / sum(df$class == posValue)
    #   fpr[i] <- sum(df$pred == posValue & df$class == negValue) / sum(df$class == negValue)
    #   if (i > 1){
    #     auc <- auc + (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i-1]) / 2
    #   }
    # }
  } else {
    auc <- NA
    tpr <- NA
    fpr <- NA
  }

  out  <- data.frame(
    TP        = TP,
    TN        = TN,
    FP        = FP,
    FN        = FN,
    sens      = sens,
    spec      = spec,
    ppv       = ppv,
    npv       = npv,
    f1        = f1,
    mcc       = mccNum / ifelse(mccDen == 0, 1, mccDen),
    accuracy  =  (TP + TN) / (TP + TN + FP + FN),
    bal.acc   = .5 * (sens + spec),
    gmean     = sqrt(sens*spec),
    auc       = auc,
    n_predPos = nPos,
    n_predNeg = nNeg,
    n_NA      = length(idx))

  if(ret.as.list){
    out     <- as.list(out)
    out$roc <- roc
  }

  return(out)
}
