#' Calculate epitope classification performance indices
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
#'        calculate AUC, if a `prob` vector is passed.)
#'
#' @return list containing several performance indicators
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

calc_performance <- function(truth, pred, prob = NULL,
                             posValue = 1, negValue = -1,
                             ret.as.list = FALSE){

  assertthat::assert_that(length(posValue) == 1,
                          length(negValue) == 1,
                          length(truth) == length(pred),
                          all(stats::na.omit(truth) %in% c(posValue, negValue)),
                          all(stats::na.omit(pred)  %in% c(posValue, negValue)),
                          is.null(prob) | all(prob >= 0),
                          is.null(prob) | all(prob <= 1),
                          is.null(prob) | length(prob) == length(pred),
                          is.logical(ret.as.list))

  idx   <- which(is.na(truth) | is.na(pred))
  nPos  <- sum(pred == posValue)
  nNeg  <- sum(pred == negValue)
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
  mccDen  <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)

  if(!is.null(prob)){
    # Implementation validated in large-samples against pROC::roc
    df  <- data.frame(prob = prob, class = truth)
    tr  <- sort(unique(prob), decreasing = TRUE)
    tpr <- fpr <- numeric(length(tr))
    auc <- 0
    for (i in seq_along(tr)){
      df$pred <- ifelse(df$prob >= tr[i],
                        yes = posValue,
                        no  = negValue)
      tpr[i] <- sum(df$pred == posValue & df$class == posValue) / sum(df$class == posValue)
      fpr[i] <- sum(df$pred == posValue & df$class == negValue) / sum(df$class == negValue)
      if (i > 1){
        auc <- auc + (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i-1]) / 2
      }
    }
  } else {
    auc = NULL
    tpr = NULL
    fpr = NULL
  }

  if (!ret.as.list){
    out <- data.frame(
      TP        = TP,
      TN        = TN,
      FP        = FP,
      FN        = FN,
      sens      = TP / ifelse(TP + FN == 0, 1, TP + FN),
      spec      = TN / ifelse(TN + FP == 0, 1, TN + FP),
      ppv       = TP / ifelse(TP + FP == 0, 1, TP + FP),
      npv       = TN / ifelse(TN + FN == 0, 1, TN + FN),
      f1        = 2 * TP / ifelse(2 * TP + FP + FN == 0, 1, 2 * TP + FP + FN),
      mcc       = mccNum / ifelse(mccDen == 0, 1, sqrt(mccDen)),
      accuracy  =  (TP + TN) / (TP + TN + FP + FN),
      auc       = auc,
      #tpr       = tpr,
      #fpr       = fpr,
      n_predPos = nPos,
      n_predNeg = nNeg,
      n_NA      = length(idx))
  } else {
    out <- list(
      TP        = TP,
      TN        = TN,
      FP        = FP,
      FN        = FN,
      sens      = TP / ifelse(TP + FN == 0, 1, TP + FN),
      spec      = TN / ifelse(TN + FP == 0, 1, TN + FP),
      ppv       = TP / ifelse(TP + FP == 0, 1, TP + FP),
      npv       = TN / ifelse(TN + FN == 0, 1, TN + FN),
      f1        = 2 * TP / ifelse(2 * TP + FP + FN == 0, 1, 2 * TP + FP + FN),
      mcc       = mccNum / ifelse(mccDen == 0, 1, sqrt(mccDen)),
      accuracy  =  (TP + TN) / (TP + TN + FP + FN),
      auc       = auc,
      tpr       = tpr,
      fpr       = fpr,
      n_predPos = nPos,
      n_predNeg = nNeg,
      n_NA      = length(idx))
  }

  return(out)
}
