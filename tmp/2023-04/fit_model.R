#' Fit Random Forest model to epitope data
#'
#' This function fits a Random Forest model to epitope data (after feature
#' calculation using [calc_features()]). Feature column names should start with
#' *feat_*, and the class attribute should be called *Class*. All other columns
#' are ignored when fitting the model.
#'
#' @param data.train data frame containing the training data (one or more
#' numerical predictors and one **Class** attribute).
#' @param data.test data frame containing the test data, in the same format as
#' `data.train`.
#' @param rnd.seed seed for random number generator
#' @param ncpus number of cores to use.
#' @param threshold probability threshold for attributing a prediction as
#' *positive*.
#' @param ... options to be passed down to the Random Forest implementation (see
#' [ranger::ranger()] for details)
#'
#' @return List containing the fitted model, predictions for  and several performance indicators.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%

fit_model <- function(data.train,
                      data.test = NULL,
                      rnd.seed = NULL,
                      ncpus = 1,
                      threshold = 0.5,
                      ...){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(data.train),
                          assertthat::is.count(rnd.seed),
                          assertthat::is.count(ncpus),
                          is.numeric(threshold),
                          length(threshold) == 1,
                          threshold >= 0, threshold <= 1)

  if(!is.null(rnd.seed)) set.seed(rnd.seed)

  # Some preliminary preprocessing
  Class = NA # just to stop a CRAN check message, no effect on the code
  data.train <- dplyr::mutate(as.data.frame(data.train),
                              dplyr::across(!dplyr::starts_with("feat_"), as.character),
                              dplyr::across(dplyr::starts_with("feat_"), function(x) {x + 1e-9}),
                              Class = as.factor(Class))

  if (!is.null(data.test)){
    data.test <- dplyr::mutate(as.data.frame(data.test),
                               dplyr::across(!dplyr::starts_with("feat_"), as.character),
                               dplyr::across(dplyr::starts_with("feat_"), function(x) {x + 1e-9}),
                               Class = as.factor(Class))
  }

  myRF <- ranger::ranger(Class ~ .,
                         data = dplyr::select(data.train,
                                              dplyr::starts_with("feat_"),
                                              Class),
                         probability = TRUE,
                         num.threads = ncpus,
                         ...)


  outlist <- list(rf_model    = myRF,
                  rf_oob.err  = myRF$prediction.error,
                  rf_preds    = NULL,
                  rf_probs    = NULL,
                  rf_class    = NULL,
                  rf_perf     = NULL)


  if(!is.null(data.test)){
    test.pred <- stats::predict(myRF,
                                data = dplyr::select(data.test,
                                                     dplyr::starts_with("feat_"),
                                                     Class))
    rf_probs <- test.pred$predictions[, "1"]
    rf_class <- ifelse(rf_probs >= threshold, 1, -1)

    if(!("Info_UID" %in% names(data.test))){
      data.test$Info_UID <- data.test$Info_protein_id
    }
    rf_preds <- data.test %>%
      dplyr::select(dplyr::all_of(c("Info_UID", "Info_center_pos", "Class"))) %>%
      dplyr::bind_cols(pred_prob = rf_probs) %>%
      dplyr::mutate(Class       = as.numeric(as.character(Class)),
                    pred_class = as.numeric(as.character(rf_class)))

    myperf <- calc_performance(truth = data.test$Class,
                               pred  = as.factor(rf_class),
                               prob  = rf_probs,
                               posValue = "1",
                               negValue = "-1")

    outlist$rf_preds <- rf_preds
    outlist$rf_probs <- rf_probs
    outlist$rf_class <- rf_class
    outlist$rf_perf  <- myperf
  }

  return(outlist)

}
