#' Fit Random Forest model to epitope data
#'
#' Fits a Random Forest model to epitope data, using the data splits and
#' features previously calculated using [make_data_splits()] and
#' [calc_features()], respectively.
#'
#' Function [make_data_splits()] defines data splits based on (protein or
#' peptide) similarity. The split identifiers are stored in
#' `peptides.list$df$Info_split`. Function [calc_features()] calculates
#' the local and/or global features for each entry. Local features are stored in
#' `peptides.list$df` as columns starting with `feat_local_`, whereas global
#' features are stored in `peptides.list$proteins`, as columns starting with
#' `feat_global_`.
#' **NOTE**: Global features should only be used if the splitting level
#' used was "protein", otherwise they may cause contamination of performance
#' assessment due to data leakage. The splitting level of the data can be
#' checked on `peptides.list$splits.attrs$split_level`.
#'
#' @section Dealing with class imbalance:
#' Parameter `sample.rebalancing` regulates whether the resulting model attempts
#' to compensate class imbalances. If `TRUE` the Random Forest model is subject
#' to cost-sensitive training, which is done internally by setting the
#' parameter `case.weights` in the call to [ranger::ranger()] to a vector where
#' each observation of class _i_ has a weight equal to `1 / K_i`, where
#' `K_i` is the total number of cases of class `i` in the training data.
#'
#' @param peptides.list data frame containing the training data (one or more
#' numerical predictors and one **Class** attribute).
#' @param assessment.mode mode of performance assessment to use. Accepts "CV"
#' (for cross-validation using all splits in `peptides.list$df`) or "holdout"
#' (in which case the holdout split must be informed in `holdout.split`).
#' Defaults to "CV".
#' @param holdout.split name of split to be used as a holdout set.Ignored if
#' `assessment.mode = "CV"`.
#' @param threshold probability threshold for attributing a prediction as
#' *positive*.
#' @param sample.rebalancing logical: should the model try to compensate class
#' imbalances by weighted sampling of examples when training the trees in the
#' random forest? See **Dealing with class imbalance**.
#' @param use.global.features logical: should global features (potentially
#' available in `peptides.list$proteins`) be used? Should be left as the default
#' unless the user knows exactly what they're doing. See **Details**.
#' @param ncpus number of cores to use.
#' @param rnd.seed seed for random number generator. **Note**: this function
#' always returns the state of the random number generator back to its original
#' value before returning the results. Running it twice in sequence with
#' `rnd.seed = NULL` should result in exactly the same results (but it is
#' safer to simply set a specific seed, or to reuse the seed that is returned
#' in the output list).
#' @param ... other options to be passed down to [ranger::ranger()].
#'
#' @return List containing the fitted model and several performance indicators.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data

fit_model <- function(peptides.list,
                      assessment.mode = c("CV", "holdout"),
                      threshold = 0.5,
                      holdout.split = NULL,
                      sample.rebalancing = TRUE,
                      use.global.features = ifelse(peptides.list$splits.attrs$split_level == "protein", TRUE, FALSE),
                      ncpus = 1,
                      rnd.seed = NULL,
                      ...){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assessment.mode = assessment.mode[1]
  assertthat::assert_that(is.list(peptides.list),
                          all(c("df", "proteins") %in% names(peptides.list)),
                          "local.features" %in% class(peptides.list),
                          "splitted.peptide.data" %in% class(peptides.list),
                          is.character(assessment.mode),
                          length(assessment.mode) == 1,
                          assessment.mode %in% c("CV", "holdout"),
                          is.null(holdout.split) || (
                            is.character(holdout.split) &&
                              length(holdout.split) == 1 &&
                              holdout.split %in% unique(peptides.list$df$Info_split)),
                          is.numeric(threshold), length(threshold) == 1,
                          threshold >= 0, threshold <= 1,
                          is.logical(sample.rebalancing),
                          length(sample.rebalancing) == 1,
                          is.logical(use.global.features),
                          length(use.global.features) == 1,
                          assertthat::is.count(ncpus),
                          is.null(rnd.seed) || is.integer(rnd.seed))

  oldseed <- .Random.seed
  if(!is.null(rnd.seed)){
    set.seed(rnd.seed)
  }

  # Merge global features into data frame if needed/available
  df <- peptides.list$df
  if(use.global.features && "global.features" %in% class(peptides.list)){
    message("Merging global features into windowed dataframe...")
    df <- dplyr::left_join(df,
                           dplyr::select(peptides.list$proteins,
                                         -dplyr::starts_with("TSeq"), -c("DB")),
                           by = c("Info_protein_id" = "UID"))
  }

  if (assessment.mode == "holdout"){
    message("Fitting model in 'holdout' mode...")
    ho.preds   <- fit_RF_holdout(df, sample.rebalancing, holdout.split, threshold, ncpus, ...)

  } else if (assessment.mode == "CV"){

  } else stop("assessment mode '", assessment.mode, "' not recognized.")



#
#   # ========================================================================== #
#   # Estimate model performance
#   if (assessment.mode == "holdout"){
#     # Set holdout set (case weights = 0)
#     df <- peptides.list$df %>%
#       dplyr::mutate("weight" = ifelse(.data$Info_split == holdout.split, 0, 1))
#
#
#
#     df <- df %>%
#       dplyr::mutate(pred_prob  = myRF$predictions[, "1"],
#                     pred_class = ifelse(myRF$predictions[, "1"] >= threshold, "1", "-1"))
#
#     myperf <- calc_performance(truth = df$Class[df$Info_split == holdout.split],
#                                pred  = df$pred_class[df$Info_split == holdout.split],
#                                prob  = df$pred_prob[df$Info_split == holdout.split],
#                                posValue = "1",
#                                negValue = "-1",
#                                ret.as.list = TRUE)
#
#   }




  .Random.seed <- oldseed
  return()

}
