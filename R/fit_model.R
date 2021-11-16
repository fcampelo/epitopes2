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
#' @section Performance assessment:
#' This function has the following modes of performance assessment:
#'
#' \itemize{
#'     \item If both _holdout.split_ and _CV.folds_ are `NULL`, no performance
#'     assessment is done. A model is fit on the full data and returned.
#'     \item If _holdout.split_ is the valid name of a data split in
#'     `peptides.list$df$Info_split` **AND** _CV.folds_ is `NULL`, a model is
#'     trained on the full data except the _holdout.split_ and then applied to
#'     predict the labels for _holdout.split_. The performance returned
#'     corresponds to the performance on _holdout.split_.
#'     \item If _CV.folds_ is a vector of valid names of data splits in
#'     `peptides.list$df$Info_split`, the splits named are used as
#'     cross-validation folds. The performance returned
#'     corresponds to the average cross-validation performance using the
#'     _CV.folds_. Notice that if _CV.folds_ is not `NULL` the value of
#'     _holdout.split_ is ignored.
#' }
#'
#' **IMPORTANT**: in all cases, the model returned by this routine is a model
#' trained on the **full data** (fit after the performance is assessed on
#' holdout splits or on cross-validation folds).
#'
#' @param peptides.list data frame containing the training data (one or more
#' numerical predictors and one **Class** attribute).
#' @param holdout.split name of split to be used as a holdout set. If `NULL`
#' then the full data is used for training.
#' @param CV.folds vector with the names of the splits to be used as CV folds.
#' If `NULL` no cross-validation is performed.
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
                      threshold = 0.5,
                      holdout.split = NULL,
                      CV.folds = NULL,
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
                          is.null(holdout.split) || (
                            is.character(holdout.split) &&
                              length(holdout.split) == 1 &&
                              holdout.split %in% unique(peptides.list$df$Info_split)),
                          is.null(CV.folds) || (
                            is.character(CV.folds) &&
                              length(CV.folds) > 1 &&
                              all(CV.folds %in% unique(peptides.list$df$Info_split))),
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
  df <- peptides.list$df %>%
    dplyr::mutate(Class = as.factor(.data$Class))
  if(use.global.features && "global.features" %in% class(peptides.list)){
    if(peptides.list$splits.attrs$split_level == "peptide"){
      warning("Using global features when split_level is 'peptide' is\n",
              "likely to result in data leakage through protein-level\n",
              "information. Performance estimates may not accurately\n",
              "represent expected generalisation behaviour.")
    }

    message("Merging global features into windowed dataframe...")
    df <- dplyr::left_join(df,
                           dplyr::select(peptides.list$proteins,
                                         -dplyr::starts_with("TSeq"), -c("DB")),
                           by = c("Info_protein_id" = "UID"))
  }

  full.model <- NULL
  if(is.null(CV.folds)){
    # hold-out OR no assessment mode
    # Fit random forest
    message("Fitting model...")
    mymodel <- fit_RF(df, sample.rebalancing, holdout.split, threshold, ncpus, ...)
    perf    <- mymodel$perf
    if(is.null(holdout.split)) full.model <- mymodel$RF.model

  } else {
    # Cross-validation

  }

  # Fit full model here
  if (is.null(full.model)){
    message("Fitting full model")
    mymodel <- fit_RF(df, sample.rebalancing, holdout.split = NULL, threshold, ncpus, ...)
    full.model <- mymodel$RF.model
  }

  .Random.seed <- oldseed

  # Assemble return list
  # return()

}
