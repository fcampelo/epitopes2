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
#' Parameter `sample.rebalancing` regulates how the modelling routine attempts
#' to compensate class imbalances. The following strategies are available:
#' \itemize{
#'     \item `"by_tree"`: changes the sampling procedure for each tree in the
#'     random forest, so that each observation has a sampling probability
#'     inversely proportional to its class prevalence (resulting in samples
#'     that are, on average, balanced for each tree).
#'     Works by setting the parameter `case.weights` of [ranger::ranger()].
#'     \item `"undersampling"`: undersamples the majority class to reach a
#'     balanced training set.
#'     \item any other value: no sample rebalancing is done.
#' }
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
#'     _CV.folds_. If _CV.folds_ is not `NULL` the value of
#'     _holdout.split_ is ignored.
#' }
#'
#' **IMPORTANT**: by default, the model returned by this routine is a model
#' trained on the **full data** (fit after the performance is assessed on
#' holdout splits or on cross-validation folds). This can be regulated by
#' parameter `return.model`.
#'
#' @param peptides.list data frame containing the training data (one or more
#' numerical predictors and one **Class** attribute).
#' @param holdout.split name of split to be used as a holdout set. If `NULL`
#' then the full data is used for training. Ignored if `CV.folds` is not `NULL`.
#' @param CV.folds vector with the names of the splits to be used for cross validation.
#' If `NULL` no cross-validation is performed.
#' @param threshold probability threshold for attributing a prediction as
#' *positive*.
#' @param sample.rebalancing character: should the model try to compensate class
#' imbalances during training? See **Dealing with class imbalance** for details.
#' @param use.global.features logical: should global features (potentially
#' available in `peptides.list$proteins`) be used? Should be left as the default
#' unless the user knows exactly what they're doing. See **Details**.
#' @param return.model model to be returned. Accepts `"full"` (return model
#' trained on the full data); `"partial"` (model trained on all data except the
#' holdout split; or on the set of data specified in all CV folds); or `"none"`
#' (does not return a model).
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
                      sample.rebalancing = 'by_tree',
                      use.global.features = ifelse(peptides.list$splits.attrs$split_level == "protein", TRUE, FALSE),
                      return.model = "full",
                      ncpus = 1,
                      rnd.seed = NULL,
                      ...){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.list(peptides.list),
                          all(c("df", "proteins") %in% names(peptides.list)),
                          "local.features" %in% class(peptides.list),
                          "splitted.peptide.data" %in% class(peptides.list),
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
                          is.character(sample.rebalancing),
                          length(sample.rebalancing) == 1,
                          is.logical(use.global.features),
                          length(use.global.features) == 1,
                          is.character(return.model),
                          length(return.model) == 1,
                          return.model %in% c("none", "partial", "full"),
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

  # prepare some parameters for the RF model
  min.node.size <- ifelse("min.node.size" %in% ...names(),
                          ...elt(which(...names() == "min.node.size")),
                          peptides.list$peptide.attrs$min_peptide)

  ntrees <- NULL
  if("num.trees" %in% ...names()) {
    ntrees <- ...elt(which(...names() == "num.trees"))
  }

  output.model <- NULL
  if(is.null(CV.folds)){
    # hold-out OR no assessment mode
    res <- fit_RF(df = df,
                  sample.rebalancing = sample.rebalancing,
                  holdout.split = holdout.split,
                  threshold = threshold,
                  ncpus = ncpus,
                  min.node.size = min.node.size,
                  ntrees = ntrees,
                  ...)

    perf     <- res$perf
    perflist <- NULL

    if (return.model == "partial"){
      output.model <- res$RF.model
    } else if (return.model == "full"){
      output.model <- fit_RF(df = df,
                             sample.rebalancing = sample.rebalancing,
                             threshold = threshold,
                             ncpus = ncpus,
                             min.node.size = min.node.size,
                             ntrees = ntrees,
                             ...)$RF.model
    }

  } else {
    # Cross-validation mode
    perflist <- vector("list", length(CV.folds))
    for (i in seq_along(CV.folds)){
      tmpdf <- dplyr::filter(df, Info_split %in% CV.folds)
      res   <- fit_RF(df = tmpdf,
                      sample.rebalancing = sample.rebalancing,
                      holdout.split = CV.folds[i],
                      threshold = threshold,
                      ncpus = ncpus,
                      min.node.size = min.node.size,
                      ntrees = ntrees,
                      ...)
      perflist[[i]] <- res$perf
    }

    perf <- lapply(perflist,
                   function(x){
                     x$roc <- NULL
                     as.data.frame(x)
                   }) %>%
      dplyr::bind_rows() %>%
      dplyr::summarise(across(everything(), mean)) %>%
      as.list()

    if (return.model == "partial"){
      output.model <- fit_RF(df = tmpdf,
                             sample.rebalancing = sample.rebalancing,
                             threshold = threshold,
                             ncpus = ncpus,
                             min.node.size = min.node.size,
                             ntrees = ntrees,
                             ...)$RF.model
    } else if (return.model == "full"){
      output.model <- fit_RF(df = df,
                             sample.rebalancing = sample.rebalancing,
                             threshold = threshold,
                             ncpus = ncpus,
                             min.node.size = min.node.size,
                             ntrees = ntrees,
                             ...)$RF.model
    }
  }

  .Random.seed <- oldseed

  # Assemble return list
  peptides.list$model       <- output.model
  peptides.list$model.perf  <- perf
  peptides.list$CV.perflist <- perflist
  peptides.list$model.attrs <- list(threshold = threshold,
                                    holdout.split = holdout.split,
                                    CV.folds = CV.folds,
                                    sample.rebalancing = sample.rebalancing,
                                    use.global.features = use.global.features,
                                    return.model = return.model,
                                    rnd.seed = rnd.seed,
                                    other.args = as.list(substitute(list(...)))[-1])

  return(peptides.list)

}
