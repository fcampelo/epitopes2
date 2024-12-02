#' Apply minimum Redundancy Maximum Relevance Feature Selection (mRMR)
#'
#' This function is adapted (almost verbatim) from package `colino` by Steven Pawley
#' (<https://github.com/stevenpawley/colino>).
#' `step_select_mrmr` creates a *specification* of a recipe step that will apply
#' minimum Redundancy Maximum Relevance Feature Selection (mRMR) to numeric
#' data. The top `top_p` scoring features, or features whose scores occur in the
#' top percentile `threshold` will be retained as new predictors.
#'
#' @param recipe 	A recipe object. The step will be added to the sequence of
#'   operations for this recipe
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See selections() for more details. For the tidy
#'   method, these are not currently used
#' @param role Not used by this step since no new variables are created
#' @param trained A logical to indicate if the quantities for preprocessing have
#'   been estimated
#' @param outcome A character string specifying the name of response variable
#'   used to evaluate mRMR.
#' @param top_p An integer that will be used to select the number of best
#'   scoring features.
#' @param threshold A numeric value between 0 and 1 representing the percentile
#'   of best scoring features to select. For example `threshold = 0.9` will
#'   retain only predictors with scores in the top 90th percentile and a smaller
#'   threshold will select more features. Note that `top_p` and `threshold` are
#'   mutually exclusive but either can be used in conjunction with `cutoff` to
#'   select the top-ranked features and those that have filter scores that are
#'   larger than the cutoff value.
#' @param cutoff A numeric value where predictors with _larger_ absolute filter
#'   scores than the cutoff will be retained. A value of `NA` implies that this
#'   criterion will be ignored.
#' @param threads An integer specifying the number of threads to use for
#'   processing. The default = 0 uses all available threads.
#' @param exclude A character vector of predictor names that will be removed
#'   from the data. This will be set when `prep()` is used on the recipe and
#'   should not be set by the user.
#' @param scores A tibble with 'variable' and 'scores' columns containing the
#'   names of the variables and their mRMR scores. This parameter is only
#'   produced after the recipe has been trained.
#' @param skip A logical. Should the step be skipped when the recipe is baked by
#'   bake.recipe()? While all operations are baked when prep.recipe() is run,
#'   some operations may not be able to be conducted on new data (e.g.
#'   processing the outcome variable(s)). Care should be taken when using skip =
#'   TRUE as it may affect the computations for subsequent operations.
#' @param id 	A character string that is unique to this step to identify it.
#' @return A step_select_mrmr object.
#' @keywords datagen
#' @concept preprocessing
#' @concept supervised_filter
#' @export
#' @details
#'
#' The recipe will stop if all of `top_p`, `threshold` and `cutoff` are left
#' unspecified.
#'

step_select_mrmr <- function(
    recipe, ...,
    outcome = NULL,
    role = NA,
    trained = FALSE,
    top_p = NA,
    threshold = NA,
    cutoff = NA,
    threads = 0,
    exclude = NULL,
    scores = NULL,
    skip = FALSE,
    id = recipes::rand_id("select_mrmr")) {

  recipes::recipes_pkg_check("praznik")

  terms <- recipes::ellipse_check(...)

  recipes::add_step(
    recipe,
    step_select_mrmr_new(
      terms = terms,
      trained = trained,
      outcome = outcome,
      role = role,
      top_p = top_p,
      threshold = threshold,
      cutoff = cutoff,
      threads = threads,
      exclude = exclude,
      scores = scores,
      skip = skip,
      id = id
    )
  )
}

step_select_mrmr_new <- function(terms, role, trained,
                                 outcome, top_p, threshold,
                                 cutoff, threads, exclude,
                                 scores, skip, id) {

  recipes::step(
    subclass = "select_mrmr",
    terms = terms,
    role = role,
    trained = trained,
    outcome = outcome,
    top_p = top_p,
    threshold = threshold,
    cutoff = cutoff,
    threads = threads,
    exclude = exclude,
    scores = scores,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_select_mrmr <- function(x, training, info = NULL, ...) {
  # extract response and predictor names
  y_name <- recipes::recipes_eval_select(x$outcome, training, info)
  y_name <- y_name[1]
  x_names <- recipes::recipes_eval_select(x$terms, training, info)

  # check criteria
  check_criteria(x$top_p, x$threshold, match.call())
  check_zero_one(x$threshold)
  x$top_p <- check_top_p(x$top_p, length(x_names))

  if (length(x_names) > 0) {

    call <- rlang::call2(
      .fn = "MRMR",
      .ns = "praznik",
      X = rlang::quo(training[, x_names]),
      Y = rlang::quo(training[[y_name]]),
      k = length(x_names),
      threads = x$threads
    )

    res <- rlang::eval_tidy(call)

    res <- dplyr::tibble(
      variable = names(res$selection),
      score = res$score
    )

    exclude <-
      dual_filter(res$score, x$top_p, x$threshold, x$cutoff, maximize = TRUE)

  } else {
    exclude <- character()
  }

  step_select_mrmr_new(
    terms = x$terms,
    trained = TRUE,
    role = x$role,
    outcome = y_name,
    top_p = x$top_p,
    threshold = x$threshold,
    cutoff = x$cutoff,
    threads = x$threads,
    exclude = exclude,
    scores = res,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_select_mrmr <- function(object, new_data, ...) {
  if (length(object$exclude) > 0) {
    new_data <- new_data[, !(colnames(new_data) %in% object$exclude)]
  }
  dplyr::as_tibble(new_data)
}

#' @export
print.step_select_mrmr <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("mRMR feature selection")

    if (recipes::is_trained(x)) {
      n <- length(x$exclude)
      cat(paste0(" (", n, " excluded)"))
    }
    cat("\n")

    invisible(x)
  }

#' @rdname step_select_mrmr
#' @param x A `step_select_mrmr` object.
#' @param type A character with either 'terms' (the default) to return a
#'   tibble containing the variables that have been removed by the filter step,
#'   or 'scores' to return the scores for each variable.
#' @export
tidy.step_select_mrmr <- function(x, type = "terms", ...) {
  tidy_filter_step(x, type)
}

#' @export
tunable.step_select_mrmr <- function(x, ...) {
  dplyr::tibble(
    name = c("top_p", "threshold", "cutoff"),
    call_info = list(
      list(pkg = "epitopes", fun = "top_p"),
      list(pkg = "dials", fun = "threshold", range = c(0, 1)),
      list(pkg = "epitopes", fun = "cutoff")
    ),
    source = "recipe",
    component = "step_select_mrmr",
    component_id = x$id
  )
}

#' @export
required_pkgs.step_select_mrmr <- function(x, ...) {
  c("praznik")
}


# For step_select_mrmr
check_zero_one <- function(x) {
  if (is.na(x)) {
    return(x)
  } else {
    if (is.numeric(x)) {
      if (x >= 1 | x <= 0) {
        rlang::abort("`threshold` should be on (0, 1).")
      }
    } else {
      rlang::abort("`threshold` should be numeric.")
    }
  }
  return(x)
}

check_top_p <- function(x, n) {
  # checks on x (top_p) and n (number of features)
  if (is.na(x)) {
    return(x)
  }

  if (!is.numeric(x)) {
    rlang::abort("`top_p` should be numeric.")
  }

  if (!is.integer(x)) {
    x <- as.integer(x)
  }

  msg <- paste0("`top_p` should be on (1, ", n, ") based on the number of features available.")

  # return top_n = all features if top_n > n
  if (x >= n) {
    rlang::warn(msg)
    x <- min(n - 1, x)

    # return a single feature if top_p < 1
  } else if (x < 1) {
    rlang::warn(msg)
    x <- 1
  }

  return(x)
}

check_criteria <- function(top_p, threshold, cl) {
  if (is.na(top_p) & is.na(threshold)) {
    msg <- paste0(
      "For `",
      cl[[1]],
      "`, `top_p` and `threshold` cannot both be missing."
    )
    rlang::abort(msg)
  }
  invisible(NULL)
}

#' Select features using `top_p` or `threshold`.
#'
#' Feature selection using either the `top_p` or `threshold` features OR
#' `cutoff` where cutoff refers to the absolute numeric value of the feature
#' importance scores.
#'
#' @details
#' `dual_filter` selects feature that are selected using either (`top_p`,
#' `threshold`) or `cutoff` or both. If top_p/threshold and cutoff are both used
#' then features are selected using OR. For example, if top_p selects features 1
#' & 2, and threshold selects features 1 & 3, then the selected features =
#' 1,2,3.
#'
#' @param x a named numeric vector of scores per feature
#' @param top_p an integer specifying the number of top-performing features to
#'   retain
#' @param threshold a numeric with percentile of top-performing features to
#'   retain. For example, `threshold = 0.9` will only retain features that are
#'   in the top 90th percentile. A smaller value of threshold will select
#'   more features.
#' @param cutoff a numeric with the value that represents the cutoff in the
#'   scores in `x` by which to retain/discard features.
#' @param maximize logical to indicate whether `top_p`, `threshold` and `cutoff`
#'   are used to keep features where high scores = 'best' (maximize = TRUE) or
#'   where low scores = 'best' (maximize = FALSE).
#'
#' @return character vector of feature names to exclude
#' @keywords internal
dual_filter <- function(x, top_p, threshold, cutoff, maximize) {
  if (!is.na(top_p) & !is.na(threshold)) {
    rlang::abort("`top_p` and `threshold` are mutually exclusive")
  }

  na_x <- x[is.na(x)]
  x <- x[!is.na(x)]
  x <- sort(x, decreasing = maximize)

  p <- length(x)

  # assign logical selection variable using top_p
  if (!is.na(top_p)) {
    top_p_lgl <- seq_along(x) <= top_p
  } else {
    top_p_lgl <- rep(FALSE, p)
  }

  # assign logical selection variable using threshold
  if (!is.na(threshold)) {
    p_to_exceed <- stats::quantile(x, threshold)

    if (maximize) {
      threshold_lgl <- x >= p_to_exceed
    } else {
      threshold_lgl <- x < p_to_exceed
    }

  } else {
    threshold_lgl <- rep(FALSE, p)
  }

  # assign logical selection variable using cutoff
  if (!is.na(cutoff)) {
    if (maximize) {
      cutoff_lgl <- x >= cutoff
    } else {
      cutoff_lgl <- x <= cutoff
    }

  } else {
    cutoff_lgl <- rep(FALSE, p)
  }

  keep_lgl <- top_p_lgl | threshold_lgl | cutoff_lgl
  excluded <- c(names(x)[!keep_lgl], names(na_x))

  return(excluded)
}

tidy_filter_step <- function(x, type = "terms") {
  if (recipes::is_trained(x)) {
    if (type == "terms") {
      res <- dplyr::tibble(terms = x$exclude)
    } else if (type == "scores") {
      res <- x$scores
      res <- res[order(res$score, decreasing = TRUE), ]
    }

  } else {
    res <- dplyr::tibble(terms = rlang::na_chr)
  }
  res$id <- x$id
  res
}

#' Parameter functions for feature selection recipes
#'
#' Feature selection recipes allow the top-performing features to be selected
#' using three parameters. `cutoff` is for selecting features using the absolute
#' value in the filter methods scores.
#'
#' @param range A two-element vector holding the _defaults_ for the smallest and
#'   largest possible values, respectively.
#' @param trans A `trans` object from the `scales` package, such as
#'   `scales::log10_trans()` or `scales::reciprocal_trans()`. If not provided,
#'   the default is used which matches the units used in `range`. If no
#'   transformation, `NULL`.
#'
#' @return A function with classes "quant_param" and "param"
#' @export
#'
#' @examples
#' cutoff(c(3.5, 15))
cutoff <- function(range = c(dials::unknown(), dials::unknown()), trans = NULL) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(FALSE, FALSE),
    trans = trans,
    label = c(cutoff = "Absolute cutoff threshold for the feature scores")
  )
}
