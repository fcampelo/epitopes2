#' Calculate classic features for epitope prediction
#'
#' This function is used to calculate several distinct families of features for
#' epitope prediction. These are 'classic' features, based on statistical
#' summaries and average physico-chemical characteristics of the local
#' neighborhood of each residue in a protein sequence, or of the protein as a
#' whole.
#'
#' The following features are calculated based on the implementations available
#' in package
#' [**protr**](https://cran.r-project.org/package=protr). As of December 2023,
#' the following groups of features are supported by the **epitopes** package
#' (see the [protr vignette](https://CRAN.R-project.org/package=protr/vignettes/protr.html#3_package_overview)
#' for details on each of these):
#' \itemize{
#'    \item "AAC" - Amino acid composition
#'    \item "DC"  - Dipeptide composition
#'    \item "TC"  - Tripeptide composition
#'    \item  CTD descriptors:
#'    \itemize{
#'        \item "CTDC" - Composition
#'        \item "CTDT" - Transition
#'        \item "CTDD" - Distribution
#'    }
#'    \item "CTriad" - Conjoint triad descriptors
#'    \item  Quasi-sequence-order descriptors:
#'    \itemize{
#'        \item "SOCN" - Sequence-order-coupling number (with maximum lag
#'        `nlag = 3`)
#'        \item "QSO" - Quasi-sequence-order descriptors (with maximum lag
#'        `nlag = 3` and weighting factor `w = 0.1`)
#'    }
#'    \item Proteochemometric Modeling descriptors:
#'    \itemize{
#'        \item "ScalesGap" - Scales-based descriptors derived by Principal
#'        Components Analysis (using all properties in the `protr::AAindex`
#'        matrix, `pc = 5` and `lag = 3`)
#'    }
#' }
#'
#' **NOTE**: in all feature groups except "ScalesGap", invalid AA codes
#' (B, J, O, U, X, Y) are removed from the strings prior to feature calculation.
#' In "ScalesGap" these codes are replaced by a gap indicator, "-".
#'
#' Besides those, the following features are also available based on native
#' implementations:
#' \itemize{
#'     \item "Entropy" - the Shannon entropy of a sequence.
#'     \item "Atoms" - the number of C, H, N, O, S atoms in the sequence
#'     \item "MolWeight" - the total molecular weight of the peptide
#'     \item "AAtypes" - the proportion of AAs of each type (acidic, aliphatic,
#'     acidic, etc.)
#'     \item "BLOSUM" - BLOSUM-derived descriptors (same as
#'     [protr::extractBLOSUM()] with `submat = "AABLOSUM62"`,
#'     `k = 5`, `lag = 3` and `scale = TRUE`)
#'     \item "LegacyFeatures" - calculates the features used in paper
#'     \doi{10.1093/bioinformatics/btab536}
#' }
#'
#' Each feature group may be used for peptides or full protein sequences
#' Note, however, the warning from the **protr** documentation:
#' "*Users need to intelligently evaluate the underlying*
#' *details of the descriptors provided, instead of using protr with their data*
#' *blindly, especially for the descriptor types with more flexibility. It*
#' *would be wise for the users to use some negative and positive control*
#' *comparisons where relevant to help guide interpretation of the results.*".
#' Users should therefore be savvy when choosing which features to use for
#' epitope prediction, and the choice should ideally be guided by domain
#' expertise. Certain feature groups may not make sense for short peptides,
#' as they may be almost completely uninformative (e.g., tripeptide composition,
#' "TC"); or may require specific overriding of standard parameters (e.g.,
#' parameter `lambda` for "PAAC" must be smaller than the length of the shortest
#' peptide).
#'
#' @section **Feature Vectors**:
#' Input vector `features` is used to define which features are calculated
#' for proteins or peptides. This input parameter must be a
#' character vector, where each element is one of the names listed in
#' **Details**. For more information on the features calculated by **protr**,
#' check `?protr::extractXYZ`, replacing `XYZ` by the group abbreviation (see
#' **Details** or the documentation of the **protr** package for the list of
#' available feature groups). Notice that, for consistency purposes, the user
#' parameters of all **protr** features are kept fixed in this routine.
#' If the user wishes to add other features (or the same features with distinct
#' parameters) they can calculate those separately
#' and bind them to `X` (or to `X$df`, if X is a peptide list).
#' All feature columns calculated in this function
#' will have names starting with "feat_", and if externally-calculated features
#' are added they should follow the same pattern (for compatibility with the
#' other functions in this package).
#'
#' @param X a peptide.list object (returned by [extract_labelled_data()] or a
#' data frame containing one column with sequences for feature
#' calculation.
#' @param seqs_column name of the column containing the sequences
#' @param features vector of names of features to be calculated.
#' See **Feature Vector** for details.
#' @param ncpus positive integer, number of cores to use
#'
#' @return Updated `X` with features appended as columns
#' (directly to `X` if it is a data.frame, or to `X$df` if it is a peptide.list object)
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
calc_features_classic <- function(X,
                                  seqs_column,
                                  features,
                                  ncpus = 1){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(X) | is.peptide.list(X),
                          is.character(seqs_column),
                          length(seqs_column) == 1,
                          is.character(features),
                          length(features) > 0,
                          assertthat::is.count(ncpus))

  if(is.data.frame(X)){
    assertthat::assert_that(nrow(X) > 0,
                            seqs_column %in% names(X))
  } else {
    assertthat::assert_that(nrow(X$df) > 0,
                            seqs_column %in% names(X$df))
  }

  if(is.peptide.list(X)) df <- X$df else df <- X
  # ========================================================================== #

  SEQs <- df %>%
    dplyr::select(dplyr::all_of(seqs_column)) %>%
    unlist() %>% unname()

  myres <- vector("list", length(features))
  for (i in seq_along(features)){
    message("Calculating features: ", features[i])

    cl <- set_mc(ncpus)
    myres[[i]] <- mypblapply(X = SEQs,
                             FUN = call_feat_functions,
                             feat.name = features[i],
                             myargs = get_feature_args(features[i]),
                             ncpus = ncpus) %>%
      dplyr::bind_rows()
    close_mc(cl)
  }

  y <- dplyr::bind_cols(myres)

  idx <- which(names(df) %in% names(y))

  if(length(idx) > 0) df <- df[, -idx]

  df <- dplyr::bind_cols(df, y)

  attr(df, "features") <- features

  if(is.peptide.list(X)) X$df <- df else X <- df

  return(X)
}
