#' Calculate features for epitope prediction
#'
#' This function is used to calculate several distinct families of features for
#' epitope prediction.
#'
#' Two major types of features can be calculated:
#' \itemize{
#'     \item _local features_, which are calculated based on the local
#'     neighbourhood of each AA position (column *Info_window* of the windowed
#'     data frame).
#'     \item _global features_, which are calculated using the full sequence
#'     of the protein (listed in column *TSeq_sequence*  of the protein
#'     data frame).
#' }
#'
#' The following features are calculated based on the implementations available
#' in package
#' [**protr**](https://cran.r-project.org/package=protr). As of August 2021,
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
#'        \item "BLOSUM" - BLOSUM -derived descriptors (with
#'        `submat = "AABLOSUM62"`, `k = 5` and `lag = 3`)
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
#'     \item "MolWeight" - the total molecular weight of the sequence
#'     \item "AAtypes" - the proportion of AAs of each type (acidic, aliphatic,
#'     acidic, etc.)
#'     \item "LegacyFeatures" - calculates the features used in paper
#'     \doi{10.1093/bioinformatics/btab536}
#' }
#'
#' Each feature group may be used either at the local or global level -
#' which does not mean it _should_ be. The **protr** documentation provides
#' the following warning: "*Users need to intelligently evaluate the underlying*
#' *details of the descriptors provided, instead of using protr with their data*
#' *blindly, especially for the descriptor types with more flexibility. It*
#' *would be wise for the users to use some negative and positive control*
#' *comparisons where relevant to help guide interpretation of the results.*".
#' Users should therefore be savvy when choosing which features to use for epitope
#' prediction, and the choice should ideally be guided by domain expertise.
#' Certain feature groups may not make sense as local features,
#' as the (usually very short) local substrings will not allow the features to
#' be informative (e.g., tripeptide composition, "TC"); or may require
#' specific overriding of standard parameters (e.g., parameter
#' `lambda` for "PAAC" must be smaller than the length of the shortest local
#' string).
#'
#' **IMPORTANT**: if the splitting level intended to be used for modelling is
#' "peptide" (this can be checked in `peptides.list$splits.attrs$split_level`)
#' then global features should be avoided, as their use can result in
#' data leakage across splits and contaminate performance calculations.
#'
#' @section **Feature Vectors**:
#' Input vectors `local.features` and `global.features` are used to define
#' which features are calculated at either level. These input parameters must be
#' character vectors, where each element is one of the names listed in
#' **Details**. For more information on the features calculated by **protr**,
#' check `?protr::extractXYZ`, replacing `XYZ` by the group abbreviation (see
#' **Details** or the documentation of the **protr** package for the list of
#' available feature groups). Notice that, for consistency purposes, the user
#' parameters of all **protr** features are kept fixed in this routine.
#' If the user wishes to add other features (or the same features with distinct
#' parameters) they can calculate those separately
#' and bind them to `peptides.list$df` (for local features) or
#' `peptides.list$proteins` (for global ones). All feature columns should have
#' names starting with "feat_local_" or "feat_global_".
#'
#' @param peptides.list list object returned by [make_data_splits()].
#' @param local.features list of features to be calculated
#' at the local neighbourhood (`peptides.list$df$Info_window`) level.
#' @param global.features lists of features to be calculated
#' at the global level (`peptides.list$proteins$TSeq_sequence`).
#' See **Feature Lists** for details.
#' @param ... any other parameters (currently unused)
#'
#' @return Updated `peptides.list` object, with local features added as columns
#' to `peptides.list$df`, and global features added as columns to
#' `peptides.list$proteins`.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
calc_features <- function(peptides.list,
                          local.features = character(),
                          global.features = character(),
                          ...){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.list(peptides.list),
                          all(c("df", "proteins") %in% names(peptides.list)),
                          is.character(local.features),
                          is.character(global.features))
  # ========================================================================== #

  # Prevent paralellisation errors in Windows
  # TODO: fix this
  # if(.Platform$OS.type == "windows") {
  #   message('Attention: multicore feature calculations are not currently supported for Windows')
  #   ncpus <- 1
  # }

  message("Calculating features:")
  # Calculate local features
  if(length(local.features) > 0) {
    for (i in seq_along(local.features)){
      y <- call_feat_functions(SEQs      = peptides.list$df$Info_window,
                               feat.name = local.features[i],
                               txt.opts  = c("local", "df"),
                               dfnames   = names(peptides.list$df))

      if(is.data.frame(y)) {
        torm <- which(names(peptides.list$df) %in% names(y))
        if(length(torm) > 0) peptides.list$df <- peptides.list$df[, -torm]
        peptides.list$df <- dplyr::bind_cols(peptides.list$df, y)
      }
    }
    class(peptides.list) <- unique(c(class(peptides.list), "local.features"))
    peptides.list$feature.attrs$local.features <- local.features
  }

  # Calculate global features
  if(length(global.features) > 0) {
    for (i in seq_along(global.features)){
      y <- call_feat_functions(SEQs      = peptides.list$proteins$TSeq_sequence,
                               feat.name = global.features[i],
                               txt.opts  = c("global", "proteins"),
                               dfnames   = names(peptides.list$proteins))

      if(is.data.frame(y)) {
        peptides.list$proteins <- cbind(peptides.list$proteins, y)
        if(is.data.frame(y)) {
          torm <- which(names(peptides.list$proteins) %in% names(y))
          if(length(torm) > 0) peptides.list$proteins <- peptides.list$proteins[, -torm]
          peptides.list$proteins <- dplyr::bind_cols(peptides.list$proteins, y)
        }
      }
    }
    class(peptides.list) <- unique(c(class(peptides.list), "global.features"))
    peptides.list$feature.attrs$global.features <- global.features

    if(peptides.list$splits.attrs$split_level == "peptide"){
      warning("Global features should not be used for classification when\n",
              "the data is split at the 'peptide' level, since it can result\n",
              "in data leakage. Proceed with care.")
    }
  }
  return(peptides.list)
}
