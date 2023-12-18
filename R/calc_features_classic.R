#' Calculate classic features for epitope prediction
#'
#' This function is used to calculate several distinct families of features for
#' epitope prediction. These are 'classic' features, based on statistical
#' summaries and average physico-chemical characteristics of the local
#' neighborhood of each residue in a protein sequence, or of the protein as a
#' whole.
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
#' **IMPORTANT**: if these features are to be used for epitope prediction,
#' the **global** features should be avoided or approached with extreme care, as
#' their use may result in accidental data leakage across splits and contaminate
#' performance calculations.
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
#' @param peptides.list list object containing a data frame `df` (with labelled
#' peptides, 1 row per residue) and a data frame `proteins` (with protein
#' information). Commonly returned by make_data_splits().
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
calc_features_classic <- function(peptides.df = NULL,
                                  proteins.df = NULL,
                                  features = character()){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.null(peptides.df) | is.data.frame(peptides.df),
                          is.null(proteins.df) | is.data.frame(proteins.df),
                          (is.null(proteins.df) + is.null(peptides.df)) == 1,
                          is.character(features),
                          length(features) > 0)
  # ========================================================================== #

  message("Calculating features:")
  if(!is.null(peptides.df)){
    assertthat::assert_that("Info_window" %in% names(peptides.df),
                            nrow(peptides.df) > 0)
    SEQs <- peptides.df$Info_window
  } else if(!is.null(proteins.df)){
    assertthat::assert_that("TSeq_sequence" %in% names(proteins.df),
                            nrow(proteins.df) > 0)
    SEQs <- proteins.df$TSeq_sequence
  }


  for (i in seq_along(features)){
    y <- call_feat_functions(SEQs = SEQs, feat.name = features[i])

    # if(is.data.frame(y)) {
    #   torm <- which(names(peptides.list$df) %in% names(y))
    #   if(length(torm) > 0) peptides.list$df <- peptides.list$df[, -torm]
    #   peptides.list$df <- dplyr::bind_cols(peptides.list$df, y)
    # }
  }
  # class(peptides.list) <- unique(c(class(peptides.list), "local.features"))
  # peptides.list$feature.attrs$local.features <- local.features




  # # Calculate global features
  # if(length(global.features) > 0) {
  #   for (i in seq_along(global.features)){
  #     y <- call_feat_functions(SEQs      = peptides.list$proteins$TSeq_sequence,
  #                              feat.name = global.features[i],
  #                              txt.opts  = c("global", "proteins"),
  #                              dfnames   = names(peptides.list$proteins))
  #
  #     if(is.data.frame(y)) {
  #       peptides.list$proteins <- cbind(peptides.list$proteins, y)
  #       if(is.data.frame(y)) {
  #         torm <- which(names(peptides.list$proteins) %in% names(y))
  #         if(length(torm) > 0) peptides.list$proteins <- peptides.list$proteins[, -torm]
  #         peptides.list$proteins <- dplyr::bind_cols(peptides.list$proteins, y)
  #       }
  #     }
  #   }
  #   class(peptides.list) <- unique(c(class(peptides.list), "global.features"))
  #   peptides.list$feature.attrs$global.features <- global.features
  #
  #   if(peptides.list$splits.attrs$split_level == "peptide"){
  #     warning("Global features should not be used for classification when\n",
  #             "the data is split at the 'peptide' level, since it can result\n",
  #             "in data leakage. Proceed with care.")
  #   }
  # }
  # return(peptides.list)
}
