#' Split epitope data.
#'
#' Split the windowed epitope data returned by [extract_peptides()] into
#' non-overlapping subsets. Proteins / peptides with similarities higher than
#' a predefined threshold are always placed in the same split to prevent data
#' leakage in machine learning. This routine tries to simultaneously approximate
#' the user-defined proportions and maintain the overall class balance within
#' each split.
#'
#' @section **Grouping strategy**:
#' The first step of this routine is to group the observations at either the
#' proteins or the peptide level (depending on input parameter `split_level`).
#' Local alignment scores for all pairs of sequences are calculated using the
#' implementation of the Smith-Waterman algorithm available in function
#' [Biostrings::pairwiseAlignment()], with default parameters and the scoring
#' matrix defined in `substitution_matrix`. These scores are returned as
#' element **SW.scores** of the output list.
#'
#' The dissimilarity matrix is calculated based on **SW.scores**, as:
#'
#' `diss[i,j] = 1 - SW.scores[i,j] / min(SW.scores[i,i], SW.scores[j,j])`
#'
#' which gives a value between 0 (perfect similarity) and 1 (maximum
#' dissimilarity). `diss[i,j]` will be 0 if and only if the shorter sequence is
#' fully and perfectly contained in the longer one; and will be 1 if and only
#' if the Smith-Waterman alignment score is zero.
#'
#' The dissimilarity matrix (returned as element **diss.matrix** of the output
#' list) is used to calculate a hierarchical clustering of
#' the sequences, and input parameter `similarity_threshold` is then used to
#' define the resulting similarity clusters. Single-linkage is used, to
#' guarantee that any pair of sequences with similarity
#' greater than `similarity_threshold` will always be contained within the same
#' cluster.
#' The resulting clusters (returned as element **clusters** of the output list)
#' represent the allocation units that are considered by the optimisation
#' to split the data.
#'
#' @section **Optimisation Problem**:
#' This routine attempts to simultaneously minimise two objectives: (i) the
#' difference between the actual proportion of data within each split and the
#' desired levels informed by `target_props`, and (ii) the difference between the
#' proportion of _positive_ observations within each split and the overall
#' proportion in the data. A simple linear aggregation strategy is used to
#' define the following optimisation problem. Let:
#' \itemize{
#'    \item nC:  number of clusters.
#'    \item K:   number of splits.
#'    \item xi:  integer variable defining the split to which cluster i is allocated.
#'    \item Ni+: number of _positive_ observations in cluster i.
#'    \item Ni:  total number of observations in cluster i.
#'    \item Gj*: desired proportion of data for split j.
#'    \item P*:  proportion of _positive_ observations in the whole data.
#' }
#'
#' We want to solve the problem:
#'
#' `minimize sum_j{ alpha x (Gj - Gj*)^2 + (1-alpha) x (pj - P*)^2 }`
#'
#' With:
#' \itemize{
#'     \item `xi \in {1, ..., K}`, for all `i = 1, ..., nC`
#'     \item `yij = ifelse(xi == j, 1, 0)`
#'     \item `Gj = sum_i{ yij * Ni } / sum_i{ Ni }`
#'     \item `pj = sum_i{ yij * Ni+ } / sum_i{ yij * Ni }`
#'}
#'
#' Input parameter `alpha` regulates the relative importance attributed to
#' each of the two objectives. At the limits, `alpha = 0` ignores the desired
#' proportions and
#' focuses only on defining splits with class balances that are as close as
#' possible to __P*__, whereas `alpha = 0` ignores the class balance and tries
#' to simply generate splits that are as faithful as possible to the desired
#' proportions. An approximation to the Pareto-optimal front can be obtained
#' by varying `alpha` between these two values.
#'
#' The search space of this optimisation problem has a cardinality of
#' `K ^ nC`. If the cardinality is under `10^6` possible
#' allocations then this routine performs enumerative search and is guaranteed
#' to return the global optimum. For larger search spaces a
#' constructive heuristic followed by Simulated Annealing (see [stats::optim()]
#' for details) is used. Input parameter `SAopts` can be used to pass control
#' parameters to the Simulated Annealing.
#'
#'
#' @param peptides.list list object returned by [extract_peptides()], containing
#'        the data frame of windowed epitope data (**df**) and the data frame of
#'        individual peptides (**peptides**).
#' @param proteins data frame of proteins, returned by [get_proteins()]
#'        containing all proteins listed in `peptides.list$df$Info_protein_id`.
#' @param split_level which sequences should be used to calculate similarity
#'        for splitting the data.
#'        Accepts "protein" (uses similarity of the full protein sequences,
#'        `proteins$TSeq_sequence`, to determine which observations should stay
#'        together in the splits) or "peptide" (uses similarity of the labelled
#'        peptides, `peptides.list$peptides$Info_peptide`).
#'        See **Grouping strategy** for details.
#' @param target_props numeric vector of target proportions for each split
#'        (i.e., a vector (p1, p2, ..., pK) such that 0 < pk < 1 for all k and
#'        sum(pk) = 1).
#' @param split_names optional, vector of names to be given to each split.
#' @param similarity_threshold similarity threshold for grouping observations.
#'        See **Grouping strategy** for details.
#' @param substitution_matrix character string indicating the substitution
#'        matrix to be used to calculate the peptide / protein alignments.
#'        Must be a matrix recognised by [Biostrings::substitution.matrices()]
#'        (e.g., "BLOSUM45", "PAM30", etc.)
#' @param alpha weight parameter to regulate focus on maximising match to desired
#'        split proportions (`target_props`) vs. on approximating the class
#'        balance of the full data set. Must be a numeric value between 0 and 1.
#'        See **Optimisation Problem** for details.
#' @param return_front should different tradeoff solutions (based on distinct
#'        values of alpha) be returned as well? Accepts either `NULL` (don't
#'        calculate front) or a positive integer greater than 2
#'        (number of equally spaced alpha values between 0 and 1 to evaluate).
#'        If it's a positive integer, then the tradeoff dataset will be exported
#'        as object `$splits.attrs$tradeoffs` in the output list. Note that
#'        only non-dominated tradeoffs are returned, so the number of tradeoff
#'        points is likely to be smaller than the value of `return_front`.
#'        Usually 11 or 21 points are enough to give a good rough idea of the
#'        tradeoffs.
#' @param SAopts list of control parameters to be used by the simulated
#'        annealing (SANN) algorithm. See [stats::optim()] for details.
#' @param save_folder path to folder for saving the results. It will save the
#' results as file *peptides_list.rds* (overwriting if necessary)
#' @param ncpus positive integer, number of cores to use.
#'
#' @return A list object containing:
#' \itemize{
#'    \item **df**: `peptides.list$df` updated to contain cluster number and
#'          split allocation of each entry
#'    \item **peptides**: `peptides.list$peptides` updated to contain cluster
#'          number and split allocation of each entry
#'    \item **proteins**: data frame containing data on all proteins listed in
#'          `df$Info_protein_id`.
#'    \item **peptide.attrs**: list inherited from `peptides.list`, containing
#'          relevant attributes used in the earlier call to [extract_peptides()].
#'    \item **splits.attrs**: list containing information about the splitting:
#'    \itemize{
#'        \item *split_level*: same as input parameter `split_level`
#'        \item *similarity_threshold*: same as input parameter `similarity_threshold`
#'        \item *substitution_matrix*: same as input parameter `substitution_matrix`
#'        \item *split_props*: proportion of data allocated to each split
#'        \item *split_balance*: proportion of positive observations in each split
#'        \item *target_props*: same as input parameter `target_props`
#'        \item *target_balance*: proportion of positive observations in full data
#'        \item *alpha*: same as input parameter `alpha`
#'        \item *SW.scores*: local alignment scores between each sequence (Smith-Waterman)
#'        \item *diss.matrix*: dissimilarity matrix (see **Grouping strategy** for details)
#'        \item *clusters*: `hclust` object with clustering structure.
#'        \item *cluster.alloc*: data frame summarising the split allocations.
#'    }
#'  }
#'
#' If the splitting is impossible (e.g., if the number of clusters is smaller
#' than the desired number of splits) the function throws a warning and returns
#' a list with only **SW.scores**, **diss.matrix** and **clusters**.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'

make_data_splits <- function(peptides.list,
                             proteins,
                             split_level  = "protein",
                             target_props = c(.75, .25),
                             split_names  = NULL,
                             similarity_threshold = .7,
                             substitution_matrix = "BLOSUM62",
                             alpha  = 0.5,
                             return_front = NULL,
                             SAopts = list(),
                             save_folder = NULL,
                             ncpus       = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  split_level <- tolower(split_level)
  if (is.null(split_names)) split_names <- sprintf("split_%02d_%02d",
                                                   seq_along(target_props),
                                                   round(100*target_props))

  assertthat::assert_that(is.list(peptides.list),
                          all(c("df", "peptides") %in% names(peptides.list)),
                          is.data.frame(proteins),
                          length(split_level) == 1,
                          split_level %in% c("peptide", "protein"),
                          is.numeric(target_props), length(target_props) > 1,
                          all(target_props > 0), all(target_props < 1),
                          sum(target_props) == 1,
                          is.character(split_names),
                          length(split_names) == length(target_props),
                          is.numeric(similarity_threshold),
                          length(similarity_threshold) == 1,
                          similarity_threshold > 0, similarity_threshold < 1,
                          is.numeric(alpha), length(alpha) == 1,
                          alpha >= 0, alpha <= 1,
                          is.null(return_front) || assertthat::is.count(return_front),
                          is.null(return_front) || return_front > 2,
                          is.list(SAopts),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          assertthat::is.count(ncpus))

  diss_t   <- 1 - similarity_threshold

  df <- peptides.list$df %>%
    dplyr::select(-starts_with("Info_cluster"),
                  -starts_with("Info_split"))
  peptides <- peptides.list$peptides %>%
    dplyr::select(-starts_with("Info_cluster"),
                  -starts_with("Info_split"))

  proteins <- proteins %>% dplyr::filter(.data$UID %in% df$Info_protein_id)

  message("Performing data split at ", split_level, " level")
  if(split_level == "peptide"){
    X <- peptides %>%
      dplyr::select(IDs  = .data$Info_PepID,
                    SEQs = .data$Info_peptide)
  } else if(split_level == "protein"){
    X <- proteins %>%
      dplyr::select(IDs  = .data$UID,
                    SEQs = .data$TSeq_sequence)
  }

  # ========================================================================== #
  # Extract similarity-based clusters

  # Run Smith-Waterman local alignment and build similarity score matrix
  message("Calculating similarities (normalized Smith-Waterman scores)")
  scores <- mypblapply(X   = seq_along(X$SEQs),
                       FUN = function(i){
                         patt <- rep(X$SEQs[i], times = 1 + length(X$SEQs) - i)
                         subj <- X$SEQs[i:length(X$SEQs)]
                         vals <- Biostrings::pairwiseAlignment(pattern = patt,
                                                               subject = subj,
                                                               substitutionMatrix = substitution_matrix,
                                                               type = "local",
                                                               scoreOnly = TRUE)
                         return(c(rep(NA, i - 1), vals))
                       },
                       ncpus = ncpus) %>%
    do.call(what = cbind)

  # Build denominator matrix: D_{ij} = min(scores_{i,i}, scores{j,j})
  denom <- matrix(pmin(rep(diag(scores), times = nrow(scores)),
                       rep(diag(scores), each = nrow(scores))),
                  nrow  = nrow(scores), byrow = FALSE)

  # Calculate normalized dissimilarity
  rownames(scores) <- colnames(scores) <- X$IDs
  diss <- 1 - scores / denom

  message("Extracting clusters (Hierarchical, single linkage)")
  clusters  <- stats::hclust(d = stats::as.dist(diss), method = "single")
  X$Cluster <- stats::cutree(clusters, h = diss_t)

  if(length(unique(X$Cluster)) < length(target_props)){
    warning("Impossible to divide data into ", length(target_props),
            " splits at similarity level ", similarity_threshold,
            ".\nTry a higher similarity threshold or a smaller number of splits.")
    return(list(SW.scores       = scores,
                diss.matrix     = diss,
                clusters        = clusters))
  }
  # Check how many positive / negative examples per group
  if (split_level == "protein") {
    Y <- df %>% dplyr::left_join(X, by = c("Info_protein_id" = "IDs"))
  } else if (split_level == "peptide"){
    Y <- df %>% dplyr::left_join(X, by = c("Info_PepID" = "IDs"))
  }
  Y <- Y %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(nPos = sum(.data$Class == 1),
                     nNeg = sum(.data$Class == -1),
                     N    = dplyr::n(),
                     P    = .data$nPos / dplyr::n())

  # ========================================================================== #
  # Optimise split alllocations

  if(!("maxit" %in% names(SAopts))) {
    SAopts$maxit <- min(1e5, 2000 * round(log10(length(target_props) ^ nrow(Y))))
  }
  y <- optimise_splits(Y = Y, Nstar = target_props, alpha = alpha,
                       SAopts = SAopts, ncpus = ncpus)

  Y$Split <- split_names[y$x]

  X <- dplyr::left_join(dplyr::select(X, -c("SEQs")),
                        dplyr::select(Y, c("Cluster", "Split")),
                        by = "Cluster") %>%
    dplyr::arrange(.data$Split, .data$Cluster)
  if(split_level == "peptide"){
    names(X)[1] <- "Info_PepID"
    df          <- dplyr::left_join(df, X, by = "Info_PepID")
    peptides    <- dplyr::left_join(peptides, X, by = "Info_PepID")
  } else if(split_level == "protein"){
    names(X)[1] <- "Info_protein_id"
    df <- dplyr::left_join(df, X, by = "Info_protein_id")
    peptides <- dplyr::left_join(peptides, X, by = "Info_protein_id")
  }

  names(y$solstats$Gj) <- split_names
  names(y$solstats$pj) <- split_names


  # Build tradeoff dataset if required
  if(!is.null(return_front)){
    alpha.vec <- seq(from = 0, to = 1, length.out = return_front)
    message("Building tradeoff dataset")
    for (i in seq_along(alpha.vec)){
      message("Processing tradeoff ", i, " of ", return_front)
      z <- optimise_splits(Y = Y, Nstar = target_props, alpha = alpha.vec[i],
                           SAopts = SAopts, ncpus = ncpus)
      f1 <- max(sort(target_props) - sort(z$solstats$Gj))
      f2 <- max(z$solstats$pj - sum(Y$nPos)/sum(Y$N))
      if (i == 1){
        tradeoffs <- data.frame(alpha = alpha.vec[i],
                                max.dev.perc = f1,
                                max.dev.prop = f2)
      } else {
        tradeoffs <- rbind(tradeoffs,
                           data.frame(alpha = alpha.vec[i],
                                      max.dev.perc = f1,
                                      max.dev.prop = f2))
      }
    }
    nondom <- rep(TRUE, length(alpha.vec))
    for (i in seq_along(alpha.vec)){
      idx <- matrix(as.matrix(tradeoffs[i, -1]),
                    nrow = length(alpha.vec) - 1,
                    ncol = 2, byrow = TRUE)
      nondom[i] <- !any(rowSums(idx - as.matrix(tradeoffs[-i, -1]) > 0) == 2)

    }
    tradeoffs <- tradeoffs[which(nondom), ]
    tradeoffs <- tradeoffs[-duplicated(tradeoffs[, -1]), ]
    tradeoffs <- tradeoffs[order(tradeoffs[, 3]), ]

  }

  # Build outlist
  outlist          <- peptides.list
  outlist$df       <- df %>%
    dplyr::rename(Info_cluster = c("Cluster"),
                  Info_split   = c("Split")) %>%
    dplyr::select(dplyr::starts_with("Info"), dplyr::everything())

  outlist$peptides <- peptides %>%
    dplyr::rename(Info_cluster = c("Cluster"),
                  Info_split   = c("Split")) %>%
    dplyr::select(dplyr::starts_with("Info"), dplyr::everything())

  outlist$proteins     <- proteins
  outlist$splits.attrs <- list(
    split_level          = split_level,
    similarity_threshold = similarity_threshold,
    substitution_matrix  = substitution_matrix,
    split_props          = y$solstats$Gj,
    split_balance        = y$solstats$pj,
    target_props         = target_props,
    target_balance       = y$solstats$Pstar,
    alpha                = alpha,
    SW.scores            = scores,
    diss.matrix          = diss,
    clusters             = clusters,
    cluster.alloc        = X)

  if(!is.null(return_front)){
    outlist$splits.attrs$tradeoffs <- tradeoffs
  }

  class(outlist) <- unique(c(class(outlist), "splitted.peptide.data"))

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    saveRDS(outlist, paste0(normalizePath(save_folder), "/peptides_list.rds"))
  }

  # Print resulting split statistics
  cl_tbl <- outlist$splits.attrs$cluster.alloc %>%
    group_by(Split) %>% summarise(Clusters = length(unique(Cluster)))

  message("============================================================")
  message("Data splitting summary")
  message("Splitting level: ", split_level)
  message("Similarity threshold: ",  similarity_threshold)
  message("Number of clusters found: ", length(unique(X$Cluster)))
  message("Target balance: ", signif(y$solstats$Pstar, 4))
  message("Target split proportions: ", paste(signif(target_props, 4), collapse = ", "))
  message("alpha: ",  alpha)
  for (i in seq_along(y$solstats$Gj)){
    message(names(y$solstats$Gj)[i], ": Split proportion = ",
            round(y$solstats$Gj[i], 4),
            " | Split balance: ",
            round(y$solstats$pj[i], 4),
            " | Number of clusters: ",
            cl_tbl$Clusters[i])
  }
  message("============================================================")

  # return results
  return(outlist)
}
