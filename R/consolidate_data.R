#' Consolidate epitope and protein data
#'
#' Label protein positions based on linear B-cell epitope information extracted
#' from IEDB.
#'
#' This routine is used to consolidate epitope information, by:
#'  \itemize{
#'    \item Filtering `proteins` such that only entries that appear at least
#'    once in `epitopes$protein_id` are retained
#'    \item Using the `epitopes` data (fields `protein_id`, `start_pos` and
#'    `end_pos`) to generate labels for each position of the selected proteins.
#'  }
#'
#' Prior to the consolidation, entries in `epitopes` input are removed if they:
#' \itemize{
#'   \item Lack a valid *protein_id* (i.e., if they have no corresponding
#'   entry in `proteins`)
#'   \item Lack a valid string in *epit_seq*
#'   \item Lack a valid definition in *epit_struc_def*
#'   \item Have a mismatch between the sequence in *epit_seq* and the
#'   corresponding sequence between *start_pos* and *end_pos* on the protein
#'   sequence.
#' }
#'
#' Entries in the `proteins` data frame with a value of *TSeq_seqtype* different
#' from "protein" are also removed.
#'
#' @param epitopes data frame of epitope data (returned by [get_LBCE()]).
#' @param proteins data frame of protein data (returned by [get_proteins()]).
#' @param save_folder path to folder for saving the results.
#' @param only_exact logical, should only sequences labelled as "Exact Epitope"
#'        in variable *epit_struc_def* (within `epitopes`) be considered?
#' @param set_positive how to decide whether an observation should be of the
#'        "Positive" (+1) or "Negative" (-1) class? Use "any" to set a sequence
#'        as positive if $n_positive > 0$, "mode" to set it if
#'        $n_positive >= n_negative$, or "all" to set it if $n_negative == 0$.
#'        Defaults to "mode".
#' @param ncpus positive integer, number of cores to use
#'
#' @return Data frame containing all proteins mentioned that appear at least
#' once in `epitopes$protein_id` (one row per aminoacid residue), containing
#' the consolidated information extracted from `epitopes`.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data

consolidate_data <- function(epitopes, proteins,
                             save_folder     = NULL,
                             only_exact      = FALSE,
                             set_positive    = c("any", "mode", "all"),
                             ncpus           = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(epitopes),
                          is.data.frame(proteins),
                          is.logical(only_exact), length(only_exact) == 1,
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          is.character(set_positive),
                          assertthat::is.count(ncpus))

  if(all(tolower(set_positive) == "any")) {
    set_positive <- "any"
  } else if(all(tolower(set_positive) == "all")) {
    set_positive <- "all"
  } else set_positive <- "mode"


  # ========================================================================== #
  # Initial preprocessing

  # Filter epitopes
  message("Cleaning epitopes dataset...")
  epits <- dplyr::filter(epitopes,
                         .data$protein_id %in% unique(proteins$UID),         # must have a corresponding protein
                         !is.na(.data$epit_seq),                             # must have an epit_seq
                         !is.na(.data$epit_struc_def)) %>%                   # must have a valid epit_struc_def
    dplyr::mutate(prot_substr = mapply(function(id, st, en){                    # filter by correct sequence in protein
      substr(proteins[which(proteins$UID == id), "TSeq_sequence"], st, en)},
      id = .data$protein_id, st = .data$start_pos, en = .data$end_pos)) %>%
    dplyr::filter(toupper(.data$epit_seq) == toupper(.data$prot_substr)) %>%
    dplyr::select(-"prot_substr")

  if (only_exact){
    epits <- dplyr::filter(epits, .data$epit_struc_def == "Exact Epitope")
  }

  # Filter proteins
  message("Cleaning proteins datset...")
  df <- dplyr::filter(proteins,
                      .data$UID %in% unique(epits$protein_id),
                      .data$TSeq_seqtype == "protein",
                      !is.na(.data$TSeq_sequence))

  idx <- which(is.na(df$TSeq_taxid))
  for(i in seq_along(idx)){
    df$TSeq_taxid[idx[i]] <- epits$sourceOrg_id[which(epits$protein_id == df$UID[idx[i]])][1]
  }

  # ========================================================================== #
  # Build long data frame

  # Build long protein data frame
  message("Building long data frame: proteins")
  df <- mypblapply(split(df, seq(nrow(df))),
                   FUN = function(x){
                     data.frame(Info_organism_id = x$TSeq_taxid,
                                Info_protein_id  = x$UID,
                                Info_pos         = 1:nchar(x$TSeq_sequence),
                                Info_AA          = strsplit(x$TSeq_sequence,
                                                            split = "")[[1]])},
                   ncpus = ncpus) %>%
    dplyr::bind_rows()


  # Build long epitope data frame
  message("Building long data frame: epitopes")
  epit_summary <- mypblapply(split(epits, seq(nrow(epits))),
                             FUN = function(x){
                               data.frame(Info_protein_id   = x$protein_id,
                                          Info_pos          = x$start_pos:x$end_pos,
                                          Info_pubmed_id    = x$pubmed_id,
                                          Info_epitope_id   = x$epitope_id,
                                          Info_sourceOrg_id = x$sourceOrg_id,
                                          Info_host_id      = x$host_id,
                                          Info_nPos         = x$n_Positive,
                                          Info_nNeg         = x$n_Negative,
                                          Info_type         = x$epit_struc_def)},
                             ncpus = ncpus) %>%
    dplyr::bind_rows() %>%
    #
    # Consolidate information by protein-position
    dplyr::group_by(.data$Info_protein_id, .data$Info_pos) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   ~ paste(.x, collapse = ",")),
                     .groups = "drop") %>%
    #
    # Remove duplicated information from specific fields
    dplyr::mutate(Info_pubmed_id    = get_uniques(.data$Info_pubmed_id),
                  Info_epitope_id   = get_uniques(.data$Info_epitope_id),
                  Info_sourceOrg_id = get_uniques(.data$Info_sourceOrg_id),
                  Info_host_id      = get_uniques(.data$Info_host_id),
                  Info_type         = get_uniques(.data$Info_type))

  # function to determine class:
  # TODO - make the aggregation strategy (by direct vote of assays?
  # majority vote of epitope_ids? other?) a user-defined parameter
  make_class <- function(pos, neg, set_positive){
    pos <- as.numeric(strsplit(pos, ",")[[1]])
    neg <- as.numeric(strsplit(neg, ",")[[1]])
    class <- switch(set_positive,
                    mode = sign(sum(sign(-1 + 2 * (pos > neg)))  - .1),
                    any  = -1 + 2 * any(-1 + 2 * (pos > 0) > 0),
                    all  = -1 + 2 * all(neg == 0 & pos > 0))
    return(class)
  }

  # Join epitope information onto long protein data frame
  message("Consolidating data...")
  df <- df %>%
    dplyr::left_join(epit_summary, by = c("Info_protein_id", "Info_pos")) %>%
    dplyr::select(-c("Info_sourceOrg_id")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Class = make_class(.data$Info_nPos,
                                     .data$Info_nNeg,
                                     set_positive)) %>%
    dplyr::ungroup()

  # Propagate original attributes from epitopes onto df
  tmp1 <- attributes(epits)
  idx <- which(!(names(tmp1) %in% names(attributes(df))))
  for (i in idx){
    attr(df, names(tmp1)[i]) <- tmp1[[i]]
  }

  # Add this function's attributes
  attr(df, "only_exact") <- only_exact
  attr(df, "set_positive") <- set_positive


  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    saveRDS(df, paste0(normalizePath(save_folder), "/labelled_prots.rds"))
  }

  return(df)
}
