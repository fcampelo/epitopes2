#' Extract peptides and windowed representation for epitope prediction.
#'
#' Extract relevant peptides from the consolidated protein-epitope data
#' (generated using [consolidate_data()]) and build column with neighbourhood
#' of each labeled aminoacid residue.
#'
#' @param df data frame of consolidated protein-epitope data, returned by
#' [consolidate_data()].
#' @param min_peptide positive integer, shortest peptide to be considered
#' (applies to both **positive** and **negative** observations)
#' @param max_epitope positive integer, longest peptide to be considered (only
#' applies to **positive** observations.)
#' @param window_size positive integer, size of the local neighbourhood to be
#' considered.
#' @param save_folder path to folder for saving the results. It will save the
#' results as file *peptides_list.rds* (overwriting if necessary)
#'
#' @return List containing two data frames:
#'
#' \itemize{
#'    \item **df**:  data frame containing the labeled positions of `df`
#'    (filtered by *min_peptide* and *max_epitope*) with a new column
#'    **Info_window** containing the local neighbourhood of each position
#'    (to be used for feature calculation)
#'    \item **peptides**: data frame containing one peptides per row, together
#'    with its calculated Class attribute value.
#' }
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data

extract_peptides <- function(df,
                             min_peptide = 8, max_epitope = 30,
                             window_size = (2 * min_peptide) - 1,
                             save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(df),
                          assertthat::is.count(min_peptide),
                          assertthat::is.count(max_epitope),
                          min_peptide <= max_epitope,
                          assertthat::is.count(window_size),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)

  my.attrs <- attributes(df)
  filter.attrs <- my.attrs[names(my.attrs) %in% c("orgIDs", "hostIDs", "removeIDs")]
  consolidate.attrs <- my.attrs[names(my.attrs) %in% c("only_exact", "set_positive")]

  # ========================================================================== #
  # Identify contiguous labelled peptides in each protein
  message("Identifying contiguous labeled regions...")
  df <- df %>%
    dplyr::group_by(.data$Info_protein_id) %>%
    dplyr::mutate(IsBreak = find_breaks(.data$Class),
                  Info_PepID   = ifelse(is.na(.data$Class), NA,
                                        paste0(.data$Info_protein_id, ":",
                                               cumsum(.data$IsBreak)))) %>%
    dplyr::group_by(.data$Info_PepID) %>%
    dplyr::mutate(Info_peptide_length = (!is.na(.data$Class)) * dplyr::n()) %>%
    dplyr::ungroup()


  # Extract individual contiguous peptides of length between min_peptide and
  # max_peptide
  message("Extracting labelled peptides...")
  peptides <- df %>%
    dplyr::filter(!is.na(.data$Class),
                  .data$Info_peptide_length >= min_peptide,
                  (.data$Class == -1) | (.data$Info_peptide_length <= max_epitope)) %>%
    dplyr::group_by(.data$Info_PepID) %>%
    dplyr::summarise(Info_organism_id    = dplyr::first(.data$Info_organism_id),
                     Info_host_id        = dplyr::first(.data$Info_host_id),
                     Info_protein_id     = dplyr::first(.data$Info_protein_id),
                     Info_start_pos      = dplyr::first(.data$Info_pos),
                     Info_end_pos        = dplyr::last(.data$Info_pos),
                     Info_peptide        = paste(.data$Info_AA, collapse = ""),
                     Info_peptide_length = dplyr::n(),
                     Class               = dplyr::first(.data$Class),
                     .groups = "drop")


  # Update df to include local neighbourhood and remove NA regions and those
  # that do not comply with min_peptide and max_epitope
  message("Extracting windows...")
  df <- df %>%
    dplyr::group_by(.data$Info_protein_id) %>%
    dplyr::mutate(Info_window = make_windows(.data$Info_AA,
                                             .data$Class,
                                             window_size)) %>%
    dplyr::filter(!is.na(.data$Class),
                  .data$Info_peptide_length >= min_peptide,
                  (.data$Class == -1) | (.data$Info_peptide_length <= max_epitope)) %>%
    dplyr::select(.data$Info_PepID, dplyr::everything(),
                  -.data$IsBreak, -.data$Info_peptide_length, .data$Class) %>%
    dplyr::ungroup()


  outlist <- list(df                = df,
                  peptides          = peptides,
                  filter.attrs      = filter.attrs,
                  consolidate.attrs = consolidate.attrs,
                  peptide.attrs     = list(min_peptide = min_peptide,
                                           max_epitope = max_epitope,
                                           window_size = window_size))

  class(outlist) <- unique(c(class(outlist), "peptide.list"))

  # Check save folder and save files
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    saveRDS(outlist, paste0(normalizePath(save_folder), "/peptides_list.rds"))
  }

  return(outlist)
}


