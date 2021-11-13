#' Save information from a peptide.list object to files .
#'
#' Saves the full information from objects of type `peptide.list` (e.g., output
#' of functions [extract_peptides()], [make_data_splits()] or [calc_features()])
#' to a set of CSV files for easier sharing and cross-platform use.
#'
#' @param peptide.list object of class `peptide.list`.
#' @param save_folder path to folder for saving the resulting CSV files.
#'
#' @return Vector of file names generated (invisibly).
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export

save_peptide_list <- function(peptide.list, save_folder){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that("peptide.list" %in% class(peptide.list),
                          is.character(save_folder),
                          length(save_folder) == 1)

  if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
  # ========================================================================== #

  saved <- character()

  # save main data frame
  message("Saving main peptides data frame")
  write.csv(peptide.list$df, paste0(save_folder, "peptides_df_main.csv"),
            row.names = FALSE)
  saved <- c(saved, "peptides_df_main.csv")

  # save summary peptides data frame
  message("Saving summary peptides data frame")
  write.csv(peptide.list$peptides, paste0(save_folder, "peptides_df_summary.csv"),
            row.names = FALSE)
  saved <- c(saved, "peptides_df_summary.csv")

  # save function call parameters
  fcpars <- rbind(
    # filter attributes
    data.frame(Function  = "filter_epitopes",
               Parameter = names(peptide.list$filter.attrs),
               Value     = unname(unlist(peptide.list$filter.attrs))),
    # consolidate attributes
    data.frame(Function  = "consolidate_data",
               Parameter = names(peptide.list$consolidate.attrs),
               Value     = unname(unlist(peptide.list$consolidate.attrs))),
    # peptide extraction attributes
    data.frame(Function  = "peptide.attrs",
               Parameter = names(peptide.list$peptide.attrs),
               Value     = unname(unlist(peptide.list$peptide.attrs))))



  message("Saving function call parameters")
  write.csv(peptide.list$peptides, paste0(save_folder, "peptides_df_summary.csv"),
            row.names = FALSE)
  saved <- c(saved, "peptides_df_summary.csv")



}
