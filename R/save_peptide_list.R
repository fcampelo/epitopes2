#' Save information from a peptide.list object to files .
#'
#' Saves the full information from objects of type `peptide.list` (e.g., output
#' of functions [extract_peptides()], [make_data_splits()] or [calc_features()])
#' to a set of CSV files for easier sharing and cross-platform use.
#'
#' @param peptide.list object of class `peptide.list`.
#' @param save_folder path to folder for saving the resulting CSV files.
#'
#' @return Vector with names of the files saved (invisibly).
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @importFrom utils write.csv
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
  write.csv(peptide.list$df, paste0(save_folder, "/peptides_df_main.csv"),
            row.names = FALSE)
  saved <- c(saved, "peptides_df_main.csv")

  # save summary peptides data frame
  message("Saving summary peptides data frame")
  write.csv(peptide.list$peptides, paste0(save_folder, "/peptides_df_summary.csv"),
            row.names = FALSE)
  saved <- c(saved, "peptides_df_summary.csv")

  # save summary proteins data frame
  if("proteins" %in% names(peptide.list)){
    message("Saving proteins data frame")
    write.csv(peptide.list$proteins, paste0(save_folder, "/proteins_df.csv"),
              row.names = FALSE)
    saved <- c(saved, "proteins_df.csv")
  }

  # save function call parameters
  message("Saving function call parameters")
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
    data.frame(Function  = "extract_peptides",
               Parameter = names(peptide.list$peptide.attrs),
               Value     = unname(unlist(peptide.list$peptide.attrs))))

  # Add split parameters, if available
  if("splits.attrs" %in% names(peptide.list)){
    tmp <- peptide.list$splits.attrs[which(names(peptide.list$splits.attrs) %in% names(formals(epitopes::make_data_splits)))]
    tmp$target_props <- paste(tmp$target_props, collapse = ",")
    fcpars <- rbind(fcpars,
                    data.frame(Function  = "make_data_splits",
                               Parameter = names(tmp),
                               Value     = unname(unlist(tmp))))
  }

  # Add feature calculation parameters, if available
  if("feature.attrs" %in% names(peptide.list)){
    tmp <- peptide.list$feature.attrs[which(names(peptide.list$feature.attrs) %in% names(formals(epitopes::calc_features)))]
    tmp$local.features <- paste(tmp$local.features, collapse = ",")
    tmp$global.features <- paste(tmp$global.features, collapse = ",")
    fcpars <- rbind(fcpars,
                    data.frame(Function  = "calc_features",
                               Parameter = names(tmp),
                               Value     = unname(unlist(tmp))))
  }

  write.csv(fcpars, paste0(save_folder, "/function_call_parameters.csv"),
            row.names = FALSE)
  saved <- c(saved, "function_call_parameters.csv")


  invisible(saved)
}
