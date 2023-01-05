#' Retrieve basic pre-computed datasets for linear B-cell epitope prediction
#'
#' This function can be used to retrieve pre-computed base data sets
#' extracted from the Immune Epitope DataBase (IEDB). The datasets are:
#' * epitopes.rds (containing consolidated LBCE data)
#' * proteins.rds (containing all proteins related to epitopes.rds)
#' * taxonomy.rds (containing all taxonomic relations from the organisms listed
#' in epitopes.rds)
#' * protein_dissimilarity.rds (normalised local dissimilarity matrix for all
#' proteins in proteins.rds, calculated using [Diamond](https://github.com/bbuchfink/diamond)).
#'
#' @param datasets character vector with the dataset(s) to be retrieved. Accepts
#' any combination of "epitopes", "proteins", "taxonomy" or "protein_dissimilarity".
#' @param timeout maximum download time per dataset.
#' @param save_folder path to folder for saving the results.
#'
#' @return This function returns `TRUE` (invisibly) upon completion.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_basic_data <- function(save_folder,
                           datasets = c("epitopes", "proteins",
                                        "taxonomy", "protein_dissimilarity"),
                           timeout = 120){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.null(save_folder) | is.character(save_folder),
                          length(save_folder) <= 1,
                          is.character(datasets),
                          all(datasets %in% c("epitopes", "proteins",
                                              "taxonomy", "protein_dissimilarity")))

  oldTO <- as.numeric(options("timeout"))
  options(timeout = timeout)

  baseURL <- "https://raw.githubusercontent.com/epitopes-dataset/ESM1b_IEDB_LBCE/main/basic_datasets/"

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
  }

  for (i in seq_along(datasets)){
    fn <- paste0(baseURL, datasets[i], ".rds")
    utils::download.file(url = fn, destfile = paste0(save_folder, "/", datasets[i], ".rds"))
  }

  options(timeout = oldTO)
  invisible(TRUE)
}

