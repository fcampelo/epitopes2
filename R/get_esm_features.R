#' Retrieve pre-calculated ESM-1b protein features.
#'
#' This function retrieves pre-calculated ESM-1b features given a
#' list of protein IDs.
#'
#' @param Info_protein_id vector of protein IDs.
#' @param path string containing the path to the folder where the feature
#' files are located. There must be a tab-separated file in this folder called
#' `filemap.tsv` with columns `Info_protein_id` (protein ids), `folder`
#' (relative subfolder containing that protein) and `Info_organism_id`
#' (comma-separated string .
#' @param save_folder path to folder for saving the results.
#'
#' @return A list vector containing the data frames with the ESM features
#' for each protein. IDs that are not available in `path` are returned as
#' empty list positions.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_esm_features <- function(Info_protein_id, path, save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  t0 <- Sys.time()
  ok_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(is.null(save_folder) | (is.character(save_folder)),
                          length(save_folder) <= 1,
                          any(class(uids) %in% ok_classes),
                          length(uids) >= 1,
                          is.character(path),
                          length(path) == 1)

  # Make sure UIDs are unique
  uids <- unique(as.character(unlist(uids)))

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    df_file <- paste0(normalizePath(save_folder), "/esm_features.rds")
  }

  reslist <- vector("list", length = length(uids))
  nerr    <- Inf

  ## STOPPED HERE




  # Save results to file
  if(!is.null(save_folder)){
    saveRDS(object = reslist, file = df_file)
  }

  message("Done!\n", length(reslist), " protein files retrieved.\n",
          length(errlist), " retrieval errors.")

  return(reslist)
}
