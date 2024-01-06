#' Retrieve pre-calculated ESM-1b protein features.
#'
#' This function retrieves pre-calculated ESM-1b features given a
#' vector of protein IDs or a vector of taxonomy IDs.
#'
#' @param Info_protein_id vector of protein IDs.
#' @param Info_organism_id vector of organism taxonomy IDs.
#' @param path string containing the path to the folder where the feature
#' files are located. There must be a tab-separated file in this folder called
#' `filemap.tsv` with columns `Info_protein_id` (protein ids), `folder`
#' (relative subfolder containing that protein) and `Info_organism_id`
#' (comma-separated string with the taxonomy IDs of pathogens associated with
#' the protein in the IEDB records).
#' @param save_folder path to folder for saving the results.
#' @param ncpus positive integer, number of cores to use
#'
#' @return A named list vector containing the data frames with the ESM features
#' for each protein. IDs that are not available in `path` are returned as
#' empty list positions. This list vector has attributes `OK` (vector of
#' ids of proteins successfully retrieved) and `FAIL` (unsuccessful ids).
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_esm_features <- function(path,
                             Info_protein_id = NULL,
                             Info_organism_id = NULL,
                             save_folder = NULL,
                             ncpus = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  ok_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(is.null(save_folder) | (is.character(save_folder)),
                          length(save_folder) <= 1,
                          any(class(Info_protein_id) %in% ok_classes),
                          is.null(Info_protein_id) || length(Info_protein_id) >= 1,
                          any(class(Info_organism_id) %in% ok_classes),
                          is.null(Info_organism_id) || length(Info_organism_id) >= 1,
                          (is.null(Info_organism_id) + is.null(Info_protein_id)) <= 1,
                          is.character(path), length(path) == 1, dir.exists(path),
                          assertthat::is.count(ncpus))

  # Read file map
  filemap <- read.table(paste0(path, "/filemap.tsv"), sep = "\t", header = TRUE)


  # ====== Prepare protein ID list
  if(is.null(Info_protein_id)) Info_protein_id <- character()

  if(!is.null(Info_organism_id)) {
    for (i in seq_along(Info_organism_id)){
      idx <- which(sapply(strsplit(filemap$Info_organism_id, split = ","),
                          function(x, id) {any(x == id)},
                          id = Info_organism_id[i]))
      if (length(idx) > 0) {
        Info_protein_id <- c(Info_protein_id, filemap$Info_protein_id[idx])
      } else {
        message("Warning: No records found for Info_organism_id ", Info_organism_id[i])
      }
    }
  }

  Info_protein_id <- unique(as.character(unlist(Info_protein_id)))


  # ====== Check save folder and create file name
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    df_file <- paste0(normalizePath(save_folder), "/esm_features.rds")
  }

  reslist <- mypblapply(Info_protein_id,
                        function(id, path, filemap){
                          idx <- which(filemap$Info_protein_id == id)
                          fp <- paste0(path, "/", filemap$folder[idx], "/", id, ".rds")
                          if(file.exists(fp)){
                            x <- readRDS(fp)
                            x$Info_protein_id <- id
                            return(cbind(x[, ncol(x), drop = FALSE], x[, -ncol(x)]))
                          } else {
                            return(NULL)
                          }
                        }, path = path, filemap = filemap, ncpus = ncpus)

  names(reslist) <- Info_protein_id

  errs <- unlist(sapply(reslist, is.null))

  attr(reslist, "OK") <- Info_protein_id[which(!errs)]
  attr(reslist, "FAIL") <- Info_protein_id[which(errs)]

  # Save results to file
  if(!is.null(save_folder)){
    saveRDS(object = reslist, file = df_file)
  }

  message("Done!\n", length(reslist) - sum(errs), " protein files retrieved.\n",
          sum(errs), " protein ids not found.")

  return(reslist)
}
