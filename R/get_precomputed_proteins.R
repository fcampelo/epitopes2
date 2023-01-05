#' Retrieve protein data with ESM-1b features
#'
#' This function can be used to retrieve pre-computed features for proteins
#' containing at least one peptide labelled as a linear B-cell epitope (or
#' non-LBCE) in IEDB. The features are calculated using the ESM-1b embedder
#' ([https://github.com/facebookresearch/esm](https://github.com/facebookresearch/esm))
#'
#' @param proteinIDs A character vector with protein IDs.
#' @param ntries, number of times the routine will try to retrieve each dataset
#' before giving up
#' @param save_folder path to folder for saving the results.
#'
#' @return Vector of retrieval status
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_precomputed_proteins <- function(proteinIDs,
                                     save_folder,
                                     ntries = 3){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(save_folder), length(save_folder) == 1,
                          is.character(proteinIDs), length(proteinIDs) >= 1,
                          assertthat::is.count(ntries))

  baseURL <- "https://github.com/epitopes-dataset/XXXREPOXXX/raw/main/data/XXXFILEXXX"

  file_map <- utils::read.csv(url("https://github.com/epitopes-dataset/ESM1b_IEDB_LBCE_proteins_1/raw/main/file_map.csv"),
                              quote = "")
  file_map$filename <- paste0(file_map$protID, ".rds")
  file_map$protID <- gsub("\\_part..$", "", file_map$protID)

  # Check save folder
  if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)

  resvec <- rep(FALSE, length(proteinIDs))
  names(resvec) <- proteinIDs

  # Try retrieving from github repo
  for (i in seq_along(proteinIDs)){
    idx <- file_map[grep(proteinIDs[i], file_map$protID, fixed = TRUE), ]
    x <- NA
    if(nrow(idx) > 0){
      idx <- idx[order(idx$filename), ]
      message("Retrieving data for protein ID ", proteinIDs[i], " ", appendLF = FALSE)
      cc <- 0
      while (cc < ntries){
        errk <- FALSE
        tryCatch({
          for (k in 1:nrow(idx)){
            furl <- gsub("XXXREPOXXX", idx$which.repo[k], baseURL)
            furl <- gsub("XXXFILEXXX", idx$filename[k], furl)
            tmp <- readRDS(url(furl))
            if (k == 1){
              x <- tmp
            } else {
              x <- rbind(x, tmp)
            }
          }
        },
        warning = function(c) {errk <<- TRUE},
        error   = function(c) {errk <<- TRUE},
        finally = NULL)

        if(errk){
          cc <- cc + 1
          if(cc == ntries) message("- unable to retrieve.")
        } else {
          saveRDS(x, paste0(save_folder, "/", proteinIDs[i], ".rds"))
          resvec[i] <- TRUE
          message(" done!")
          break
        }
      }
    } else {
      message("No precomputed data found for protein ID ", proteinIDs[i], ". Skipping.")
    }
  }

  invisible(resvec)
}
