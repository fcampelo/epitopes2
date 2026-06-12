#' Download and unzip IEDB database (XML export)
#'
#' This function is used to retrieve the full IEDB database from
#' [IEDB](https://www.iedb.org) and extract it to a target folder.
#'
#' @param url URL of the *.zip* file containing the full IEDB
#'        export (XML).
#' @param save_folder Path to folder for extracting the results. Defaults to
#'        "iedb_yyyymmdd", where *yyyymmdd* is replaced by the current date.
#' @param remove_zip logical flag: should the *.zip* file be deleted after
#'        extraction?
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_IEDB <- function(url = "https://www.iedb.org/downloader.php?file_name=doc/iedb_export.zip",
                     save_folder = NULL,
                     remove_zip  = TRUE){
  # ========================================================================== #
  # Sanity checks and initial definitions
  if(is.null(save_folder)) save_folder <- paste0("iedb_",
                                                 gsub("-", "", Sys.Date()))

  assertthat::assert_that(is.character(url), length(url) == 1,
                          is.character(save_folder), length(save_folder) == 1,
                          is.logical(remove_zip), length(remove_zip) == 1)

  # Download the file into save_folder and unzip it.
  if(!dir.exists(save_folder)) dir.create(save_folder)
  cat("\nDownloading file:\n")
  utils::download.file(url, destfile = paste0(save_folder, "/iedb.zip"),
                       quiet = FALSE)
  cat("Unzipping file into folder:", save_folder,
      "\n(This may take a while)")
  utils::unzip(paste0(save_folder, "/iedb.zip"), exdir = save_folder)

  if(remove_zip){
    cat("\nRemoving ZIP file")
    file.remove(paste0(save_folder, "/iedb.zip"))
  }


  invisible(TRUE)
}
