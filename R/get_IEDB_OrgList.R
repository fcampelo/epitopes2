#' Download and unzip IEDB Organism List (XML export)
#'
#' This function is used to retrieve the full IEDB Organism List from
#' [IEDB](https://www.iedb.org) and extract it to a target folder.
#'
#' @param url URL of the *.zip* file containing the full IEDB
#'        export (XML).
#' @param save_folder Path to folder for extracting the results. Defaults to
#'        "IEDBOrganismList_yyyymmdd", where *yyyymmdd* is replaced by the current date.
#' @param remove_zip logical flag: should the *.zip* file be deleted after
#'        extraction?
#' @param timeout Timeout (in seconds) to be applied when downloading the IEDB
#'        export. Increase for slow connections (the export is about 500Mb).
#'
#' @author Felipe Campelo (\email{fcampelo@@gmail.com})
#'
#' @return The function returns `TRUE` upon completion.
#'
#' @export
#'

get_IEDB_OrgList <- function(url = "https://iedb.org/downloader.php?file_name=doc/OrganismList.zip",
                     save_folder = NULL,
                     remove_zip  = TRUE,
                     timeout = 3600){
  # ========================================================================== #
  # Sanity checks and initial definitions
  if(is.null(save_folder)) save_folder <- paste0("IEDBOrganismList_",
                                                 gsub("-", "", Sys.Date()))

  assertthat::assert_that(is.character(url), length(url) == 1,
                          is.character(save_folder), length(save_folder) == 1,
                          is.logical(remove_zip), length(remove_zip) == 1,
                          assertthat::is.count(timeout))

  oldTO <- getOption("timeout")
  options(timeout = timeout)

  # Download the file into save_folder and unzip it.
  if(!dir.exists(save_folder)) dir.create(save_folder)
  message("Downloading file:\n")
  utils::download.file(url, destfile = paste0(save_folder, "/IEDBOrganismList.zip"),
                       quiet = FALSE)
  message("Unzipping file into folder: ", save_folder,
          "\n(This may take a while)")
  utils::unzip(paste0(save_folder, "/IEDBOrganismList.zip"), exdir = save_folder)

  if(remove_zip){
    message("Removing ZIP file")
    file.remove(paste0(save_folder, "/IEDBOrganismList.zip"))
  }

  options(timeout = oldTO)

  invisible(TRUE)
}
