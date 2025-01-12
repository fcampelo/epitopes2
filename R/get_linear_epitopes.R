#' Extract linear epitopes from XML files retrieved from IEDB.
#'
#' This function is used to extract information for *linear epitopes*
#' from the XML files exported using the functionality provided by
#' [IEDB](https://www.iedb.org/).
#' It assumes that the user has downloaded the *Complete Database Export* from
#' the *XML Database Export* field in IEDB's
#' [Database Export](https://www.iedb.org/database_export_v3.php) and extracted
#' it in a given folder, which is passed as an argument to the function. This
#' can be easily done with [get_IEDB()].
#'
#' @param data_folder path (either relative or absolute) to the directory
#'        containing the XML files
#' @param epitope_type type of linear epitope to retrieve ("T" or "B")
#' @param ncpus positive integer, number of cores to use
#' @param save_folder path to folder for saving the output.
#'
#' @return A data frame containing the epitope data.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @examples
#' my.dir   <- system.file("extdata/xml_examples", package="epitopes")
#' epitopes <- get_linear_epitopes(my.dir)
#'

get_linear_epitopes <- function(data_folder,
                                epitope_type = "B",
                                ncpus = 1,
                                save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(data_folder),
                          length(data_folder) == 1,
                          dir.exists(data_folder),
                          is.character(epitope_type),
                          length(epitope_type) == 1,
                          epitope_type %in% c("T", "B"),
                          assertthat::is.count(ncpus),
                          is.null(save_folder) | is.character(save_folder),
                          length(save_folder) <= 1)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    df_file <- paste0(normalizePath(save_folder), "/L", epitope_type, "CE.rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/epitopes_retrieval_errlist.rds")
  }

  # Get file list and initialise variables
  filelist    <- dir(normalizePath(data_folder), pattern = ".xml",
                     full.names = TRUE)

  # ==================================================
  t <- Sys.time()
  message("Processing ", length(filelist), " files using ", ncpus, " cores",
          "\nStarted at ", as.character(t), "\n")

  df <- mypblapply(ncpus = ncpus,
                   X     = filelist,
                   FUN   = process_xml_file,
                   type  = epitope_type,
                   toexport = list("process_individual_epitope_T",
                                   "process_individual_epitope_B"))

  td <- Sys.time() - t
  message("Ended at ", as.character(Sys.time()),
          "\nElapsed time: ", signif(as.numeric(td), 3), " ", attr(td, "units"))

  erridx  <- which(sapply(df, function(x) is.character(x) && x == "Error"))
  errlist <- basename(filelist[erridx])
  if(length(erridx) > 0) df <- df[-erridx]

  emptidx <- which(sapply(df, function(x) {is.null(x) || nrow(x) == 0}))
  if(length(emptidx) > 0) df <- df[-emptidx]

  if (length(df) > 0) {
    df <- dplyr::bind_rows(df)
  } else {
    df <- data.frame()
  }

  if(!is.null(save_folder)){
    saveRDS(object = df, file = df_file)
    if(length(erridx) > 0) saveRDS(object = errlist, file = errfile)
  }

  message("Done!\n", nrow(df), " epitopes retrieved.\n",
          length(errlist), " processing errors.")

  return(dplyr::as_tibble(df))
}

#' Exgtract LBCE data from IEDB export
#'
#' Alias function for [get_linear_epitopes()] with type = "B"
#' @inheritParams get_linear_epitopes
#' @export
get_LBCE <- function(data_folder,
                     ncpus = 1,
                     save_folder = NULL){

  get_linear_epitopes(data_folder,
                      epitope_type = "B",
                      ncpus,
                      save_folder)
}

#' Exgtract LTCE data from IEDB export
#'
#' Alias function for [get_linear_epitopes()] with type = "T"
#' @inheritParams get_linear_epitopes
#' @export
get_LTCE <- function(data_folder,
                     ncpus = 1,
                     save_folder = NULL){

  get_linear_epitopes(data_folder,
                      epitope_type = "T",
                      ncpus,
                      save_folder)
}
