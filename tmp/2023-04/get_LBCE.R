#' Extract linear B Cell epitopes from XML files retrieved from IEDB.
#'
#' This function is used to extract information for *linear B-cell epitopes*
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
#' @param ncpus positive integer, number of cores to use
#' @param save_folder path to folder saving the output.
#'
#' @return A *data.table* containing the epitope data is returned invisibly.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @examples
#' my.dir   <- system.file("extdata/xml_examples", package="epitopes")
#' epitopes <- get_LBCE(my.dir)
#'

get_LBCE <- function(data_folder,
                     ncpus = 1,
                     save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(data_folder), length(data_folder) == 1,
                          dir.exists(data_folder),
                          assertthat::is.count(ncpus),
                          is.null(save_folder) | (is.character(save_folder) && length(save_folder) == 1))

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    df_file <- paste0(normalizePath(save_folder), "/epitopes.rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/epit_errlist.rds")
  }

  # Get file list and initialise variables
  filelist    <- dir(normalizePath(data_folder), pattern = ".xml",
                     full.names = TRUE)

  # ==================================================
  t <- Sys.time()
  cat("Processing", length(filelist), "files using", ncpus, "cores",
      "\nStarted at", as.character(t), "\n")

  df <- mypblapply(ncpus = ncpus,
                   X     = filelist,
                   FUN   = process_xml_file,
                   type  = "B")

  td <- Sys.time() - t
  cat("Ended at", as.character(Sys.time()),
      "\nElapsed time:", signif(as.numeric(td), 3), attr(td, "units"))

  erridx  <- which(sapply(df, function(x) is.character(x) && x == "Error"))
  errlist <- basename(filelist[erridx])
  if(length(erridx) > 0) df <- df[-erridx]

  emptidx <- which(sapply(df, function(x) nrow(x) == 0))
  df      <- data.table::rbindlist(df[-emptidx])

  class(df) <- c(class(df), "LBCE_dt")

  if(!is.null(save_folder)){
    saveRDS(object = df,      file = df_file)
    saveRDS(object = errlist, file = errfile)
  }

  cat("\nDone!\n", nrow(df), "epitopes retrieved.\n",
      length(errlist), "processing errors.")

  invisible(df)
}
