#' # Deprecated functions
#'
#' #' Extract linear B Cell epitopes from XML files retrieved from IEDB.
#' #'
#' #' This function is deprecated. Please use get_LBCE() instead.
#' #'
#' #' @param data_folder path (either relative or absolute) to the directory
#' #'        containing the XML files
#' #' @param save_folder path to folder for saving the results.
#' #'
#' #'
#' get_linear_bcell_epitopes <- function(data_folder = "./",
#'                                       save_folder = NULL){
#'
#'   .Deprecated("get_LBCE()")
#'   # ========================================================================== #
#'   # Sanity checks and initial definitions
#'   assertthat::assert_that(is.character(data_folder),
#'                           length(data_folder) == 1,
#'                           dir.exists(data_folder),
#'                           is.null(save_folder) | (is.character(save_folder)),
#'                           is.null(save_folder) | length(save_folder) == 1)
#'
#'
#'   # =======================================
#'   # Get file list and initialise variables
#'   filelist    <- dir(data_folder, pattern = ".xml", full.names = TRUE)
#'   filelist    <- gsub("//", "/", filelist, fixed = TRUE)
#'
#'   # Check save folder and create file names
#'   if(!is.null(save_folder)) {
#'     if(!dir.exists(save_folder)) dir.create(save_folder)
#'     df_file <- normalizePath(paste0(save_folder, "/epitopes.rds"))
#'     # errfile <- normalizePath(paste0(save_folder, "/xml_load_errors.rds"))
#'   }
#'
#'   # ==================================================
#'   cat("Processing files:\n")
#'   df <- pbmcapply::pbmclapply(X        = filelist,
#'                               FUN      = process_xml_file,
#'                               mc.cores = 1)
#'
#'
#'   cat("\nProcessing resulting list:\n")
#'   df <- pbmcapply::pbmclapply(X        = df,
#'                               FUN      = data.table::rbindlist,
#'                               mc.cores = 1)
#'   df           <- data.frame(data.table::rbindlist(df))
#'
#'   if(!is.null(save_folder)){
#'     saveRDS(object = df,      file = df_file)
#'     # saveRDS(object = errlist, file = errfile)
#'   }
#'
#'   invisible(df)
#' }
#'
#'
#' #' Retrieve protein sequences and data from GenBank
#' #'
#' #' This function is deprecated. Please use get_proteins() instead.
#' #'
#' #' @param uids A list of potein IDs provided a character vector.
#' #' @param save_folder path to folder for saving the results.
#' #'
#' #'
#'
#' retrieve_protein_data <- function(uids, save_folder = NULL){
#'
#'   .Deprecated("get_proteins()")
#'   # ========================================================================== #
#'   # Sanity checks and initial definitions
#'   assertthat::assert_that(is.character(uids),
#'                           length(uids) >= 1,
#'                           is.null(save_folder) | (is.character(save_folder)),
#'                           is.null(save_folder) | length(save_folder) == 1)
#'
#'   # Check save folder and create file names
#'   if(!is.null(save_folder)) {
#'     if(!dir.exists(save_folder)) dir.create(save_folder)
#'     df_file <- normalizePath(paste0(save_folder, "/proteins.rds"))
#'     errfile <- normalizePath(paste0(save_folder, "/proteins_retrieval_errors.rds"))
#'   }
#'
#'   # Retrieving proteins using individual requests rather than (more efficient)
#'   # batch requests, to catch and treat efetch() or parsing errors more easily.
#'   df <- pbmcapply::pbmclapply(X       = uids,
#'                               FUN      = retrieve_single_protein,
#'                               mc.cores = 1)
#'
#'   errlist <- uids[sapply(df, is.null)]
#'   df      <- data.frame(data.table::rbindlist(df,
#'                                               use.names = TRUE,
#'                                               fill      = TRUE))
#'
#'   # Save proteins to file and clean up outfile
#'   if(!is.null(save_folder)){
#'     saveRDS(object = df,      file = df_file)
#'     saveRDS(object = errlist, file = errfile)
#'   }
#'
#'   invisible(list(proteins = df,
#'                  errlist  = errlist))
#' }
