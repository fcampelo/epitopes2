#' Process taxonomic information from IEDB OrganismList XML.
#'
#'
#' @param file_path path (either relative or absolute) to the OrganismList.xml
#' file retrieved from the IEDB database export (e.g., via [get_IEDB_OrgList()]).
#' @param cl a SOCK cluster object created using [epitopes::set_mc_cluster()], or
#'        `NULL` if parallel processing is not desired.
#' @param save_folder path to folder for saving the output.
#'
#' @return A data frame containing the taxonomic data.
#'
#' @author Felipe Campelo (\email{fcampelo@@gmail.com})
#'
#' @export
#'

process_IEDB_OrgList <- function(file_path,
                                 cl = NULL,
                                 save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(file_path),
                          length(file_path) == 1,
                          file.exists(file_path),
                          is.null(cl) | "SOCKcluster" %in% class(cl),
                          is.null(save_folder) | is.character(save_folder),
                          length(save_folder) <= 1)

  # Check save folder and create file names
  if(!is.null(save_folder) && !dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }

  message("Importing XML file (", gsub("\\..+$", "", Sys.time()), ")")

  # Initialise error flag
  errk <- FALSE

  # Load XML file
  tryCatch({
    invisible(utils::capture.output(
      xmlfile <- XML::xmlParse(file_path)))},
    error   = function(c) {errk <<- TRUE},
    warning = function(c) {errk <<- TRUE},
    finally = NULL)

  if (errk) return("Error importing XML file")

  df <- XML::xmlToDataFrame(xmlfile,
                            colClasses = rep("character", 4),
                            nodes = XML::getNodeSet(xmlfile, "/OrganismList/Organism"))

  message("Reading nodes (", gsub("\\..+$", "", Sys.time()), ")")
  nodes <- XML::getNodeSet(xmlfile, "/OrganismList/Organism")
  df    <- pbapply::pblapply(nodes,
                             function(c){
                               c <- XML::getChildrenStrings(c)
                               data.frame(as.list(c))},
                             cl = cl)

  # ========stopped here


  # td <- Sys.time() - t
  # message("Ended at ", as.character(Sys.time()),
  #         "\nElapsed time: ", signif(as.numeric(td), 3), " ", attr(td, "units"))
  #
  # erridx  <- which(sapply(df, function(x) is.character(x) && x == "Error"))
  # errlist <- basename(filelist[erridx])
  # if(length(erridx) > 0) df <- df[-erridx]
  #
  # emptidx <- which(sapply(df, function(x) {is.null(x) || nrow(x) == 0}))
  # if(length(emptidx) > 0) df <- df[-emptidx]
  #
  # if (length(df) > 0) {
  #   df <- dplyr::bind_rows(df)
  # } else {
  #   df <- data.frame()
  # }
  #
  # if(!is.null(save_folder)){
  #   saveRDS(object = df, file = df_file)
  #   if(length(erridx) > 0) saveRDS(object = errlist, file = errfile)
  # }
  #
  # message("Done!\n", nrow(df), " epitopes retrieved.\n",
  #         length(errlist), " processing errors.")
  #
  # return(dplyr::as_tibble(df))
}
