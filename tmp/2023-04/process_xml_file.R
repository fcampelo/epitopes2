process_xml_file <- function(filename, type = "B", ...){

  if(!file.exists(filename)) return(NULL)

  # Initialise error flag
  errk <- FALSE

  # Load XML file as a list
  tryCatch({
    invisible(utils::capture.output(
      list_data <- XML::xmlToList(XML::xmlParse(filename))))},
    error   = function(c) {errk <<- TRUE},
    warning = function(c) {errk <<- TRUE},
    finally = NULL)

  if (errk) return("Error")

  if (length(list_data$Reference$Epitopes) < 1) return(data.table::data.table())

  # ============= ONLY VALID LIST_DATA OBJECTS CROSS THIS LINE ============= #
  list_data$filename <- basename(filename)
  if (type == "B") {
    outlist <- lapply(seq_along(list_data$Reference$Epitopes),
                      FUN  = process_individual_epitope_B,
                      list_data = list_data)
  } else if (type == "T"){
    stop("T-cell epitope extraction still under development.")
    #outlist <- lapply(list_data$Reference$Epitopes,
    #                  FUN = process_individual_epitope_T,
    #                  file_id = nullcheck(list_data$Reference$ReferenceId))
  } else stop("epitope type not recognised")

  return(data.table::rbindlist(outlist))
}
