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

  if (length(list_data$Reference$Epitopes) < 1) return(data.frame())

  # ============= ONLY VALID LIST_DATA OBJECTS CROSS THIS LINE ============= #
  list_data$filename <- basename(filename)
  if (type == "B") {
    outlist <- lapply(seq_along(list_data$Reference$Epitopes),
                      FUN  = process_individual_epitope_B,
                      list_data = list_data)
  } else if (type == "T"){
    # TODO: implement T cell epitope processor
    stop("T-cell epitope extraction still under development.")
  } else stop("epitope type not recognised")

  df <- do.call(rbind, outlist)
  return(df)
}
