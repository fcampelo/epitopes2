#' Read saved output from SVMtrip
#'
#' This function is used to read the saved HTML output from SVMtrip
#'
#' @param res.file path to the results file. Must be an HTML file containing the
#' SVMtrip output.
#' @param protID ID of the protein represented in `res.file`
#' @param proteins data frame containing protein data related to the predictions
#' in `res.file`. Each position in a protein must be represented by a row in
#' this data frame. Must have at least the columns *Info_UID* (with protein
#' IDs) and *Info_center_pos* (with position on the protein).
#' @param ... Currently unused.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

read_svmtrip <- function(res.file, protID, proteins, ...){
  # ========================================================================== #
  # Sanity checks and initial definitions
  protCols <- c("Info_UID", "Info_center_pos")
  assertthat::assert_that(file.exists(res.file),
                          is.character(protID), length(protID) == 1,
                          is.data.frame(proteins),
                          all(protCols %in% names(proteins)))

  # ========================================================================== #

  proteins <- dplyr::select(as.data.frame(proteins), dplyr::all_of(protCols))
  proteins <- proteins[proteins$Info_UID == protID, ]

  proteins$SVMtrip_prob  <- 0
  proteins$SVMtrip_class <- -1

  preds <- XML::readHTMLTable(res.file, header = TRUE)[[1]]

  if(!is.null(preds)){
    pos <- strsplit(preds$Location, split = " - ")
    for (j in seq_along(pos)){
      proteins$SVMtrip_prob[pos[[j]][1]:pos[[j]][2]]  <- as.numeric(preds$Score[j])
      proteins$SVMtrip_class[pos[[j]][1]:pos[[j]][2]] <- 1
    }
  }

  return(proteins)
}
