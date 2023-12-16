#' Read saved output from ABCPred
#'
#' This function is used to read the saved HTML output from ABCPred.
#'
#' @param res.file path to the results file. Must be an HTML file containing the
#' ABCPred output.
#' @param protID ID of the protein represented in `res.file`
#' @param proteins data frame containing protein data related to the predictions
#' in `res.file`. Each position in a protein must be represented by a row in
#' this data frame. Must have at least the columns *Info_UID* (with protein
#' IDs) and *Info_center_pos* (with position on the protein).
#' @param threshold probability threshold for setting the class of a peptide as
#' positive. If `NULL` the routine employs the threshold value used in the
#' ABCPred call that generated the data in `res_file`.
#' @param ... Currently unused.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

read_abcpred <- function(res.file, protID, proteins, threshold = NULL, ...){
  # ========================================================================== #
  # Sanity checks and initial definitions
  protCols <- c("Info_UID", "Info_pos")
  assertthat::assert_that(file.exists(res.file),
                          is.character(protID), length(protID) == 1,
                          is.data.frame(proteins),
                          all(protCols %in% names(proteins)),
                          is.null(threshold) || is.numeric(threshold),
                          is.null(threshold) || length(threshold) == 1)

  # ========================================================================== #

  proteins <- dplyr::select(as.data.frame(proteins), dplyr::all_of(protCols))
  proteins <- proteins[proteins$Info_UID == protID, ]

  preds        <- XML::readHTMLTable(res.file, header = FALSE)[[2]]
  names(preds) <- preds[1, ]
  preds        <- preds[-c(1, which(is.na(preds[, 1]))), -5]
  preds[, -2]  <- lapply(preds[, -2], as.numeric)

  proteins$ABCpred_class <- -1
  proteins$ABCpred_prob  <- 0

  for (j in 1:nrow(preds)){
    stpos <- preds$`Start position`[j]
    idx   <- stpos:(stpos + nchar(preds$Sequence[j]) - 1)
    proteins$ABCpred_prob[idx]  <- pmax(proteins$ABCpred_prob[idx],
                                        preds$Score[j])
    if(is.null(threshold)) {
      proteins$ABCpred_class[idx] <- 1
    } else {
      proteins$ABCpred_class[idx] <- -1 + 2 * as.numeric(proteins$ABCpred_prob[idx] > threshold)
    }
  }

  return(proteins)

}
