#' Read saved output from Bepipred2
#'
#' This function is used to read the saved CSV output from Bepipred2.
#'
#' @param res.file path to the results file. Must be a CSV file saved from the
#' Bepipred 2.0 output.
#' @param threshold threshold for converting prediction probabilities to
#' positive class predictions.
#' @param ... Currently unused.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

read_bepipred2 <- function(res.file, threshold = 0.5, ...){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(file.exists(res.file))

  # ========================================================================== #

  preds <- utils::read.csv(res.file, header = TRUE, stringsAsFactors = FALSE)
  preds <- preds[, c("Entry", "Position", "EpitopeProbability")]

  names(preds) <- c("Info_UID", "Info_center_pos", "Bepipred2_prob")
  preds$Bepipred2_class <- -1 + 2 * (preds$Bepipred2_prob > threshold)

  return(preds)

}
