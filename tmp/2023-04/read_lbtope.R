#' Read saved output from LBtope
#'
#' This function is used to read the saved TXT output from LBtope
#'
#' @param res.file path to the results file. Must be a TXT file saved from the
#' LBtope output.
#' @param threshold threshold for converting prediction probabilities to
#' positive class predictions.
#' @param ... Currently unused.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

read_lbtope <- function(res.file, threshold = 0.5, ...){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(file.exists(res.file))

  # ========================================================================== #

  preds  <- utils::read.csv(res.file, header = FALSE, sep = "\t")
  pos    <- c(which(is.na(preds$V2)), nrow(preds) + 1)
  ids <- sapply(preds$V1[pos],
                function(x){
                  gsub(">", "", strsplit(x, split = " ")[[1]][4])},
                USE.NAMES = FALSE)

  preds$Info_UID <- ""
  preds$Info_center_pos <- 0
  for (i in 1:(length(pos) - 1)){
    idx <- (pos[i] + 1):(pos[i + 1] - 1)
    preds$Info_UID[idx] <- ids[i]
    preds$Info_center_pos[idx] <- seq_along(idx)
  }
  preds <- preds[-pos[1:(length(pos) - 1)], c(4, 5, 3)]
  names(preds)[3] <- "LBtope_prob"
  preds$LBtope_prob <- preds$LBtope_prob / 100
  preds$LBtope_class <- -1 + 2 * (preds$LBtope_prob > threshold)

  return(preds)

}
