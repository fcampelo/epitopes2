#' Read saved output from iBCE-EL
#'
#' This function is used to read the saved HTML output from ABCPred.
#'
#' @param res.file path to the results file. Must be an HTML file containing the
#' iBCE-EL output.
#' @param split_string short string used to separate protein ID from position in
#' the iBCE-EL output. This assumes that the iBCE-EL results were generated
#' from a FASTA file in which each peptide had the ID string set as
#' *proteinID<split_string>position*.
#' @param ... Currently unused.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

read_ibceel <- function(res.file, split_string = "pp", ...){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(file.exists(res.file))

  # ========================================================================== #

  preds  <- XML::readHTMLTable(res.file, header = TRUE)[[1]]
  idvars <- strsplit(preds$`FASTA ID`, split = split_string)

  preds$Info_UID        <- sapply(idvars, function(x) x[1])
  preds$Info_UID        <- gsub("-", "_", preds$Info_UID, fixed = TRUE)
  preds$Info_center_pos <- sapply(idvars, function(x) x[2])
  preds$`iBCE-EL_prob`  <- pmax(0, pmin(1, as.numeric(preds$Prob)))
  preds$`iBCE-EL_class` <- -1 + 2 * as.numeric(preds$`PIP or Non-PIP` == "BCE")
  preds$Info_center_pos <- as.numeric(preds$Info_center_pos)
  preds <- preds[, c("Info_UID", "Info_center_pos",
                     "iBCE-EL_prob", "iBCE-EL_class")]

  return(preds)

}
