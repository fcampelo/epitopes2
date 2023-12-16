.onAttach <- function(...) {
  op2 <- options()
  op2$pboptions_oldchar   <- op2$pboptions$char
  op2$pboptions_olduse_lb <- op2$pboptions$use_lb
  op2$pboptions$char   <- ">"
  op2$pboptions$use_lb <- TRUE
  options(op2)
}


.onDetach <- function(...) {
  op2 <- options()

  op2$pboptions$char      <- op2$pboptions_oldchar
  op2$pboptions$use_lb    <- op2$pboptions_olduse_lb
  op2$pboptions_oldchar   <- NULL
  op2$pboptions_olduse_lb <- NULL

  options(op2)
}
