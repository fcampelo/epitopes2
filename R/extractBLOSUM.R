extractBLOSUM <- function (x, eig, AABLOSUM62, ...)
{
  ### adapted from package protr to use internal BLOSUM62 matrix (preventing
  ### errors emerging from matrix protr::AABLOSUM62 not being exported to
  ### the search path of package epitopes.)
  if (protr::protcheck(x) == FALSE) {
    stop("x has unrecognized amino acid type")
  }

  k      <- 5
  lag    <- 3
  accmat <- matrix(0, k, nchar(x))

  A           <- eig$vectors
  B           <- eig$values
  rownames(A) <- rownames(AABLOSUM62)

  x.split <- strsplit(x, "")[[1]]
  for (i in 1:nchar(x)) accmat[, i] <- A[x.split[i], 1:k]

  return(protr::acc(accmat, lag))
}
