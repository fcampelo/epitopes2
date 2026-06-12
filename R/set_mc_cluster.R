#' Set multicore cluster
#'
#' This function is a wrapper around [parallel::makePSOCKcluster()] followed
#' by [parallel::clusterEvalQ()], which loads all required packages in the
#' individual cluster workers.
#'
#' @param ncpus number of cores to use
#'
#' @return a SOCK cluster object
#'
#' @export
#'
set_mc_cluster <- function(ncpus){
  cl <- parallel::makePSOCKcluster(names = ncpus, setup_timeout = 2)
  ignore <- parallel::clusterEvalQ(cl = cl,
                                   {
                                     require(dplyr)
                                     require(assertthat)
                                     require(XML)
                                     require(rlang)
                                     require(rentrez)
                                     require(pbapply)
                                     require(protr)
                                     require(stringr)
                                     require(BiocManager)
                                     require(R.utils)
                                     require(moses)
                                     require(seqinr)
                                     require(parallel)
                                     require(utils)
                                     require(Biostrings)
                                     invisible(NULL)
                                   })
  return(cl)
}
