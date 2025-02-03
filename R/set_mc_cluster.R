#' Set multicore cluster
#'
#' @param ncpus number of ores to use
#'
#'
set_mc_cluster <- function(ncpus){
  cl <- parallel::makePSOCKcluster(names = ncpus, setup_timeout = 2)
  ignore <- parallel::clusterEvalQ(cl = cl,
                                   {
                                     require(dplyr)
                                     require(assertthat)
                                     require(XML)
                                     require(rlang)
                                     require(reutils)
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
