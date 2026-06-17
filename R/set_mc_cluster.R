#' Set multicore cluster
#'
#' This function is a wrapper around [parallel::makePSOCKcluster()] followed
#' by [parallel::clusterEvalQ()], which loads all required packages in the
#' individual cluster workers.
#'
#' @param ncpus number of cores to use
#' @param setup_timeout parameter passed to [parallel::makePSOCKcluster](parallel::makePSOCKcluster)
#' @param pkgs_to_load packages to load in the workers. Defaults to the packages imported by `epitopes`.
#'
#' @return a SOCK cluster object
#'
#' @export
#'
set_mc_cluster <- function(ncpus,
                           setup_timeout = 2,
                           pkgs_to_load = c("dplyr", "assertthat", "XML",
                                            "rlang", "rentrez", "pbapply",
                                            "protr", "stringr", "BiocManager",
                                            "R.utils", "moses", "seqinr",
                                            "parallel", "utils", "Biostrings")){

  assertthat::assert_that(is.character(pkgs_to_load),
                          is.numeric(setup_timeout),
                          setup_timeout > 0,
                          assertthat::is.count(ncpus))
  cl <- parallel::makePSOCKcluster(names = ncpus, setup_timeout = setup_timeout)
  parallel::clusterExport(cl, varlist = "pkgs_to_load")
  ignore <- parallel::clusterEvalQ(cl = cl,
                                   {
                                     lapply(pkgs_to_load,
                                            \(p) require(p, character.only = TRUE))
                                     invisible(TRUE)
                                   })
  return(cl)
}
