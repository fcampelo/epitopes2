# Little field check function
nullcheck <- function(x) { ifelse(is.null(x), yes = NA, no = x) }


# ======================================================================
# Progress bar function
mypb <- function(i, max_i, t0, npos){
  nb <- max(1, ceiling(max_i / npos))
  if (i == 0){
    pbstr <- paste0("\n  |", paste(rep("_", npos), collapse = ""), "|")
    cat(pbstr, "0% processed. Elapsed time: 0s")
  } else if (i >= (max_i - 1)) {
    pbstr <- paste(rep(">", times = npos), collapse = "")
    td <- Sys.time() - t0
    perc_done <- 100
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  } else if (!(i %% nb)) {
    nn <- i / nb
    pbstr <- paste(rep(c(">", "_"), times = c(nn, npos - nn)),
                   collapse = "")
    td <- Sys.time() - t0
    perc_done <- round(100 * i / max_i, digits = 0)
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  }
  invisible(NULL)
}

set_mc <- function(ncpus){
  cl <- max(1, min(ncpus, parallel::detectCores() - 1))
  if (cl > 1){
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(ncpus, setup_timeout = 1)
    }
  }
  return(cl)
}

close_mc <- function(cl){
  # Stop cluster
  if("cluster" %in% class(cl)) parallel::stopCluster(cl)
}


mypblapply <- function(X, FUN, ncpus, ...){
  if (ncpus == 1){
    t0 <- Sys.time()
    res <- vector("list", length(X))
    for (i in seq_along(X)){
      mypb(i, length(X), t0, 40)
      res[[i]] <- FUN(X[[i]], ...)
    }
    # res <- lapply(X   = X,
    #               FUN = FUN,
    #               ...)
  } else {

    cl <- set_mc(ncpus)
    if(.Platform$OS.type == "windows"){
      res <- pbapply::pblapply(cl  = cl,
                               X   = X,
                               FUN = FUN,
                               ...)
    } else {
      res <- pbapply::pblapply(cl  = cl,
                               X   = X,
                               FUN = FUN,
                               mc.preschedule = FALSE,
                               ...)
    }
    close_mc(cl)

  }

  return(res)
}
