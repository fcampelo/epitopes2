# Little field check function
nullcheck <- function(x) { ifelse(is.null(x), yes = NA, no = x) }

# function to extract unique entries from comma-separated character strings
# and return them as comma-separated character strings
get_uniques <- function(x){
  sapply(x,
         function(y) {
           paste(unname(unique(unlist(strsplit(y, split = ",")))), collapse = ",")
         })
}

# Function to find points where a vector changes values
find_breaks <- function(x){
  x[is.na(x)] <- Inf
  xl <- dplyr::lag(x, default = -Inf)
  return(x != xl)
}

# Function to build local neighbourhoods for non-NA peptides
make_windows <- function(x, Class, window_size){
  imax    <- length(x)
  noNA    <- which(!is.na(Class))
  windows <- rep(NA, length(x))
  windows[noNA] <- sapply(noNA,
                          function(y){
                            idx <- (y - floor(window_size/2)):(y + floor(window_size/2))
                            idx[which(idx <= 0)] <- 2 - idx[which(idx <= 0)]
                            idx[which(idx > imax)] <- 2 * imax - idx[which(idx > length(x))]
                            paste(x[idx], collapse = "")
                          })
  return(windows)
}

set_mc <- function(ncpus){
  if (ncpus > 1 && .Platform$OS.type == "windows"){
    cl <- parallel::makeCluster(ncpus, setup_timeout = 2)
  } else {
    cl <- max(1, min(ncpus, parallel::detectCores() - 1))
  }
  return(cl)
}

close_mc <- function(cl){
  # Stop cluster
  if("cluster" %in% class(cl)) parallel::stopCluster(cl)
  invisible(TRUE)
}


mypblapply <- function(X, FUN, ncpus, toexport = list(), ...){
  cl  <- set_mc(ncpus)

  if(ncpus > 1 && length(toexport) > 0 && .Platform$OS.type == "windows"){
    parallel::clusterExport(cl = cl,
                            varlist = toexport)
  }
  res <- pbapply::pblapply(cl = cl, X = X, FUN = FUN, ...)

  close_mc(cl)
  return(res)
}



# ======================================================================
# Progress bar function
mypb <- function(i, max_i, t0, npos){
  nb <- max(1, ceiling(max_i / npos))
  if (i == 0){
    pbstr <- paste0("  |", paste(rep("_", npos), collapse = ""), "|")
    cat(pbstr, "0% processed. Elapsed time: 0s")
  } else if (i >= (max_i - 1)) {
    pbstr <- paste(rep(">", times = npos), collapse = "")
    td <- Sys.time() - t0
    perc_done <- 100
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  } else if (!(i %% nb)) {
    nn <- i / nb
    pbstr <- paste(rep(c("+", "_"), times = c(nn, npos - nn)),
                   collapse = "")
    td <- Sys.time() - t0
    perc_done <- round(100 * i / max_i, digits = 0)
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  }
  invisible(NULL)
}





# ======================================================================
# auxiliary functions for optimise_splits()

# Objective function
objfun <- function(x, alpha, Y, Nstar, ...){
  tmp <- getstats(x, Y, Nstar)
  sum(alpha * (tmp$Gj - Nstar)^2 + (1 - alpha) * (tmp$pj - tmp$Pstar)^2)
}

# Auxiliary function for OF
getstats <- function(x, Y, Nstar){
  Pstar <- sum(Y$nPos) / sum(Y$N)
  ymatr <- matrix(round(x), ncol = length(Nstar), nrow = length(x), byrow = FALSE)
  ymatr <- ymatr == matrix(seq_along(Nstar), ncol = length(Nstar), nrow = length(x), byrow = TRUE)
  Gj    <- colSums(ymatr * Y$N) / sum(Y$N)
  pj    <- colSums(ymatr * Y$nPos) / (colSums(ymatr * Y$N) + 1e-12)
  return(list(Gj = Gj, pj = pj, Pstar = Pstar))
}

# Movement function
neighbour <- function(x, Nstar, Y, ...){
  # Cast x as an allocation list
  xl <- lapply(seq_along(Nstar), function(i){seq_along(x)[x == i]})
  xl_movable <- lapply(xl,
                       function(x){x[x %in% Y$Cluster[is.na(Y$split)]]})

  # randomize which neighbourhood to use:
  neighs <- c("taskmove", "swap")
  move   <- sample(neighs, 1)

  if (move == "taskmove"){
    if (length(which(sapply(xl_movable, length) > 1)) == 0){
      move <- "swap"
    } else {
      tryCatch({
        from <- sample(which(sapply(xl_movable, length) > 1), 1)
        from <- c(from, xl_movable[[from]][sample.int(length(xl_movable[[from]]), 1)])
        possible <- (1:length(xl))[-from[1]]
        if(length(possible) > 1) to <- sample(possible, 1) else to <- possible

        xl[[to]]   <- c(xl[[to]], from[2])
        xl[[from[1]]] <- xl[[from[1]]][-which(xl[[from[1]]] == from[2])]
      },
      error   = function(c) {},
      finally = NULL)
    }
  }

  if (move == "swap"){
    groups <- sample(which(sapply(xl_movable, length) > 0), 2)
    if(length(groups) == 2 && !any(sapply(xl_movable, length)==0)){
      ids    <- c(xl_movable[[groups[1]]][sample.int(length(xl_movable[[groups[1]]]), 1)],
                  xl_movable[[groups[2]]][sample.int(length(xl_movable[[groups[2]]]), 1)])
      xl[[groups[1]]] <- c(xl[[groups[1]]], ids[2])
      xl[[groups[1]]] <- xl[[groups[1]]][-which(xl[[groups[1]]] == ids[1])]
      xl[[groups[2]]] <- c(xl[[groups[2]]], ids[1])
      xl[[groups[2]]] <- xl[[groups[2]]][-which(xl[[groups[2]]] == ids[2])]
    }
  }

  # Cast xl back to vector format
  xnew <- unlist(xl)
  names(xnew) <- unlist(mapply(rep, seq_along(xl), sapply(xl, length)))
  xnew <- as.numeric(names(xnew)[order(xnew)])
  names(xnew) <- names(x)
  tmp <- names(xnew)
  if(any(is.na(tmp))) xnew <- xnew[-which(is.na(tmp))]

  # if(length(xnew) != length(x)) {
  #   errl <- list(x=x, Nstar=Nstar, Y=Y,
  #                xl = xl, xl_movable = xl_movable,
  #                xnew = xnew, move = move)
  #   if (move == "taskmove") {
  #     errl$from = from
  #     errl$to = to
  #   } else {
  #     errl$groups = groups
  #     errl$ids = ids
  #   }
  #   saveRDS(errl, "tmp.rds")
  # }

  return(xnew)
}

# Constructive Heuristic
makesol <- function(alpha, Y, Nstar){
  P        <- sum(Y$nPos) / sum(Y$N)
  x        <- ifelse(is.na(Y$split), 0, Y$split)
  names(x) <- Y$Cluster
  Y2       <- Y[is.na(Y$split), ]

  Cap   <- (Nstar * sum(Y$N))
  for (i in seq_along(Cap)){
    Cap[i] <- Cap[i] - sum(Y$N[Y$split == i], na.rm = TRUE)
  }

  while(nrow(Y2) > 0){
    # Get split with largest capacity:
    split.idx <- which.max(Cap)
    # Check which allocation would result in largest objfun reduction
    tmpal <- c(0, 0, Inf)
    for (i in seq_along(Y2$Cluster)){
      tmpx <- x
      tmpx[which(names(tmpx)==Y2$Cluster[i])] <- split.idx
      tmpy <- objfun(tmpx, alpha, Y, Nstar)
      if (tmpy < tmpal[3]){
        tmpal <- c(i, Y2$Cluster[i], tmpy)
      }
    }
    x[names(x)==tmpal[2]]   <- split.idx
    Cap[split.idx] <- Cap[split.idx] - Y2$N[tmpal[1]]
    Y2 <- Y2[-tmpal[1], ]
  }

  tmp <- names(x)
  if(any(is.na(tmp))) x <- x[-which(is.na(tmp))]
  return(x)
}
