#  Solve the following optimization problem:
#  Let:
#  \itemize{
#     \item nC:  number of clusters.
#     \item K:   number of splits.
#     \item xi:  integer variable defining the split to which cluster i is allocated.
#     \item Ni+: number of _positive_ observations in cluster i.
#     \item Ni:  total number of observations in cluster i.
#     \item Gj*: desired proportion of data for split j.
#     \item P*:  proportion of _positive_ observations in the whole data.
#  }
#
#  The problem is:
#
#  `minimize sum_j{ alpha x (Gj - Gj*)^2 + (1-alpha) x (pj - P*)^2 }`
#
#  With:
#  \itemize{
#      \item `xi \in {1, ..., K}`, for all `i = 1, ..., nC`
#      \item `yij = ifelse(xi == j, 1, 0)`
#      \item `Gj = sum_i{ yij * Ni } / sum_i{ Ni }`
#      \item `pj = sum_i{ yij * Ni+ } / sum_i{ yij * Ni }`
# }

optimise_splits <- function(Y, Nstar, alpha, SAopts, ncpus, id_force_splitting){
  # TODO: generalise to more than 2 classes.

  # === Initial definitions === #
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
  neighbour <- function(x, Nstar, ...){
    # Cast x as an allocation list
    xl <- lapply(seq_along(Nstar), function(i){seq_along(x)[x == i]})

    # randomize which neighbourhood to use:
    neighs <- c("taskmove", "swap", "chain")
    move    <- sample(neighs, 1)

    split.ord <- cbind(seq_along(Nstar),
                       sapply(xl, function(y){sample(seq_along(y), 1)}))
    split.ord <- split.ord[sample(seq_along(Nstar)), ]
    if (move == "taskmove"){
      if(length(xl[[split.ord[1, 1]]]) < 2){
        # If move would result in an empty split, change to "swap"
        move <- "swap"
      } else {
        xl[[split.ord[2, 1]]] <- c(xl[[split.ord[2, 1]]], xl[[split.ord[1, 1]]][split.ord[1, 2]])
        xl[[split.ord[1, 1]]] <- xl[[split.ord[1, 1]]][-split.ord[1, 2]]
      }
    }

    if (move == "swap"){
      xl[[split.ord[2, 1]]] <- c(xl[[split.ord[2, 1]]], xl[[split.ord[1, 1]]][split.ord[1, 2]])
      xl[[split.ord[1, 1]]] <- c(xl[[split.ord[1, 1]]], xl[[split.ord[2, 1]]][split.ord[2, 2]])
      xl[[split.ord[1, 1]]] <- xl[[split.ord[1, 1]]][-split.ord[1, 2]]
      xl[[split.ord[2, 1]]] <- xl[[split.ord[2, 1]]][-split.ord[2, 2]]
    }

    if (move == "chain"){
      split.ord <- rbind(split.ord, split.ord[1, ])
      for(i in 1:(nrow(split.ord) - 1)){
        xl[[split.ord[i+1, 1]]] <- c(xl[[split.ord[i+1, 1]]], xl[[split.ord[i, 1]]][split.ord[i, 2]])
        xl[[split.ord[i, 1]]] <- xl[[split.ord[i, 1]]][-split.ord[i, 2]]
      }
    }

    # Cast xl back to vector format
    xnew <- unlist(xl)
    names(xnew) <- unlist(mapply(rep, seq_along(xl), sapply(xl, length)))
    xnew <- as.numeric(names(xnew)[order(xnew)])

    return(xnew)
  }

  # Constructive Heuristic
  makesol <- function(alpha, Y, Nstar){
    P     <- sum(Y$nPos) / sum(Y$N)
    Cap   <- Nstar * sum(Y$N)
    x0    <- numeric(nrow(Y))
    Y2    <- Y
    while(nrow(Y2) > 0){
      # Get split with largest capacity:
      split.idx <- which.max(Cap)
      # Check which allocation would result in largest objfun reduction
      tmpal <- c(0, 0, Inf)
      for (i in seq_along(Y2$Cluster)){
        tmpx <- x0
        tmpx[Y2$Cluster[i]] <- split.idx
        tmpy <- objfun(tmpx, alpha, Y, Nstar)
        if (tmpy < tmpal[3]){
          tmpal <- c(i, Y2$Cluster[i], tmpy)
        }
      }
      x0[tmpal[2]]   <- split.idx
      Cap[split.idx] <- Cap[split.idx] - Y2$N[tmpal[1]]
      Y2 <- Y2[-tmpal[1], ]
    }
    return(x0)
  }

  # === Run optimisation === #

  if(!is.null(id_force_splitting)){
    idx <- unique(unlist(sapply(id_force_splitting,
                                function(x,ids){grep(x,ids)},
                                ids = Y$txids)))

    Yp <- Y[idx, ]
    Np <- Nstar[1:min(length(Nstar), length(idx))]
    Np <- Np / sum(Np)

    x <- optimise_splits(Yp, Np, alpha, SAopts, ncpus, NULL)
  }

  # TODO: stopped here

  message("Optimising splits with alpha = ", alpha,
          "\n(Number of possibilities: ~", signif(length(Nstar) ^ nrow(Y), 3), ")")
  if(length(Nstar) ^ nrow(Y) < 1e6){
    message("Running exhaustive search...")
    # If the search space is small enough, exhaustive search suffices
    states <- do.call(expand.grid,
                      lapply(1:nrow(Y), function(x){1:length(Nstar)}))
    states <- as.list(as.data.frame(t(states)))
    y <- unlist(mypblapply(states, objfun, ncpus = ncpus,
                           alpha = alpha, Y = Y, Nstar = Nstar))

    assignment <- states[[which.min(y)]]
    cost       <- min(y)
  } else {
    x0 <- makesol(alpha, Y, Nstar)
    message("Initial solution built: Cost = ", signif(objfun(x0, alpha, Y, Nstar), 3))
    message("Running Simulated Annealing (maxit = ", SAopts$maxit, ")")
    y  <- stats::optim(par = x0, fn = objfun, gr = neighbour, method = "SANN",
                       alpha = alpha, Y = Y, Nstar = Nstar,
                       control = SAopts)
    assignment <- y$par
    cost       <- y$value
  }
  message("Done! Final cost = ", signif(cost, 3))

  # Cast x as an allocation list
  xl <- lapply(seq_along(Nstar), function(i){seq_along(assignment)[assignment == i]})
  solstats <- getstats(assignment, Y, Nstar)

  return(list(x = assignment, cost = cost, solstats = solstats, xl = xl))
}
