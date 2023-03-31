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

  # === Run optimisation === #

  # Pre-optimise splitting of id_force_splitting
  if(!is.null(id_force_splitting)){
    idx <- unique(unlist(sapply(id_force_splitting,
                                function(x,ids){grep(x,ids)},
                                ids = Y$txids)))

    Yp <- Y[idx, ]
    Yp$split <- NA
    Yp$Realcluster <- Yp$Cluster
    Yp$Cluster <- 1:nrow(Yp)
    Np <- Nstar[1:min(length(Nstar), length(idx))]
    Np <- Np / sum(Np)
    optp <- list(maxit = min(1e5,
                             2000 * round(log10(length(Np) ^ nrow(Yp)))))

    x0 <- makesol(alpha, Yp, Np)

    Yp <- cbind(Yp[, 1], split = x0)
    Y <- dplyr::left_join(Y, Yp, by = c("Cluster"))
  } else {
    Y$split <- NA
  }

  nstates <- length(Nstar) ^ sum(is.na(Y$split))
  message("Optimising splits with alpha = ", alpha,
          "\n(Number of possibilities: ~", signif(nstates, 3), ")")
  if(nstates < SAopts$maxit){
    message("Running exhaustive search...")
    # If the search space is small enough, exhaustive search suffices
    Yn <- Y[is.na(Y$split), ]
    states <- do.call(expand.grid,
                      lapply(Yn$Cluster, function(x){1:length(Nstar)}))
    states <- as.list(as.data.frame(t(states)))
    y <- unlist(mypblapply(states, objfun, ncpus = ncpus,
                           alpha = alpha, Y = Y, Nstar = Nstar))

    assignment <- states[[which.min(y)]]
    cost       <- min(y)
  } else {
    message("Running heuristic search...")
    y <- makesol(alpha, Y, Nstar)
    assignment <- y
    if(SAopts$torun == TRUE){
      SAopts$torun <- NULL
      y  <- stats::optim(par = y, fn = objfun, gr = neighbour, method = "SANN",
                         alpha = alpha, Y = Y, Nstar = Nstar,
                         control = SAopts)
      assignment <- y$par
    }
  }

  # Cast x as an allocation list
  xl <- lapply(seq_along(Nstar), function(i){seq_along(assignment)[assignment == i]})
  solstats <- getstats(assignment, Y, Nstar)

  return(list(x = assignment, solstats = solstats, xl = xl))
}
