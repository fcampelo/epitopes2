#' Smoothing of prediction vectors
#'
#' This function runs a sliding window smoother on prediction vectors.
#'
#' @param x vector of predictions (must be numeric)
#' @param window_size size of the sliding window. Only used if `type` is
#'        either `"mean"` or `"mode"`.
#' @param type type of smoothing function to use. Accepts "mode" (simple
#'        majority vote), "mean" or "minsize".
#' @param minsize smallest sequence length for positive predictions.
#'        Stretches of positive-labelled positions shorter than `minsize` are
#'        flipped to negative (only works if `type == "minsize`).
#'
#' @return vector of smoothed prediction values.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
smooth_predictions <- function(x, window_size = 15,
                               type = "mode",
                               minsize = 8){

  valid_types <- c("mean", "mode", "minsize")
  assertthat::assert_that(is.numeric(x),
                          assertthat::is.count(window_size),
                          length(x) >= window_size,
                          is.character(type),
                          length(type) == 1,
                          type %in% valid_types,
                          assertthat::is.count(minsize))

  if(type %in% c("mean", "mode")){

    xlen <- length(x)
    y    <- x
    for (i in 1:xlen){
      # Define range covered by the window
      i1 <- max(1, min(i - floor(window_size / 2),
                       xlen - window_size + 1))
      i2 <- min(xlen, i1 + window_size - 1)

      # Perform smoothing
      if (type == "mean"){
        y[i] <- mean(x[i1:i2])
      } else {
        uniqx <- unique(x[i1:i2])
        y[i] <-uniqx[which.max(tabulate(match(x[i1:i2], uniqx)))]
      }
    }
  } else if (type == "minsize"){
    lx <- c(1, x[1:(length(x) - 1)])
    IsBreak = (x == 1 & lx == -1)
    if(x[1] == 1) IsBreak[1] <- TRUE
    df       <- data.frame(ind = seq_along(x), x = x, IsBreak = IsBreak)
    df       <- df[df$x == 1, ]
    df$EpNum <- cumsum(df$IsBreak)
    EpSize   <- table(df$EpNum)
    ToRM     <- as.numeric(names(EpSize)[EpSize < minsize])
    ToRM     <- df$ind[df$EpNum %in% ToRM]
    y        <- x
    y[ToRM]  <- -1
  }

  return(y)
}
