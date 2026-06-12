# Function to extract windows from a given row of df
# (local scope only - not exported)
extract_windows <- function(x, ws){

  # Pad string at both ends
  pseq  <- x$TSeq_sequence
  nw    <- nchar(pseq)
  char1 <- paste(rep(substr(pseq, 1, 1), (ws - 1) / 2),
                 collapse = "")
  charN <- paste(rep(substr(pseq, nw, nw), (ws - 1) / 2),
                 collapse = "")
  pseq  <- paste(char1, pseq, charN, collapse = "", sep = "")

  if(x$df_type == "prot"){
    start_pos <- 1 + (ws - 1) / 2
    stop_pos  <- nchar(pseq) - (ws - 1) / 2
  } else {
    start_pos <- x$start_pos + (ws - 1) / 2
    stop_pos  <- x$end_pos + (ws - 1) / 2
    nw        <- stop_pos - start_pos + 1
  }

  # Initialise data frame
  x <- data.table::as.data.table(x)[rep(1, nw), ]

  torm <- c("TSeq_sequence", "df_type", "start_pos", "end_pos")
  torm <- torm[torm %in% names(x)]
  if(length(torm) > 0) x <- x [, (torm) := NULL]

  x$center_pos <- (start_pos - (ws - 1) / 2):(stop_pos - (ws - 1) / 2)
  x$window_seq <- sapply(start_pos:stop_pos,
                         function(i, pseq, wl){ substr(pseq, i - wl, i + wl) },
                         pseq = pseq,
                         wl   = (ws - 1) / 2)

  return(x)
}
