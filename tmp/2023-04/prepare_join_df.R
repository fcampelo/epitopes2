#' Merge and filter epitope and protein data
#'
#' Merges previously extracted epitopes and proteins, and performs some
#' sanity checks and filtering.
#'
#' Entries in the `epitopes` input are removed if they:
#' \itemize{
#'   \item Lack a valid protein ID (i.e., one that has a corresponding entry in
#'   `proteins`)
#'   \item Lack a valid string in *epit_seq*
#'   \item Lack a valid definition in *epit_struc_def*
#'   \item Have sequences shorter than `min_epit` or longer than `max_epit`.
#'   \item Have a mismatch between the sequence in *epit_seq* and the
#'   corresponding sequence between *start_pos* and *end_pos* on the protein
#'   sequence (see description of parameter `pos.mismatch.rm`).
#' }
#'
#' @param epitopes data frame of epitope data (returned by [get_LBCE()]).
#' @param proteins data frame of protein data (returned by [get_proteins()]).
#' @param save_folder path to folder for saving the results.
#' @param min_epit positive integer, shortest epitope to be considered
#' @param max_epit positive integer, longest epitope to be considered
#' @param only_exact logical, should only sequences labelled as "Exact Epitope"
#'        in variable *epit_struc_def* (within `epitopes`) be considered?
#' @param pos.mismatch.rm should epitopes with position mismatches be removed?
#'        Use "all" (default) for removing any position mismatch or "align"
#'        if the routine should attempt to
#'        search the epitope sequence in the protein sequence.
#' @param set.positive how to decide whether an observation should be of the
#'        "Positive" (+1) class? Use "any" to set a sequence as positive if
#'        $n_positive > 0$, "mode" to set it if $n_positive >= n_negative$,
#'        or "all" to set it if $n_negative == 0$. Defaults to "mode".
#'
#' @return A data.table object containing the merged data frame
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%

prepare_join_df <- function(epitopes, proteins,
                            save_folder     = NULL,
                            min_epit        = 8,
                            max_epit        = 25,
                            only_exact      = FALSE,
                            pos.mismatch.rm = c("all", "align"),
                            set.positive    = c("any", "mode", "all")){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(epitopes),
                          is.data.frame(proteins),
                          assertthat::is.count(min_epit),
                          assertthat::is.count(max_epit),
                          min_epit <= max_epit,
                          is.logical(only_exact), length(only_exact) == 1,
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          is.character(pos.mismatch.rm),
                          is.character(set.positive))

  if(all(tolower(pos.mismatch.rm) == "align")) {
    pos.mismatch.rm <- "align"
  } else pos.mismatch.rm <- "all"

  if(all(tolower(set.positive) == "any")) {
    set.positive <- "any"
  } else if(all(tolower(set.positive) == "all")) {
    set.positive <- "all"
  } else set.positive <- "mode"

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_file <- paste0(normalizePath(save_folder), "/joined_df_", ymd, ".rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/joined_df_errlist_", ymd, ".rds")
  }

  # ========================================================================== #
  # Initial preprocessing

  # Join epitopes and proteins dataframes, preliminary feature transformation
  df <- dplyr::left_join(epitopes, proteins,
                         by = c("protein_id" = "UID"))

  # Filter by missing info:
  noseq  <- is.na(df$epit_seq)
  noprot <- is.na(df$TSeq_sequence)
  nodef  <- is.na(df$epit_struc_def)
  rm.idx <- which(noseq | noprot | nodef)

  rm.df  <- data.frame(epitope_id = df$epitope_id[rm.idx],
                       excl_crit  = "missing data",
                       stringsAsFactors = FALSE)
  df <- df[-rm.idx, ]

  # filter by sequence mismatch
  # Check if epitopes are in their indicated positions
  refstr <- mapply(FUN   = substr,
                   x     = df$TSeq_sequence,
                   start = df$start_pos,
                   stop  = df$end_pos,
                   USE.NAMES = FALSE)
  isok <- (df$epit_seq == refstr)
  rm.idx <- which(is.na(isok) | !isok)
  if(pos.mismatch.rm == "align"){
    # Attempt to fix some inconsistencies before removing
    pos    <- stringr::str_locate_all(string  = df$TSeq_sequence[rm.idx],
                                      pattern = df$epit_seq[rm.idx])
    tofix  <- which(sapply(pos, function(l){(nrow(l) == 1) && !any(is.na(l))}))
    pos    <- do.call(rbind, pos[tofix])
    df.idx <- rm.idx[tofix]
    rm.idx <- rm.idx[-tofix]

    df$start_pos[df.idx] <- pos[, 1]
    df$end_pos[df.idx]   <- pos[, 2]
  }

  rm.df  <- rbind(rm.df,
                  data.frame(epitope_id = df$epitope_id[rm.idx],
                             excl_crit  = "position mismatch",
                             stringsAsFactors = FALSE))
  df <- df[-rm.idx, ]


  # Filter by epitope length
  el <- nchar(df$epit_seq)
  rm.idx <- which((el < min_epit) | (el > max_epit))
  rm.df  <- rbind(rm.df,
                  data.frame(epitope_id = df$epitope_id[rm.idx],
                             excl_crit  = "epitope size",
                             stringsAsFactors = FALSE))
  df <- df[-rm.idx, ]


  # Filter by exact epitopes
  if (only_exact){
    rm.idx <- which(df$epit_struc_def != "Exact Epitope")
    rm.df  <- rbind(rm.df,
                    data.frame(epitope_id = df$epitope_id[rm.idx],
                               excl_crit  = "non-Exact Epitope",
                               stringsAsFactors = FALSE))
    df <- df[-rm.idx, ]
  }

  # TODO: Combine redundant entries
  #multidx <- table(df$epitope_id)
  #multidx <- multidx[multidx > 1]
  # df <- as.data.frame(df) %>%
  #   dplyr::group_by("epitope_id")


  # Set class
  if (set.positive == "any"){
    df$Class <- -1 + 2 * as.numeric(df$n_Positive > 0)
  } else if (set.positive == "all") {
    df$Class <- -1 + 2 * (df$n_Negative == 0)
  } else {
    df$Class <- -1 + 2 * (df$n_Positive >= df$n_Negative)
  }

  class(df) <- c("data.table", "data.frame", "joined_epit_dt")
  attr(df, "min_epit") <- min_epit

  if(!is.null(save_folder)){
    saveRDS(object = df,    file = df_file)
    saveRDS(object = rm.df, file = errfile)
  }

  invisible(df)
}
