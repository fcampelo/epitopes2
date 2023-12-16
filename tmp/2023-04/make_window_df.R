#' Assemble a sliding window representation of epitope or protein data
#'
#' Casts epitope or protein data into a dataframe where rows correspond to a
#' fixed-length window centred on consecutive positions of each sequence.
#'
#' The sliding window runs from the first to the last positions of each sequence
#' in the **df** input (which in turn comes from running [prepare_join_df()]).
#' If a window extends beyond the limits of the protein it is padded with the
#' first (or last) letter of the sequence.
#'
#' @param df data frame of epitope data (returned by [prepare_join_df()]) or
#'        protein data (returned by [get_proteins()]).
#' @param save_folder path to folder for saving the results.
#' @param window_size positive integer, size of window to use. If **df** is a
#'        *data.table* of class *joined_epitope_dt* (returned from
#'        [prepare_join_df()]) then a standard value for this parameter is
#'        used automatically if **window_size** = _NULL_ (this standard value is
#'        calculated as $(2 x min_epitope - 1)$, based on the value of
#'        **min_epitope** used in the call to [prepare_join_df()].
#' @param ncpus positive integer, number of cores to use
#'
#' @return A *data.table* object containing the sliding window representation of
#' the input data.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom data.table :=

make_window_df <- function(df,
                           save_folder     = NULL,
                           window_size     = NULL,
                           ncpus           = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.null(window_size) | assertthat::is.count(window_size),
                          is.data.frame(df),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          assertthat::is.count(ncpus))

  if(is.null(window_size)){
    if(!("min_epit" %in% names(attributes(df)))){
      stop("df must have a min_epit attribute when window_size == NULL")
    } else {
      assertthat::assert_that(assertthat::is.count(attr(df, "min_epit")))
      window_size <- 2 * attr(df, "min_epit") - 1
    }
  }

  type <- "prot"
  if("joined_epit_dt" %in% class(df)){
    type <- "epit"
  } else if(!("protein_dt" %in% class(df))){
    cat("\nNote: Cannot determine if 'df' is an epitope or protein dataframe.",
        "\nTrying to treat as 'protein' (default). It may result in errors.")
  }

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_file <- paste0(normalizePath(save_folder), "/windowed_", type, "_df_",
                      ymd, ".rds")
  }

  # ========================================================================== #
  # Generate dataframe by sliding windows
  # extract_windows() is an internal function defined in "extract_windows.R"
  cat("\nExtracting windows:\n")

  # Columns to remove:
  torm <- c("pubmed_id", "year", "epit_name", "evid_code",
            "epit_struc_def", "epit_seq", "n_assays", "bcell_id", "assay_type",
            "assay_class", "TSeq_seqtype", "TSeq_defline", "DB", "TSeq_sid",
            "TSeq_length")

  torm <- torm[torm %in% names(df)]
  if(length(torm) > 0) {
    df  <- df[, (torm) := NULL]
  }

  # Convert df into a list of lists (for lapply)
  X <- lapply(purrr::pmap(as.list(df), list),
              function(x, t){
                x$df_type <- t
                return(x)},
              t = type)

  wdf <- mypblapply(ncpus = ncpus,
                    X     = X,
                    FUN   = extract_windows,
                    ws    = window_size)

  wdf <- data.table::rbindlist(wdf)
  class(wdf) <- c(class(wdf), paste0("windowed_", type, "_dt"))

  nm <- names(wdf)
  idx <- which(nm != "Class")
  nm[idx] <- paste0("Info_", nm[idx])
  names(wdf) <- nm
  if("Class" %in% names(wdf)){
    wdf <- data.table::setcolorder(wdf , c(idx, which(names(wdf) == "Class")))
  }

  # Save resulting dataframe to file
  if(!is.null(save_folder)) {
    saveRDS(wdf, file = df_file)
  }
  cat("Done!\n")

  invisible(wdf)
}
