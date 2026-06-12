#' Label protein positions based on epitope data
#'
#' Takes a data.table of data of class *protein_dt* (returned by
#' [get_proteins()]) or *windowed_prot_dt* (returned by [make_window_df()]) and
#' labels each protein position based on known epitope/non-epitope regions
#' documented in an **epitopes** data.table of class *joined_epit_dt* (returned
#' by [prepare_join_df()] or by [filter_epitopes()]).
#'
#' @param epitopes data frame of epitope data, returned by [prepare_join_df()]
#'        or [filter_epitopes()].
#' @param proteins data frame of protein data, returned by [get_proteins()] or
#'        [make_window_df()].
#' @param set.positive how to decide whether a position should be labeled as
#'        "Positive" (+1). Use "any" to set a position as positive if
#'        it is labeled as $+1$ in at least one entry of **epitopes**; "mode" to
#'        set it by majority voting; or "all" to only label a position as
#'        Positive if it has at least one occurrence as $+1$ and none as $-1$.
#'        Unlabelled positions receive **NA** in their *Class* column.
#' @param ncpus number of cores to use
#' @param save_folder path to folder for saving the results.
#'
#' @return A data table of class *windowed_prot_dt* with columns containing the
#' number of positive and negative labels found for each position of each
#' protein, plus a *Class* column calculated according to *set.positive*.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

label_proteins <- function(proteins, epitopes,
                           set.positive = c("any", "mode", "all"),
                           ncpus = 1,
                           save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  prot_classes <- c("windowed_prot_dt", "protein_dt")
  epit_classes <- c("joined_epit_dt")
  assertthat::assert_that(is.data.frame(epitopes),
                          is.data.frame(proteins),
                          any(class(proteins) %in% prot_classes),
                          any(class(epitopes) %in% epit_classes),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          is.character(set.positive),
                          assertthat::is.count(ncpus))


  if(all(tolower(set.positive) == "any")) {
    set.positive <- "any"
  } else if(all(tolower(set.positive) == "all")) {
    set.positive <- "all"
  } else set.positive <- "mode"

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_file <- paste0(normalizePath(save_folder), "/", ymd,
                      "_df_labelled_prots.rds")
  }

  # Prepare protein data.table for labeling
  if(!("windowed_prot_dt" %in% class(proteins))){
    # Convert proteins into windowed data format
    proteins <- make_window_df(proteins, window_size = 1, ncpus = ncpus)
  }

  # names(proteins)[which(names(proteins) == "Info_window_seq")] <- "Info_AA"

  cat("\nExtracting label positions from 'epitopes'\n")
  myf <- function(x){
    st <- as.numeric(x$start_pos)
    en <- as.numeric(x$end_pos)
    data.frame(UID     = x$protein_id,
               pos     = st:en,
               epit_id = x$epitope_id,
               nPos    = x$n_Positive,
               nNeg    = x$n_Negative)
  }

  tmp <- mypblapply(ncpus = ncpus,
                    X     = purrr::pmap(as.list(epitopes), list),
                    FUN   = myf)

  cat(" - Done!\n")

  # Aggregate multiply-labelled entries. The variable names are initialised
  # below just to prevent NOTEs on CRAN. The dplyr functions use references to
  # variables internal to df)
  UID <- pos <- nPos <- nNeg <- epit_id <- NULL
  tmp <- data.table::rbindlist(tmp, use.names = TRUE)
  tmp <- tmp[, list(Info_epit_id = paste(unique(epit_id), collapse = ","),
                    Info_nPos    = sum(nPos),
                    Info_nNeg    = sum(nNeg)),
             by = list(UID, pos)]

  # Join epitope labels into proteins dataset
  proteins <- dplyr::left_join(proteins, tmp,
                               by = c("Info_UID" = "UID",
                                      "Info_center_pos" = "pos"))

  # Set class
  if (set.positive == "any"){
    proteins$Class <- -1 + 2 * as.numeric(proteins$Info_nPos > 0)
  } else if (set.positive == "all") {
    proteins$Class <- -1 + 2 * (proteins$Info_nNeg == 0)
  } else {
    proteins$Class <- -1 + 2 * (proteins$Info_nPos >= proteins$Info_nNeg)
  }

  if(!is.null(save_folder)) {
    saveRDS(proteins, file = df_file)
  }

  return(proteins)

}




#' Label protein positions based on epitope data
#'
#' Takes a protein data frame and
#' labels each protein position based on known epitope/non-epitope regions
#' documented in an **epitopes** dataframe.
#'
#' @param epitopes data frame of epitope data.
#' @param proteins data frame of protein data.
#' @param ncpus number of cores to use.
#'
#' @return An updated protein dataframe with a Class attribute calculated for
#' each position.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

label_proteins2 <- function(proteins, epitopes, ncpus = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(epitopes),
                          is.data.frame(proteins),
                          assertthat::is.count(ncpus))

  # Prepare protein data for labeling
  if(!("Info_center_pos" %in% names(proteins))){
    # Convert proteins into windowed data format
    proteins <- make_window_df(as.data.table(proteins),
                               window_size = 1, ncpus = ncpus)
  }

  cat("\nExtracting label positions from 'epitopes'\n")
  myf <- function(start_pos, end_pos, protein_id, epitope_id, nPos, nNeg){
    st <- as.numeric(start_pos)
    en <- as.numeric(end_pos)
    if (!any(is.na(c(st, en)))){
      data.frame(UID     = protein_id,
                 pos     = st:en,
                 epit_id = epitope_id,
                 nPos    = nPos,
                 nNeg    = nNeg)
    } else {
      data.frame(UID = character(), pos = numeric(), epit_id = character(),
                 nPos = numeric(), nNeg = numeric())
    }
  }

  tmp <- parallel::mcmapply(myf,
                            epitopes$start_pos,
                            epitopes$end_pos,
                            epitopes$protein_id,
                            epitopes$epitope_id,
                            epitopes$n_Positive,
                            epitopes$n_Negative,
                            mc.cores = ncpus,
                            SIMPLIFY = FALSE)
  tmp <- data.table::rbindlist(tmp, use.names = TRUE)
  cat(" - Done!\n")

  # Aggregate multiply-labelled entries.
  tmp$Class <- -1 + 2 * (tmp$nPos >= tmp$nNeg)
  tmp <- tmp[, list(Info_epit_id = paste(unique(epit_id), collapse = ","),
                    Info_nPos    = sum(Class == 1),
                    Info_nNeg    = sum(Class == -1)),
             by = list(UID, pos)]

  # Join epitope labels into proteins dataset
  proteins <- dplyr::left_join(proteins, tmp,
                               by = c("Info_UID" = "UID",
                                      "Info_center_pos" = "pos"))

  # Set class
  if (set.positive == "any"){
    proteins$Class <- -1 + 2 * as.numeric(proteins$Info_nPos > 0)
  } else if (set.positive == "all") {
    proteins$Class <- -1 + 2 * (proteins$Info_nNeg == 0)
  } else {
    proteins$Class <- -1 + 2 * (proteins$Info_nPos >= proteins$Info_nNeg)
  }

  return(proteins)

}

