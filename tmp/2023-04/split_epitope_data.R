#' Split epitope data based on epitope, protein or organism IDs.
#'
#' Takes a data.table of data of class *windowed_epit_dt* (returned by
#' [make_window_df()]) and split it into mutually exclusive sets of
#' observations, based on columns *Info_sourceOrg_id*, *Info_protein_id* or
#' *Info_epitope_id*.
#'
#' If the sum of **split_perc** is less than 100 an extra split is generated
#' with the remaining observations - e.g., `split_perc = c(50, 30)` results in
#' three sets with an approximately 50/30/20% split *of the total observations.*
#' If the sum is greater than 100 the splits are linearly scaled down so that
#' the sum becomes 100. Note that the split percents correspond to the number of
#' observations, not the number of unique IDs.
#'
#' This function will attempt to approximate the desired split levels, but
#' depending on the size of **wdf** set and the desired **split_level** it may
#' not be possible (e.g., if `split_level = "org"` and a single organism
#' corresponds to 90% of the data, one of the splits will necessarily correspond
#' to at least 90% of the data, regardless of the values informed in
#' `split_perc`.
#'
#' If a BLASTp file is provided the routine will keep any pairs of proteins
#' having (coverage >= **coverage_threshold** AND
#' identity >=  **identity_threshold**) under the same split. This is useful to
#' prevent accidental data leakage due to quasi-identical proteins with
#' different UIDs. **NOTE**: this only works if `split_level == "prot`.
#'
#' @param wdf data table of class *windowed_epit_dt* (returned by
#'        [make_window_df()])
#' @param split_level which level should be used for splitting? Use "org" for
#'        splitting by source organism ID, "prot" by protein ID or "epit" by
#'        epitope ID. When "prot" is used the routine attempts to identify
#'        different protein versions and treat them as a single unit for
#'        splitting purposes.
#' @param split_perc numeric vector of desired splitting percentages. See
#'        Details.
#' @param split_names optional character vector with short names for each split.
#' @param blast_file path to file containing all-vs-all BLASTp alignment results
#'        for all proteins in **wdf**. See Details.
#' @param coverage_threshold coverage threshold for grouping proteins by
#' similarity, see Details.
#' @param identity_threshold identity threshold for grouping proteins by
#' similarity, see Details.
#' @param save_folder path to folder for saving the results.
#'
#' @return A list object containing the split data tables.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

split_epitope_data <- function(wdf,
                               split_level = "prot",
                               split_perc = c(70, 30),
                               split_names = NULL,
                               save_folder = NULL,
                               blast_file  = NULL,
                               coverage_threshold = 80,
                               identity_threshold = 80){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that("windowed_epit_dt" %in% class(wdf),
                          is.character(split_level), length(split_level) == 1,
                          split_level %in% c("org", "prot", "epit"),
                          is.numeric(split_perc),
                          all(sapply(split_perc, assertthat::is.count)),
                          is.null(blast_file) | is.character(blast_file),
                          is.null(split_names) | is.character(split_names))

  # Check and adjust split sizes if necessary
  ssp <- sum(split_perc)
  if (ssp < 100){
    split_perc <-c(split_perc, 100 - ssp)
  } else if (ssp > 100){
    split_perc <- 100 * split_perc / ssp
  }
  if (ssp >= 100 & length(split_perc) == 1){
    cat("\nNo splitting performed (split_perc = 100).")
    invisible(wdf)
  }
  nsplits <- length(split_perc)

  # Set up split names
  if (is.null(split_names)){
    split_names <- 1:nsplits
  } else if (length(split_names) < nsplits){
    cat("\nLength of split names smaller than number of splits.",
    "\nUsing split numbers instead.")
    split_names <- 1:nsplits
  }


  # Determine splits by splitting column
  # (Note that protein IDs ignore the trailing version number)
  id_var <- switch(split_level,
                   org  = wdf$Info_sourceOrg_id,
                   prot = gsub("\\.[1-9]+$", "", wdf$Info_protein_id),
                   epit = wdf$Info_epitope_id)

  # Check similarity based on the blast file (if)
  if (split_level == "prot" && !is.null(blast_file)){
    blast <- utils::read.csv(blast_file, sep = "\t",
                             header = FALSE,
                             stringsAsFactors = FALSE)
    names(blast) <- c("QueryID", "SubjectID", "Alignment_length",
                      "Query_length", "Subject_length", "Num_matches",
                      "Perc_identity", "Query_coverage", "Num_mismatches",
                      "Num_gaps", "Query_match_start", "Query_match_end",
                      "Subject_match_start", "Subject_match_end",
                      "E_value", "Score")

    blast$QueryID <- gsub("\\.[1-9]+$", "", blast$QueryID)
    blast$SubjectID <- gsub("\\.[1-9]+$", "", blast$SubjectID)

    blast <- data.table::as.data.table(blast)

    # Initialise internal data.table variable names to prevent CRAN notes.
    Query_coverage <- Perc_identity <- QueryID <- SubjectID <- NULL

    # Filter relevant blast entries and variables
    blast <- blast[(Query_coverage >= coverage_threshold) & (Perc_identity >= identity_threshold), ]
    blast <- blast[!duplicated(t(apply(blast[, 1:2], 1, sort))), list(QueryID, SubjectID)]
    blast <- apply(blast, 1, paste, collapse=" ")
  } else {
    blast <- character()
  }

  # Get the frequencies of occurrence of each ID
  divs <- as.data.frame(sort(table(id_var), decreasing = TRUE),
                        stringsAsFactors = FALSE)
  divs$Freq <- 100 * divs$Freq / length(id_var)

  # Determine which IDs go into which splits
  split_ids <- vector("list", nsplits)
  names(split_ids) <- split_names
  for(i in 1:nsplits) split_ids[[i]] <- list(id_type = split_level,
                                             IDs = character(),
                                             IDgroups = vector("list"),
                                             Perc = 0)
  sc <- 1   # split counter
  nc <- 0   # number of attempts to attribute
  while (nrow(divs) > 0){
    IDgroup <- divs$id_var[1]
    go <- TRUE
    cc <- 1
    while(go){
      idx <- grep(IDgroup[cc], blast)
      if(length(idx) > 0){
        x <- blast[idx] # get entries that have cands[cc]
        blast <- blast[-idx]
        x <- gsub(paste0("\\s*", IDgroup[cc], "\\s*"), "", x)
        IDgroup <- c(IDgroup, x)
        cc <- cc + 1
      } else {
        go <- FALSE
      }
    }

    Ptot <- sum(divs$Freq[divs$id_var %in% IDgroup])

    # If the current split (sc) can accommodate all the data in IDgroup
    if(Ptot <= split_perc[sc] - split_ids[[sc]]$Perc){
      split_ids[[sc]]$IDs      <- c(split_ids[[sc]]$IDs, IDgroup)
      split_ids[[sc]]$Perc     <- split_ids[[sc]]$Perc + Ptot
      split_ids[[sc]]$IDgroups <- c(split_ids[[sc]]$IDgroups, list(IDgroup))
      divs <- divs[-which(divs$id_var %in% IDgroup), ]
      nc <- 0
    } else {
      nc <- nc + 1
    }
    if (nc > nsplits) break

    sc <- 1 + sc %% nsplits
  }

  # Allocate any remaining ones
  if (length(divs) > 0){
    sidx <- which.max(split_perc - sapply(split_ids, function(x) x$Perc))
    split_ids[[sidx]]$IDs  <- c(split_ids[[sidx]]$IDs, divs$id_var)
    split_ids[[sidx]]$Perc <- split_ids[[sidx]]$Perc + sum(divs$Freq)
  }

  for (i in seq_along(split_ids)){
    split_ids[[i]]$wdf <- wdf[which(id_var %in% split_ids[[i]]$IDs), ]
  }

  return(split_ids)
}


