#' Filter peptides by taxonomy and host
#'
#' Filters a peptides dataframe by source organism and/or host IDs.
#'
#' Taxonomic filtering works hierarchically - e.g., if a genus id is passed in
#' `orgIDs`, all species under that genus are removed.
#' `removeIDs` is generally used to censor out one or more subgroups under those
#' listed in `orgIDs`. For instance, if we wanted to get data for all
#' orthopoxvirus (txID:10242) except Variola (txID:10255) we would make
#' `orgIDs = 10242` and `removeIDs = 10255`.
#'
#' @inheritParams taxonomy_filter
#'
#' @return data frame with peptide data filtered by the criteria in **orgIDs**,
#' **hostIDs** and **removeIDs**.
#'
#' @author Felipe Campelo (\email{fcampelo@@gmail.com})
#'

taxonomy_filter_df <- function(df            = NULL,
                               tax_list      = NULL,
                               orgIDs        = NULL,
                               hostIDs       = NULL,
                               removeIDs     = NULL,
                               orgID_column  = "sourceOrg_id",
                               hostID_column = "host_id") {


  # Standardise relevant variables:
  ids <- data.frame(org  = as.character(df[, which(names(df) == orgID_column), drop = TRUE]),
                    host = as.character(df[, which(names(df) == hostID_column), drop = TRUE]))

  # Function to extract the relevant organism IDs using the taxonomy data
  fextr <- function(x, tid){
    ul <- unique(c(x$Taxonomy$UID, x$UID))
    if(any(tid %in% ul)) return(x$UID)
  }

  # Function to get the filtering indices
  fmatch <- function(pat, trg){
    str <- unlist(strsplit(pat, split = ",", fixed = TRUE))
    any(str %in% trg)
  }

  # Filter by OrgIDs
  if (!is.null(orgIDs)){
    target_org  <- unlist(sapply(tax_list, FUN = fextr, tid = orgIDs))
    target_org  <- unique(c(target_org, orgIDs))
    idx1        <- which(sapply(ids$org, FUN = fmatch, trg = target_org))
    if(length(idx1) > 0){
      ids  <- ids[idx1, ]
      df   <- df[idx1, ]
    } else {
      ids <- numeric()
      df  <- numeric()
    }
  }

  # Filter by hostIDs
  if (!is.null(hostIDs) && nrow(df) > 0){
    target_host <- unlist(sapply(tax_list, FUN = fextr, tid = hostIDs))
    target_host <- unique(c(target_host, hostIDs))
    idx2        <- which(sapply(ids$host, FUN = fmatch, trg = target_host))
    if(length(idx2) > 0){
      ids  <- ids[idx2, ]
      df   <- df[idx2, ]
    } else {
      ids <- numeric()
      df  <- numeric()
    }
  }

  # Filter by removeIDs
  if (!is.null(removeIDs) && nrow(df) > 0){
    target_rm   <- unlist(sapply(tax_list, FUN = fextr, tid = removeIDs))
    target_rm   <- unique(c(target_rm, removeIDs))
    idx3        <- which(sapply(ids$org, FUN = fmatch, trg = target_rm))
    if (length(idx3) > 0){
      ids  <- ids[-idx3, ]
      df   <- df[-idx3, ]
    }
  }

  df <- dplyr::as_tibble(df)

  # Set attributes to output data frame:
  attr(df, "orgIDs")    <- nullcheck(orgIDs)
  attr(df, "hostIDs")   <- nullcheck(hostIDs)
  attr(df, "removeIDs") <- nullcheck(removeIDs)

  return(df)

}
