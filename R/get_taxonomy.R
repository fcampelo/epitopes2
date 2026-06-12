#' Retrieve taxonomic classification tables from NCBI.
#'
#' This function retrieves the taxonomic information of a vector of organism
#' IDs, from the NCBI Taxonomy data base.
#'
#' @param uids vector of organism IDs to retrieve taxonomical information.
#' @param save_folder path to folder for saving the results.
#' @param consolidate logical, should the results of each element of uid be
#' consolidated into a single data frame? (defaults to `FALSE` for compatibility
#' with older versions)
#' @param IEDBOrgFile path to the IEDB  OrganismList.XML file, if available.
#' This can be retrieved using [get_IEDB_OrgList()] and is used to fix the
#' taxonomy of some entries in the IEDB database that are not present in the
#' NCBI taxononomy database (e.g., 10000253 - a lineage under 10090,
#' _Mus musculus_).
#'
#' @return A list containing the information for each element of uids
#'
#' @author Felipe Campelo (\email{fcampelo@@gmail.com})
#'
#' @export
#'
#'
#' @examples
#' uids <- c("6282", # O. volvulus,
#'           "9606") # H. sapiens
#' get_taxonomy(uids, consolidate = TRUE)
#'

get_taxonomy <- function(uids,
                         save_folder = NULL,
                         consolidate = FALSE,
                         IEDBOrgFile = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  t0 <- Sys.time()
  ok_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(is.null(save_folder) | (is.character(save_folder)),
                          length(save_folder) <= 1,
                          any(class(uids) %in% ok_classes),
                          length(uids) >= 1,
                          is.logical(consolidate), length(consolidate) == 1,
                          is.null(IEDBOrgFile) | (is.character(IEDBOrgFile) &&
                                                    length(IEDBOrgFile) == 1 &&
                                                    file.exists(IEDBOrgFile)))

  # Extract unique Taxonomy IDs for retrieval
  if(is.character(uids)){
    uids <- lapply(uids,
                   function(x){strsplit(x, split = ",")[[1]]})
  } else if (!is.null(uids)){
    uids <- format(unlist(uids), scientific = FALSE, trim = TRUE)
  }

  uids <- unique(uids)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    df_file <- paste0(normalizePath(save_folder), "/taxonomy.rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/taxonomy_retrieval_errlist.rds")
  }

  errlist <- seq_along(uids)
  reslist <- vector("list", length = length(uids))
  nerr    <- Inf

  # Try to retrieve taxonomies from NCBI (one by one)
  while(length(errlist) < nerr && length(errlist) > 0){
    nerr <- length(errlist)
    message("Trying to retrieve ", length(errlist), " entries from NCBI (db = taxonomy)")
    cc <- 0
    for (idx in errlist){
      errk <- FALSE
      suppressMessages(
        tryCatch({
          tt <- entrez_fetch(db = "taxonomy",
                             id = as.numeric(uids[idx]),
                             retmode = "xml", rettype = "xml")
          ttp <- XML::xmlTreeParse(tt, useInternalNodes = TRUE)
          reslist[[idx]]$Taxonomy <- data.frame(
            ScientificName = XML::xpathSApply(ttp,
                                              "//TaxaSet/Taxon/LineageEx/Taxon/ScientificName",
                                              XML::xmlValue),
            Rank = XML::xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/Rank",
                                    XML::xmlValue),
            UID  = XML::xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/TaxId",
                                    XML::xmlValue),
            stringsAsFactors = FALSE)

          reslist[[idx]]$Target <-  data.frame(
            ScientificName = XML::xpathSApply(ttp,
                                              "//TaxaSet/Taxon/ScientificName",
                                              XML::xmlValue),
            Rank = XML::xpathSApply(ttp, "//TaxaSet/Taxon/Rank",
                                    XML::xmlValue),
            UID  = XML::xpathSApply(ttp, "//TaxaSet/Taxon/TaxId",
                                    XML::xmlValue),
            stringsAsFactors = FALSE)
        },
        warning = function(c) {errk <<- TRUE},
        error   = function(c) {errk <<- TRUE},
        finally = NULL))

      if(!errk){
        reslist[[idx]]$UID <- uids[idx]
      }

      # Print progress bar
      mypb(i = cc, max_i = length(errlist), t0 = t0, npos = 30)
      cc <- cc + 1

      # NCBI limits requests to three per second
      Sys.sleep(0.3)
    }
    errlist <- which(sapply(reslist, function(x) {is.null(x$UID)}))
  }

  if(length(errlist) > 0) reslist <- reslist[-errlist]
  errlist <- uids[errlist]

  reslist <- lapply(reslist,
                    \(x, consolidate){
                      if(consolidate){
                        x$Taxonomy <- rbind(x$Taxonomy, x$Target)
                        x$Target   <- NULL
                      }
                      x$consolidated <- consolidate
                      return(x)
                    }, consolidate = consolidate)


  # Fix taxonomy for entries with no taxonomic information if IEDBOrgFile is provided
  # These are mostly sub-species lineages that are recorded in the IEDB database
  idx   <- which(sapply(reslist, \(x) nrow(x$Taxonomy) == 0))
  rlids <- sapply(reslist, \(x) x$UID)
  searchids <- rlids[idx]

  if(length(idx) > 0 & !is.null(IEDBOrgFile)){
    message("Trying to fix taxonomy for ", length(idx), " entries using IEDB OrganismList.XML")

    # Read IEDB OrganismList.XML
    message("Reading and extracting data from OrganismList.XML...")
    xmldata <- XML::xmlTreeParse(IEDBOrgFile, useInternalNodes = TRUE)
    IEDB.orgdf <- data.frame(
      OrganismId = as.character(XML::getNodeSet(xmldata, "//Organism/OrganismId", fun = XML::xmlValue)),
      ParentTaxId = as.character(XML::getNodeSet(xmldata, "//Organism/ParentTaxId", fun = XML::xmlValue)),
      OrganismName = as.character(XML::getNodeSet(xmldata, "//Organism/OrganismName", fun = XML::xmlValue)))

    for(k in seq_along(searchids)){
      searchstack <- searchids[k]
      go <- TRUE
      while(go){
        tr <- IEDB.orgdf[sapply(searchstack, \(x) which(IEDB.orgdf$OrganismId == x)) , ]
        rlpos <- which(rlids == tr$ParentTaxId[nrow(tr)])
        if(length(rlpos) == 1 && nrow(reslist[[rlpos]]$Taxonomy) > 0){
          go <- FALSE
        } else {
          searchstack <- c(searchstack, tr$ParentTaxId[nrow(tr)])
        }
      }

      tr <- tr[nrow(tr):1, ]
      reslist[[idx[k]]]$Taxonomy <- rbind(reslist[[rlpos]]$Taxonomy,
                                          data.frame(ScientificName = tr$OrganismName,
                                                     Rank = "no rank",
                                                     UID = tr$OrganismId,
                                                     stringsAsFactors = FALSE))
    }
  }

  # Save results to file
  if(!is.null(save_folder)){
    saveRDS(object = reslist, file = df_file)
    saveRDS(object = errlist, file = errfile)
  }


  message("Done!\n", length(reslist), " taxonomies retrieved.\n",
          length(errlist), " retrieval errors.")

  return(reslist)
}
