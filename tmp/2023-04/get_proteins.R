#' Retrieve protein sequences and data from GenBank and Uniprot
#'
#' This function is used to retrieve data from
#' Genbank's protein database for
#' given protein IDs. If an ID is not available from Genbank the function will
#' try to retrieve it from Uniprot.
#'
#' Queries are processed one by one (rather than in batch) to enable treatment
#' of individual inconsistencies (e.g., wrong UIDs, queries that return a
#' different identifier, etc.). This makes this routine substantially slower,
#' but considerably more robust to errors.
#'
#' @param uids A character vector with protein IDs.
#' @param save_folder path to folder for saving the results.
#'
#' @return A data frame with the extracted proteins.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_proteins <- function(uids, save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  t0 <- Sys.time()
  assertthat::assert_that(is.null(save_folder) | (is.character(save_folder) && length(save_folder) == 1),
                          is.character(uids),
                          length(uids) >= 1)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    df_file <- paste0(normalizePath(save_folder), "/00_proteins.rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/00_prots_not_retrieved.rds")
    tmpf    <- tempfile(fileext = ".rds", tmpdir = save_folder)
  }

  errlist <- seq_along(uids)
  reslist <- vector("list", length = length(uids))
  nerr    <- Inf

  # First try retrieving from NCBI/protein
  cat("Trying to retrieve", length(reslist), "proteins",
      "\nStarted at", as.character(t0),
      "\nThis may take a while...")
  while(length(errlist) < nerr && length(errlist) > 0){
    nerr <- length(errlist)
    cat("\nTrying to retrieve", length(errlist),
        "entries from NCBI (db = protein)")
    cc <- 0
    for (idx in errlist){
      # Try fetching data
      tryCatch({
        x <- reutils::efetch(uid = uids[idx],
                             db      = "protein",
                             rettype = "fasta",
                             retmode = "xml")

        reslist[idx] <- XML::xmlToList(x$get_content())
      },
      warning = function(c) {errk <<- TRUE},
      error   = function(c) {errk <<- TRUE},
      finally = NULL)

      if(!is.null(reslist[[idx]]) &
         !("ERROR" %in% names(reslist[[idx]]))){
        reslist[[idx]]$UID <- uids[idx]
        reslist[[idx]]$DB  <- "NCBI protein"
      }
      # Print progress bar
      mypb(i = cc, max_i = length(errlist), t0 = t0, npos = 40)
      cc <- cc + 1

      # save tmp results (if needed)
      if(!is.null(save_folder) && !(cc %% 100)){
        saveRDS(object = list(reslist = reslist, errlist = errlist,
                              idx = idx, uids = uids),
                file = tmpf)
      }
    }
    errlist <- which(sapply(reslist, function(x) {is.null(x$UID)}))
  }

  # save tmp results (if needed)
  if(!is.null(save_folder)){
    saveRDS(object = list(reslist = reslist, errlist = errlist,
                          idx = idx, uids = uids),
            file = tmpf)
  }

  if (length(errlist) > 0){
    # Try retrieving remaining ids from Uniprot
    nerr <- Inf
    while(length(errlist) < nerr && length(errlist) > 0){
      nerr <- length(errlist)
      cat("\nTrying to retrieve", length(errlist), "entries from Uniprot")
      cc <- 0
      for (idx in errlist){
        # Try fetching data
        errk <- FALSE
        tryCatch({
          myurl <- paste0("https://www.uniprot.org/uniprot/",
                          uids[idx], ".fasta")
          x     <- utils::read.csv(myurl, header = FALSE, sep = ";")
        },
        warning = function(c) {errk <<- TRUE},
        error   = function(c) {errk <<- TRUE},
        finally = NULL)

        if(!errk){
          seq   <- paste0(x[-1, 1], collapse = "")
          #>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
          mystr <- strsplit(paste(x[1, ], collapse = ";"), "|", fixed = TRUE)[[1]][3]
          reslist[[idx]]$TSeq_seqtype  <- "protein"
          reslist[[idx]]$TSeq_accver   <- NA
          reslist[[idx]]$TSeq_taxid    <- strsplit(strsplit(mystr, " OX=")[[1]][2],
                                                   " ")[[1]][1]
          reslist[[idx]]$TSeq_orgname  <- strsplit(strsplit(mystr, " OS=")[[1]][2],
                                                   " OX=")[[1]][1]
          reslist[[idx]]$TSeq_defline  <- strsplit(mystr, " OS=")[[1]][1]
          reslist[[idx]]$TSeq_length   <- nchar(seq)
          reslist[[idx]]$TSeq_sequence <- seq

          reslist[[idx]]$UID <- uids[idx]

          dbstr <- strsplit(paste(x[1, ], collapse = ";"), "|", fixed = TRUE)[[1]][1]
          reslist[[idx]]$DB  <- ifelse(dbstr == ">tr",
                                       yes = "UniProtKB/TrEMBL",
                                       no = ifelse(dbstr == ">sp",
                                                   yes = "UniProtKB/Swiss-Prot",
                                                   no  = "Uniprot/other"))

        }
        # Print progress bar
        mypb(i = cc, max_i = length(errlist), t0 = t0, npos = 30)
        cc <- cc + 1
      }
      errlist <- which(sapply(reslist, function(x) {is.null(x$UID)}))
    }
  }

  if(length(errlist) > 0) reslist <- reslist[-errlist]
  errlist <- uids[errlist]

  df <- data.table::rbindlist(reslist,
                              use.names = TRUE,
                              fill      = TRUE)

  class(df) <- c(class(df), "protein_dt")

  # Save results to file
  if(!is.null(save_folder)){
    saveRDS(object = df,      file = df_file)
    saveRDS(object = errlist, file = errfile)
    if(file.exists(tmpf)) file.remove(tmpf)
  }

  cat("\nDone!\n", nrow(df), "proteins retrieved.\n",
      length(errlist), "retrieval errors.")

  invisible(df)
}
