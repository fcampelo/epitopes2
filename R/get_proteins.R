#' Retrieve protein sequences and data from GenBank and Uniprot
#'
#' This function is used to retrieve data from Genbank's protein database for
#' given protein IDs. If an ID is not available from Genbank the function will
#' try to retrieve it from Uniprot.
#'
#' Queries are processed one by one (rather than in batch) to enable treatment
#' of individual inconsistencies (e.g., wrong UIDs, queries that return a
#' different identifier, etc.). This makes this routine substantially slower,
#' but considerably more robust to errors.
#'
#' @param uids A character vector with protein IDs.
#' @param DBs A character vector listing the databases to be searched.
#' Valid entries are "ncbi", "uniprot" and "uniprot-archived". Anything else is
#' ignored.
#' @param blocksize size of the retrieval blocks. Needs to be kept below 500.
#' @param save_folder path to folder for saving the results.
#' @param block.timeout positive integer: timeout for trying to retrieve each
#' block
#'
#' @return A data frame with the extracted proteins.
#'
#' @author Felipe Campelo (\email{fcampelo@@gmail.com})
#'
#' @export
#'

get_proteins <- function(uids,
                         DBs = c("ncbi", "uniprot", "uniprot-archived"),
                         blocksize = 250, save_folder = NULL,
                         block.timeout = max(60, blocksize)){

  # ========================================================================== #
  # Sanity checks and initial definitions
  DBs <- tolower(DBs)
  assertthat::assert_that(is.null(save_folder) | is.character(save_folder),
                          length(save_folder) <= 1,
                          is.character(uids),
                          length(uids) >= 1,
                          is.character(DBs), length(DBs) >= 1,
                          any(DBs %in% c("ncbi", "uniprot", "uniprot-archived")),
                          assertthat::is.count(blocksize),
                          assertthat::is.count(block.timeout),
                          blocksize < 500)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    df_file <- paste0(normalizePath(save_folder), "/proteins_", as.character(Sys.Date()), ".rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/protein_retrieval_errors_", as.character(Sys.Date()), ".rds")
    tmpf    <- tempfile(pattern = "get_proteins_tmpfile_", fileext = ".rds", tmpdir = save_folder)

    if(file.exists(tmpf)) file.remove(tmpf)
  }

  if(any(is.na(uids))) uids <- uids[-which(is.na(uids))]
  uids <- unique(uids)
  prots <- data.frame(Info_protein_id           = uids,
                      Info_protein_id_clean     = gsub("\\.[0-9]+$", "", uids),
                      Info_protein_version      = NA,
                      Info_protein_all_ids      = NA,
                      Info_protein_sequence     = NA,
                      Info_protein_database     = NA)

  queries <- unique(prots$Info_protein_id_clean)
  nq <- length(queries) + 1
  idx <- lapply((0:floor(length(queries) / blocksize)),
                function(i){
                  unique(pmin(length(queries), (i*blocksize + 1):((i+1)*blocksize)))})
  if(!(length(queries) %% blocksize)) idx <- idx[-length(idx)]

  # ========================================================================== #

  ## Try retrieving from NCBI/protein
  if("ncbi" %in% DBs && length(queries) > 0){
    while(nq > length(queries)){
      reslist <- vector("list", length(idx))
      nq <- length(queries)
      message("\nNCBI-Protein: Trying to retrieve ", nq, " proteins in ", length(idx), " blocks\n")

      for (i in seq_along(idx)){
        t0 <- Sys.time()
        message(sprintf("\rBlock %03d of %03d: Started on %s", i, length(idx), as.character(t0)))
        # Try fetching data
        tryCatch({
          R.utils::withTimeout(
            {
              x <- rentrez::entrez_fetch(db = "protein",
                                         id = queries[idx[[i]]],
                                         retmode = "xml",
                                         rettype = "xml")
              x <- XML::xmlToList(x)

              if(length(x) == 1){
                reslist[[i]] <- data.frame(
                  Info_protein_id_clean = nullcheck(x$`GBSeq_primary-accession`),
                  Info_protein_version  = nullcheck(x$`GBSeq_accession-version`),
                  Info_protein_all_ids  = nullcheck(paste(x$`GBSeq_other-seqids`, collapse = ";")),
                  Info_protein_sequence = nullcheck(toupper(x$GBSeq_sequence)))
              } else {
                reslist[[i]] <- data.frame(
                  Info_protein_id_clean = sapply(x, function(c) {nullcheck(c$`GBSeq_primary-accession`)}),
                  Info_protein_version  = sapply(x, function(c) {nullcheck(c$`GBSeq_accession-version`)}),
                  Info_protein_all_ids  = sapply(x, function(c) {nullcheck(paste(c$`GBSeq_other-seqids`, collapse = ";"))}),
                  Info_protein_sequence = sapply(x, function(c) {nullcheck(toupper(c$GBSeq_sequence))})
                )
              }
            }, timeout = block.timeout)
        },
        TimeoutException  = function(c) message("\n\nTimeout - consider increasing block.timeout"),
        warning = function(c) message("Warning happened"),
        error   = function(c) message("Error happened"),
        finally = NULL)
      }

      retrieved <- dplyr::bind_rows(reslist)
      if(nrow(retrieved) == 0) break

      # Remove duplicates
      ii <- which(duplicated(retrieved))
      if(length(ii) > 0) retrieved <- retrieved[-ii, ]

      # Remove NAs
      ii <- which(is.na(retrieved$Info_protein_sequence))
      if(length(ii) > 0) retrieved <- retrieved[-ii, ]

      if(nrow(retrieved) == 0) break

      # Get indices of each query on the retrieved dataframe
      ii <- sapply(queries,
                   function(id) {
                     k <- grep(id, retrieved$Info_protein_all_ids, fixed = TRUE)
                     ifelse(length(k) == 1, k, NA)
                   },
                   simplify = TRUE)

      # Remove NAs
      if(any(is.na(ii))) ii <- ii[-which(is.na(ii))]

      # Build index map and merge retrieved dataframe
      retdf <- dplyr::left_join(
        data.frame(Info_protein_id_clean = names(ii),
                   ret.id = retrieved$Info_protein_id_clean[ii]),
        retrieved,
        by = c("ret.id" = "Info_protein_id_clean"))

      # Get indices for merging into final prots dataframe
      ii <- sapply(prots$Info_protein_id_clean,
                   function(id) {
                     k <- grep(id, retdf$Info_protein_all_ids, fixed = TRUE)
                     ifelse(length(k) == 1, k, NA)
                   },
                   simplify = TRUE)
      ii <- data.frame(prots.idx = 1:nrow(prots),
                       retdf.idx = ii)
      if(any(is.na(ii$retdf.idx))) ii <- ii[!is.na(ii$retdf.idx), ]

      # Update final prots dataframe
      prots$Info_protein_version[ii$prots.idx]  <- retdf$Info_protein_version[ii$retdf.idx]
      prots$Info_protein_all_ids[ii$prots.idx]  <- retdf$Info_protein_all_ids[ii$retdf.idx]
      prots$Info_protein_sequence[ii$prots.idx] <- retdf$Info_protein_sequence[ii$retdf.idx]
      prots$Info_protein_database[ii$prots.idx] <- "NCBI-Protein"

      # Extract remaining (not retrieved) ids
      queries <- unique(prots$Info_protein_id_clean[which(is.na(prots$Info_protein_sequence))])
      if(length(queries) == 0) break

      # Update blocksizes
      blocksize <- ceiling(blocksize * length(queries) / nq)
      idx <- lapply((0:floor(length(queries) / (blocksize))),
                    function(i){
                      unique(pmin(length(queries), (i*blocksize + 1):((i+1)*blocksize)))})
      if(!(length(queries) %% blocksize)) idx <- idx[-length(idx)]

      if(!is.null(save_folder)) saveRDS(object = prots, file = tmpf)
    }
    message("\rNCBI-Protein: Finished!\t\t\t\t\t\t\t\t")
  }



  ## ============ Try retrieving ids from Uniprot
  if("uniprot" %in% DBs && length(queries) > 0){

    blocksize <- min(25, blocksize)
    idx <- lapply((0:floor(length(queries) / blocksize)),
                  function(i){
                    unique(pmin(length(queries), (i*blocksize + 1):((i+1)*blocksize)))})
    if(!(length(queries) %% blocksize)) idx <- idx[-length(idx)]

    nq <- length(queries) + 1

    while(nq > length(queries)){
      reslist <- vector("list", length(idx))
      nq <- length(queries)
      message("\nUniprotKB: Trying to retrieve ", nq, " proteins in ", length(idx), " blocks\n")

      for (i in seq_along(idx)){
        t0 <- Sys.time()
        message(sprintf("\rBlock %03d of %03d: Started on %s", i, length(idx), as.character(t0)))
        # Try fetching data
        tryCatch({
          R.utils::withTimeout(
            {
              tmp <- lapply(queries[idx[[i]]],
                            function(id){
                              tryCatch({
                                myurl <- paste0("https://rest.uniprot.org/uniprotkb/",
                                                id, ".fasta")
                                seqs <- protr::readFASTA(myurl, seqonly = FALSE)
                                data.frame(Info_protein_id_clean = id,
                                           Info_protein_all_ids  = names(seqs),
                                           Info_protein_sequence = unname(seqs)[[1]])},
                                warning = function(c) cat("!"),
                                error   = function(c) data.frame(Info_protein_id_clean = character(),
                                                                 Info_protein_all_ids  = character(),
                                                                 Info_protein_sequence = character()),
                                finally = NULL)})

              reslist[[i]] <- dplyr::bind_rows(tmp)
            }, timeout = block.timeout)
        },
        TimeoutException  = function(c) message("\n\nTimeout - consider increasing block.timeout"),
        warning = function(c) message("Warning(s) happened"),
        error   = function(c) message("Error happened"),
        finally = NULL)

      }

      retrieved <- dplyr::bind_rows(reslist)
      if(nrow(retrieved) == 0) break

      # Remove duplicates
      ii <- which(duplicated(retrieved))
      if(length(ii) > 0) retrieved <- retrieved[-ii, ]

      # Remove NAs
      ii <- which(is.na(retrieved$Info_protein_sequence))
      if(length(ii) > 0) retrieved <- retrieved[-ii, ]

      if(nrow(retrieved) == 0) break

      # Get indices of each query on the retrieved dataframe
      ii <- sapply(queries,
                   function(id) {
                     k <- grep(id, retrieved$Info_protein_all_ids, fixed = TRUE)
                     ifelse(length(k) == 1, k, NA)
                   },
                   simplify = TRUE)

      # Remove NAs
      if(any(is.na(ii))) ii <- ii[-which(is.na(ii))]

      # Build index map and merge retrieved dataframe
      retdf <- dplyr::left_join(
        data.frame(Info_protein_id_clean = names(ii),
                   ret.id = retrieved$Info_protein_id_clean[ii]),
        retrieved,
        by = c("ret.id" = "Info_protein_id_clean"))

      # Get indices for merging into final prots dataframe
      ii <- sapply(prots$Info_protein_id_clean,
                   function(id) {
                     k <- grep(id, retdf$Info_protein_all_ids, fixed = TRUE)
                     ifelse(length(k) == 1, k, NA)
                   },
                   simplify = TRUE)
      ii <- data.frame(prots.idx = 1:nrow(prots),
                       retdf.idx = ii)
      if(any(is.na(ii$retdf.idx))) ii <- ii[!is.na(ii$retdf.idx), ]

      # Update final prots dataframe
      prots$Info_protein_version[ii$prots.idx]  <- retdf$Info_protein_id_clean[ii$retdf.idx]
      prots$Info_protein_all_ids[ii$prots.idx]  <- retdf$Info_protein_all_ids[ii$retdf.idx]
      prots$Info_protein_sequence[ii$prots.idx] <- retdf$Info_protein_sequence[ii$retdf.idx]
      prots$Info_protein_database[ii$prots.idx] <- "UniprotKB"

      # Extract remaining (not retrieved) ids
      queries <- unique(prots$Info_protein_id_clean[which(is.na(prots$Info_protein_sequence))])

      # Update blocksizes
      blocksize <- ceiling(blocksize * length(queries) / nq)
      idx <- lapply((0:floor(length(queries) / blocksize)),
                    function(i){
                      unique(pmin(length(queries), (i*blocksize + 1):((i+1)*blocksize)))})
      if(!(length(queries) %% blocksize)) idx <- idx[-length(idx)]

      if(!is.null(save_folder)) saveRDS(object = prots, file = tmpf)
    }
    message("\rUniprotKB: Finished!\t\t\t\t\t\t\t\t")
  }

  ## ============ Try retrieving ids from Uniprot-Archived
  if("uniprot-archived" %in% DBs && length(queries) > 0){

    blocksize <- min(25, blocksize)
    idx <- lapply((0:floor(length(queries) / blocksize)),
                  function(i){
                    unique(pmin(length(queries), (i*blocksize + 1):((i+1)*blocksize)))})
    if(!(length(queries) %% blocksize)) idx <- idx[-length(idx)]

    nq <- length(queries) + 1

    while(nq > length(queries)){
      reslist <- vector("list", length(idx))
      nq <- length(queries)
      message("\nUniprotKB-Archived: Trying to retrieve ", nq, " proteins in ", length(idx), " blocks\n")

      for (i in seq_along(idx)){
        t0 <- Sys.time()
        message(sprintf("\rBlock %03d of %03d: Started on %s", i, length(idx), as.character(t0)))
        # Try fetching data
        tryCatch({
          R.utils::withTimeout(
            {
              tmp <- lapply(queries[idx[[i]]],
                            function(id){
                              tryCatch({
                                myurl <- paste0("https://rest.uniprot.org/unisave/",
                                                id, "?format=fasta&versions=1")
                                seqs <- protr::readFASTA(myurl, seqonly = FALSE)
                                data.frame(Info_protein_id_clean = id,
                                           Info_protein_all_ids  = names(seqs),
                                           Info_protein_sequence = unname(seqs)[[1]])},
                                warning = function(c) cat("!"),
                                error   = function(c) data.frame(Info_protein_id_clean = character(),
                                                                 Info_protein_all_ids  = character(),
                                                                 Info_protein_sequence = character()),
                                finally = NULL)})

              reslist[[i]] <- dplyr::bind_rows(tmp)
            }, timeout = block.timeout)
        },
        TimeoutException  = function(c) message("\n\nTimeout - consider increasing block.timeout"),
        warning = function(c) message("Warning(s) happened"),
        error   = function(c) message("Error happened"),
        finally = NULL)

      }

      retrieved <- dplyr::bind_rows(reslist)
      if(nrow(retrieved) == 0) break

      # Remove duplicates
      ii <- which(duplicated(retrieved))
      if(length(ii) > 0) retrieved <- retrieved[-ii, ]

      # Remove NAs
      ii <- which(is.na(retrieved$Info_protein_sequence))
      if(length(ii) > 0) retrieved <- retrieved[-ii, ]

      if(nrow(retrieved) == 0) break

      # Get indices of each query on the retrieved dataframe
      ii <- sapply(queries,
                   function(id) {
                     k <- grep(id, retrieved$Info_protein_all_ids, fixed = TRUE)
                     ifelse(length(k) == 1, k, NA)
                   },
                   simplify = TRUE)

      # Remove NAs
      if(any(is.na(ii))) ii <- ii[-which(is.na(ii))]

      # Build index map and merge retrieved dataframe
      retdf <- dplyr::left_join(
        data.frame(Info_protein_id_clean = names(ii),
                   ret.id = retrieved$Info_protein_id_clean[ii]),
        retrieved,
        by = c("ret.id" = "Info_protein_id_clean"))

      # Get indices for merging into final prots dataframe
      ii <- sapply(prots$Info_protein_id_clean,
                   function(id) {
                     k <- grep(id, retdf$Info_protein_all_ids, fixed = TRUE)
                     ifelse(length(k) == 1, k, NA)
                   },
                   simplify = TRUE)
      ii <- data.frame(prots.idx = 1:nrow(prots),
                       retdf.idx = ii)
      if(any(is.na(ii$retdf.idx))) ii <- ii[!is.na(ii$retdf.idx), ]

      # Update final prots dataframe
      prots$Info_protein_version[ii$prots.idx]  <- retdf$Info_protein_id_clean[ii$retdf.idx]
      prots$Info_protein_all_ids[ii$prots.idx]  <- retdf$Info_protein_all_ids[ii$retdf.idx]
      prots$Info_protein_sequence[ii$prots.idx] <- retdf$Info_protein_sequence[ii$retdf.idx]
      prots$Info_protein_database[ii$prots.idx] <- "UniprotKB-Archived/Deleted"

      # Extract remaining (not retrieved) ids
      queries <- unique(prots$Info_protein_id_clean[which(is.na(prots$Info_protein_sequence))])

      # Update blocksizes
      blocksize <- ceiling(blocksize * length(queries) / nq)
      idx <- lapply((0:floor(length(queries) / blocksize)),
                    function(i){
                      unique(pmin(length(queries), (i*blocksize + 1):((i+1)*blocksize)))})
      if(!(length(queries) %% blocksize)) idx <- idx[-length(idx)]

      if(!is.null(save_folder)) saveRDS(object = prots, file = tmpf)
    }

    message("\rUniprotKB-Archived: Finished!\t\t\t\t\t\t\t\t")
  }

  # Give it one last try based on a common naming inconsistency
  ii <- which(is.na(prots$Info_protein_sequence))
  tmp <- prots[ii, ]
  tmp$Info_protein_id_clean <- gsub("\\_.$", "", tmp$Info_protein_id)

  if(any(tmp$Info_protein_id_clean != tmp$Info_protein_id)){
    lasttry <- get_proteins(tmp$Info_protein_id_clean, blocksize = 1)
    ii <- ii[which(!is.na(lasttry$Info_protein_sequence))]
    lasttry <- lasttry[which(!is.na(lasttry$Info_protein_sequence)), ]

    prots$Info_protein_version[ii]  <- lasttry$Info_protein_id_clean
    prots$Info_protein_all_ids[ii]  <- lasttry$Info_protein_all_ids
    prots$Info_protein_sequence[ii] <- lasttry$Info_protein_sequence
    prots$Info_protein_database[ii] <- lasttry$Info_protein_database
  }

  errlist <- prots$Info_protein_id[which(is.na(prots$Info_protein_sequence))]

  if(any(is.na(prots$Info_protein_sequence))) prots <- prots[-which(is.na(prots$Info_protein_sequence)), ]

  if(!is.null(save_folder)) {
    saveRDS(object = prots, file = df_file)
    if(length(errlist) > 0) saveRDS(object = errlist, file = errfile)
    if(file.exists(tmpf)) file.remove(tmpf)
  }

  return(prots)
}

