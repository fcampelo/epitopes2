# Function to extract FASTA file for ESM2 feature calculation
make_esm_fasta <- function(prot.df, save_folder,
                           break_pieces = TRUE,
                           chunk_size = 1022,
                           step_size = 512){

  if(any(duplicated(prot.df))) prot.df <- prot.df[-which(duplicated(prot.df)), ]

  prot.df$SEQs <- gsub("[^ABCDEFGHIKLMNOPQRSTUVWXYZ]", "<unk>", prot.df$SEQs)

  if (break_pieces){
    # Break longer sequences for feature calculation and later mean-aggregation.
    prot1 <- prot.df[nchar(prot.df$SEQs) <= chunk_size, ]
    prot2 <- prot.df[nchar(prot.df$SEQs) > chunk_size, ]

    if (nrow(prot2) > 0){
      prot2_exp <- prot2[-(1:nrow(prot2)), ]

      prot2_exp <- lapply(1:nrow(prot2),
                          function(i, ids, seqs, chunk_size, step_size){
                            ll      <- nchar(seqs[i])
                            nblocks <- 1 + ceiling((ll - chunk_size) / step_size)
                            tmp     <- data.frame(id = character(nblocks),
                                                  seq = character(nblocks))
                            for (k in 1:nblocks){
                              st <- 1 + (k - 1) * step_size
                              en <- min(ll, chunk_size + (k - 1) * step_size)
                              if(en == ll) st <- ll - chunk_size + 1
                              tmp$seq[k] <- substr(seqs[i],
                                                   start = st,
                                                   stop = en)
                              tmp$id[k]  <- paste0(ids[i],
                                                   "__", st,
                                                   "_", en)
                            }
                            return(tmp)
                          },
                          ids        = prot2$IDs,
                          seqs       = prot2$SEQs,
                          chunk_size = chunk_size,
                          step_size  = step_size) %>%
        dplyr::bind_rows() %>%
        dplyr::as_tibble()

      names(prot2_exp) <- names(prot1)
      prot1 <- dplyr::bind_rows(prot1, prot2_exp)
    }
  } else {
    prot1 <- prot.df
  }

  seqinr::write.fasta(as.list(prot1$SEQs),
                      names = prot1$IDs,
                      file.out = paste0(save_folder, "/proteins_masked_blocked.fa"))

  invisible(prot1)

}
