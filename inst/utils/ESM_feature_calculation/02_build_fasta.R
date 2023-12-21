require(dplyr)
require(seqinr)

# Sanity checks and initial definitions
assertthat::assert_that(is.data.frame(prot.df),
                        is.character(id_column), length(id_column) == 1,
                        is.character(seqs_column), length(seqs_column) == 1,
                        is.character(save_folder), length(save_folder) == 1,
                        assertthat::is.count(chunk_size),
                        assertthat::is.count(step_size),
                        all(c(seqs_column, id_column) %in% names(prot.df)))


if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
# ========================================================================== #

# Write FASTA files to pass to ESM model
SEQs <- prot.df %>%
  dplyr::select(dplyr::all_of(seqs_column)) %>%
  unlist() %>% unname()

IDs <- prot.df %>%
  dplyr::select(dplyr::all_of(id_column)) %>%
  unlist() %>% unname()

seqinr::write.fasta(as.list(SEQs),
                    names = IDs,
                    file.out = paste0(save_folder, "/proteins.fa"))


# Breat longer sequences for feature calculation and later mean-aggregation.
prot1 <- prot.df[nchar(SEQs) <= chunk_size, c(id_column, seqs_column)]
prot2 <- prot.df[nchar(SEQs) > chunk_size, c(id_column, seqs_column)]

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
                      ids        = prot2$Info_protein_id,
                      seqs       = prot2$Info_protein_sequence,
                      chunk_size = chunk_size,
                      step_size  = step_size) %>%
    dplyr::bind_rows() %>%
    dplyr::as_tibble()

  names(prot2_exp) <- names(prot1)
  prot1 <- dplyr::bind_rows(prot1, prot2_exp)

}

SEQs <- prot1 %>%
  dplyr::select(dplyr::all_of(seqs_column)) %>%
  unlist() %>% unname()

IDs <- prot1 %>%
  dplyr::select(dplyr::all_of(id_column)) %>%
  unlist() %>% unname()

SEQs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "<mask>", SEQs)


seqinr::write.fasta(as.list(SEQs),
                    names = IDs,
                    file.out = paste0(save_folder, "/proteins_masked.fa"))
