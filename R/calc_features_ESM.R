#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
calc_ESM_features(proteins,
                  py_script_path,
                  model = "esm1b_t33_650M_UR50S",
                  save_folder = "./",
                  max_seq = 1022,
                  step_size = 256){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(proteins),
                          "TSeq_sequence" %in% names(proteins),
                          "UID" %in% names(proteins),
                          is.character(model), length(model) == 1,
                          is.character(save_folder), length(save_folder) == 1)

  if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
  # ========================================================================== #

  ### PART 1: write FASTA files to pass to ESM model

  seqinr::write.fasta(as.list(proteins$TSeq_sequence),
                      names = proteins$UID,
                      file.out = paste0(save_folder, "/proteins.fa"))

  # ESM-1b returns errors when trying to predict embeddings for proteins longer
  # than 1022 residues (see https://github.com/facebookresearch/esm/issues/49).
  # A possible solution (suggested in https://github.com/facebookresearch/esm/issues/21#issuecomment-763217386)
  # is to break longer sequences in smaller chunks, predict independently and then
  # aggregate in post-processing, averaging where needed. We follow this strategy
  # here.

  prot1 <- proteins[proteins$TSeq_length <= max_seq, ]
  prot2 <- proteins[proteins$TSeq_length > max_seq, ]

  prot2_exp <- prot2[-(1:nrow(prot2)), ]

  prot2_exp <- mypblapply(seq_along(prot2$UID),
                          function(i){
                            ll <- prot2$TSeq_length[i]
                            nblocks <- 1 + ceiling((ll - max_seq) / stepsize)
                            tmp <- prot2[rep(i, nblocks), ]
                            for (k in 1:nblocks){
                              st <- 1 + (k - 1) * stepsize
                              en <- min(ll, max_seq + (k - 1) * stepsize)
                              if(en == ll) st <- ll - max_seq + 1
                              tmp$TSeq_length[k] <- en - st + 1
                              tmp$TSeq_sequence[k] <- substr(prot2$TSeq_sequence[i],
                                                             start = st,
                                                             stop = en)
                              tmp$UID[k] <- paste0(prot2$UID[i],
                                                   "__", st,
                                                   "_", en)
                              tmp$TSeq_accver[k] <- tmp$UID[k]
                            }
                            return(tmp)
                          }) %>%
    dplyr::bind_rows()

  proteins2 <- dplyr::bind_rows(prot1, prot2_exp)
  proteins2$TSeq_sequence <- gsub("B|Z|X|J", "<mask>", proteins2$TSeq_sequence)

  seqinr::write.fasta(as.list(proteins2$TSeq_sequence),
                      names = proteins2$UID,
                      file.out = paste0(save_folder, "/proteins_masked.fa"))

  # ========================================================================== #

  ### PART 2: invoke ESM model

  cmdline <- paste0("python ", py_script_path, "/extract.py ", model, " ",
                    save_folder, "/proteins_masked.fa ",
                    save_folder, "/proteins_features_", model,
                    "--include per_tok --nogpu --repr_layers 33")

}
