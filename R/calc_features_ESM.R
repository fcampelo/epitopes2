#'
#' ESM-1b returns errors when trying to predict embeddings for proteins longer
#' than 1022 residues (see https://github.com/facebookresearch/esm/issues/49).
#' A possible solution (suggested in
#' https://github.com/facebookresearch/esm/issues/21#issuecomment-763217386)
#' is to break longer sequences in smaller chunks, predict independently and then
#' aggregate in post-processing, averaging where needed. We follow this strategy
#' here.
#'
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'



calc_ESM_features(prot.df,
                  py_script_path,
                  id_column  = "Info_protein_id",
                  seqs_column = "Info_protein_sequence",
                  save_folder = "./esm_features",
                  model = "esm1b_t33_650M_UR50S",
                  max_seq = 1022,
                  step_size = 50){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(prot.df),
                          is.character(id_column), length(id_column) == 1,
                          is.character(seqs_column), length(seqs_column) == 1,
                          is.character(py_script_path),
                          length(py_script_path) == 1,
                          file.exists(py_script_path),
                          is.character(save_folder),
                          length(save_folder) == 1,
                          is.character(model), length(model) == 1,
                          assertthat::is.count(max_seq),
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
  prot1 <- prot.df[nchar(SEQs) <= max_seq, c(id_column, seqs_column)]
  prot2 <- prot.df[nchar(SEQs) > max_seq, c(id_column, seqs_column)]

  if (nrow(prot2) > 0){
    prot2_exp <- prot2[-(1:nrow(prot2)), ]

    prot2_exp <- lapply(1:nrow(prot2),
                        function(i, ids, seqs, max_seq, step_size){
                          ll <- nchar(seqs[i])
                          nblocks <- 1 + ceiling((ll - max_seq) / step_size)
                          tmp <- data.frame(id = character(nblocks),
                                            seq = character(nblocks))
                          for (k in 1:nblocks){
                            st <- 1 + (k - 1) * step_size
                            en <- min(ll, max_seq + (k - 1) * step_size)
                            if(en == ll) st <- ll - max_seq + 1
                            tmp$seq[k] <- substr(seqs[i],
                                                 start = st,
                                                 stop = en)
                            tmp$id[k] <- paste0(ids[i],
                                                "__", st,
                                                "_", en)
                          }
                          return(tmp)
                        },
                        ids = prot2$Info_protein_id,
                        seqs = prot2$Info_protein_sequence,
                        max_seq = max_seq,
                        step_size = step_size) %>%
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

  SEQs <- gsub("B|J|O|U|X|Z", "<mask>", SEQs)


  seqinr::write.fasta(as.list(SEQs),
                      names = IDs,
                      file.out = paste0(save_folder, "/proteins_masked.fa"))

  # ========================================================================== #

  ### PART 2: invoke ESM model

  cmdline <- paste0("python ", py_script_path, " ", model, " ",
                    save_folder, "/proteins_masked.fa ",
                    save_folder, "/proteins_features_", model,
                    "--include per_tok --nogpu --repr_layers 33")

}
