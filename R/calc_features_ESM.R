#' Calculate classic features for epitope prediction
#'
#' This function is used to calculate features using the ESM embedders,
#' particularly ESM-1b. It is a **very** computationally demanding model,
#' and usually cannot run in regular computers (this function is intended
#' for use in high-performance workstations or clusters). The function also
#' needs a good internet connection to run, as the ESM model requires the
#' download of certain elements to be run.
#' To check the details and install instructions of ESM, please see
#' <https://github.com/facebookresearch/esm> and examples therein. A clone
#' of that repository is also available at <https://github.com/fcampelo/esm>.
#'
#' @section Sequence length and non-standard AAs:
#' ESM-1b returns errors when trying to calculate features for proteins longer
#' than 1024 residues (see https://github.com/facebookresearch/esm/issues/49).
#' A possible solution (suggested in, e.g.,
#' https://github.com/brianhie/evolocity/issues/2)
#' is to break longer sequences in smaller overlapping windows, predict
#' independently and then aggregate in post-processing, averaging where needed.
#' We follow this strategy here.
#'
#' For the feature calculation, any non-standard AA character is replaced by
#' the `<mask>` placeholder.
#'
#' @param prot.df dataframe containing a column with proteins
#' sequences
#' @param py_script_path path to the ESM Python script (NOTE: this is
#' provided as inst/utils/extract.py in this package's folder structure. The
#' script is a slightly modified version of the one available in
#' <https://github.com/facebookresearch/esm/blob/main/scripts/extract.py>).
#' @param id_column name of column in `prot.df` containing the unique protein
#' IDs.
#' @param seqs_column name of column in `prot.df` containing the protein
#' sequences
#' @param save_folder path to folder for saving the output.
#' @param model string with the full model name
#' (see <https://github.com/facebookresearch/esm>)
#' @param model.params string with model options to be used
#' @param chunk_size size of chunk to be used when proteins are
#' longer than 1024 residues (usually as large as possible and < 1024)
#' @param step_size step size to use for processing long proteins. Smaller
#' values are better, but more computationally intensive.
#' @param ncpus positive integer, number of cores to use
#' @param feat_prefix prefix to be added to the feature names.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'

calc_ESM_features <- function(prot.df,
                              py_script_path,
                              id_column  = "Info_protein_id",
                              seqs_column = "Info_protein_sequence",
                              save_folder = "./esm1b_features",
                              model = "esm1b_t33_650M_UR50S",
                              model.params = "--include per_tok --nogpu --repr_layers 33",
                              chunk_size = 1000,
                              step_size = 50,
                              ncpus = 1,
                              feat_prefix = "feat_esm1b_"){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(prot.df),
                          is.character(id_column), length(id_column) == 1,
                          is.character(seqs_column), length(seqs_column) == 1,
                          is.character(py_script_path),
                          length(py_script_path) == 1,
                          file.exists(py_script_path),
                          is.character(save_folder), length(save_folder) == 1,
                          is.character(model), length(model) == 1,
                          is.character(model.params), length(model.params) == 1,
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

  # ========================================================================== #

  # Invoke ESM model

  cmdline <- paste0("python ", py_script_path, " ", model, " ",
                    save_folder, "/proteins_masked.fa ",
                    save_folder, "/proteins_features_", model,
                    " ", model.params)

  system(cmdline)


  # Concatenate output
  concatenate_ESM_output(save_folder, feat_prefix, ncpus)

  return(TRUE)

}
