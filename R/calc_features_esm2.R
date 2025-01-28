#'
#'
#'
#'
#'
#' @author Felipe Campelo (\email{fcampelo@@gmail.com})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#'

# TODO: Check https://cran.r-project.org/web/packages/reticulate/vignettes/python_dependencies.html
# to add Python environment configurations to package epitopes

calc_features_esm2 <- function(X,
                               mode = "run",
                               seqs_column = "Info_protein_sequence",
                               ids_column  = "Info_protein_id",
                               ncpus = 1,
                               save_folder = "./data",
                               model_spec  = "esm2_t33_650M_UR50D",
                               model_opts  = "--include per_tok --repr_layers 33",
                               delete_protein_csv_files = FALSE){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(X) | is.peptide.list(X),
                          is.character(seqs_column),
                          length(seqs_column) == 1,
                          assertthat::is.count(ncpus))

  if(is.data.frame(X)){
    assertthat::assert_that(nrow(X) > 0,
                            seqs_column %in% names(X))
  } else {
    assertthat::assert_that(nrow(X$proteins) > 0,
                            seqs_column %in% names(X$proteins))
  }

  # Create save folder if needed
  if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)


  if(is.peptide.list(X)) prot.df <- X$proteins else prot.df <- X

  prot.df <- data.frame(IDs  = as.character(prot.df[, which(names(prot.df) == ids_column), drop = TRUE]),
                        SEQs = as.character(prot.df[, which(names(prot.df) == seqs_column), drop = TRUE]))

  # ========================================================================== #

  # Build FASTA file for ESM2 calculations. This breaks up
  # proteins longer than 1022 AA into chunks of 1022 with a step size of 512.
  # The results are later aggregated using concatenate_esm_outputs()
  prot1 <- make_esm_fasta(prot.df,
                          save_folder = save_folder,
                          chunk_size  = 1022,
                          step_size   = 512)

  scriptpath <- system.file("utils/Python/extract_ESM_pertoken.py", package = "epitopes")

  cmdline1 <- paste0("python3 ", scriptpath, " ", model_spec, " ", save_folder,
                     "/proteins_masked_blocked.fa ", save_folder, "/csv ",
                     model_opts)

  if(mode == "run"){
    message(sprintf("\nCalling ESM2 model (%s) in Python.\nThis may take a while...", model_spec))
    system(cmdline1, invisible = FALSE)

    message(sprintf("\nConcatenating ESM2 output.\nThis may take a while..."))
    protfeats <- concatenate_esm_outputs(csv_folder = paste0(save_folder, "/csv"),
                                         filenames  = paste0(prot1$IDs, ".csv"),
                                         save_folder = save_folder,
                                         ncpus = ncpus,
                                         delete_originals = delete_protein_csv_files)

    if(is.data.frame(X)){

    } else {
      X$df <- X$df %>%
        dplyr::select(-dplyr::starts_with("feat_esm2_")) %>%
        dplyr::left_join(protfeats, by = c("Info_protein_id", "Info_pos"))

      return(X)
    }

  } else {
    message("FASTA file built for ESM processing. To generate the features:")
    message("\n 1) Run the following on the command line:\n> ", cmdline1)
    message("\n 2) Run the following in R:\n> protfeats <- concatenate_esm_outputs(csv_folder = '<save_folder>/csv', save_folder = '<save_folder>', ...)")
    message("\n 3) Left-join protfeats onto your data frame using protein ID and position as joining variables")
    return(NULL)
  }
}
