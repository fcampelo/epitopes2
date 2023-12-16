#' Builds a partially labelled dataset from protein and epitope sequences
#'
#' This function extracts proteins sequences and labels them based on IEDB
#' epitope data.
#'
#' @inheritParams make_OrgSpec_datasets
#' @param prot_IDs IDs of proteins to be extracted for the data set
#'
#' @return Data frame containing the resulting dataset.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
make_proteins_dataset <- function(epitopes, proteins, taxonomy_list, prot_IDs,
                                  orgIDs          = NULL,
                                  removeIDs       = NULL,
                                  hostIDs         = NULL,
                                  min_epit        = 8,
                                  max_epit        = 25,
                                  only_exact      = FALSE,
                                  pos.mismatch.rm = "all",
                                  set.positive    = "mode",
                                  window_size     = 2 * min_epit - 1,
                                  max.N           = 2,
                                  save_folder     = "./",
                                  ncpus           = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  id_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(class(orgIDs)    %in% id_classes,
                          class(hostIDs)   %in% id_classes,
                          class(removeIDs) %in% id_classes,
                          assertthat::is.count(min_epit),
                          assertthat::is.count(max_epit),
                          min_epit <= max_epit,
                          pos.mismatch.rm %in% c("all", "align"),
                          set.positive    %in% c("any", "mode", "all"),
                          is.logical(only_exact) & length(only_exact) == 1,
                          assertthat::is.count(ncpus),
                          is.character(save_folder),
                          length(save_folder) == 1)

  if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
  # ========================================================================== #

  # Join and filter epitope/protein data
  jdf <- prepare_join_df(epitopes = epitopes, proteins = proteins,
                         min_epit = min_epit, max_epit = max_epit,
                         only_exact = only_exact,
                         pos.mismatch.rm = pos.mismatch.rm,
                         set.positive = set.positive)

  jdf <- filter_epitopes(jdf,
                         orgIDs = orgIDs,
                         removeIDs = removeIDs,
                         hostIDs   = hostIDs,
                         tax_list  = taxonomy_list)

  res_prot  <- label_proteins(proteins[proteins$UID %in% prot_IDs, ],
                              epitopes = jdf,
                              set.positive = set.positive,
                              ncpus = ncpus)
  res_prot  <- res_prot[, c("Info_UID", "Info_center_pos", "Class")]

  wres_prot <- make_window_df(proteins[proteins$UID %in% prot_IDs, ],
                              window_size = window_size, ncpus = ncpus)

  wres_prot <- dplyr::left_join(wres_prot, res_prot,
                                by = c("Info_UID", "Info_center_pos"))
  wres_prot <- calc_features(wres_prot, max.N = 2, ncpus = ncpus)

  saveRDS(wres_prot, paste0(save_folder, "/prots_df.rds"))
  seqinr::write.fasta(as.list(gsub("[BJXZ]", "", wres_prot$Info_window_seq)),
                      names = paste(wres_prot$Info_UID, wres_prot$Info_center_pos, sep = "pp"),
                      file.out = "./data/splits/prots_windows.fasta",
                      as.string = TRUE, nbchar = 10000)
  seqinr::write.fasta(as.list(res_prot$TSeq_sequence), names = res_prot$UID,
                      file.out = paste0(save_folder, "/prots.fasta"),
                      as.string = TRUE, nbchar = 10000)

  return(wres_prot)
}

