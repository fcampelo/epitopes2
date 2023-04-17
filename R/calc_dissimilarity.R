#' Invoke DIAMOND to calculate alignment-based similarity scores
#'
#' This function calculates the dissimilarity score between protein sequences
#' in data frame `proteins`. NOTE: this function requires
#' DIAMOND to be installed. Depending on the OS, this may mean a proper
#' system-wide installation or simply the presence of the executable in the
#' appropriate folder. Please check
#' [https://github.com/bbuchfink/diamond/wiki/2.-Installation](https://github.com/bbuchfink/diamond/wiki/2.-Installation)
#' for details.
#'
#' @section Dissimilarity calculation:
#' Dissimilarity is calculated based on the local alignment scores returned by
#' DIAMOND. Alignments shorter than `min_align` are ignored. The dissimilarity
#' between a pair of protein sequences is calculated as $1 - max(pident)/100$,
#' where `pident` is the vector of local alignment scores returned for that
#' pair of proteins. When no alignment is found for a given pair, the
#' dissimilarity is set to 1.
#'
#' @param proteins data frame of protein sequence data. Must have at least
#' two columns, `UID` (with the protein IDs) and `TSeq_sequence` (with the
#' protein sequences). This data frame is usually created using
#' [get_proteins()].
#' @param min_align smallest alignment size to consider.
#' @param save_folder path to folder for saving the results. Use `NULL`
#' if saving to file is not desired. Defaults to
#' the current working directory.
#' @param diamond_folder path to folder where DIAMOND can be run. Defaults to
#' the current working directory (which works, e.g., when a system-wide install
#' is present).
#' @param mop_up logical: should the DIAMOND files be deleted upon completion?
#'
#' @return Dissimilarity matrix for the proteins queried.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'
calc_dissimilarity <- function(proteins,
                               min_align = 6,
                               diamond_folder = "./",
                               save_folder = "./",
                               verbose = TRUE,
                               mop_up = FALSE){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(proteins),
                          "TSeq_sequence" %in% names(proteins),
                          "UID" %in% names(proteins),
                          is.character(diamond_folder),
                          length(diamond_folder) <= 1,
                          is.null(save_folder) | is.character(save_folder),
                          length(save_folder) <= 1,
                          assertthat::is.count(min_align),
                          is.logical(verbose), length(verbose) == 1)
  vrb <- verbose
  if(!dir.exists(diamond_folder)) dir.create(diamond_folder, recursive = TRUE)
  # ========================================================================== #
  mymessage(vrb, "Calculating similarities using DIAMOND...")

  seqinr::write.fasta(sequences = as.list(proteins$TSeq_sequence),
                      names     = proteins$UID,
                      file.out  = paste0(diamond_folder, "/proteins.fa"))

  cmdline <- paste0(diamond_folder, "/diamond makedb --in ",
                    diamond_folder, "/proteins.fa -d ",
                    diamond_folder, "/proteins-reference ",
                    ifelse(vrb, "--verbose", "--quiet"))
  system(cmdline)

  cmdline <- paste0(diamond_folder, "/diamond blastp -d ",
                    diamond_folder, "/proteins-reference -q ",
                    diamond_folder, "/proteins.fa --ultra-sensitive -b 1 -o ",
                    diamond_folder, "/protein-matches.tsv ",
                    ifelse(vrb, "--verbose", "--quiet"))

  system(cmdline)

  # read results
  scores <- utils::read.csv(paste0(diamond_folder, "/protein-matches.tsv"),
                            sep = "\t",
                            header = FALSE,
                            stringsAsFactors = FALSE)
  names(scores) <- c("qseqid", "sseqid", "pident", "length",
                     "mismatch", "gapopen", "qstart", "qend",
                     "sstart", "send", "evalue", "bitscore")

  scores <- scores %>%
    dplyr::filter(.data$length >= min_align) %>%
    dplyr::group_by(.data$qseqid, .data$sseqid) %>%
    dplyr::arrange(dplyr::desc(.data$pident)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::first),
                     .groups = "drop") %>%
    dplyr::mutate(diss = 1 - pident / 100) %>%
    dplyr::select("qseqid", "sseqid", "diss") %>%
    tidyr::pivot_wider(names_from = "sseqid",
                       values_from = "diss",
                       values_fill = 1) %>%
    dplyr::arrange(qseqid) %>%
    as.data.frame()

  rownames(scores) <- scores$qseqid

  # # Make the dissimilarity scores matrix
  scores <- scores %>%
    dplyr::select(order(colnames(scores)),
                  -"qseqid")

  pnames <- rownames(scores)
  scores <- as.matrix(scores)

  protIDs <- names(seqinr::read.fasta(paste0(diamond_folder, "/proteins.fa"),
                                      as.string = TRUE))
  missing <- protIDs[which(!(protIDs %in% rownames(scores)))]

  if(length(missing) > 0){
    scores[(nrow(scores) + 1):(nrow(scores) + length(missing)), ] <- 1
    scores[, (ncol(scores) + 1):(ncol(scores) + length(missing))] <- 1
    colnames(scores) <- c(pnames, missing)
    rownames(scores) <- c(pnames, missing)
  }

  if(!is.null(save_folder)){
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    saveRDS(scores, paste0(save_folder, "/protein_dissimilarity.rds"))
  }

  if(mop_up){
    file.remove(paste0(diamond_folder,
                       c("/protein-matches.tsv",
                         "/proteins-reference.dmnd",
                         "/proteins.fa")))
  }

  return(scores)
}
