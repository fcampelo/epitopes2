# retrieve_single_protein <- function(uid, wait = 0){
#
#   # Try retrieving from GenBank
#   tmpfile <- tempfile("tmp_", tmpdir = ".", fileext = ".xml")
#   errk    <- FALSE
#   tryCatch(
#     {reutils::efetch(uid     = uid,
#                      db      = "protein",
#                      rettype = "fasta",
#                      outfile = tmpfile);
#       prot_row <- XML::xmlToDataFrame(XML::xmlParse(tmpfile),
#                                       stringsAsFactors = FALSE)},
#     warning = function(c) {errk <<- TRUE},
#     error   = function(c) {errk <<- TRUE},
#     finally = NULL)
#
#   # Clean up tmpfile
#   if(file.exists(tmpfile)) file.remove(tmpfile)
#
#
#   # If error try on uniprot:
#   if (!errk) {
#     prot_row$source_DB <- "Genbank"
#   } else {
#     errk <- FALSE
#     tryCatch({
#       myurl <- paste0("https://www.uniprot.org/uniprot/", uid, ".fasta")
#       x     <- read.csv(myurl, header = TRUE)
#       x     <- paste0(x[, 1], collapse = "")
#       prot_row <- data.frame(TSeq_seqtype  = NA,
#                              TSeq_accver   = uid,
#                              TSeq_taxid    = NA,
#                              TSeq_orgname  = NA,
#                              TSeq_defline  = NA,
#                              TSeq_length   = nchar(x),
#                              TSeq_sequence = x,
#                              source_DB = "Uniprot",
#                              stringsAsFactors = FALSE)},
#       warning = function(c) {errk <<- TRUE},
#       error   = function(c) {errk <<- TRUE},
#       finally = NULL)
#   }
#
#   # If error in retrieval return NULL
#   if (errk) return(NULL)
#
#   prot_row$molecule_id <- uid
#   if (wait > 0) Sys.sleep(wait * (1 + stats::runif(1)))
#
#   return(prot_row)
# }




#
#
# get_proteins <- function(uids, save_folder = NULL){
#
#   # ========================================================================== #
#   # Sanity checks and initial definitions
#   assertthat::assert_that(is.character(uids),
#                           length(uids) >= 1,
#                           is.null(save_folder) | (is.character(save_folder)),
#                           is.null(save_folder) | length(save_folder) == 1)
#
#
#   # Check save folder and create file names
#   if(!is.null(save_folder)) {
#     if(!dir.exists(save_folder)) dir.create(save_folder)
#     df_file <- paste0(normalizePath(save_folder), "/proteins.rds")
#     errfile <- paste0(normalizePath(save_folder), "/proteins_retrieval_errors.rds")
#   }
#
#   # Retrieving proteins using individual requests rather than (more efficient)
#   # batch requests, to catch and treat efetch() or parsing errors more easily.
#   cat("\nRetrieving proteins:\n")
#   df <- pbmcapply::pbmclapply(X        = uids,
#                               FUN      = retrieve_single_protein,
#                               wait     = 0.1,
#                               mc.cores = 1)
#
#   errlist <- uids[sapply(df, is.null)]
#   df      <- data.frame(data.table::rbindlist(df,
#                                               use.names = TRUE,
#                                               fill      = TRUE))
#
#   # Save proteins to file and clean up outfile
#   if(!is.null(save_folder)){
#     saveRDS(object = df,      file = df_file)
#     saveRDS(object = errlist, file = errfile)
#   }
#
#   invisible(list(proteins = df,
#                  errlist  = errlist))
# }




# prepare_protein_for_prediction <- function(proteins,
#                                            window_size,
#                                            save_folder  = NULL,
#                                            ncpus        = 1){
#
#   # ========================================================================== #
#   # Sanity checks and initial definitions
#   assertthat::assert_that(assertthat::is.count(window_size),
#                           is.data.frame(proteins),
#                           is.null(save_folder) | (is.character(save_folder)),
#                           is.null(save_folder) | length(save_folder) == 1,
#                           assertthat::is.count(ncpus))
#
#   # Set up parallel processing
#   if ((.Platform$OS.type == "windows") & (ncpus > 1)){
#     cat("\nAttention: multicore not currently available for Windows.\n
#         Forcing ncpus = 1.")
#     ncpus <- 1
#   } else {
#     available.cores <- parallel::detectCores()
#     if (ncpus >= available.cores){
#       cat("\nAttention: ncpus too large, we only have ", available.cores,
#           " cores.\nUsing ", available.cores - 1,
#           " cores for run_experiment().")
#       ncpus <- available.cores - 1
#     }
#   }
#
#   # Check save folder and create file names
#   if(!is.null(save_folder)) {
#     if(!dir.exists(save_folder)) dir.create(save_folder)
#     df_file <- paste0(normalizePath(save_folder), "/df_windowed.rds")
#   }
#
#   # ========================================================================== #
#   # Initial preprocessing
#
#   # Join epitopes and proteins dataframes, preliminary feature transformation
#   df <- proteins %>%
#     dplyr::transmute(protein_id    = as.character(.data$molecule_id),
#                      protein_seq   = as.character(.data$TSeq_sequence),
#                      protein_len   = nchar(.data$TSeq_sequence),
#                      protein_taxid = as.character(.data$TSeq_taxid),
#                      org_name      = as.character(.data$TSeq_orgname))
#
#   # ========================================================================== #
#   # Generate dataframe by sliding windows
#   # get_prot_window() is defined at the bottom of this file.
#   cat("\nExtracting windows:\n")
#   mydf <- pbmcapply::pbmclapply(X = purrr::pmap(as.list(df), list),
#                                 FUN            = get_prot_window,
#                                 ws             = window_size,
#                                 mc.cores       = ncpus,
#                                 mc.preschedule = FALSE)
#
#   cat("\nAssembling windowed dataframe...")
#   mydf <- data.frame(data.table::rbindlist(mydf))
#
#   # Save resulting dataframe and error IDs to file
#   if(!is.null(save_folder)) {
#     saveRDS(mydf, file = df_file)
#   }
#
#   invisible(mydf)
# }
#
# get_prot_window <- function(x, ws){
#   wdf <- data.frame(window_seq = rep(NA_character_, x$protein_len),
#                     taxid      = x$protein_taxid,
#                     protein_id = x$protein_id,
#                     center_pos = rep(NA_integer_, x$protein_len),
#                     stringsAsFactors = FALSE)
#   for (i in 1:x$protein_len){
#     i1   <- max(1, min(i - floor(ws / 2), x$protein_len - ws + 1))
#     i2   <- min(x$protein_len, i1 + ws - 1)
#     wdf$window_seq[i] <- substr(x$protein_seq, i1, i2)
#     wdf$center_pos[i] <- i
#   }
#   return(wdf)
# }



# process_individual_epitope_T <- function(ep, file_id){
#
#   not_valid <- FALSE
#   # If it is not a linear T-Cell epitope, then ignore
#   # -->>> ASSUMPTION: All assays of same type
#   # -->>> ASSUMPTION: All linear T-Cell epitopes appear as
#   # -->>> "FragmentOfANaturalSequenceMolecule" and contain a
#   # -->>> field named "LinearSequence"
#   not_valid <- not_valid | is.null(ep$Assays$TCell)
#   not_valid <- not_valid | is.null(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence)
#
#   if(not_valid) return(NULL)
#
#   # ============= ONLY LINEAR B-CELL EPITOPES CROSS THIS LINE ============= #
#
#   # Extract relevant fields.
#   host_id        <- nullcheck(ep$Assays$TCell$Immunization$HostOrganism$OrganismId)
#   sourceOrg_id   <- nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceOrganismId)
#   epitope_id     <- nullcheck(ep$EpitopeId)
#   molecule_id    <- nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceMolecule$GenBankId)
#   seq            <- nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence)
#   evid_code      <- nullcheck(ep$EpitopeEvidenceCode)
#   epit_struc_def <- nullcheck(ep$EpitopeStructureDefines)
#   qual_measure   <- nullcheck(ep$Assays$TCell$AssayInformation$QualitativeMeasurement)
#
#
#   # Extract start and end positions, checking against sequence length
#   stpos  <- as.numeric(nullcheck(ep$ReferenceStartingPosition))
#   endpos <- as.numeric(nullcheck(ep$ReferenceEndingPosition))
#   sqlen  <- ifelse(is.na(seq),
#                    yes = endpos - stpos + 1,
#                    no  = nchar(seq))
#   if (any(is.na(c(stpos, endpos))) || (endpos - stpos + 1 != sqlen)){
#     stpos  <- as.numeric(nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$StartingPosition))
#     endpos <- as.numeric(nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$EndingPosition))
#   }
#
#   sqlen  <- ifelse(is.na(seq),
#                    yes = endpos - stpos + 1,
#                    no  = nchar(seq))
#   if (all(!is.na(c(stpos, endpos))) && (endpos - stpos + 1 == sqlen)){
#     start_pos <- stpos
#     end_pos   <- endpos
#   }
#
#   if(!exists("start_pos")) start_pos <- NA
#   if(!exists("end_pos"))   end_pos   <- NA
#
#   return(data.frame(file_id        = file_id,
#                     host_id        = host_id,
#                     sourceOrg_id   = sourceOrg_id,
#                     epitope_id     = epitope_id,
#                     molecule_id    = molecule_id,
#                     start_pos      = start_pos,
#                     end_pos        = end_pos,
#                     epit_len       = sqlen,
#                     seq            = seq,
#                     evid_code      = evid_code,
#                     epit_struc_def = epit_struc_def,
#                     qual_measure   = qual_measure,
#                     stringsAsFactors = FALSE))
# }



# get_LTCE <- function(data_folder = "./",
#                      ncpus = 1,
#                      save_folder = NULL){
#
#
#   # ========================================================================== #
#   # Sanity checks and initial definitions
#   assertthat::assert_that(is.character(data_folder),
#                           length(data_folder) == 1,
#                           dir.exists(data_folder),
#                           assertthat::is.count(ncpus),
#                           is.null(save_folder) | (is.character(save_folder)),
#                           is.null(save_folder) | length(save_folder) == 1)
#
#
#   # Set up parallel processing
#   if ((.Platform$OS.type == "windows") & (ncpus > 1)){
#     cat("\nAttention: multicore not currently available for Windows.\n
#         Forcing ncpus = 1.")
#     ncpus <- 1
#   } else {
#     available.cores <- parallel::detectCores()
#     if (ncpus >= available.cores){
#       cat("\nAttention: ncpus too large, we only have ", available.cores,
#           " cores.\nUsing ", available.cores - 1,
#           " cores for run_experiment().")
#       ncpus <- available.cores - 1
#     }
#   }
#
#
#   # =======================================
#   # Get file list and initialise variables
#   filelist    <- dir(data_folder, pattern = ".xml", full.names = TRUE)
#   filelist    <- gsub("//", "/", filelist, fixed = TRUE)
#
#   # Check save folder and create file names
#   if(!is.null(save_folder)) {
#     if(!dir.exists(save_folder)) dir.create(save_folder)
#     df_file <- paste0(normalizePath(save_folder), "/epitopes.rds")
#     #errfile <- paste0(normalizePath(save_folder), "/xml_load_errors.rds")
#   }
#
#   # ==================================================
#   cat("Processing files:\n")
#
#   df <- pbmcapply::pbmclapply(X = filelist, FUN = process_xml_file, type = "T",
#                               mc.cores = ncpus, mc.preschedule = FALSE)
#
#   #errlist <- filelist[which(sapply(df, function(x) length(x$Epitope) == 0))]
#
#   cat("\nProcessing resulting list:\n")
#   df <- pbmcapply::pbmclapply(df, data.table::rbindlist, mc.cores = 1L)
#   df <- data.frame(data.table::rbindlist(df))
#
#   if(!is.null(save_folder)){
#     saveRDS(object = df,      file = df_file)
#     #saveRDS(object = errlist, file = errfile)
#   }
#
#   invisible(df)
# }




#
# # Sequence Entropy
# feat_SeqEntr <- function(x){
#   cnts <- stringr::str_count(x, LETTERS)
#   cnts <- cnts[cnts > 0] / sum(cnts)
#   -sum(cnts * log2(cnts))
# }
#
# # Number of Atoms
# feat_NumAtom <- function(x, aa_prop){
#   cnts <- stringr::str_count(x, aa_prop$One_letter_code)
#   as.data.frame(t(colSums(cnts * aa_prop[, -1])))
# }
#
# # AA composition
# feat_AAComp <- function(x, AAtypes){
#   as.data.frame(
#     lapply(AAtypes,
#            function(pat, x){
#              sum(stringr::str_count(x, pat)) / nchar(x)},
#            x = x))
# }
#
# # AA Descriptors
# feat_AADesc <- function(x){
#   x <- strsplit(x, split = "")[[1]]
#   x <- x[-which(x == "X")]
#   f <- sapply(x, Peptides::aaDescriptors)
#   f <- data.table::as.data.table(t(rowMeans(f, na.rm = TRUE)))
# }
#
#
# # Conjoint Triads
# feat_CT <- function(x, aa_prop){
#   idx <- stringr::str_locate(paste(aa_prop$One_letter_code,
#                                    collapse = ""),
#                              strsplit(x, "")[[1]])[, 1]
#   x <- paste(aa_prop$CT_group[idx], collapse = "")
#   ct <- stringr::str_sub(x,
#                          start = 1:(nchar(x) - 2),
#                          end   = 3:nchar(x))
#
#   # Return percent occurrence of CTs
#   as.data.frame(t(as.matrix(table(ct)))) / length(ct)
# }
#
# # Molecular Weight
# feat_MolWei <- function(x, aa_prop){
#   cnts <- stringr::str_count(x, aa_prop$One_letter_code)
#   sum(cnts * aa_prop$Amino_acid_molecular_weight)
# }
#
# # N-peptides
# feat_NPep <- function(x, N){
#   i1 <- 1:(nchar(x) - N + 1)
#   i2 <- i1 + N - 1
#   y <- stringr::str_sub(x, i1, i2)
#   as.data.frame(t(as.matrix(table(y)))) / length(y)
# }



#'
#'
#' #' Calculate qualitative aminoacid composition of peptides
#' #'
#' #' This function is used to calculate the percent composition of peptides, in
#' #' terms of nine AA types: "Tiny", "Small", "Aliphatic", "Aromatic",
#' #' "NonPolar", "Polar", "Charged", "Basic" and "Acidic".
#' #'
#' #' This function is a direct adaptation of [Peptides::aaComp()]
#' #'
#' #' @param df a *data.table* of class *windowed_epit_dt* or *windowed_prot_dt*,
#' #'        generated by [make_window_df()]
#' #' @param ncpus positive integer, number of cores to use
#' #'
#' #' @return updated **df** data.table containing the AA percentage features
#' #' appended as columns.
#' #'
#' #' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#' #'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#' #'
#' calc_aa_composition <- function(df, ncpus){
#'
#'   cat("\nCalculating AA composition\n")
#'
#'   AAtypes <- list(Tiny      = c("A", "C", "G", "S", "T"),
#'                   Small     = c("A", "B", "C", "D", "G", "N", "P", "S", "T", "V"),
#'                   Aliphatic = c("A", "I", "L", "V"),
#'                   Aromatic  = c("F", "H", "W", "Y"),
#'                   NonPolar  = c("A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y"),
#'                   Polar     = c("D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z"),
#'                   Charged   = c("B", "D", "E", "H", "K", "R", "Z"),
#'                   Basic     = c("H", "K", "R"),
#'                   Acidic    = c("B", "D", "E", "Z"))
#'
#'   tmp <- mypblapply(ncpus   = ncpus,
#'                     X       = df$Info_window_seq,
#'                     FUN     = feat_AAComp,
#'                     AAtypes = AAtypes)
#'
#'   tmp        <- data.table::rbindlist(tmp, use.names = TRUE)
#'   names(tmp) <- paste0("feat_Perc_", names(tmp))
#'   newvars    <- names(tmp)
#'
#'   return(cbind(df, tmp))
#' }
#'
#'
#'
#'
#' #' Calculate a set of mean AA descriptors
#' #'
#' #' This function is used to calculate 66 aminoacid descriptors based on
#' #' function [Peptides::aaDescriptors()]. These descriptors apply to each
#' #' individual aminoacid in a peptide, and are summarised for the peptide using a
#' #' simple mean. The features returned are provided in the `Details`.
#' #'
#' #'   \itemize{
#' #'      \item\code{crucianiProperties}
#' #'      \itemize{
#' #'          \item PP1: Polarity
#' #'          \item PP2: Hydrophobicity
#' #'          \item PP3: H-bonding
#' #'      }
#' #'      \item\code{kideraFactors} (the four first features are essentially
#' #'      pure physical properties; the remaining six are linear combinations of
#' #'      several properties):
#' #'      \itemize{
#' #'          \item KF1: Helix/bend preference
#' #'          \item KF2: Side-chain size
#' #'          \item KF3: Extended structure preference
#' #'          \item KF4: Hydrophobicity
#' #'          \item KF5 - KF10: linear combinations of other characteristics
#' #'      }
#' #'      \item\code{zScales} (based on physicochemical properties of the AAs
#' #'      including NMR data and thin-layer chromatography data):
#' #'      \itemize{
#' #'          \item Z1: Lipophilicity
#' #'          \item Z2: Steric properties (Steric bulk/Polarizability)
#' #'          \item Z3: Electronic properties (Polarity / Charge)
#' #'          \item Z4-Z5: relate electronegativity, heat of formation,
#' #'          electrophilicity and hardness.
#' #'      }
#' #'      \item\code{FASGAI} (based on physicochemical properties of the AAs
#' #'      including NMR data and thin-layer chromatography data):
#' #'      \itemize{
#' #'          \item F1: Hydrophobicity index
#' #'          \item F2: Alpha and turn propensities
#' #'          \item F3: Bulky properties
#' #'          \item F4: Compositional characteristic index
#' #'          \item F5: Local flexibility
#' #'          \item F6: Electronic properties
#' #'      }
#' #'      \item\code{tScales} (based on 67 common topological descriptors of
#' #'      amino acids. These topological descriptors are based on the connectivity
#' #'      table of amino acids alone, and to not explicitly consider 3D properties
#' #'      of each structure):
#' #'      \itemize{
#' #'          \item T1 - T5
#' #'      }
#' #'      \item\code{VHSE scales} (principal component score Vectors of
#' #'      Hydrophobic, Steric, and Electronic properties. Derived from PCA on
#' #'      independent families of 18 hydrophobic properties, 17 steric properties,
#' #'      and 15 electronic properties):
#' #'      \itemize{
#' #'          \item VHSE1 - VHSE2: Hydrophobic properties
#' #'          \item VHSE3 - VHSE4: Steric properties
#' #'          \item VHSE5 - VHSE8: Electronic properties
#' #'      }
#' #'      \item\code{ProtFP descriptors}: (constructed from a large initial
#' #'      selection of indices obtained from the AAindex database for all 20
#' #'      naturally occurring amino acids.):
#' #'      \itemize{
#' #'          \item ProtFP1-ProtFP8
#' #'      }
#' #'      \item\code{stScales} (proposed by Yang et al.'2010, taking 827
#' #'      properties into account which are mainly constitutional, topological,
#' #'      geometrical, hydrophobic, electronic, and steric properties of AAs):
#' #'      \itemize{
#' #'          \item ST1 - ST8
#' #'      }
#' #'      \item\code{BLOSUM indices} (derived of physicochemical properties that
#' #'      have been subjected to a VARIMAX analyses and an alignment matrix of the
#' #'      20 natural AAs using the BLOSUM62 matrix):
#' #'      \itemize{
#' #'          \item BLOSUM1 - BLOSUM10
#' #'      }
#' #'      \item\code{MS-WHIM scores} (derived from 36 electrostatic potential
#' #'      properties derived from the three-dimensional structure of the 20
#' #'      natural amino acids):
#' #'      \itemize{
#' #'          \item MSWHIM1 - MSWHIM3
#' #'      }
#' #' }
#' #'
#' #' See the documentation for [Peptides::aaDescriptors()] for further information
#' #' and references.
#' #'
#' #' @param df a *data.table* of class *windowed_epit_dt* or *windowed_prot_dt*,
#' #'        generated by [make_window_df()]
#' #' @param ncpus positive integer, number of cores to use
#' #'
#' #' @return updated **df** data.table containing the AA descriptors features
#' #' appended as columns.
#' #'
#' #' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#' #'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#' #'
#' #'
#' calc_aa_descriptors <- function(df, ncpus){
#'
#'   cat("\nCalculating AA descriptors\n")
#'
#'   tmp <- mypblapply(ncpus = ncpus,
#'                     X     = df$Info_window_seq,
#'                     FUN   = feat_AADesc)
#'
#'   tmp <- data.table::rbindlist(tmp)
#'
#'   # Prepare feature names
#'   col_names  <- paste0("feat_", colnames(Peptides::aaDescriptors("K")))
#'   names(tmp) <- gsub(pattern = "\\.1$", replacement = "", x = col_names)
#'
#'   return(cbind(df, tmp))
#' }
#'
#'
#'
#'
#' #' Calculate conjoint triads
#' #'
#' #' This function is used to calculate the conjoint triad frequencies for
#' #' each peptide
#' #'
#' #' @param df a *data.table* of class *windowed_epit_dt* or *windowed_prot_dt*,
#' #'        generated by [make_window_df()]
#' #' @param ncpus positive integer, number of cores to use
#' #'
#' #' @return updated **df** data.table containing the conjoint triad frequencies
#' #'
#' #' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#' #'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#' #'
#' calc_conjoint_triads <- function(df, ncpus){
#'
#'   cat("\nCalculating conjoint triads\n")
#'
#'   # Load AA properties
#'   aa_prop <- readRDS(system.file("extdata", "amino_acid_propensity.rds",
#'                                  package = "epitopes"))
#'   aa_prop <- aa_prop[, c("One_letter_code", "CT_group")]
#'
#'   tmp <- mypblapply(ncpus   = ncpus,
#'                     X       = df$Info_window_seq,
#'                     FUN     = feat_CT,
#'                     aa_prop = aa_prop)
#'
#'   # Add a dummy element with all possible triads
#'   dummy <- as.data.frame(matrix(0, ncol = 7 ^ 3, nrow = 1))
#'   names(dummy) <- apply(X   = expand.grid(0:6, 0:6, 0:6),
#'                         FUN = paste, MARGIN = 1, collapse = "")
#'   tmp[[length(tmp) + 1]] <- dummy
#'
#'   # Binds results
#'   tmp <- data.table::rbindlist(tmp, use.names = TRUE, fill = TRUE)
#'   tmp[is.na(tmp)] <- 0
#'
#'   # Removes dummy observation
#'   tmp <- tmp[-nrow(tmp), ]
#'
#'   names(tmp) <- paste0("feat_CT", names(tmp))
#'
#'   return(cbind(df, tmp))
#' }
#'
#'
#'
#'
#' #' Calculate molecular weight of peptide sequences
#' #'
#' #' This function is used to calculate the molecular weights for each peptide
#' #'
#' #' @param df a *data.table* of class *windowed_epit_dt* or *windowed_prot_dt*,
#' #'        generated by [make_window_df()]
#' #' @param ncpus positive integer, number of cores to use
#' #'
#' #' @return updated **df** data.table containing the molecular weight features
#' #' appended as columns.
#' #'
#' #' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#' #'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#' #'
#'
#' calc_molecular_weight <- function(df, ncpus){
#'
#'   cat("\nCalculating molecular weight\n")
#'
#'   # Load AA properties
#'   aa_prop <- readRDS(system.file("extdata", "amino_acid_propensity.rds",
#'                                  package = "epitopes"))
#'   aa_prop <- aa_prop[, c("One_letter_code", "Amino_acid_molecular_weight")]
#'
#'   tmp <- mypblapply(ncpus   = ncpus,
#'                     X       = df$Info_window_seq,
#'                     FUN     = feat_MolWei,
#'                     aa_prop = aa_prop)
#'
#'   return(cbind(df, feat_molecular_weight = unlist(tmp)))
#' }
#'
#'
#'
#'
#' #' Calculate frequencies of N-peptides
#' #'
#' #' This function is used to calculate the frequencies of N-peptide subsequences
#' #' for each peptide
#' #'
#' #' @param df a *data.table* of class *windowed_epit_dt* or *windowed_prot_dt*,
#' #'        generated by [make_window_df()]
#' #' @param N the length of the subsequences to consider.
#' #' @param ncpus positive integer, number of cores to use
#' #'
#' #' @return updated **df** data.table containing the N-peptide frequencies
#' #'
#' #' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#' #'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#' #'
#' #' @export
#' #'
#' calc_Npeptide_composition <- function(df, N, ncpus){
#'
#'   # Creates a list of tables contaning the dipeptide percentages for each window sequence in df
#'   cat(paste0("\nCalculating ", N, "-mer composition\n"))
#'
#'   tmp <- mypblapply(ncpus = ncpus,
#'                     X     = df$Info_window_seq,
#'                     FUN   = feat_NPep,
#'                     N     = N)
#'
#'   # Add a dummy element with all possible n-mers
#'   dummy <- as.data.frame(matrix(0, ncol = 20 ^ N, nrow = 1))
#'   args <- vector(mode = "list", length = N)
#'   for (i in seq_along(args)) args[[i]] <- get_aa_codes()
#'   names(dummy) <- apply(X   = do.call(expand.grid, args),
#'                         FUN = paste, MARGIN = 1, collapse = "")
#'   tmp[[length(tmp) + 1]] <- dummy
#'
#'   # Binds results
#'   tmp <- data.table::rbindlist(tmp, use.names = TRUE, fill = TRUE)
#'   tmp[is.na(tmp)] <- 0
#'
#'   # Removes dummy observation
#'   tmp <- tmp[-nrow(tmp), ]
#'
#'   names(tmp) <- paste0("feat_Perc_", names(tmp))
#'
#'   return(cbind(df, tmp))
#' }
#'
#'
#'
#'
#'
#' #' Calculate number of each type of atom in peptide sequences
#' #'
#' #' This function is used to calculate the number of specific atoms in each
#' #' peptide
#' #'
#' #' @param df a *data.table* of class *windowed_epit_dt* or *windowed_prot_dt*,
#' #'        generated by [make_window_df()]
#' #' @param ncpus positive integer, number of cores to use
#' #'
#' #' @return updated **df** data.table containing the number of atoms of elements
#' #'         C, H, O, N and S.
#' #'
#' #' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#' #'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#' #'
#' calc_number_of_atoms <- function(df, ncpus){
#'
#'   cat("\nCalculating atomic composition\n")
#'
#'   # Load AA properties
#'   aa_prop <- readRDS(system.file("extdata", "amino_acid_propensity.rds",
#'                                  package = "epitopes"))
#'   aa_prop <- aa_prop[, c(2, 4:8)]
#'   names(aa_prop)[-1] <- c("C", "H", "N", "O", "S")
#'
#'   tmp <- mypblapply(ncpus   = ncpus,
#'                     X       = df$Info_window_seq,
#'                     FUN     = feat_NumAtom,
#'                     aa_prop = aa_prop)
#'
#'   tmp <- data.table::rbindlist(tmp, use.names = TRUE)
#'   names(tmp) <- paste0("feat_", names(tmp), "_atoms")
#'
#'   return(cbind(df, tmp))
#'
#' }
#'
#'
#'
#'
#' #' Calculate sequence entropy of peptides
#' #'
#' #' This function is used to calculate the sequence entropies,
#' #' Entropy = -sum(p_k * log2(p_k))
#' #'
#' #' @param df a *data.table* of class *windowed_epit_dt* or *windowed_prot_dt*,
#' #'        generated by [make_window_df()]
#' #' @param ncpus positive integer, number of cores to use
#' #'
#' #' @return updated **df** data.table containing the sequence entropy of each
#' #'         peptide
#' #'
#' #' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#' #'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#' #'
#'
#' calc_sequence_entropy <- function(df, ncpus){
#'
#'   # Entropy = -SUM(p(e)log2[p(e)])
#'   cat("\nCalculating sequence entropy\n")
#'
#'   tmp <- mypblapply(ncpus = ncpus,
#'                     X     = df$Info_window_seq,
#'                     FUN   = feat_SeqEntr)
#'
#'   return(cbind(df, feat_seq_entropy = unlist(tmp)))
#'
#' }
#'
