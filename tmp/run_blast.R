# Not currently used in the package.


# run_blast <- function(BLAST_path, proteins, ncpus){
#
#   # Check if BLAST is installed
#   errk <- FALSE
#   tryCatch({
#     invisible(utils::capture.output(
#       blast_version <- system("blastp -version", intern = TRUE)[1]))},
#     error   = function(c) {errk <<- TRUE},
#     warning = function(c) {errk <<- TRUE},
#     finally = NULL)
#
#   if (errk){
#     stop(paste0("\nBLAST+ not found.",
#                 "\nPlease follow the instructions in",
#                 "\nhttps://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download",
#                 "\nto set up BLAST+ on your machine.",
#                 "\n**************************************"))
#   }
#
#   # Run BLASTp to determine protein similarities
#   if(!dir.exists(BLAST_path)) dir.create(BLAST_path, recursive = TRUE)
#   fn <- gsub("//", "/", paste0(BLAST_path, "/proteins.fasta"), fixed = TRUE)
#
#   x        <- proteins$TSeq_sequence
#   names(x) <- proteins$UID
#   x        <- Biostrings::AAStringSet(x)
#   Biostrings::writeXStringSet(x, filepath = fn, format="fasta")
#
#   message("===========================================\nBuilding BLASTp database")
#   system(paste0("makeblastdb -in ", fn, " -dbtype prot"))
#
#   message("===========================================\nRunning BLASTp\n")
#   system(paste0("blastp -query ", fn, " -db ", fn, " -seg no ",
#                 "-outfmt '6 qseqid sseqid pident qcovhsp' > ",
#                 fn, "-BLAST"))
#
#   blast <- utils::read.csv(paste0(fn, "-BLAST"), sep = "\t",
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
#
#   names(blast) <- c("QueryID", "SubjectID", "Perc_identity", "Query_coverage")
#
#   # blast$QueryID   <- gsub("\\.[1-9]+$", "", blast$QueryID)
#   # blast$SubjectID <- gsub("\\.[1-9]+$", "", blast$SubjectID)
#
#   return(blast)
#
# }
