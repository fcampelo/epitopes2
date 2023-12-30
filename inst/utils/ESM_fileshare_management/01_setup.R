# library(dplyr)
# library(git2r)
#
# # IMPORTANT: should existing files be overwritten?
# overwrite <- TRUE
#
# # Github info:
# my_username <- "fcampelo"
# my_usermail <- "fcampelo@gmail.com"
#
# # Repo and file parameters
# max_repo_size   <- 5e9  # 5GB
# max_file_size   <- 75e6 # 75MB
#
# # Folder containing the subfolders for each github repo
# containing_folder <- "/users/c/campelof/esm1b/output"
#
# # Local folder containing ESM-1b output:
# esm_folder <- "/users/c/campelof/epitopes-secret-sauce/inst/utils/ESM_feature_calculation/output/esm1b_features/out/proteins_rds"
#
# # Proteins file
# proteins <- readRDS("/users/c/campelof/epitopes-secret-sauce/data/proteins.rds")
#
# # Download file map
# file_map <- read.csv("https://raw.githubusercontent.com/epitopes-dataset/ESM1b_IEDB_LBCE_proteins_1/main/file_map.csv")
#
# # Extract repo folders from enclosing folder:
# repo_folders <- dir(containing_folder, pattern = "ESM1b_IEDB_LBCE_proteins_")
# repo_folders_path <- dir(containing_folder, pattern = "ESM1b_IEDB_LBCE_proteins_", full.names = TRUE)
#
#
#
# ### Define functions
# copy_newfile_to_repo <- function(fp, fn, local_path, max_file_size, prot, repo){
#
#   repo_files <- dir(paste0(local_path, "/data"))
#
#   # Read new file
#   X <- readRDS(fp)
#
#   # check that the sizes are consistent and append protein info if
#   # all is OK
#   if (nrow(X) != nchar(prot$TSeq_sequence)) {
#     return("Error: sizes not consistent")
#   } else {
#     X <- X %>%
#       mutate(Info_pos = 1:nrow(X),
#              Info_AA  = strsplit(prot$TSeq_sequence, split = "")) %>%
#       select(starts_with("Info"), everything())
#     names(X) <- gsub("[\\ ]+$", "", names(X))
#   }
#
#   toRM <- which(gsub("\\_part[0-9]+\\.rds", ".rds", repo_files) == fn)
#   git2r::rm_file(repo, paste0(local_path, "/data/", repo_files[toRM]))
#   git2r::commit(repo, paste0("Removed protein ", prot$UID))
#   git2r::push(repo, credentials = cred_token())
#
#   toADD <- character()
#
#   if (file.size(fp) < 0.95 * max_file_size){
#     saveRDS(X, paste0(local_path, "/data/", fn))
#     toADD <- c(toADD, paste0(local_path, "/data/", fn))
#
#     git2r::add(repo, toADD)
#
#     errk <- FALSE
#     tryCatch({
#       git2r::commit(repo, paste0("Updated protein ", prot$UID))
#       git2r::push(repo, credentials = cred_token())},
#       warning = function(c) {errk <<- TRUE},
#       error   = function(c) {errk <<- TRUE},
#       finally = NULL)
#
#     if (errk){
#
#       if(file.exists("upload_errors.csv")){
#         upload_errors <- read.csv("upload_errors.csv")
#       } else {
#         upload_errors <- data.frame(repo = character(), file = character())
#       }
#
#       upload_errors <- rbind(upload_errors,
#                              data.frame(repo = local_path,
#                                         file = fn))
#       write.csv(upload_errors, "upload_errors.csv", row.names = FALSE)
#       cat("\tFAAAAAAAILURE")
#     } else {
#       if(file.exists("processed_prots.csv")){
#         processed_prots <- read.csv("processed_prots.csv")
#       } else {
#         processed_prots <- data.frame(repo = character(), file = character())
#       }
#       processed_prots <- rbind(processed_prots,
#                                data.frame(repo = local_path,
#                                           file = fn))
#       write.csv(processed_prots, "processed_prots.csv", row.names = FALSE)
#       cat("\tFUYOOOOOH!")
#     }
#   }
# }
#
# # } else {
# #   nchunks   <- ceiling(file.size(fp) / max_file_size)
# #   nr        <- floor(nrow(X) / nchunks)
# #   for (k in 1:nchunks){
# #     st <- (k - 1) * nr + 1
# #     en <- min(k * nr, nrow(X))
# #     partfn <- paste0(local_path, "/data/",
# #                      gsub("\\.rds", "", fn),
# #                      sprintf("_part%02d.rds", k))
# #     saveRDS(X[st:en, ], partfn)
# #     toADD <- c(toADD, partfn)
# #   }
# # }
