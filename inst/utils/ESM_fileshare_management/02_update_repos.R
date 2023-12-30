### NOT WORKING - IGNORE

# newfiles <- dir(esm_folder, pattern = ".rds")
# newfiles_path <- dir(esm_folder, pattern = ".rds", full.names = TRUE)
#
# for (i in seq_along(repo_folders_path)){
#   local_path <- repo_folders_path[i]
#   repo <- repository(local_path)
#   tmp <- git2r::fetch(repo, "origin")
#
#   files_in_repo <- dir(paste0(local_path, "/data"))
#   prots_in_repo <- files_in_repo %>%
#     gsub(pattern = "\\.rds|\\_part[0-9]+\\.rds$", replacement = "") %>%
#     unique()
#
#   for (j in seq_along(prots_in_repo)){
#     cat(sprintf("\nRepo: %s, Prot: %s :: ", repo_folders[i], prots_in_repo[j]))
#     if ((paste0(prots_in_repo[j], ".rds") %in% newfiles) && overwrite){
#       cat("trying... ")
#       # Path to new file:
#       idx <- which(newfiles == paste0(prots_in_repo[j], ".rds"))
#       fp  <- newfiles_path[idx]
#       fn <- newfiles[idx]
#       prot <- proteins %>% dplyr::filter(UID == prots_in_repo[j])
#       copy_newfile_to_repo(fp, fn, local_path, max_file_size, prot, repo)
#     } else {
#       cat("nothing done.")
#     }
#   }
# }

## STOPPED HERE
##
## fs <- file.size(newfiles_path)
## bigfiles <- newfiles[which(fs >= 0.95 * max_file_size)]
## bigfiles_path <- newfiles_path[which(fs >= 0.95 * max_file_size)]





#
#
#
#       ignore <- file.copy(from = repo_files$filename[k],
#                           to   = paste0(mypath, "/data"),
#                           overwrite = TRUE)
#
#       git2r::add(myrepo, gsub("../output/proteins_esm1b_R/", "data/",
#                               repo_files$filename[k], fixed = TRUE))
#
#       errk <- FALSE
#       tryCatch({
#         git2r::commit(myrepo, paste0("Added data file ",
#                                      gsub("../output/proteins_esm1b_R/", "",
#                                           repo_files$filename[k], fixed = TRUE)))
#         git2r::push(myrepo, credentials = cred_token())},
#         warning = function(c) {errk <<- TRUE},
#         error   = function(c) {errk <<- TRUE},
#         finally = NULL)
#
#       if (errk){
#         upload_errors <- rbind(upload_errors,
#                                data.frame(repo = repos[i],
#                                           file = repo_files$filename[k]))
#         write.csv(upload_errors, "../output/upload_errors.csv", row.names = FALSE)
#         cat("\tFAIL")
#       } else {
#         cat("\tSUCCESS")
#       }
#     }
#
#   }
# }
