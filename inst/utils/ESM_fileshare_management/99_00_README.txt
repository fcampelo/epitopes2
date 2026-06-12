NOT WORKING - IGNORE


Manage github repositories

This series of functions/ routines is used to update the Github repositories containing the precalculated ESM-1b features.

Note: this assumes that the repositories already exist in a folder and that the machine is registered with a Github token for the
"epitopes-dataset" Github organization. If not, please clone the repos in the local machine and set up githu before you start.

As a note, the code below was used in the first setting up of the repositories. I'm leaving it here in case it may come in handy in the future:

=====================================================

# set up token for the epitopes-dataset organisation
# and register it:
# Sys.setenv(GITHUB_PAT = "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")

# install.packages(c("git2r"))
library(git2r)
library(dplyr)

max_repo_size   <- 5e9  # 5Gb

my_username <- "fcampelo"
my_usermail <- "fcampelo@gmail.com"

# Get total size of files
sizes <- sapply(dir("../output/proteins_esm1b_R", full.names = TRUE),
                function(x) file.info(x)$size)

repo_split <- data.frame(filename = names(sizes),
                         protID   = gsub("\\.\\./output/proteins\\_esm1b\\_R/|\\.rds", "", names(sizes)),
                         filesize = unname(sizes),
                         which.repo = sprintf("ESM1b_IEDB_LBCE_proteins_%d",
                                              ceiling(cumsum(unname(sizes)) / max_repo_size)))

repos <- unique(repo_split$which.repo)
upload_errors <- data.frame(repo = character(),
                            file = character())

for (i in seq_along(repos)){
  cat("\n\n REPO", i)
  mypath <- paste0("../output/", repos[i])

  if(dir.exists(mypath)) unlink(mypath, recursive = TRUE)

  # Clone repo to local
  myrepo <- git2r::clone(
    url = sprintf("https://github.com/epitopes-dataset/ESM1b_IEDB_LBCE_proteins_%d.git", i),
    local_path = mypath,
    credentials = git2r::cred_token())

  git2r::config(myrepo, user.name = my_username, user.email = my_usermail)

  ## ===========================================================
  ## TODO: add loop to remove / commit/ push all existing files
  ## ===========================================================

  # Add README
  writeLines(paste0("## LBCE Protein dataset - part ", i),
             file.path(mypath, "README.md"))

  git2r::add(myrepo, file.path(mypath, "README.md"))

  tryCatch({
    git2r::commit(myrepo, "Added README")
    git2r::push(myrepo, credentials = cred_token())},
    warning = function(c) {invisible(FALSE)},
    error   = function(c) {invisible(FALSE)},
    finally = NULL)


  write.csv(repo_split[, -1], paste0(mypath, "/file_map.csv"),
            row.names = FALSE, quote = FALSE)
  git2r::add(myrepo, "file_map.csv")
  tryCatch({
    git2r::commit(myrepo, "Added file list of protein repos")
    git2r::push(myrepo, credentials = cred_token())},
    warning = function(c) {invisible(FALSE)},
    error   = function(c) {invisible(FALSE)},
    finally = NULL)


  repo_files <- repo_split %>%
    filter(which.repo == repos[i])

  if(!dir.exists(paste0(mypath, "/data"))) dir.create(paste0(mypath, "/data"))
  for (k in 1:nrow(repo_files)){
    cat("\nRepo", i, "\tFile", k, "/", nrow(repo_files))
    ignore <- file.copy(from = repo_files$filename[k],
                        to   = paste0(mypath, "/data"),
                        overwrite = TRUE)

    git2r::add(myrepo, gsub("../output/proteins_esm1b_R/", "data/",
                            repo_files$filename[k], fixed = TRUE))

    errk <- FALSE
    tryCatch({
      git2r::commit(myrepo, paste0("Added data file ",
                                   gsub("../output/proteins_esm1b_R/", "",
                                        repo_files$filename[k], fixed = TRUE)))
      git2r::push(myrepo, credentials = cred_token())},
      warning = function(c) {errk <<- TRUE},
      error   = function(c) {errk <<- TRUE},
      finally = NULL)

    if (errk){
      upload_errors <- rbind(upload_errors,
                             data.frame(repo = repos[i],
                                        file = repo_files$filename[k]))
      write.csv(upload_errors, "../output/upload_errors.csv", row.names = FALSE)
      cat("\tFAIL")
    } else {
      cat("\tSUCCESS")
    }
  }
}

# for (i in 1:length(repos)){
#   cat("\n\n REPO", i)
#   mypath <- paste0("../output/", repos[i])
#   myrepo <- repository(mypath)
#   git2r::rm_file(myrepo, paste0(mypath, "/file_map.csv"))
#   git2r::commit(myrepo, "Removed list of protein repos")
#   write.csv(repo_split[, -1], paste0(mypath, "/file_map.csv"),
#             row.names = FALSE, quote = FALSE)
#   git2r::add(myrepo, "file_map.csv")
#   git2r::commit(myrepo, "Updated file list of protein repos")
#   git2r::push(myrepo, credentials = cred_token())
# }

