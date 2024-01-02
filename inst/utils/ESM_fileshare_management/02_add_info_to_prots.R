library(dplyr)

proteins <- readRDS("../../../data/proteins.rds")

feat_path <- "d:/IEDB_proteins_ESM1b_features/features/"

dirs <- dir(feat_path, full.names = TRUE)
nfiles <- sapply(dirs, function(d) length(dir(d)))

# Check number of files in each folder
# (should be 250 except for the last one)
which(nfiles != 250)

# Test file integrity
for (i in seq_along(dirs)){
  cat("\n")
  filelist <- dir(dirs[i], full.names = TRUE)
  fn       <- dir(dirs[i], full.names = FALSE)
  for (j in seq_along(filelist)){
    cat(sprintf("\nFolder %02d: protein %03d: %s", i, j, fn[j]))
    X <- readRDS(filelist[j])
    prot <- proteins %>% filter(UID == gsub("\\.rds", "", fn[j]))
    if (nrow(prot) != 1 || nrow(X) != nchar(prot$TSeq_sequence)){
      myerror <- data.frame(UID = gsub("\\.rds", "", fn[j]),
                            prot.entries = nrow(prot),
                            nrow.X = nrow(X),
                            nchar.prot = ifelse(
                              nrow(prot) == 0, NA,
                              nchar(prot$TSeq_sequence[1])))
      write.table(myerror, file = "add_info_errors.csv",
                  append = file.exists("add_info_errors.csv"),
                  quote = FALSE, sep = ",",
                  row.names = FALSE,
                  col.names = !file.exists("add_info_errors.csv"))
      cat(" - ERROR! <<<<-----")
    } else {
      X <- X %>%
        mutate(Info_pos = 1:nrow(X),
               Info_AA  = strsplit(prot$TSeq_sequence, split = "")) %>%
        select(starts_with("Info"), everything())
      names(X) <- gsub("[\\ ]+$", "", names(X))

      saveRDS(X, file = filelist[j])
      cat(" - OK!")
    }
  }
}

myerrors <- read.csv("add_info_errors.csv")
myerrors

#
# for (i in seq_along(dirs)){
#   cat("\n")
#   filelist <- dir(dirs[i], full.names = TRUE)
#   fn       <- dir(dirs[i], full.names = FALSE)
#   idx <- which(gsub("\\.rds", "", fn) %in% myerrors$UID)
#   filelist <- filelist[idx]
#   fn <- fn[idx]
#   for (j in seq_along(filelist)){
#     cat(sprintf("\nFolder %02d: protein %03d: %s", i, j, fn[j]))
#     X <- readRDS(filelist[j])
#     prot <- proteins %>% filter(UID == gsub("\\.rds", "", fn[j]))
#     if (nrow(prot) != 1 || nrow(X) != nchar(prot$TSeq_sequence)){
#       myerror <- data.frame(UID = gsub("\\.rds", "", fn[j]),
#                             prot.entries = nrow(prot),
#                             nrow.X = nrow(X),
#                             nchar.prot = ifelse(
#                               nrow(prot) == 0, NA,
#                               nchar(prot$TSeq_sequence[1])))
#       write.table(myerror, file = "add_info_errors2.csv",
#                   append = file.exists("add_info_errors2.csv"),
#                   quote = FALSE, sep = ",",
#                   row.names = FALSE,
#                   col.names = !file.exists("add_info_errors2.csv"))
#       cat(" - ERROR! <<<<-----")
#     } else {
#       X <- X %>%
#         mutate(Info_pos = 1:nrow(X),
#                Info_AA  = strsplit(prot$TSeq_sequence, split = "")) %>%
#         select(starts_with("Info"), everything())
#       names(X) <- gsub("[\\ ]+$", "", names(X))
#
#       saveRDS(X, file = filelist[j])
#       cat(" - OK!")
#     }
#   }
# }
