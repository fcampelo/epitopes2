library(dplyr)
library(pbapply)

ncpus <- parallel::detectCores() - 1

feat_path <- "d:/IEDB_proteins_ESM1b_features/features/"

dirs <- dir(feat_path, full.names = TRUE)

errlist  <- data.frame(folder = character(),
                       file   = character())

filelist <- data.frame(protein = character(),
                       folder  = character())
for (j in seq_along(dirs)){
  cat("\n")
  t0 <- Sys.time()
  fl <- dir(dirs[j], full.names = TRUE)
  fn <- dir(dirs[j], full.names = FALSE)
  filelist <- rbind(filelist,
                    data.frame(protein = gsub("\\.rds", "", fn),
                               folder  = sprintf("Folder%02d", j)))
  write.csv(filelist, "filelist.csv", row.names = FALSE, quote = FALSE)
  for (i in seq_along(fl)){
    cat(sprintf("\rFolder %02d/%02d file %03d/%03d",
                j, length(dirs), i, length(fl)))
    X <- readRDS(fl[i])
    if(!all(c("Info_pos", "Info_AA") %in% names(X))){
      cat("\t *")
      errlist <- rbind(errlist,
                       data.frame(folder = sprintf("Folder%02d", j),
                                  file   = fn[i]))
      write.csv(errlist, "errlist.csv", row.names = FALSE, quote = FALSE)
    }
  }
  a <- difftime(Sys.time(), t0)
  cat(":", as.numeric(a), attr(a, "units"))
}

head(filelist)

