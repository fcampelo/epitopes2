library(dplyr)
library(pbapply)

ncpus <- parallel::detectCores() - 1

feat_path <- "d:/IEDB_proteins_ESM1b_features/features/"

dirs <- dir(feat_path, full.names = TRUE)

if(Sys.info()["sysname"] == "Windows"){
  cl <- parallel::makeCluster(ncpus, setup_timeout = 2)
  parallel::clusterExport(cl, c("dirs"))
} else {
  cl <- ncpus
}

errs <- pblapply(seq_along(dirs),
                 function(i){
                   fl <- dir(dirs[i], full.names = TRUE)
                   fn <- dir(dirs[i], full.names = FALSE)
                   errlist <- data.frame(folder = character(),
                                         file   = character())
                   for (i in seq_along(fl)){
                     X <- readRDS(fl[i])
                     if(!all(c("Info_pos", "Info_AA") %in% names(X))){
                       errlist <- rbind(errlist,
                                        data.frame(folder = i,
                                                   file   = fn[i]))
                     }
                   }
                 }, cl = cl)
