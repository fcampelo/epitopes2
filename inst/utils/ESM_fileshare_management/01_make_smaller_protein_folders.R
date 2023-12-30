filedir <- dir("../ESM_feature_calculation/output/esm1b_features/out/proteins_rds/", full.names=TRUE)

nfiles <- length(filedir)
fpf <- 250

nfolders <- ceiling(nfiles/fpf)

for (i in 1:nfolders){
  folder <- sprintf("folder%02d/", i)
  if (!dir.exists(folder)) dir.create(folder)
  cat("\n", folder, ": ")

  st <- 1 + (i - 1)*fpf
  en <- min(i*fpf, nfiles)
  for (k in st:en){
    file.copy(filedir[k], folder, overwrite = TRUE)
    if(!(k%%10)) cat(".")
  }
}
