require(reticulate)
require(pbapply)
require(dplyr)

outdir <- paste0(SAVE_FOLDER, "/proteins_rds")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

if(!reticulate::py_module_available("pandas"))
  reticulate::py_install("pandas")

pd <- reticulate::import("pandas")

fl    <- dir(SAVE_FOLDER, pattern = ".pkl", full.names = TRUE)
fn    <- paste0("/",
                unname(sapply(dir(SAVE_FOLDER,
                                  pattern = ".pkl"),
                              function(x)
                                strsplit(x,
                                         split = "\\/|\\_\\_|\\.pkl")[[1]][1])))
fn_un <- unique(fn)

for (i in seq_along(fn_un)){
  cat(sprintf("\rConcatenating %04d of %04d", i, length(fn_un)))
  idx <- grep(paste0(gsub(".", "\\.", fn_un[i], fixed = TRUE), "$"), fn)
  X   <- lapply(fl[idx], X <- reticulate::py_load_object)
  if (length(X) == 1){
    X <- X[[1]]
  } else {
    for(k in seq_along(idx)){
      cat(".")
      rng <- strsplit(fl[idx][k], split = "\\_\\_")[[1]][2]
      rng <- as.numeric(strsplit(rng, split = "\\_|\\.")[[1]][1:2])
      X[[k]]$pos <- seq(rng[1], rng[2])
    }
    X <- dplyr::bind_rows(X) %>%
      dplyr::group_by(.data$pos) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), mean)) %>%
      dplyr::select(-pos)
  }
  names(X) <- paste0(feat_prefix, names(X))
  saveRDS(X, paste0(outdir, fn[idx[1]], ".rds"))
}

