#' @importFrom dplyr %>%
#' @importFrom rlang .data

concatenate_ESM_output <- function(save_folder, feat_prefix, ncpus){

  if(!reticulate::py_module_available("pandas")) reticulate::py_install("pandas")

  pd <- reticulate::import("pandas")

  outdir <- paste0(save_folder, "/proteins_esm1b_rds")
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  fl    <- dir(save_folder, pattern = ".pkl", full.names = TRUE)
  fn    <- paste0("/",
                  unname(sapply(dir(save_folder,
                                    pattern = ".pkl"),
                                function(x)
                                  strsplit(x,
                                           split = "\\/|\\_\\_|\\.pkl")[[1]][1])))
  fn_un <- unique(fn)

  .ignore <- mypblapply(
    X   = seq_along(fn_un),
    FUN = function(i, fn, fl, fn_un){
      idx <- grep(paste0(gsub(".", "\\.", fn_un[i], fixed = TRUE), "$"), fn)
      X <- lapply(fl[idx], X <- reticulate::py_load_object)
      if (length(X) == 1){
        X <- X[[1]]
      } else {
        for(k in seq_along(idx)){
          rng <- strsplit(fl[idx][k], split = "\\_\\_")[[1]][2]
          rng <- as.numeric(strsplit(rng, split = "\\_|\\.")[[1]][1:2])
          X[[k]]$pos <- seq(rng[1], rng[2])
        }
        X <- dplyr::bind_rows(X) %>%
          dplyr::group_by(.data$pos) %>%
          dplyr::summarise(dplyr::across(dplyr::everything(), mean))
        X$pos <- NULL
      }
      names(X) <- paste0(feat_prefix, names(X))
      saveRDS(X, paste0(outdir, fn[idx[1]], ".rds"))
      invisible(TRUE)
    },
    fn = fn, fl = fl, fn_un = fn_un,
    ncpus = ncpus)

  invisible(TRUE)
}
