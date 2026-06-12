#' Consolidate the outputs from calc_features_esm2
#'
#' @param cl a SOCK cluster object created using [epitopes::set_mc_cluster()], or
#'        `NULL` if parallel processing is not desired.
#'
#'
#' @author Felipe Campelo (\email{fcampelo@@gmail.com})
#'
#' @export
#'
#' @importFrom dplyr %>%
#'

concatenate_esm_outputs <- function(csv_folder,
                                    save_folder = NULL,
                                    filenames = NULL,
                                    save_file = NULL,
                                    feat_prefix = "feat_esm2_",
                                    cl = NULL,
                                    delete_originals = FALSE){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(csv_folder), length(csv_folder) == 1,
                          dir.exists(csv_folder),
                          is.null(save_folder) | (is.character(save_folder) & length(save_folder) == 1),
                          is.null(save_file) | (is.character(save_file) & length(save_file) == 1),
                          is.null(cl) | "SOCKcluster" %in% class(cl),
                          is.logical(delete_originals), length(delete_originals) == 1)



  if(!is.null(save_folder) && !dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)

  if(is.null(filenames)) {
    fl    <- dir(csv_folder, pattern = ".csv", full.names = TRUE)
    fn    <- paste0("/",
                    unname(sapply(dir(csv_folder, pattern = ".csv"),
                                  function(x)
                                    strsplit(x,
                                             split = "\\/|\\_\\_|\\.csv")[[1]][1])))
  } else {
    fl <- paste0(csv_folder, "/", filenames)
    fn    <- paste0("/",
                    unname(sapply(filenames,
                                  function(x)
                                    strsplit(x,
                                             split = "\\/|\\_\\_|\\.csv")[[1]][1])))
  }

  fn_un <- unique(fn)

  reslist <- pbapply::pblapply(seq_along(fn_un),
                               function(i, fn_un, fn, fl, feat_prefix){
                                 idx <- paste0(gsub(".", "\\.", fn_un[i], fixed = TRUE), "$")
                                 idx <- gsub("+", "\\+", idx, fixed = TRUE)
                                 idx <- grep(idx, fn)
                                 xx   <- lapply(fl[idx], utils::read.csv)
                                 if (length(xx) == 1){
                                   xx <- xx[[1]]
                                   if(ncol(xx) == 1) {
                                     xx <- as.data.frame(t(xx))
                                   }
                                 } else {
                                   for(k in seq_along(idx)){
                                     rng <- strsplit(fl[idx][k], split = "\\_\\_")[[1]][2]
                                     rng <- as.numeric(strsplit(rng, split = "\\_|\\.")[[1]][1:2])
                                     xx[[k]]$Info_pos <- seq(rng[1], rng[2])
                                   }
                                   xx <- dplyr::bind_rows(xx) %>%
                                     dplyr::group_by(dplyr::across(dplyr::all_of("Info_pos"))) %>%
                                     dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") %>%
                                     dplyr::select(-c("Info_pos"))
                                 }

                                 names(xx) <- paste0(feat_prefix, 1:ncol(xx))
                                 xx$Info_protein_id <- gsub("/", "", fn_un[i], fixed = TRUE)
                                 return(xx)
                               },
                               fn_un = fn_un, fn = fn, fl  = fl,
                               feat_prefix = feat_prefix,
                               cl = cl)

  # Concatenate results
  X <- reslist %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(dplyr::across(dplyr::all_of("Info_protein_id"))) %>%
    dplyr::mutate(Info_pos = 1:dplyr::n()) %>%
    dplyr::select(dplyr::starts_with("Info_"), dplyr::everything())


  if(!is.null(save_folder)) {
    if(is.null(save_file)) {
      save_file <- paste0(save_folder, "/esm_features_proteins.rds")
    } else {
      save_file <- paste0(save_folder, "/", save_file)
    }
    saveRDS(X, save_file)
  }

  if(delete_originals) file.remove(fl)
  return(X)

}
