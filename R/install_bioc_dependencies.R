#' Install Bioconductor dependencies
#'
#' This function installs the latest versions of all Bioconductor dependencies
#' of the `epitopes` package, namely:
#' \itemize{
#'    \item Biostrings
#' }
#'
#' It is essential that the package(s) above be installed for `epitopes` to
#' work properly.
#'
#' @param bioc.args list containing further arguments to be passed
#' down to `BiocManager::install()`.
#' @param force logical: reinstall already-installed packages?
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   install_bioc_dependencies()
#' }
#'
#' @return No return value, called for side effects.

install_bioc_dependencies <- function(bioc.args = list(),
                                      force = FALSE){
  # ================== Sanity checks ==================
  assertthat::assert_that(is.list(bioc.args),
                          is.logical(force), length(force) == 1)

  pkgs <- c("Biostrings")
  makeInst <- FALSE

  if (force){
    makeInst <- TRUE
    bioc.args$force <- TRUE
  } else {
    x <- rownames(utils::installed.packages())
    idx  <- which(!pkgs %in% x)
    if (length(idx) > 0){
      pkgs <- pkgs[idx]
      makeInst <- TRUE
    }
  }

  if (makeInst){
    message(paste0("\nInstalling package(s): ",
                   paste(pkgs, collapse = ", "),
                   "\nfrom BioConductor version ",
                   BiocManager::version()))
    bioc.args$pkgs <- pkgs
    bioc.args$ask  <- TRUE

    do.call(BiocManager::install, bioc.args)
  }
  invisible(TRUE)
}
