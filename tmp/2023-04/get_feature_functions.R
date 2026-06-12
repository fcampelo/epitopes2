#' Print available feature calculation functions
#'
#' Prints the feature calculation functions available in the **epitopes**
#' package
#'
#' This routine prints the names of the feature calculation functions available
#' in the **epitopes** package, which are used within **calc_features()**.
#'
#' @return vector of function names.
#'
#' @examples
#' get_feature_functions()
#'
#' @export

get_feature_functions <- function(){

  # Get only functions with "calc_" in the name
  con.list <- ls("package:epitopes")
  con.list <- con.list[grep(pattern = "calc_", con.list)]
  con.list <- con.list[-which(con.list == "calc_features")]
  return(con.list)
}
