call_feat_functions <- function(feat.name, SEQ, myargs){

  internal.functions <- c("extractAAtypes", "extractAtoms", "extractEntropy",
                          "extractLegacyFeatures", "extractMolWeight")

  fn <- paste0("extract", feat.name)

  if (!(fn %in% internal.functions)) fn <- paste0("protr::", fn)

  # Remove or replace invalid AA codes, depending on feature:
  if (grepl("Gap$", fn)) {
    SEQ <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "-", toupper(SEQ))
  } else {
    SEQ <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", SEQ)
  }

  #myargs   <- get_feature_args(feat.name)
  myargs$x <- SEQ

  if(feat.name == "BLOSUM") {
    require(protr)
    AABLOSUM62 <- protr::AABLOSUM62
  }

  y <- do.call(what = eval(parse(text = fn)), args = myargs, quote = TRUE)

  if(!is.data.frame(y)) y <- as.data.frame(t(y))

  y <- dplyr::rename_with(y,
                          ~ ifelse(feat.name == .x,
                                   paste0("feat_", .x),
                                   paste0("feat_", feat.name, "_", .x)))

  return(y)
}


# call_feat_functions <- function(feat.name, SEQs){
#
#   internal.functions <- c("extractAAtypes", "extractAtoms", "extractEntropy",
#                           "extractLegacyFeatures", "extractMolWeight")
#
#   fn <- paste0("extract", feat.name)
#
#   # Check if function exists
#   # if(!(fn %in% ls('package:protr') | fn %in% ls('package:epitopes'))){
#   #   warning("Function ", fn, "() not found.\nSkipping...")
#   #   return(FALSE)
#   # }
#
#   if (!(fn %in% internal.functions)) fn <- paste0("protr::", fn)
#
#   # Remove or replace invalid AA codes, depending on feature:
#   if (grepl("Gap$", fn)) {
#     SEQs <- sapply(SEQs, function(x){gsub("[^ACDEFGHIKLMNPQRSTVWY]", "-", toupper(x))})
#   } else {
#     SEQs <- sapply(SEQs, function(x){gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", x)})
#   }
#
#   myargs <- get_feature_args(feat.name)
#
#   if(feat.name == "BLOSUM") {
#     AABLOSUM62 <- protr::AABLOSUM62
#   }
#
#   y <- lapply(SEQs,
#               function(x, fn, myargs){
#                 myargs$x <- x
#                 do.call(eval(parse(text = fn)),
#                         args = myargs)},
#               fn = fn, myargs = myargs)
#
#   y <- dplyr::bind_rows(unname(y))
#
#   y <- dplyr::rename_with(y,
#                           ~ ifelse(feat.name == .x,
#                                    paste0("feat_", .x),
#                                    paste0("feat_", feat.name, "_", .x)))
#
#   return(y)
# }
