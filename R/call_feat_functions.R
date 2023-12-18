call_feat_functions <- function(SEQs, feat.name){

  internal.functions <- c("extractAAtypes", "extractAtoms", "extractEntropy",
                          "extractLegacyFeatures", "extractMolWeight",
                          "extractBLOSUM")

  fn <- paste0("extract", feat.name)

  # Check if function exists
  # if(!(fn %in% ls('package:protr') | fn %in% ls('package:epitopes'))){
  #   warning("Function ", fn, "() not found.\nSkipping...")
  #   return(FALSE)
  # }

  if (!(fn %in% internal.functions)) fn <- paste0("protr::", fn)

  # Remove or replace invalid AA codes, depending on feature:
  if (grepl("Gap$", fn)) {
    SEQs <- sapply(SEQs, function(x){gsub("[^ACDEFGHIKLMNPQRSTVWY]", "-", toupper(x))})
  } else {
    SEQs <- sapply(SEQs, function(x){gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", x)})
  }

  myargs <- get_feature_args(feat.name)

  AABLOSUM62 <- protr::AABLOSUM62 # just to load the matrix into the search path


  y <- lapply(SEQs,
              function(x, fn, myargs){
                myargs$x <- x
                do.call(eval(parse(text = fn)),
                        args = myargs)},
              fn = fn, myargs = myargs) %>%
    unname() %>%
    dplyr::bind_rows()
    # dplyr::rename_with(~ ifelse(feat.name == .x,
    #                             paste0("feat_", txt.opts[1], "_", .x),
    #                             paste0("feat_", txt.opts[1], "_", feat.name, "_", .x)))

  return(y)
}
