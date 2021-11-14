call_protr <- function(SEQs, feat.name, txt.opts, dfnames, ncpus){

  fn <- paste0("extract", feat.name)

  # Check if function exists
  if(!(fn %in% ls('package:protr') | fn %in% ls('package:epitopes'))){
    warning("Function ", fn, "() not found.\nSkipping...")
    return(FALSE)
  }

  if (fn %in% ls('package:protr')) fn <- paste0("protr::", fn)

  # Remove or replace invalid AA codes, depending on feature:
  if (grepl("Gap$", fn)) {
    SEQs <- sapply(SEQs, function(x){gsub("B|J|O|U|X|Y", "-", x)})
  } else {
    SEQs <- sapply(SEQs, function(x){gsub("B|J|O|U|X|Y", "", x)})
  }

  message("   ", txt.opts[1], ":", feat.name)
  myargs <- get_feature_args(feat.name)
  y <- mypblapply(SEQs,
                  function(x, fn, myargs){
                    myargs$x <- x
                    do.call(eval(parse(text = fn)),
                            args = myargs)},
                  fn = fn, myargs = myargs,
                  ncpus = ncpus) %>%
    dplyr::bind_rows() %>%
    dplyr::rename_with(~ ifelse(feat.name == .x,
                                paste0("feat_", txt.opts[1], "_", .x),
                                paste0("feat_", txt.opts[1], "_", feat.name, "_", .x)))

  return(y)
}
