get_aa_codes <- function(nletters = 1){
  if (nletters == 1){
    return(c("A","C","D","E","F",
             "G","H","I","K","L",
             "M","N","P","Q","R",
             "S","T","V","W","Y"))
  } else if (nletters == 3){
    return(c("ALA", "CYS", "ASP", "GLU", "PHE",
             "GLY", "HIS", "ILE", "LYS", "LEU",
             "MET", "ASN", "PRO", "GLN", "ARG",
             "SER", "THR", "VAL", "TRP", "TYR"))
  } else stop("Only 1 or 3 letter codes are available")
}
