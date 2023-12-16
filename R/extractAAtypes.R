extractAAtypes <- function(x){
  AAtypes <- list(Tiny      = c("A", "C", "G", "S", "T"),
                  Small     = c("A", "B", "C", "D", "G", "N", "P", "S", "T", "V"),
                  Aliphatic = c("A", "I", "L", "V"),
                  Aromatic  = c("F", "H", "W", "Y"),
                  NonPolar  = c("A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y"),
                  Polar     = c("D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z"),
                  Charged   = c("B", "D", "E", "H", "K", "R", "Z"),
                  Basic     = c("H", "K", "R"),
                  Acidic    = c("B", "D", "E", "Z"))

  f_AAc <- as.data.frame(
    lapply(AAtypes,
           function(pat, x){
             sum(strsplit(x, split = "")[[1]] %in% pat) / nchar(x)},
           x = x))

  return(f_AAc)
}
