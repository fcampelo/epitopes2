extractEntropy <- function(x){
  cnts <- as.data.frame(t(as.matrix(table(strsplit(x, split = "")[[1]]))))
  cnts <- cnts / sum(cnts)
  return(data.frame(Entropy = -sum(cnts * log2(cnts))))
}
