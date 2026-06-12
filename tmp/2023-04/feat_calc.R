feat_calc <- function(x, max.N){

  # Prepare auxiliary information
  AAtypes <- list(Tiny      = c("A", "C", "G", "S", "T"),
                  Small     = c("A", "B", "C", "D", "G", "N", "P", "S", "T", "V"),
                  Aliphatic = c("A", "I", "L", "V"),
                  Aromatic  = c("F", "H", "W", "Y"),
                  NonPolar  = c("A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y"),
                  Polar     = c("D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z"),
                  Charged   = c("B", "D", "E", "H", "K", "R", "Z"),
                  Basic     = c("H", "K", "R"),
                  Acidic    = c("B", "D", "E", "Z"))

  aa_prop <- data.frame(
    Code1    = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
    CT_group = c(0, 1, 6, 6, 2, 0, 4, 2, 5, 2, 3, 4, 2, 4, 5, 3, 3, 0, 4, 3),
    nC = c(3, 3, 4, 5, 9, 2, 6, 6, 6, 6, 5, 4, 5, 5, 6, 3, 4, 5, 11, 9),
    nH = c(7, 7, 7, 9, 11, 5, 9, 13, 14, 13, 11, 8, 9, 10, 14, 7, 9, 11, 12, 11),
    nN = c(1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 2, 1, 2, 4, 1, 1, 1, 2, 1),
    nO = c(2, 2, 4, 4, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 3, 3, 2, 2, 3),
    nS = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    AA_MW = c(89.09, 121.16, 133.10, 147.13, 165.19,
              75.07, 155.16, 131.17, 146.19, 131.17,
              149.21, 132.12, 115.13, 146.15, 174.20,
              105.09, 119.12, 117.15, 204.22, 181.19))

  exception_AA <- c("B", "Z", "X", "J")

  # ==================== Calculate features ==================== #

  # ----- Sequence Entropy -----
  cnts <- stringr::str_count(x, LETTERS)
  cnts <- cnts[cnts > 0] / sum(cnts)
  f_SE   <- data.frame(feat_seq_entropy = -sum(cnts * log2(cnts)))


  # ----- Number of Atoms -----
  cnts <- stringr::str_count(x, aa_prop$Code1)
  f_NAt <- as.data.frame(t(colSums(cnts * aa_prop[, c("nC", "nH", "nN", "nO", "nS")])))
  names(f_NAt) <- paste0("feat_", gsub("n", "", names(f_NAt)), "_atoms")


  # ----- Molecular weight -----
  f_MW <- data.frame(feat_molecular_weight = sum(cnts * aa_prop$AA_MW))


  # ----- AA composition -----
  f_AAc <- as.data.frame(
    lapply(AAtypes,
           function(pat, x){
             sum(stringr::str_count(x, pat)) / nchar(x)},
           x = x))
  names(f_AAc) <- paste0("feat_Perc_", names(f_AAc))


  # ----- AA Descriptors -----
  y <- strsplit(x, split = "")[[1]]
  if (any(y %in% exception_AA)) y <- y[-which(y %in% exception_AA)]
  f_AAd <- as.data.frame(t(colMeans(Peptides::aaDescriptors(y),
                                    na.rm = TRUE)))
  names(f_AAd) <- paste0("feat_", gsub("\\.[1-9]$", "", names(f_AAd)))


  # ----- Conjoint triads -----
  f_CT <- as.data.frame(matrix(0, ncol = 7 ^ 3, nrow = 1))
  names(f_CT) <- apply(X   = expand.grid(0:6, 0:6, 0:6),
                       FUN = paste, MARGIN = 1, collapse = "")
  f_CT <- f_CT[, order(names(f_CT))]

  idx <- stringr::str_locate(paste(aa_prop$Code1, collapse = ""),
                             strsplit(x, "")[[1]])[, 1]
  y <- paste(aa_prop$CT_group[idx], collapse = "")
  y <- gsub("NA", "*", y)
  ct <- stringr::str_sub(y,
                         start = 1:(nchar(y) - 2),
                         end   = 3:nchar(y))
  tmp  <- as.data.frame(t(as.matrix(table(ct)))) / length(ct)
  torm <- grep("\\*", names(tmp))
  if(length(torm) > 0) tmp <- tmp[, -torm]

  f_CT[, sapply(names(tmp), function(k) which(names(f_CT) == k))] <- tmp
  names(f_CT) <- paste0("feat_CT", names(f_CT))


  # ----- N-peptide composition -----
  # Prepare dataframe with all possible columns
  f_Kmer <- as.data.frame(matrix(0, ncol = sum(20 ^ (1:max.N)), nrow = 1))
  fns <- sapply(1:max.N,
                function(i){
                  sort(apply(expand.grid(data.frame(aa_prop$Code1)[, rep(1,i)]),
                             MARGIN = 1, FUN = paste, collapse = ""))})

  names(f_Kmer) <- unlist(fns)

  for (i in 1:max.N){
    i1  <- 1:(nchar(x) - i + 1)
    i2  <- i1 + i - 1
    y   <- stringr::str_sub(x, i1, i2)
    tmp <- as.data.frame(t(as.matrix(table(y)))) / length(y)
    torm <- grep("[BJXZ]", names(tmp))
    if(length(torm) > 0) tmp <- tmp[, -torm]

    f_Kmer[, sapply(names(tmp), function(k) which(names(f_Kmer) == k))] <- tmp
  }

  names(f_Kmer) <- paste0("feat_Perc_", names(f_Kmer))

  return(cbind(f_SE, f_NAt, f_MW, f_AAc, f_AAd, f_CT, f_Kmer))
}
