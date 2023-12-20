get_feature_args <- function(feat.name){

  myargs <- switch(EXPR = feat.name,
                   AAC = list(),
                   DC = list(),
                   TC = list(),
                   CTDC = list(),
                   CTDT = list(),
                   CTDD = list(),
                   CTriad = list(),
                   SOCN = list(nlag = 3),
                   QSO = list(nlag = 3, w = 0.1),
                   ScalesGap = list(propmat = t(stats::na.omit(as.matrix(protr::AAindex[, 7:26]))),
                                    pc      = 5,
                                    lag     = 3,
                                    scale   = TRUE,
                                    silent  = TRUE),
                   Entropy = list(),
                   Atoms = list(),
                   MolWeight = list(),
                   AAtypes = list(),
                   BLOSUM = list(blnames = rownames(protr::AABLOSUM62),
                                 eig = eigen(protr::AABLOSUM62)))

  return(myargs)
}
