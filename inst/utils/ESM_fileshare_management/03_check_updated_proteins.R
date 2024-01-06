library(dplyr)

feat_path <- "d:/IEDB_proteins_ESM1b_features/features/"
epitopes_path <- "../../../data/epitopes.rds"


#=====

dirs <- dir(feat_path, full.names = TRUE)

errlist  <- data.frame(folder = character(),
                       file   = character())

filelist <- data.frame(Info_protein_id = character(),
                       folder  = character())
for (j in seq_along(dirs)){
  cat("\n")
  t0 <- Sys.time()
  fl <- dir(dirs[j], full.names = TRUE)
  fn <- dir(dirs[j], full.names = FALSE)
  filelist <- rbind(filelist,
                    data.frame(Info_protein_id = gsub("\\.rds", "", fn),
                               folder  = sprintf("Folder%02d", j)))
  write.table(filelist, "filemap.tsv",
              row.names = FALSE, quote = FALSE, sep = "\t")
  for (i in seq_along(fl)){
    cat(sprintf("\rFolder %02d/%02d file %03d/%03d",
                j, length(dirs), i, length(fl)))
    X <- readRDS(fl[i])
    if(!all(c("Info_pos", "Info_AA") %in% names(X))){
      cat("\t *")
      errlist <- rbind(errlist,
                       data.frame(folder = sprintf("Folder%02d", j),
                                  file   = fn[i]))
      write.table(errlist, "errlist.tsv",
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
  a <- difftime(Sys.time(), t0)
  cat(":", as.numeric(a), attr(a, "units"))
}

#=====

epitopes <- readRDS(epitopes_path)

get_uniques <- function(x){
  x <- strsplit(x, split = ",")
  sapply(x, function(y) paste(unique(y), collapse = ","))
}

X <- epitopes %>%
  select(protein_id, sourceOrg_id) %>%
  group_by(protein_id) %>%
  summarise(sourceOrg_id = paste(sourceOrg_id, collapse = ",")) %>%
  rename(Info_protein_id = protein_id,
         Info_organism_id = sourceOrg_id) %>%
  mutate(Info_organism_id = get_uniques(Info_organism_id))

filelist <- filelist %>%
  left_join(X, by = "Info_protein_id")

write.table(filelist, "filemap.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")
