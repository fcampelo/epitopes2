extractMolWeight <- function(x){
  y <- data.frame(
  Code1    = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
  AA_MW = c(89.09, 121.16, 133.10, 147.13, 165.19,
            75.07, 155.16, 131.17, 146.19, 131.17,
            149.21, 132.12, 115.13, 146.15, 174.20,
            105.09, 119.12, 117.15, 204.22, 181.19))

  cnts <- as.data.frame(as.matrix(table(strsplit(x, split = "")[[1]])))
  cnts$Code1 <- rownames(cnts)
  cnts <-  dplyr::left_join(y, cnts, by = "Code1")
  cnts$AA_MW <- cnts$AA_MW * cnts$V1
  return(data.frame(MolWeight = sum(cnts$AA_MW, na.rm = TRUE)))

  # cnts <- y %>%
  #   dplyr::left_join(cnts, by = "Code1") %>%
  #   dplyr::mutate("AA_MW" = .data$AA_MW * .data$V1) %>%
  #   dplyr::select(-c("Code1", "V1")) %>%
  #   dplyr::summarise(MolWeight = sum(.data$AA_MW, na.rm = TRUE))

}
