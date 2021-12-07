.onAttach <- function(...) {
  # Check if Biostrings is installed
  bsp <- suppressWarnings(utils::packageDescription("Biostrings"))

  if(all(is.na(bsp))) {
    packageStartupMessage("\nPackage 'Biostrings' not detected.\n",
                          "Please run install_bioc_dependencies()\n",
                          "before using the epitopes package.")
  } else if(utils::packageVersion("Biostrings") < '2.60.0') {
    packageStartupMessage("\nPackage 'Biostrings' version 2.60.0 or later is required\n",
                          "Please run install_bioc_dependencies(force = TRUE)\n",
                          "before using the epitopes package.")
  }

}
