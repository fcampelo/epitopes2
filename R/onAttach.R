.onAttach <- function(...) {
  # Check if Biostrings is installed
  pkgs <- utils::installed.packages()[, 3]
  bsp <- grep("Biostrings", names(pkgs))

  if(length(bsp) == 0) {
    packageStartupMessage("\nPackage 'Biostrings' not detected.\n",
                          "Please run install_bioc_dependencies()\n",
                          "before using the epitopes package.")
  } else if(utils::packageVersion("Biostrings") < '2.60.0') {
    packageStartupMessage("\nPackage 'Biostrings' version 2.60.0 or later is required\n",
                          "Please run install_bioc_dependencies(force = TRUE)\n",
                          "before using the epitopes package.")
  }

}
