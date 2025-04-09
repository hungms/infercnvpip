.onLoad <- function(libname, pkgname) {
  ### load default packages
  packages <- c(
    "Seurat", "infercnv", "dplyr", "magrittr")
  invisible(lapply(packages, library, character.only = TRUE))

  ### start up settings
  set.seed(123)
  errorMessage <- NULL
}
