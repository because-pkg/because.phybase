# This file is traditionally used in R packages for startup functions.

.onAttach <- function(libname, pkgname) {
  # Print the current R package version dynamically
  pkg_version <- utils::packageVersion(pkgname)
  packageStartupMessage(sprintf("This is %s v%s", pkgname, pkg_version))
}
