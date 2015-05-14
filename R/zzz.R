.onLoad <- function(libname, pkgname) {
  options(error=traceback)
}

.onAttach <- function(libname, pkgname) {
  wellDatabaseCoverage(verbose=FALSE, is.startup=TRUE)
  invisible(NULL)
}

.onUnload <- function(libpath) {
  options(error=NULL)
  invisible(NULL)
}