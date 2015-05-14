.onLoad <- function(libname, pkgname) {
  options(error=traceback)
  wellDatabaseCoverage()
  invisible(NULL)
}

.onUnload <- function(libpath) {
  options(error=NULL)
  invisible(NULL)
}