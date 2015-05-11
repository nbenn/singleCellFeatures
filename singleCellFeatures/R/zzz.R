.onLoad <- function(libname, pkgname) {
  options(error=traceback)
  invisible(NULL)
}

.onUnload <- function(libpath) {
  options(error=NULL)
  invisible(NULL)
}