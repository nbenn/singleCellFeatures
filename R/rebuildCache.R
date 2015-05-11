#' Rebuild cache files for certain Data objects
#'
#' Given an existing MatData cache file and some feature files that were not
#' imported correctly, update and save it to the corresponding cache location.
#'
#' @param x MatData object
#' 
#' @return The updated MatData object. The updated cache file is saved to its
#'         corresponding location.
#' 
#' @examples
#' data <- MatData(PlateLocation("J101-2C"))
#' rebuildCache(data)
#' 
#' @export
rebuildCache <- function(x) {
  UseMethod("rebuildCache", x)
}

#' @export
rebuildCache.MatData <- function(x) {
  loc  <- convertToPlateLocation(x)
  errdir <- paste0(getLocalPath(loc), "/", "Errors")
  if(!dir.exists(errdir)) stop("directory ", errdir, " not found")
  data <- readMatFeatureHelper(errdir)
  x$data <- c(x$data, data)
  x$data <- x$data[order(names(x$data))]
  message("added features:\n  ", paste0(names(data), collapse="\n  "))
  saveToCache(x, force.write=TRUE)
  return(x)
}

#' @export
rebuildCache.default <- function(x) {
  stop("can only deal with MatData objects.")
}