#' Build cache files for certain Data objects
#'
#' Given a MatData/WellData object, save it to the corresponding cache location
#'
#' @param x           MatData/WellData object
#' @param force.write Logical value indicating whether the cache file is
#'                    overwritten if already present
#' 
#' @return NULL (invisibly). The cache file is saves to its corresponding
#'         location.
#' 
#' @examples
#' plate <- PlateData(PlateLocation("J101-2C"))
#' well  <- extractWell(plate, WellLocation(plate, "A", 5), FALSE)
#' saveToCache(well)
#' 
#' @export
saveToCache <- function(x, force.write=FALSE) {
  UseMethod("saveToCache", x)
}

#' @export
saveToCache.MatData <- function(x, force.write=FALSE) {
  loc  <- convertToPlateLocation(x)
  file <- getCacheFilenameData(loc)
  dir  <- dirname(file)
  if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  if(!file.exists(file) | force.write) {
    saveRDSMC(x$data, file)
  }
  # remove outdated well caches
  dir <- dirname(getCacheFilenameData(WellLocation(getBarcode(loc), "A", 1)))
  if(dir.exists(dir)) unlink(dir, recursive=TRUE)
  invisible(NULL)
}

#' @export
saveToCache.WellData <- function(x, force.write=FALSE) {
  if(checkCompletenessFeature(x)$result & checkCompletenessImage(x)) {
    loc  <- convertToWellLocation(x)
    file <- getCacheFilenameData(loc)
    dir  <- dirname(file)
    if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
    if(!file.exists(file) | force.write) {
      saveRDS(x, file=file, compress="xz")
    }
  } else {
    message("well data is incomplete: not writing cache file.")
  }
  invisible(NULL)
}

#' @export
saveToCache.PlateAggregate <- function(x, force.write=FALSE) {
  loc    <- convertToPlateLocation(x)
  file   <- getCacheFilenameMeta(loc)
  dir    <- dirname(file)
  if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  if(!file.exists(file) | force.write) {
    saveRDS(x$data, file=file)
  }
  invisible(NULL)
}

#' @export
saveToCache.default <- function(x, force.write=FALSE) {
  stop("can only deal with MatData/WellData/PlateAggregate objects.")
}