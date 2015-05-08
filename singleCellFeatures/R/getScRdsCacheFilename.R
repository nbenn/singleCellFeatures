#' Get path/filename of a PlateLocation/WellLocation object single cell data
#' cache
#'
#' Upon downloading from openBIS and importin to R, all single cell data
#' associated with a PlateLocation/WellLocation is cached in an .rds file. This
#' function gets the filepath of this cache file.
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The path/filename of the single cell data cache of
#'         PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getScRdsCacheFilename(plate)
#' 
#' @export
getScRdsCacheFilename <- function(x) {
  UseMethod("getScRdsCacheFilename", x)
}

#' @export
getScRdsCacheFilename.WellLocation <- function(x) {
  # load file containing some needed paths
  data(settingsDatabase, envir = environment())
  name <- paste0(x$plate, "_", x$row, x$col, "_sc_all.rds")
  path <- paste(settings.database$openBIS.data, x$space, x$group, x$experiment, 
                x$plate, "WellData", name, sep="/")
  return(path)
}

#' @export
getScRdsCacheFilename.PlateLocation <- function(x) {
  # load file containing some needed paths
  data(settingsDatabase, envir = environment())
  name <- paste0(x$plate, "_sc_all.rds")
  path <- paste(settings.database$openBIS.data, x$space, x$group, x$experiment, 
                x$plate, name, sep="/")
  return(path)
}

#' @export
getScRdsCacheFilename.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}