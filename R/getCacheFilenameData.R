#' Get path/filename of a PlateLocation/WellLocation object single cell data
#' cache
#'
#' Upon downloading from openBIS and importing to R, all single cell data
#' associated with a PlateLocation/WellLocation is cached in an .rds file. This
#' function retrieves the filepath of this cache file.
#'
#' @param x PlateLocation/WellLocation object
#' 
#' @return The path/filename of the single cell data cache of
#'         PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getCacheFilenameData(plate)
#' 
#' @export
getCacheFilenameData <- function(x) {
  UseMethod("getCacheFilenameData", x)
}

#' @export
getCacheFilenameData.WellLocation <- function(x) {
  # load file containing some needed paths
  config <- configGet()
  name <- paste0(x$plate, "_", x$row, x$col, "_sc_all.rds")
  path <- paste(config$dataStorage$dataDir, x$space, x$group, x$experiment,
                x$plate, "WellData", name, sep="/")
  return(path)
}

#' @export
getCacheFilenameData.PlateLocation <- function(x) {
  # load file containing some needed paths
  config <- configGet()
  name <- paste0(x$plate, "_sc_all.rds")
  path <- paste(config$dataStorage$dataDir, x$space, x$group, x$experiment,
                x$plate, name, sep="/")
  return(path)
}

#' @export
getCacheFilenameData.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}