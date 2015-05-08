#' Get the path/filename of the metadata chache for a
#' PlateLocation/WellLocation object
#'
#' Given a PlateLocation/WellLocation object, return the path/filename of the
#' corresponding metadata chache file.
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The path/filename of the metadata chache for a 
#'         PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getMetadataCacheFilename(plate)
#' 
#' @export
getMetadataCacheFilename <- function(x) {
  UseMethod("getMetadataCacheFilename", x)
}

#' @export
getMetadataCacheFilename.DataLocation <- function(x) {
  data(settingsDatabase, envir = environment())
  name <- paste0(x$plate, "_metadata.rds")
  path <- paste(settings.database$openBIS.data, x$space, x$group, x$experiment, 
                x$plate, name, sep="/")
  return(path)
}

#' @export
getMetadataCacheFilename.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}