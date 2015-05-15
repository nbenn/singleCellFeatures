#' Get the local path of a PlateLocation/WellLocation object 
#'
#' Given a PlateLocation/WellLocation object, return the local path
#'
#' @param x A PlateLocation/WellLocation object
#' 
#' @return The local path of the PlateLocation/WellLocation object (as a
#'         string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getLocalPath(plate)
#' 
#' @export
getLocalPath <- function(x) {
  UseMethod("getLocalPath", x)
}

#' @export
getLocalPath.DataLocation <- function(x) {
  # load file containing some needed paths
  config <- configGet()
  path <- paste(config$dataStorage$dataDir, x$space, x$group, x$experiment,
                x$plate, sep="/")
  return(path)
}

#' @export
getLocalPath.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}