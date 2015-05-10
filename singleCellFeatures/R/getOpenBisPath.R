#' Get the path of a PlateLocation/WellLocation object on openBIS
#'
#' Given a PlateLocation/WellLocation object, return the openBIS location path
#'
#' @param x A PlateLocation/WellLocation object
#' 
#' @return The openBIS location path of the PlateLocation/WellLocation object
#'         (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getOpenBisPath(plate)
#' 
#' @export
getOpenBisPath <- function(x) {
  UseMethod("getOpenBisPath", x)
}

#' @export
getOpenBisPath.DataLocation <- function(x) {
  path <- paste0("/", x$space, "/", x$group, "/", x$experiment, "/", x$plate)
  return(path)
}

#' @export
getOpenBisPath.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}