#' Get the pathogen name of a PlateLocation/WellLocation object
#'
#' Give a PlateLocation/WellLocation object, return the pathogen name
#'
#' @param x A PlateLocation/WellLocation object
#' 
#' @return The pathogen name of the PlateLocation/WellLocation object (as a
#'         string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' pathogen <- getPathogen(plate)
#' 
#' @export
getPathogen <- function(x) {
  UseMethod("getPathogen", x)
}

#' @export
getPathogen.DataLocation <- function(x) {
  name1 <- unlist(strsplit(x$group, "_"))[1]
  name2 <- unlist(strsplit(x$experiment, "-"))[1]
  if(name1 != name2) stop("unrecognized naming scheme")
  return(tolower(name1))
}

#' @export
getPathogen.WellData <- function(x) {
  return(getPathogen(x$meta))
}

#' @export
getPathogen.WellMetadata <- function(x) {
  return(x$experiment.pathogen)
}

#' @export
getPathogen.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}