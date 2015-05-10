#' Get the well name of a WellLocation object
#'
#' Given a WellLocation object, return the well name
#'
#' @param x A WellLocation object
#' 
#' @return The well name of the WellLocation object (as a string)
#' 
#' @examples
#' plate <- WellLocation("J101-2C", "B", 15)
#' well.name <- getWellName(plate)
#' 
#' @export
getWellName <- function(x) {
  UseMethod("getWellName", x)
}

#' @export
getWellName.WellLocation <- function(x) {
  return(paste0(x$row, x$column))
}

#' @export
getWellName.default <- function(x) {
  stop("can only deal with WellLocation objects.")
}