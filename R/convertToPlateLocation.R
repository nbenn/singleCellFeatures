#' Covert compatible objects to PlateLocation
#'
#' Build PlateLocation object from other objects, such as PlateData,
#' MatData or WellLocation
#'
#' @param The object to be used as basis for new PlateLocation object
#' 
#' @return The new PlateLocation object
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- convertToPlateLocation(plate)
#' 
#' @export
convertToPlateLocation <- function(x) {
  UseMethod("convertToPlateLocation", x)
}

#' @export
convertToPlateLocation.DataLocation <- function(x) {
  barcode <- getBarcode(x)
  return(PlateLocation(barcode))
}

#' @export
convertToPlateLocation.Data <- function(x) {
  barcode <- getBarcode(x)
  return(PlateLocation(barcode))
}

#' @export
convertToPlateLocation.PlateAggregate <- function(x) {
  return(x$plate)
}

#' @export
convertToPlateLocation.default <- function(x) {
  stop("can only deal with PlateData/MatData/WellLocation objects.")
}