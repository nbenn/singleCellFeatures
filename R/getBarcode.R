#' Get the barcode of a PlateLocation/WellLocation object
#'
#' Given a PlateLocation/WellLocation object, return the barcode
#'
#' @param x A Data/DataLocation object
#' 
#' @return The barcode of the PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' barcode <- getBarcode(plate)
#' 
#' @export
getBarcode <- function(x) {
  UseMethod("getBarcode", x)
}
#' @export
getBarcode.DataLocation <- function(x) {
  return(x$plate)
}

#' @export
getBarcode.MatData <- function(x) {
  return(x$meta$plate.barcode)
}
#' @export
getBarcode.PlateData <- function(x) {
  return(x$meta$plate.barcode)
}
#' @export
getBarcode.WellData <- function(x) {
  return(x$meta$plate.barcode)
}
#' @export
getBarcode.ImageData <- function(x) {
  return(x$plate)
}

#' @export
getBarcode.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) or ",
       "Data (ImageData/WellData/PlateData) objects.")
}