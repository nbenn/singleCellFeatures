#' Covert compatible objects to WellLocation
#'
#' Build WellLocation object from other objects, such as WellData.
#'
#' @param The object to be used as basis for new WellLocation object
#' 
#' @return The new WellLocation object
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- convertToWellLocation(plate)
#' 
#' @export
convertToWellLocation <- function(x) {
  UseMethod("convertToWellLocation", x)
}

#' @export
convertToWellLocation.WellData <- function(x) {
  barcode <- getBarcode(x)
  row <- x$meta$well.row
  col <- x$meta$well.col
  return(WellLocation(barcode, row, col))
}

#' @export
convertToWellLocation.default <- function(x) {
  stop("can only deal with WellData objects.")
}