#' Extract WellData objects form a PlateData object
#' 
#' Given a PlateData object and a single WellLocation object of a list of
#' WellLocation objects, return either a PlateData object containing only the 
#' specified subset of wells or a single WellData/a list of WellData objects.
#'
#' @param x          The original PlateData object, holding the superset of
#'                   wells
#' @param wells      A single WellLocation object/a list of WellLocation objects
#'                   specifying the targets
#' @param keep.plate A boolean indicating whether to keep the PlateData object
#'                   intact or return WellData object(s)
#'
#' @return Either a PlateData or a WellData object holding only the data
#'         corresponding to the specified subset of wells
#'
#' @examples
#' plate  <- PlateData(PlateLocation("J101-2C"))
#' a23 <- extractWells(plate, WellLocation("J101-2C", "A", 23), FALSE)
#' 
#' @export
extractWells <- function(x, wells, keep.plate=TRUE) {
  UseMethod("extractWells", x)
}
#' @export
extractWells.PlateData <- function(x, wells, keep.plate=TRUE) {
  if(any(class(wells) == "WellLocation")) wells <- list(wells)
  if(!all(sapply(
    wells,
    function(well) {
      return(any(class(well) == "WellLocation"))
    }
  ))) stop("can only work with a list of WellLocation objects")
  barcodes <- sapply(wells, getBarcode)
  if(length(unique(barcodes)) != 1)
    stop("can only deal with WellLocations on the same plate")
  if(unique(barcodes) != getBarcode(x))
    stop("WellLocations have to be on the same plate as the data")
  well.names <- sapply(wells, getWellName)
  if(!keep.plate) {
    if(length(well.names) == 1) return(x$data[[well.names]])
    else return(x$data[well.names])
  }
  else {
    x$data <- x$data[well.names]
    return(x)
  }
}
#' @export
extractWells.default <- function(x, wells, keep.plate=TRUE) {
  stop("can only deal with PlateData objects.")
}