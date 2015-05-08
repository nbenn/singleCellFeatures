#' @export
extractWells <- function(x, wells, keep.plate=TRUE) {
  UseMethod("extractWells", x)
}
#' @export
extractWells.PlateData <- function(x, wells, keep.plate=TRUE) {
  if(any(class(wells) == "WellLocation")) wells <- list(wells)
  if(!all(sapply(wells, function(well) {
    return(any(class(well) == "WellLocation"))
  }))) stop("can only work with a list of WellLocation objects")
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