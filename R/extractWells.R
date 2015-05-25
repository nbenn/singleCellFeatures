#' Extract WellData objects form a PlateData object
#' 
#' Given a PlateData object and a single WellLocation object of a list of
#' WellLocation objects, return either a PlateData object containing only the 
#' specified subset of wells or a single WellData/a list of WellData objects.
#'
#' @param x          The original PlateData object, holding the superset of
#'                   wells
#' @param wells      A single WellLocation object/a list of WellLocation
#'                   objects, a single integer/vector of integers (in 1:384),
#'                   or a single string/vector of strings (e.g. A5) specifying
#'                   the targets.
#' @param keep.plate A boolean indicating whether to keep the PlateData object
#'                   intact or return WellData object(s)
#'
#' @return Either a PlateData or a WellData object holding only the data
#'         corresponding to the specified subset of wells
#'
#' @examples
#' plate  <- PlateData(PlateLocation("J101-2C"))
#' a23 <- extractWells(plate, "A23", FALSE)
#' wells <- extractWells(plate, c(135, 155))
#'
#' @export
extractWells <- function(x, ...) {
  UseMethod("extractWells", x)
}
#' @export
extractWells.PlateData <- function(x, wells, keep.plate=TRUE) {
  # turn input into a list of WellLocation objects
  if(any(class(wells) == "WellLocation")) {
    wells <- list(wells)
  } else if(is.numeric(wells)) {
    wells <- unique(wells)
    wells <- lapply(wells, function(well) {
      bc <- getBarcode(x)
      ind <- getWellIndex2D(well)
      return(WellLocation(bc, ind$wel.row, ind$wel.col))
    })
  } else if(is.character(wells)) {
    wells <- unique(wells)
    wells <- lapply(wells, function(well) {
      bc <- getBarcode(x)
      row <- substr(well, 1, 1)
      col <- substring(well, 2)
      return(WellLocation(bc, row, col))
    })
  }
  # input validation
  if(!all(sapply(wells, function(well) {
      return(any(class(well) == "WellLocation"))
  }))) stop("can only work with a list of WellLocation objects")
  barcodes <- sapply(wells, getBarcode)
  if(length(unique(barcodes)) != 1) {
    stop("can only deal with WellLocations on the same plate")
  }
  if(unique(barcodes) != getBarcode(x)) {
    stop("WellLocations have to be on the same plate as the data")
  }
  # process wells list: get names of requested wells
  well.names <- sapply(wells, getWellName)
  # save well caches
  l_ply(well.names, function(name, data) {
    saveToCache(data$data[[name]])
  }, x)
  # return results
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
extractWells.default <- function(x, ...) {
  stop("can only deal with PlateData objects.")
}