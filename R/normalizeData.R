#' Combine all data in a given Data object
#'
#' Given a PlateData/WellData/ImageData object, rbind all data contained in the
#' hierarchy below.
#'
#' @param x The MatData/PlateData/WellData/ImageData object of interest.
#' 
#' @return A list of data frames for the combined objects.
#' 
#' @examples
#' plate    <- PlateData(PlateLocation("J101-2C"))
#' plate.df <- normalizeData(plate)
#' 
#' @export
normalizeData <- function(x, ...) {
  UseMethod("normalizeData", x)
}

#' @export
normalizeData.PlateData <- function(x, features="intensity", method="Bscore") {

}

#' @export
normalizeData.WellData <- function(x) {

}

#' @export
normalizeData.ImageData <- function(x) {

}

#' @export
normalizeData.default <- function(x, ...) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}