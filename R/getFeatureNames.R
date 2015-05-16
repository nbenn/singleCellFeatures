#' List feature names of a Data object
#'
#' Given a MatData/PlateData/WellData/ImageData object, list all currently
#' stored features.
#'
#' @param x The MatData/PlateData/WellData/ImageData object of interest
#' 
#' @return A list of all currently available features
#' 
#' @examples
#' plate    <- PlateData(PlateLocation("J101-2C"))
#' features <- getFeatureNames(plate)
#' 
#' @export
getFeatureNames <- function(x) {
  UseMethod("getFeatureNames", x)
}

#' @export
getFeatureNames.MatData <- function(x) {
  return(names(x$data))
}

#' @export
getFeatureNames.PlateData <- function(x) {
  return(getFeatureNames(x$data[[1]]))
}

#' @export
getFeatureNames.WellData <- function(x) {
  return(getFeatureNames(x$data[[1]]))
}

#' @export
getFeatureNames.ImageData <- function(x) {
  vec <- unlist(lapply(x$data.vec, names))
  mat <- unlist(lapply(x$data.mat, colnames))
  lst <- unlist(lapply(x$data.lst, names))
  result <- sort(c(vec, mat, lst))
  names(result) <- NULL
  return(result)
}

#' @export
getFeatureNames.default <- function(x) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}