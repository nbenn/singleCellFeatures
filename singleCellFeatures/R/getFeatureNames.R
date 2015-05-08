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
  mat <- unlist(lapply(x$data.mat, function(group) {
    res <- colnames(group)
    if (is.null(res)) res <- names(group)
    return(res)
  }))
  lst <- unlist(lapply(x$data.lst, names))
  result <- sort(c(vec, mat, lst))
  names(result) <- NULL
  return(result)
}
#' @export
getFeatureNames.default <- function(x) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}