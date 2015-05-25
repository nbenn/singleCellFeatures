#' Remove some bad image data
#'
#' Given a PlateData/WellData object, remove all images where the cell count
#' does not lie within the 0.05 and 0.95 quantiles.
#'
#' @param x The PlateData/WellData object of interest.
#' 
#' @return A cleaned up version of the input object.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' cleaned <- cleanData(data)
#' 
#' @export
cleanData <- function(x) {
  UseMethod("cleanData", x)
}

#' @export
cleanData.PlateData <- function(x) {
  x$data <- lapply(x$data, cleanData)
  return(x)
}

#' @export
cleanData.WellData <- function(x) {
  quant1 <- x$data[[1]]$counts.quantiles[1]
  quant2 <- x$data[[1]]$counts.quantiles[2]
  discard <- sapply(x$data, function(img) {
    if(img$counts.cells > quant1 & img$counts.cells < quant2) return(TRUE)
    else return(FALSE)
  })
  if(sum(!discard) > 0) {
    message("well ", stri_pad_right(getWellName(x), 3), ": discarding ",
            sum(!discard), " images because count.cells not in [", quant1, ", ",
            quant2, "]")
    x$data <- x$data[discard]
  }
  return(x)
}

#' @export
cleanData.default <- function(x) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}