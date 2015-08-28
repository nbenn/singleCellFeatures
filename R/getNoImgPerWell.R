#' Get the number of images per Well
#'
#' Given a Data (MatData/ImageData/WellData/PlateData) object, return the
#' number of images per well (either 6 or 9).
#'
#' @param data A Data (MatData/ImageData/WellData/PlateData) object
#' 
#' @return The number of images per well (either 6 or 9).
#' 
#' @examples
#' plate <- PlateData()PlateLocation("J101-2C"))
#' n.img <- getNoImgPerWell(plate)
#' 
#' @export
getNoImgPerWell <- function(data) {
  UseMethod("getNoImgPerWell", data)
}

#' @export
getNoImgPerWell.MatData <- function(data) {
  lengths <- sapply(data$data, length)
  max.length <- max(lengths)
  ful.length <- lengths[lengths == max.length]
  full.fract <- length(ful.length) / length(lengths)
  if (max.length > 2304) n.imgs <- 9
  else if (full.fract > 0.5) n.imgs <- 6
  else {
    n.imgs <- 9
    warning("unsure about the number of images per well.")
  }
  message("assuming ", n.imgs, " images per well:\nmax legnth: ", max.length,
          ", fraction of full length features: ", format(full.fract, digits=3))

  return(n.imgs)
}

#' @export
getNoImgPerWell.PlateData <- function(data) {
  return(getNoImgPerWell(data$data[[1]]))
}

#' @export
getNoImgPerWell.WellData <- function(data) {
  return(getNoImgPerWell(data$data[[1]]))
}

#' @export
getNoImgPerWell.ImageData <- function(data) {
  res <- data$image.total
  if(!res %in% c(6, 9)) stop("no img is ", res, ", should be in c(6, 9)")
  else return(res)
}

#' @export
getNoImgPerWell.default <- function(data) {
  stop("can only deal with Data (MatData/ImageData/WellData/PlateData) ",
       "objects, not with ", class(data))
}