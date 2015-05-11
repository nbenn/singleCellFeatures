#' Extract data for images from Data objects
#' 
#' Given a Data object and a vector of image indices (either in 1:6 or in 1:9)
#' return only the data belonging to the selected images. If this extraction is
#' performed on well level, either a well with only a subset of images is
#' returned or a single ImageData object (or a list of ImageData) objects is/are
#' returned.
#'
#' @param x         The original Data object, holding the superset of images
#' @param images    A vector of integers, specifying image indices of the images
#'                  to keep
#' @param keep.well (Only at WellData level) A boolean indicating whether to keep
#'                  the WellData object intact or return ImageData object(s)
#'
#' @return Either a PlateData, WellData or ImageData object holding only the
#'         data corresponding to the specified subset of images.
#'
#' @examples
#' plate  <- PlateData(PlateLocation("J101-2C"))
#' img_11 <- extractImages(plate, 1)
#' 
#' @export
extractImages <- function(x, images, ...) {
  UseMethod("extractImages", x)
}
#' @export
extractImages.PlateData <- function(x, images) {
  x$data <- lapply(x$data, extractImages, images, TRUE)
  return(x)
}
#' @export
extractImages.WellData <- function(x, images, keep.well=TRUE) {
  n.img <- x$data[[1]]$image.total
  if(!n.img %in% c(6,9)) stop("number of images per well has to be 6 or 9.")
  if(is.null(images)) stop("images have to be 1:", n.img)
  images <- as.vector(images, mode="integer")
  if(any(images < 1) | any(images > n.img)) stop("images have to be 1:", n.img)
  images <- unique(images)
  if(!keep.well) {
    if(length(images) == 1) return(x$data[[images]])
    else return(x$data[images])
  } else {
    x$data <- x$data[images]
    return(x)
  }
}
#' @export
extractImages.default <- function(x, images) {
  stop("can only deal with WellData/PlateData objects.")
}