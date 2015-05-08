#' @export
extractImages <- function(x, images, keep.well=TRUE) {
  UseMethod("extractImages", x)
}
#' @export
extractImages.PlateData <- function(x, images, keep.well=TRUE) {
  if(!keep.well) warning("keep.well will be ignored")
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