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
#' plate.df <- meltData(plate)
#' 
#' @export
meltData <- function(x) {
  UseMethod("meltData", x)
}

#' @export
meltData.MatData <- function(x) {
  return(names(x$data))
}

#' @export
meltData.PlateData <- function(x) {
  groups <- names(x$data[[1]]$data[[1]]$data.mat)
  result <- unlist(lapply(x$data, meltData), recursive=FALSE)
  result <- lapply(groups, function(group, data) {
    regexp <- paste0("^[A-P]([1-9]|1[0-9]|2[0-4])\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    res <- do.call(rbind, res)
    return(res)
  }, result)
  names(result) <- groups
  return(result)
}

#' @export
meltData.WellData <- function(x) {
  groups <- names(x$data[[1]]$data.mat)
  result <- unlist(lapply(x$data, meltData), recursive=FALSE)
  result <- lapply(groups, function(group, data) {
    regexp <- paste0("^img_[1-3]{2}\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    res <- do.call(rbind, res)
    return(res)
  }, result)
  names(result) <- groups
  return(result)
}

#' @export
meltData.ImageData <- function(x) {
  result <- lapply(x$data.mat, function(group) {
    n.rows <- nrow(group)
    if(n.rows > 0) {
      image.group <- x$data.vec$Image[["Image.Group"]]

      image.ind      <- x$image.index
      image.ind      <- rep(image.ind, n.rows)
      dim(image.ind) <- c(n.rows, 1)
      well.ind      <- getWellIndex1D(x$well.row, x$well.col, NULL,
                                      x$image.total)
      well.ind      <- rep(well.ind, n.rows)
      dim(well.ind) <- c(n.rows, 1)
      plate      <- x$plate
      plate      <- rep(plate, n.rows)
      dim(plate) <- c(n.rows, 1)

      if(!is.null(image.group)) {
        image.group      <- rep(image.group, n.rows)
        dim(image.group) <- c(n.rows, 1)
        info      <- data.frame(image.ind, well.ind, plate, image.group,
                                stringsAsFactors=FALSE)
        colnames(info) <- c("Image.Index", "Well.Index", "Plate.Barcode",
                            "Image.Group")
      } else {
        info      <- data.frame(image.ind, well.ind, plate,
                                stringsAsFactors=FALSE)
        colnames(info) <- c("Image.Index", "Well.Index", "Plate.Barcode")
      }

      return(data.frame(group, info))
    } else {
      return(NULL)
    }
  })
  return(result)
}

#' @export
meltData.default <- function(x) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}