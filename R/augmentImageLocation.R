#' Augment the image location within a well
#'
#' Given a PlateData/WellData/ImageData object, add the "Image.Group" feature
#' which encodes how many well borders are image borders and for every
#' "Location_Center_X/Y" feature, add "Location_Shifted_X/Y" feature for where
#' the object is located within the well.
#'
#' @param x The PlateData/WellData/ImageData object of interest.
#' 
#' @return An augmented version of the input object.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' aug <- augmentImageLocation(data)
#' 
#' @export
augmentImageLocation <- function(x) {
  UseMethod("augmentImageLocation", x)
}

#' @export
augmentImageLocation.PlateData <- function(x) {
  progress.bar <- getOption("singleCellFeatures.progressBars")
  if(progress.bar != "none") {
    message("augmenting image location features:")
  }
  x$data <- llply(x$data, augmentImageLocation,
                  .progress=getOption("singleCellFeatures.progressBars"))
  return(x)
}

#' @export
augmentImageLocation.WellData <- function(x) {
  x$data <- lapply(x$data, augmentImageLocation)
  return(x)
}

#' @export
augmentImageLocation.ImageData <- function(x) {
  n.img  <- getNoImgPerWell(x)
  img.no <- x$image.index
  if(img.no < 1 | img.no > n.img) stop("impossible image number")
  img.group <- c(2, 1, 2, 1, 0, 1, 2, 1, 2)
  if(n.img == 6) {
    img.group <- rep(img.group[1:3], 2)
  }
  img.group <- img.group[img.no]
  dim(img.group) <- c(1, 1)
  colnames(img.group) <- "Image.Group"
  x$data.vec$Image <- data.frame(cbind(x$data.vec$Image, img.group))

  x$data.mat <- lapply(x$data.mat, function(group) {
    if(!is.null(group)) {
      index.x <- grep("Location_Center_X", colnames(group))
      index.y <- grep("Location_Center_Y", colnames(group))
      if(length(index.x) != length(index.y)) {
        stop("expecting equal number of center features for x & y")
      }
      length.original <- ncol(group)
      length.new      <- length(index.x)
      if(nrow(group) == 0 & length.new > 0) {
        # in case no rows are present, add features to colnames anyway
        nmes <- c(colnames(group),
                  gsub("Center_X$", "Shifted_X", colnames(group)[index.x]),
                  gsub("Center_Y$", "Shifted_Y", colnames(group)[index.y]))
        dims <- c(0, length.original + 2 * length.new)
        dim(group)      <- dims
        colnames(group) <- nmes
        return(group)
      } else if(nrow(group) > 0 & length.new > 0) {
        feat.x <- group[,index.x, drop=FALSE]
        feat.y <- group[,index.y, drop=FALSE]
        shifts <- cbind(rep(c(0, 1402, 2804), 3),
                        rep(c(2100, 1050, 0), each=3))
        if(n.img == 6) {
          shifts <- shifts[4:9,]
        }
        shifts <- shifts[img.no,]
        shifted.x <- sapply(1:ncol(feat.x), function(ind, x, shift.x) {
          return(x[,ind] + shift.x)
        }, feat.x, shifts[1])
        shifted.y <- sapply(1:ncol(feat.y), function(ind, y, shift.y) {
          return(y[,ind] + shift.y)
        }, feat.y, shifts[2])
        dim(shifted.x) <- dim(feat.x)
        dim(shifted.y) <- dim(feat.y)
        res.shift <- cbind(shifted.x, shifted.y)
        colnames(res.shift) <- c(gsub("Center_X$", "Shifted_X",
                                      colnames(group)[index.x]),
                                 gsub("Center_Y$", "Shifted_Y",
                                      colnames(group)[index.y]))
        return(cbind(group, res.shift))
      } else {
        return(group)
      }
    } else return(NULL)
  })
  return(x)
}

#' @export
augmentImageLocation.default <- function(x) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}