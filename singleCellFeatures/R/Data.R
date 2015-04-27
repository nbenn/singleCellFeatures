#' Constructor for PlateData objects
#' 
#' @export
PlateData <- function(plate, data) {
  result <- list(meta = PlateMetadata(plate),
                 data = data)
  return(structure(result, class = c("PlateData", "Data")))
}

#' Constructor for WellData objects
#' 
#' @export
WellData <- function(well, data) {
  result <- list(meta = WellMetadata(well),
                 data = data)
  return(structure(result, class = c("WellData", "Data")))
}

#' Constructor for ImageData objects
#' 
#' @export
ImageData <- function(barcode, index, vec, mat, lst) {
  ind <- getWellIndex2D(index, image=TRUE)
  result <- list(plate       = barcode,
                 well.index  = index,
                 image.index = ind$img,
                 well.row    = ind$row,
                 well.col    = ind$col,
                 data.vec    = vec,
                 data.mat    = mat,
                 data.lst    = lst)
  return(structure(result, class = c("ImageData", "Data")))
}

#' @export
getBarcode.PlateData <- function(x) {
  return(x$meta$barcode)
}
#' @export
getBarcode.WellData <- function(x) {
  return(x$meta$barcode)
}
#' @export
getBarcode.ImageData <- function(x) {
  return(x$plate)
}

#' @export
getFeatureNames <- function(x) {
  UseMethod("getFeatureNames", x)
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

#' @export
extractFeatures <- function(x, select, drop, features=NULL) {
  UseMethod("extractFeatures", x, select, drop, features)
}
#' @export
extractFeatures.PlateData <- function(x, select, drop, features=NULL) {
  if(is.null(features)) {
    features <- extractFeaturesMatchingHelper(x, select, drop)
  } else if(!is.null(select) | !is.null(drop)) {
    warning("select/drop arguments will be ignored")
    select <- NULL
    drop <- NULL
  }
  x$data <- lapply(x$data, extractFeatures, select, drop, features)
  return(x)
}
#' @export
extractFeatures.WellData <- function(x, select, drop, features=NULL) {
  if(is.null(features)) {
    features <- extractFeaturesMatchingHelper(x, select, drop)
  } else if(!is.null(select) | !is.null(drop)) {
    warning("select/drop arguments will be ignored")
    select <- NULL
    drop <- NULL
  }
  x$data <- lapply(x$data, extractFeatures, select, drop, features)
  return(x)
}
#' @export
extractFeatures.ImageData <- function(x, select, drop, features=NULL) {
  if(is.null(features)) {
    features <- extractFeaturesMatchingHelper(x, select, drop)
  } else if(!is.null(select) | !is.null(drop)) {
    warning("select/drop arguments will be ignored")
    select <- NULL
    drop <- NULL
  }
  if(length(features) == 0 | !is.vector(features, mode = "character")) {
    message("no features are removed (zero length or not a character vector).")
  } else {
    x$data.vec <- lapply(x$data.vec, function(object, feat) {
      match <- names(object) %in% features
      if(sum(match) == 0) return(NULL)
      else return(object[match])
    }, features)
    x$data.mat <- lapply(x$data.mat, function(object, feat) {
      match <- colnames(object) %in% features
      if(sum(match) == 0) return(NULL)
      else return(object[,match])
    }, features)
    x$data.lst <- lapply(x$data.lst, function(object, feat) {
      match <- names(object) %in% features
      if(sum(match) == 0) return(NULL)
      else return(object[match])
    }, features)
  }
  return(x)
}
#' @export
extractFeatures.default <- function(x, select, drop, features=NULL) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}
#' @export
extractFeaturesMatchingHelper <- function(x, select, drop) {
  features <- getFeatureNames(x)
  n.feat <- length(features)
  # if features specified, drop other files
  if (!is.null(select)) {
    # for each feature entry, get all partial matches
    keep <- lapply(select, grep, x=features)
    # remove duplicates from match index
    keep <- unique(unlist(keep))
    # drop unmatched features
    features <- features[keep]    
  }
  if(n.feat-length(features) > 0) 
    message("removing ", n.feat-length(keep), " unmatched features.")
  n.feat2 <- length(features)
  # if drop specified, drop features
  if (!is.null(drop)) {
    # for each drop entry, get all partial matches
    remove <- lapply(drop, grep, x=features)
    # remove duplicates from match index
    remove <- unique(unlist(remove))
    # drop unmatched features
    features <- features[-remove]    
  }
  if(n.feat2-length(features) > 0) 
    message("removing ", length(remove), " matched features.")
  return(features)
}

#' @export
extractImages <- function(x, images) {
  UseMethod("extractImages", x, images)
}
#' @export
extractImages.PlateData <- function(x, images) {
  x$data <- lapply(x$data, extractImages, images)
  return(x)
}
#' @export
extractImages.WellData <- function(x, images) {
  if(is.null(images)) stop("images have to be 1:9")
  images <- as.vector(images, mode="integer")
  if(any(images < 1) | any(images > 9)) stop("images have to be 1:9")
  x$data <- x$data[unique(images)]
  return(x)
}
#' @export
extractImages.default <- function(x, images) {
  stop("can only deal with WellData/PlateData objects.")
}

#' @export
extractWells <- function(x, wells, keep.plate=TRUE) {
  UseMethod("extractWells", x)
}
#' @export
extractWells.PlateData <- function(x, wells, keep.plate=TRUE) {
  if(any(class(wells) == "WellLocation")) wells <- list(wells)
  if(!all(sapply(wells, function(well) {
    return(any(class(well) == "WellLocation"))
  }))) stop("can only work with a list of WellLocation objects")
  barcodes <- sapply(wells, getBarcode)
  if(length(unique(barcodes)) != 1)
    stop("can only deal with WellLocations on the same plate")
  if(unique(barcodes) != getBarcode(x))
    stop("WellLocations have to be on the same plate as the data")
  well.names <- sapply(wells, getWellName)
  if(!keep.plate) return(x$data[well.names])
  else {
    x$data <- x$data[well.names]
    return(x)
  }
}
#' @export
extractWells.default <- function(x, wells, keep.plate=TRUE) {
  stop("can only deal with PlateData objects.")
}