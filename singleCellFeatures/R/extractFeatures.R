#' Extract features from Data objects
#' 
#' Given a Data object and either a set of regular expressions to select or
#' discard features or a set of features, return the same object type as the
#' original, but only with a subset of features.
#'
#' @param x        The original Data object, containing the superset of features
#' @param select   A vector of strings (regular expressions), specifying
#'                 features to keep
#' @param drop     A vector of strings (regular expressions), specifying
#'                 features to drop
#' @param features For performance reasons, the complete set of features to keep
#'                 can be determined in advance and specified here (in this
#'                 case, arguments select/drop will be ignored if not NULL)
#'
#' @return A Data object of the same type as x, but only storing a subset of all
#'         features.
#'
#' @examples
#' features <- c("^Cells.Location_Center_X$",
#'               "^Cells.Location_Center_Y$")
#' data <- PlateData(PlateLocation("J101-2C"))
#' data <- extractFeatures(data, features)
#' 
#' @export
extractFeatures <- function(x, select=NULL, drop=NULL, features=NULL) {
  UseMethod("extractFeatures", x)
}

#' @export
extractFeatures.MatData <- function(x, select=NULL, drop=NULL, features=NULL) {
  if(is.null(features)) {
    features <- extractFeaturesMatchingHelper(x, select, drop)
  } else if(!is.null(select) | !is.null(drop)) {
    warning("select/drop arguments will be ignored")
  }
  x$data <- x$data[which(getFeatureNames(x) %in% features)]
  return(x)
}

#' @export
extractFeatures.PlateData <- function(x, select=NULL, drop=NULL,
                                      features=NULL) {
  if(is.null(features)) {
    features <- extractFeaturesMatchingHelper(x, select, drop)
    select <- NULL
    drop <- NULL
  } else if(!is.null(select) | !is.null(drop)) {
    warning("select/drop arguments will be ignored")
    select <- NULL
    drop <- NULL
  }
  x$data <- lapply(x$data, extractFeatures, select, drop, features)
  return(x)
}

#' @export
extractFeatures.WellData <- function(x, select=NULL, drop=NULL, features=NULL) {
  if(is.null(features)) {
    features <- extractFeaturesMatchingHelper(x, select, drop)
    select <- NULL
    drop <- NULL
  } else if(!is.null(select) | !is.null(drop)) {
    warning("select/drop arguments will be ignored")
    select <- NULL
    drop <- NULL
  }
  x$data <- lapply(x$data, extractFeatures, select, drop, features)
  return(x)
}

#' @export
extractFeatures.ImageData <- function(x, select=NULL, drop=NULL,
                                      features=NULL) {
  if(is.null(features)) {
    features <- extractFeaturesMatchingHelper(x, select, drop)
  } else if(!is.null(select) | !is.null(drop)) {
    warning("select/drop arguments will be ignored")
  }
  if(length(features) == 0 | !is.vector(features, mode = "character")) {
    message("no features are removed (zero length or not a character vector).")
  } else {
    x$data.vec <- lapply(
      x$data.vec,
      function(object, feat) {
        match <- names(object) %in% features
        if(sum(match) == 0) return(NULL)
        else return(object[match])
      },
      features
    )
    mat <- unlist(lapply(x$data.mat, function(group) {
      res <- colnames(group)
      if (is.null(res)) res <- names(group)
      return(res)
    }))
    x$data.mat <- lapply(
      x$data.mat,
      function(object, feat) {
        nms <- colnames(object)
        if(!is.null(nms)) {
          match <- nms %in% features
          if(sum(match) == 0) return(NULL)
          else return(object[,match])
        } else {
          nms <- names(object)
          match <- nms %in% features
          if(sum(match) == 0) return(NULL)
          else return(object[match])
        }
      },
      features
    )
    x$data.lst <- lapply(
      x$data.lst,
      function(object, feat) {
        match <- names(object) %in% features
        if(sum(match) == 0) return(NULL)
        else return(object[match])
      },
      features
    )
  }
  return(x)
}

#' @export
extractFeatures.default <- function(x, select=NULL, drop=NULL, features=NULL) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}

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
  if(n.feat - length(features) > 0)
    message("removing ", n.feat - length(keep), " unmatched features.")
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
  if(n.feat2 - length(features) > 0)
    message("removing ", length(remove), " matched features.")
  return(features)
}