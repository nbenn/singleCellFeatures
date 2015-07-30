#' Normalize Datasets
#'
#' After running aggregate, bscore or MARS augmentation functions, this
#' function can be used for centering and scaling data.
#'
#' @param x        The Data object of interest.
#' @param select   A vector of regular expressions for selecting features to
#'                 normalize.
#' @param drop     A vector of regular expressions for dropping some of the
#'                 previously matched features.
#' @param values   The type of features to normalize.
#' @param center   NULL of the type of features to use for centering.
#' @param scale    NULL of the type of features to use for scaling
#' @param features For internal use only.
#' 
#' @return A Data object of input type with normalized features added.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' data <- augmentImageLocation(data)
#' data <- augmentCordinateFeatures(data, ellipse=1, facet=NULL,
#'                                  center.dist=FALSE, density=FALSE)
#' data <- augmentBscore(data)
#' data <- augmentMars(data)
#' 
#' data <- augmentAggregate(data, features="_MARSed$", level="plate",
#'                          func.aggr="median")
#' data <- augmentAggregate(data, features="_MARSed$", level="well",
#'                          neighbors=TRUE,  func.aggr="mad")
#' data <- normalizeData(data)
#'
#' @export
normalizeData <- function(x, ...) {
  UseMethod("normalizeData", x)
}

#' @export
normalizeData.PlateData <- function(x,
                                    select=c(".AreaShape_", ".Intensity_",
                                             ".Texture_"),
                                    drop=c("^Bacteria.", "^BlobBacteria."),
                                    values="MARSed$",
                                    center="MARSed_Aggreg_P_median$",
                                    scale="MARSed_Aggreg_N_mad$") {
  progress.bar <- getOption("singleCellFeatures.progressBars")
  features <- findFeatureHelper(x, select, drop, values, center, scale)
  if(progress.bar != "none") {
    message("normalizing plate:")
  }
  x$data <- plyr::llply(x$data, normalizeData, features, NULL, NULL, values,
                        center, scale, .progress=progress.bar)
  return(x)
}

#' @export
normalizeData.WellData <- function(x,
                                   features=NULL,
                                   select=c(".AreaShape_", ".Intensity_",
                                            ".Texture_"),
                                   drop=c("^Bacteria.", "^BlobBacteria."),
                                   values="MARSed$",
                                   center="MARSed_Aggreg_P_median$",
                                   scale="MARSed_Aggreg_N_mad$") {
  if(is.null(features)) {
    features <- findFeatureHelper(x, select, drop, values, center, scale)
  }
  x$data <- lapply(x$data, normalizeData, features, NULL, NULL, values, center,
                   scale)
  return(x)
}

#' @export
normalizeData.ImageData <- function(x,
                                    features=NULL,
                                    select=c(".AreaShape_", ".Intensity_",
                                             ".Texture_"),
                                    drop=c("^Bacteria.", "^BlobBacteria."),
                                    values="MARSed$",
                                    center="MARSed_Aggreg_P_median$",
                                    scale="MARSed_Aggreg_N_mad$") {
  if(is.null(features)) {
    features <- findFeatureHelper(x, select, drop, values, center, scale)
  }
  if(!is.list(features)) stop("features has to be a list.")
  if(!all(names(features) == c("names", "values", "center", "scale"))) {
    stop("features is incorrectly formatted.")
  }

  if(!is.null(features$center) & !is.null(features$scale)) {
    indices <- cbind(features$values, features$center, features$scale)
    res <- apply(indices, 1, function(i, dat) {
      if(length(x[[i[1]]][[i[2]]]) > 0 &
         length(x[[i[4]]][[i[5]]]) > 0 &
         length(x[[i[7]]][[i[8]]]) > 0) {
        return((x[[i[1]]][[i[2]]][,i[3]] -
                x[[i[4]]][[i[5]]][[i[6]]]) /
                x[[i[7]]][[i[8]]][[i[9]]])
      } else return(NULL)
    }, x)
  } else if(!is.null(features$center)) {
    indices <- cbind(features$values, features$center)
    res <- apply(indices, 1, function(i, dat) {
      if(length(x[[i[1]]][[i[2]]]) > 0 &
         length(x[[i[4]]][[i[5]]]) > 0) {
        return(x[[i[1]]][[i[2]]][,i[3]] -
               x[[i[4]]][[i[5]]][[i[6]]])
      } else return(NULL)
    }, x)
  } else if(!is.null(features$scale)) {
    indices <- cbind(features$values, features$scale)
    res <- apply(indices, 1, function(i, dat) {
      if(length(x[[i[1]]][[i[2]]]) > 0 &
         length(x[[i[4]]][[i[5]]]) > 0) {
        return(x[[i[1]]][[i[2]]][,i[3]] /
               x[[i[4]]][[i[5]]][[i[6]]])
      } else return(NULL)
    }, x)
  }
  ncol <- nrow(features$values)
  if(!is.null(res)) {
    nrow <- length(res) / ncol
    dim(res) <- c(nrow, ncol)
  } else {
    res <- numeric()
    dim(res) <- c(0, ncol)
  }
  colnames(res) <- paste0(features$names, "_Normed")
  x$data.mat$Norm <- res
  return(x)
}

#' @export
normalizeData.default <- function(x, ...) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}

findFeatureHelper <- function(data, select, drop, values, center, scale) {
  features <- unique(unlist(lapply(select, grep, getFeatureNames(data),
                                   value=TRUE)))
  drop.ext <- c(drop, "_Bsco", "_Aggreg_", "MARSed")
  drop.ind <- unique(unlist(lapply(drop.ext, grep, features)))
  if(length(drop.ind) > 0) features <- features[-drop.ind]

  all.val <- grep(values, getFeatureNames(data), value=TRUE)
  match.val <- unlist(lapply(paste0("^", features), function(feat, x) {
    return(any(grepl(feat, x)))
  }, all.val))
  match.val <- features[match.val]

  if(!is.null(center)) {
    all.cen <- grep(center, getFeatureNames(data), value=TRUE)
    match.cen <- unlist(lapply(paste0("^", features), function(feat, x) {
      return(any(grepl(feat, x)))
    }, all.cen))
    match.cen <- features[match.cen]
    intersection <- match.cen
  }
  if(!is.null(scale)) {
    all.sca <- grep(scale, getFeatureNames(data), value=TRUE)
    match.sca <- unlist(lapply(paste0("^", features), function(feat, x) {
      return(any(grepl(feat, x)))
    }, all.sca))
    match.sca <- features[match.sca]
    intersection <- match.sca
  }
  if(!is.null(center) & !is.null(scale)) {
    intersection <- intersect(match.cen, match.sca)
  }

  features <- intersect(match.val, intersection)
  if(length(features) == 0) stop("no features left.")
  res.val <- unlist(lapply(paste0("^", features), grep, all.val, value=TRUE))
  if(!is.null(center)) {
    res.cen <- unlist(lapply(paste0("^", features), grep, all.cen, value=TRUE))
  } else res.cen <- NULL
  if(!is.null(scale)) {
    res.sca <- unlist(lapply(paste0("^", features), grep, all.sca, value=TRUE))
  } else res.sca <- NULL

  if("PlateData" %in% class(data)) {
    struct <- getImageStructure(data$data$A1$data$img_11)
  } else if("WellData" %in% class(data)) {
    struct <- getImageStructure(data$data$img_11)
  } else if("ImageData" %in% class(data)) {
    struct <- getImageStructure(data)
  }

  loc.val <- lapply(res.val, getFeatureLocationInStructure, struct)
  loc.val <- do.call(rbind, loc.val)
  if(!is.null(center)) {
    loc.cen <- lapply(res.cen, getFeatureLocationInStructure, struct)
    loc.cen <- do.call(rbind, loc.cen)
  } else loc.cen <- NULL
  if(!is.null(scale)) {
    loc.sca <- lapply(res.sca, getFeatureLocationInStructure, struct)
    loc.sca <- do.call(rbind, loc.sca)
  } else loc.sca <- NULL

  return(list(names=features, values=loc.val, center=loc.cen, scale=loc.sca))
}

getImageStructure <- function(dat) {
  top.level <- names(dat)
  groups <- lapply(top.level, function(top, dat) {
    group.names <- names(dat[[top]])
    if(!is.null(group.names)) {
      feats <- lapply(group.names, function(group, dat) {
        res <- names(dat[[group]])
        if(is.null(res)) res <- colnames(dat[[group]])
        return(res)
      }, dat[[top]])
      feats <- feats[which(!sapply(feats, is.null))]
      if(length(feats) > 0) {
        names(feats) <- group.names
        group.names  <- feats
      }
    }
    return(group.names)
  }, dat)
  names(groups) <- top.level
  return(groups)
}

getFeatureLocationInStructure <- function(feature, struct) {
  feat <- lapply(struct, function(top, feat) {
    if(!is.null(names(top))) {
      lapply(top, function(group, feat) {
        if(any(feat %in% group)) return(match(feat, group))
        else return(NULL)
      }, feat)
    } else return(NULL)
  }, feature)
  group <- lapply(feat, function(top) {
    if(!is.null(names(top))) {
      group <- !(sapply(top, is.null))
      if(any(group)) {
        return(which(group))
      } else return(NULL)
    } else return(NULL)
  })
  top <- which(!sapply(group, is.null))
  group <- group[[top]]
  feat <- feat[[top]][[group]]
  if(struct[[top]][[group]][feat] != feature) {
    stop("could not find ", feature, " in struct.")
  }
  return(c(top, group, feat))
}