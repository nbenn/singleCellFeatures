#' Augment all features involving locations
#'
#' All features containing the strings "Location_Center_X" or
#' "Location_Center_Y" are extended to include their membership to an ellipse (
#' concentric w.r.t. image center), their membership to a tile and the distance
#' to the image center.
#'
#' @param x           The PlateData/WellData/ImageData object of interest.
#' @param ellipse     An integer corresponding to the number of ellipses.
#' @param facet       A vector of length 2 e.g. c(14, 10) corresponding to the
#'                    number of bins in x-direction (14) and y-direction (10).
#' @param center.dist NULL or a boolean whether to include distances.
#' 
#' @return The augmented version of the original object.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' aug <- augmentCordinateFeatures(data, 5, c(14, 10), TRUE)
#' 
#' @export
augmentCordinateFeatures <- function(x, ...) {
  UseMethod("augmentCordinateFeatures", x)
}

#' @export
augmentCordinateFeatures.PlateData <- function(x, ellipse=NULL, facet=NULL,
                                               center.dist=NULL) {
  output <- augmentCordinateFeatures(x$data[[1]], ellipse, facet, center.dist)
  x$data <- llply(x$data, function(dat, ell, fac, cen) {
    suppressMessages(res <- augmentCordinateFeatures(dat, ell, fac, cen))
    return(res)
  }, ellipse, facet, center.dist,
  .progress=getOption("singleCellFeatures.progressBars"))
  return(x)
}

#' @export
augmentCordinateFeatures.WellData <- function(x, ellipse=NULL, facet=NULL,
                                              center.dist=NULL) {
  output <- augmentCordinateFeatures(x$data[[1]], ellipse, facet, center.dist)
  suppressMessages(x$data <- lapply(x$data, augmentCordinateFeatures, ellipse,
                   facet, center.dist))
  return(x)
}

#' @export
augmentCordinateFeatures.ImageData <- function(x, ellipse=NULL, facet=NULL,
                                               center.dist=NULL) {
  if(is.null(ellipse) & is.null(facet) & is.null(center.dist)) {
    warning("nothing to do.")
    return(x)
  }
  image.dimensions <- c(1392, 1040)
  image.center     <- image.dimensions / 2
  if(!is.null(ellipse)) {
    ellipse <- as.integer(ellipse)
    if(ellipse < 1 | ellipse > 20) stop("expecting ellipse to be in 1:20")
    if(ellipse != 1) {
      rect <- image.dimensions[1] * image.dimensions[2]
      step <- 1:(ellipse - 1) * rect / ellipse
      message("ellipse sizes: ", paste(step, collapse=", "))
      elli <- sapply(step, function(x) {
        b <- sqrt(x / pi * image.dimensions[2] / image.dimensions[1])
        a <- image.dimensions[1] * b / image.dimensions[2]
        return(c(a, b))
      })
    } else {
      elli <- image.center - 100
      message("using a single ellipse, 100px dist from borders.")
    }
    x$data.mat <- lapply(x$data.mat, function(group) {
      if(!is.null(group)) {
        index.x <- grep("Location_Center_X", colnames(group))
        index.y <- grep("Location_Center_Y", colnames(group))
        if(length(index.x) != length(index.y)) {
          stop("expecting equal number of center features for x & y")
        }
        length.original <- ncol(group)
        length.new      <- length(index.x)
        if(length.new > 0 & nrow(group) > 0) {
          feat.x <- group[,index.x, drop=FALSE]
          feat.y <- group[,index.y, drop=FALSE]
          if (ellipse == 1) {
            # only one ellipse
            res.ell <- sapply(1:ncol(feat.x), function(ind, x, y, cent) {
              sum <- (x[,ind] - cent[1])^2 / (elli[1])^2 +
                     (y[,ind] - cent[2])^2 / (elli[2])^2
              return(1 - (sum <= 1))
            }, feat.x, feat.y, image.center)
          } else {
            # multiple ellipses
            res.ell <- sapply(1:ncol(feat.x),
                              function(ind, x, y, cent, n.ell) {
              sum <- apply(elli, 2, function(e, x, y, cent) {
                return((x[,ind] - cent[1])^2 / (e[1])^2 +
                       (y[,ind] - cent[2])^2 / (e[2])^2)
              }, x, y, cent)
              if(is.null(nrow(sum))) {
                subtr <- sum(sapply(sum, function(row) row <= 1))
              } else {
                subtr <- rowSums(apply(sum, 2, function(row) row <= 1))
              }
              res <- n.ell - subtr
              return(res)
            }, feat.x, feat.y, image.center, ellipse)
          }
          dim(res.ell) <- dim(feat.x)
          res <- cbind(group, res.ell)
          range <- (length.original + 1):(length.original + length.new)
          colnames(res)[range] <- gsub("Center_X$", "In_Ellipse",
                                       colnames(group)[index.x])
          return(res)
        } else {
          return(group)
        }
      } else return(NULL)
    })
  }
  if(!is.null(facet)) {
    facet <- as.integer(facet)
    if(!(is.vector(facet, mode="integer") & length(facet) == 2)) {
      stop("expecting an integer vector of length 2 for facet.")
    }
    facet.size <- image.dimensions / facet
    message("facet size: ", round(facet.size[1]), " (x), ",
            round(facet.size[2]), " (y)")
    x$data.mat <- lapply(x$data.mat, function(group) {
      if(!is.null(group)) {
        index.x <- grep("Location_Center_X", colnames(group))
        index.y <- grep("Location_Center_Y", colnames(group))
        if(length(index.x) != length(index.y)) {
          stop("expecting equal number of center features for x & y")
        }
        length.original <- ncol(group)
        length.new      <- length(index.x)
        if(length.new > 0 & nrow(group) > 0) {
          feat.x <- group[,index.x, drop=FALSE]
          feat.y <- group[,index.y, drop=FALSE]
          facet.x <- sapply(1:ncol(feat.x), function(ind, x, fac) {
            return(ceiling(x[,ind] / fac[1]))
          }, feat.x, facet.size)
          facet.y <- sapply(1:ncol(feat.y), function(ind, y, fac) {
            return(ceiling(y[,ind] / fac[2]))
          }, feat.y, facet.size)
          dim(facet.x) <- dim(feat.x)
          dim(facet.y) <- dim(feat.y)
          res.fac <- cbind(facet.x, facet.y)
          colnames(res.fac) <- c(gsub("Center_X$", "Facet_X",
                                      colnames(group)[index.x]),
                                 gsub("Center_Y$", "Facet_Y",
                                      colnames(group)[index.y]))
          return(cbind(group, res.fac))
        } else {
          return(group)
        }
      } else return(NULL)
    })
  }
  if(!is.null(center.dist)) {
    if(!is.logical(center.dist)) {
      stop("expecting a boolean for center.dist.")
    }
    if(center.dist) {
      message("including distance to image center.")
      x$data.mat <- lapply(x$data.mat, function(group) {
        if(!is.null(group)) {
          index.x <- grep("Location_Center_X", colnames(group))
          index.y <- grep("Location_Center_Y", colnames(group))
          if(length(index.x) != length(index.y)) {
            stop("expecting equal number of center features for x & y")
          }
          length.original <- ncol(group)
          length.new      <- length(index.x)
          if(length.new > 0 & nrow(group) > 0) {
            feat.x <- group[,index.x, drop=FALSE]
            feat.y <- group[,index.y, drop=FALSE]
            dist <- sapply(1:ncol(feat.x), function(ind, x, y, cent) {
              return(sqrt((x[,ind] - cent[1])^2 +
                          (y[,ind] - cent[2])^2))
            }, feat.x, feat.y, image.center)
            dim(dist) <- dim(feat.x)
            res <- cbind(group, dist)
            range <- (length.original + 1):(length.original + length.new)
            colnames(res)[range] <- gsub("Center_X$", "Dist_Center",
                                         colnames(group)[index.x])
            return(res)
          } else {
            return(group)
          }
        } else return(NULL)
      })
    }
  }
  return(x)
}

#' @export
augmentCordinateFeatures.default <- function(x, ...) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}