#' Augment PlateData with MARS residuals
#'
#' Given a PlateData, normalize the selected features with a MARS model
#' consisting of terms for membership to image/well centered ellipses and
#' medpolish estimated row/column effects.
#'
#' @param x         The PlateData object of interest.
#' @param features  A (list) of regular expressions used for selecting the
#'                  features to include.
#' @param drop      A (list) of regular expressions used for removing features
#'                  from the previously compiled list.
#' 
#' @return An augmented version of the input object.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' aug  <- augmentMars(data)
#' 
#' @export
augmentMars <- function(x, ...) {
  UseMethod("augmentMars", x)
}

#' @export
augmentMars.PlateData <- function(x,
                                  features=c(".AreaShape_", ".Intensity_",
                                             ".Texture_"),
                                  drop=c("^Bacteria.", "^BlobBacteria.")) {

  progress.bar <- getOption("singleCellFeatures.progressBars")

  matched.feats <- unique(unlist(lapply(features, grep, getFeatureNames(x),
                                        value=TRUE)))
  drop.ext      <- c(drop, "_Bsco", "_Aggreg_", "MARSed")
  drop.ind      <- unique(unlist(lapply(drop.ext, grep, matched.feats)))
  if(length(drop.ind) > 0) matched.feats <- matched.feats[-drop.ind]

  if(length(matched.feats) == 0) stop("no features found.")

  aug.loc <- c("Nuclei.Location_In_Ellipse_Well",
               "Nuclei.Location_In_Ellipse_Image")
  aug.bsc <- c(paste0(matched.feats, "_BscoAll"),
               paste0(matched.feats, "_BscoRow"),
               paste0(matched.feats, "_BscoCol"))
  if(!all(aug.loc %in% getFeatureNames(x))) {
    stop("could not find all needed location features. Please run\n",
         "augmentImageLocation/augmentCordinateFeatures first.")
  }
  if(!all(aug.bsc %in% getFeatureNames(x))) {
    stop("could not find all needed bscore features. Please run\n",
         "augmentBscore first.")
  }

  molten <- extractFeatures(x, features=c(matched.feats, aug.loc, aug.bsc))
  molten <- meltData(molten)
  molten <- moveFeatures(molten, from="_Bsco", to="Cells")

  na.col <- plyr::alply(molten, 2, function(col) return(any(is.na(col))))
  na.col <- unlist(na.col)
  if(sum(na.col) > 0) {
    message("removing ", sum(na.col), " columns due to many NAs:")
    plyr::l_ply(names(molten)[na.col], function(name) message("  ", name))
    matched.feats <- matched.feats[-match(names(molten)[na.col],
                                          matched.feats, nomatch=0)]
    molten <- molten[, !na.col]
  }

  zero.var <- sapply(matched.feats, function(feat, dat) {
    return(var(dat[[feat]]) != 0)
  }, molten)
  message("dropping ", sum(!zero.var), " features due to zero variance:")
  if(sum(!zero.var) > 0) {
    plyr::l_ply(matched.feats[!zero.var], function(x) message("  ", x))
    matched.feats <- matched.feats[zero.var]
  }

  message("normalizing ", length(matched.feats), " features:")
  newdat <- plyr::llply(matched.feats, function(feat, data) {
    bscore   <- paste0(feat, c("_BscoRow", "_BscoCol", "_BscoAll"))
    form <- formula(paste0(feat, " ~ ",
                           "Nuclei.Location_In_Ellipse_Well + ",
                           "Nuclei.Location_In_Ellipse_Image + ",
                           paste0(bscore, collapse=" + ")))
    model <- earth::earth(formula=form, data=data)
    return(model$residuals)
  }, molten, .progress=progress.bar)

  newdat <- do.call(cbind, newdat)
  newdat <- as.data.frame(newdat)

  if(nrow(newdat) != nrow(molten)) {
    stop("expecting MARSed features to be of same length as orig.")
  }

  newdat  <- split(newdat, molten$Well.Name)
  img.ind <- split(molten$Image.Index, molten$Well.Name)

  for(well in sapply(x$data, getWellName)) {
    ncol <- length(matched.feats)
    if(!is.null(newdat[[well]])) {
      dat.well <- as.matrix(newdat[[well]])
      nrow <- length(dat.well) / ncol
    } else {
      dat.well <- numeric()
      nrow <- 0
    }
    dim(dat.well) <- c(nrow, ncol)
    colnames(dat.well) <- paste0(matched.feats, "_MARSed")
    index <- img.ind[[well]]
    for(img in 1:length(x$data[[well]]$data)) {
      if(!is.null(index)) {
        x$data[[well]]$data[[img]]$data.mat$MARS <- dat.well[index == img,,
                                                             drop=FALSE]
      } else {
        x$data[[well]]$data[[img]]$data.mat$MARS <- dat.well
      }
    }
  }
  return(x)
}

#' @export
augmentMars.default <- function(x, ...) {
  stop("can only deal with PlateData objects.")
}