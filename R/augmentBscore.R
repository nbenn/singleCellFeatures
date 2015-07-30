#' Augment PlateData with B-scores
#'
#' Given a PlateData, aggregate the selected features within wells and
#' calculate their plate B-scores. 
#'
#' @param x         The PlateData object of interest.
#' @param features  A (list) of regular expressions used for selecting the
#'                  features to include.
#' @param drop      A (list) of regular expressions used for removing features
#'                  from the previously compiled list.
#' @param func.aggr The function used to aggregate data within wells.
#' 
#' @return An augmented version of the input object.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' aug  <- augmentBscore(data)
#' 
#' @export
augmentBscore <- function(x, ...) {
  UseMethod("augmentBscore", x)
}

#' @export
augmentBscore.PlateData <- function(x,
                                    features=c(".AreaShape_", ".Intensity_",
                                               ".Texture_"),
                                    drop=c("^Bacteria.", "^BlobBacteria."),
                                    func.aggr="mean") {

  matched.feats <- unique(unlist(lapply(features, grep, getFeatureNames(x),
                                        value=TRUE)))
  drop.ext      <- c(drop, "_Bsco", "_Aggreg_", "MARSed")
  drop.ind      <- unique(unlist(lapply(drop.ext, grep, matched.feats)))
  if(length(drop.ind) > 0) matched.feats <- matched.feats[-drop.ind]

  if(length(matched.feats) == 0) stop("no features found.")
  aggr.fun <- get(func.aggr, mode="function")

  progress.bar <- getOption("singleCellFeatures.progressBars")
  if(progress.bar != "none") {
    message("aggregating features to well level:")
  }

  dat.aggr <- extractFeatures(x, features=matched.feats)
  dat.aggr <- llply(dat.aggr$data, function(well, feat, fun) {
    res <- lapply(feat, function(feat, dat, fun) {
      res <- lapply(dat, function(type, f) {
        res <- lapply(type, function(grp, f) {
          return(grp[[f]])
        }, f)
        not.null <- which(!sapply(res, is.null))[1]
        return(res[[not.null]])
      }, feat)
      not.null <- which(!sapply(res, is.null))[1]
      res <- res[[not.null]]
      if(is.null(res)) return(NA)
      if(is.numeric(res)) return(fun(res))
      else return(unique(res))
    }, meltData(well), fun)
    names(res) <- feat
    return(res)
  }, matched.feats, aggr.fun, .progress=progress.bar)
  res <- lapply(matched.feats, function(feat, data) {
    dat <- sapply(data, function(dat, fet) {
      return(dat[[fet]])
    }, feat)
    mat.all <- mat.smp <- matrix(dat, nrow=16, ncol=24, byrow=TRUE)
    medp <- tryCatch({
      medpolish(mat.smp, eps=1e-05, maxiter=200, trace.iter=FALSE, na.rm=TRUE)
    }, warning = function(w) {
      message("  warning for ", feat)
      medpolish(mat.smp, eps=1e-05, maxiter=200, trace.iter=FALSE, na.rm=TRUE)
    })
    medp$row[is.na(medp$row)] <- 0
    medp$col[is.na(medp$col)] <- 0
    all <- rep(medp$overall, 384)
    row <- rep(medp$row, each=24)
    col <- rep(medp$col, 16)
    res <- data.frame(cbind(all, row, col))
    names(res) <- paste0(feat, c("_BscoAll", "_BscoRow", "_BscoCol"))
    return(res)
  }, dat.aggr)
  res <- do.call(cbind, res)
  for(i in 1:384) {
    add <- res[i,]
    for(j in 1:length(x$data[[i]]$data)) {
      x$data[[i]]$data[[j]]$data.vec$Bscore <- add
    }
  }
  return(x)
}

#' @export
augmentBscore.default <- function(x, ...) {
  stop("can only deal with PlateData objects.")
}