#' Augment PlateData with B-scores
#'
#' Given a PlateData, aggregate the selected features within wells and
#' calculate their plate B-scores. 
#'
#' @param x         The PlateData object of interest.
#' @param features  A regular expression used for selecting the features to
#'                  include.
#' @param aggregate The function used to aggregate data within wells.
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
augmentBscore.PlateData <- function(x, features="intensity",
                                    aggregate="mean") {
  matched.feats <- grep(features, getFeatureNames(x), ignore.case=TRUE,
                        value=TRUE)
  if(length(matched.feats) == 0) stop("no features found.")
  aggr.fun <- get(aggregate, mode="function")

  progress.bar <- getOption("singleCellFeatures.progressBars")
  if(progress.bar != "none") {
    message("aggregating features to well level:")
  }
  dat.aggr <- llply(x$data, function(well, feat, fun) {
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
    medp <- medpolish(mat.smp, eps=1e-05, maxiter=200,
                      trace.iter=FALSE, na.rm=TRUE)
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