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
#' aug  <- augmentNeighborAggregate(data)
#' 
#' @export
augmentNeighborAggregate <- function(x, ...) {
  UseMethod("augmentNeighborAggregate", x)
}

#' @export
augmentNeighborAggregate.PlateData <- function(x, features="intensity",
                                               aggregate="mad") {
  matched.feats <- grep(features, getFeatureNames(x), ignore.case=TRUE,
                        value=TRUE)
  matched.feats <- grep("_Bsco", matched.feats, value=TRUE, invert=TRUE)
  matched.feats <- grep("_NAggr_", matched.feats, value=TRUE, invert=TRUE)
  
  if(length(matched.feats) == 0) stop("no features found.")

  progress.bar <- getOption("singleCellFeatures.progressBars")
  if(progress.bar != "none") {
    message("aggregating features to well level:")
  }
  molten <- llply(x$data, function(well, feat) {
    res <- lapply(feat, function(feat, dat) {
      res <- lapply(dat, function(type, f) {
        res <- lapply(type, function(grp, f) {
          return(grp[[f]])
        }, f)
        not.null <- which(!sapply(res, is.null))[1]
        return(res[[not.null]])
      }, feat)
      not.null <- which(!sapply(res, is.null))[1]
      return(res[[not.null]])
    }, meltData(well))
    names(res) <- feat
    return(res)
  }, matched.feats, .progress=progress.bar)

  grid <- matrix(1:384, nrow=16, byrow=TRUE)
  grid <- cbind(grid[,3], grid, grid[,22])
  grid <- rbind(grid[3,], grid, grid[14,])
  col <- rep(rep(1:24, 16)) + 1
  row <- rep(1:16, each=24) + 1
  stencil <- cbind(grid[(col - 1) * 18 + row],     # center
                   grid[(col - 1) * 18 + row - 1], # north
                   grid[(col - 2) * 18 + row],     # west
                   grid[(col - 1) * 18 + row + 1], # south
                   grid[col * 18 + row])           # east

  res <- apply(stencil, 1, function(wells, feats, data, fun.name) {
    fun <- get(fun.name, mode="function")
    res <- lapply(feats, function(feat, wells, data, fun) {
      all <- sapply(wells, function(well, feat, data) {
        return(data[[well]][[feat]])
      }, feat, data)
      return(fun(unlist(all)))
    }, wells, data, fun)
    names(res) <- paste0(feats, "_NAggr_", fun.name)
    return(res)
  }, matched.feats, molten, aggregate)

  for(i in 1:384) {
    newdat <- data.frame(res[[i]])
    for(j in 1:length(x$data[[i]]$data)) {
      x$data[[i]]$data[[j]]$data.vec$NAggr <- newdat
    }
  }
  return(x)
}

#' @export
augmentNeighborAggregate.default <- function(x, ...) {
  stop("can only deal with PlateData objects.")
}