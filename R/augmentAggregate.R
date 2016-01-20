#' Augment PlateData with Aggregates
#'
#' Given a PlateData, aggregate the selected features within wells or the
#' whole plate with the supplied function (examples include mean, median, sum,
#' var, mad, etc.)
#'
#' @param x         The PlateData object of interest.
#' @param level     A string, either "well" or "plate", indicaing the level
#'                  at which to aggregate.
#' @param neighbors Logical, whether to include direct neighbors for well level
#'                  aggregation.
#' @param features  A (list) of regular expressions used for selecting the
#'                  features to include.
#' @param drop      A (list) of regular expressions used for removing features
#'                  from the previously compiled list.
#' @param aggregate The function used to aggregate data within wells.
#' 
#' @return An augmented version of the input object.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' aug  <- augmentAggregate(data)
#' 
#' @export
augmentAggregate <- function(x, ...) {
  UseMethod("augmentAggregate", x)
}

#' @export
augmentAggregate.PlateData <- function(x, level="well", neighbors=FALSE,
                                       features=c(".AreaShape_", ".Intensity_",
                                                  ".Texture_"),
                                       drop=c("^Bacteria.", "^BlobBacteria."),
                                       func.aggr="median") {

  if(!level %in% c("well", "plate")) {
    stop("expecting level to be either \"well\" or \"plate\".")
  }
  if(level != "well" & neighbors) {
    warning("ignoring parameter \"neighbors\".")
  }

  if(level == "well" & neighbors) newname <- paste0("_Aggreg_N_", func.aggr)
  else if(level == "well") newname <- paste0("_Aggreg_W_", func.aggr)
  else newname <- paste0("_Aggreg_P_", func.aggr)

  fun <- get(func.aggr, mode="function")
  na.rm <- "na.rm" %in% names(formals(fun))

  matched.feats <- unique(unlist(lapply(features, grep, getFeatureNames(x),
                                        value=TRUE)))
  drop.ext      <- c(drop, "_Bsco", "_Aggreg_")
  drop.ind      <- unique(unlist(lapply(drop.ext, grep, matched.feats)))
  if(length(drop.ind) > 0) matched.feats <- matched.feats[-drop.ind]

  if(length(matched.feats) == 0) stop("no features found.")
  molten <- extractFeatures(x, features=matched.feats)

  progress.bar <- getOption("singleCellFeatures.progressBars")

  if(level == "well") {
    if(progress.bar != "none") {
      message("aggregating features to well level:")
    }

    molten <- llply(molten$data, function(well, feat) {
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

    if(neighbors) {
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
      if(na.rm) {
        res <- apply(stencil, 1, function(wells, feats, data, fun, name) {
          res <- lapply(feats, function(feat, wells, data, fun) {
            all <- sapply(wells, function(well, feat, data) {
              return(data[[well]][[feat]])
            }, feat, data)
            return(fun(unlist(all), na.rm=TRUE))
          }, wells, data, fun)
          names(res) <- paste0(feats, newname)
          return(res)
        }, matched.feats, molten, fun, newname)
      } else {
        res <- apply(stencil, 1, function(wells, feats, data, fun, name) {
          res <- lapply(feats, function(feat, wells, data, fun) {
            all <- sapply(wells, function(well, feat, data) {
              return(data[[well]][[feat]])
            }, feat, data)
            return(fun(unlist(all)))
          }, wells, data, fun)
          names(res) <- paste0(feats, newname)
          return(res)
        }, matched.feats, molten, fun, newname)
      }
    } else {
      if(na.rm) {
        res <- lapply(molten, function(dat, fun) {
          res <- lapply(dat, function(feat, fun) {
            return(fun(feat, na.rm=TRUE))
          }, fun)
          names(res) <- paste0(names(dat), newname)
          return(res)
        }, fun)
      } else {
        res <- lapply(molten, function(dat, fun) {
          res <- lapply(dat, function(feat, fun) {
            return(fun(feat))
          }, fun)
          names(res) <- paste0(names(dat), newname)
          return(res)
        }, fun)
      }
    }

    for(i in 1:length(x$data)) {
      newdat <- data.frame(res[[i]])
      for(j in 1:length(x$data[[i]]$data)) {
        old <- x$data[[i]]$data[[j]]$data.vec$Aggregate
        if(is.null(old)) {
          x$data[[i]]$data[[j]]$data.vec$Aggregate <- newdat
        } else {
          x$data[[i]]$data[[j]]$data.vec$Aggregate <- cbind(old, newdat)
        }
      }
    }
  } else if(level == "plate") {
    molten <- meltData(molten)
    res <- sapply(matched.feats, function(feat, dat, fun) {
      res <- lapply(dat, function(type, feat) {
        res <- lapply(type, function(grp, feat) {
          return(grp[[feat]])
        }, feat)
        not.null <- which(!sapply(res, is.null))[1]
        return(res[[not.null]])
      }, feat)
      not.null <- which(!sapply(res, is.null))[1]
      if(na.rm) {
        return(fun(res[[not.null]], na.rm=TRUE))
      } else {
        return(fun(res[[not.null]]))
      }
    }, molten, fun)
    resnames <- paste0(names(res), newname)
    dim(res) <- c(1, length(res))
    res <- data.frame(res)
    names(res) <- resnames
    for(i in 1:length(x$data)) {
      for(j in 1:length(x$data[[i]]$data)) {
        old <- x$data[[i]]$data[[j]]$data.vec$Aggregate
        if(is.null(old)) {
          x$data[[i]]$data[[j]]$data.vec$Aggregate <- res
        } else {
          x$data[[i]]$data[[j]]$data.vec$Aggregate <- cbind(old, res)
        }
      }
    }
  }
  return(x)
}

#' @export
augmentAggregate.WellData <- function(x, features=c(".AreaShape_",
                                                    ".Intensity_",
                                                    ".Texture_"),
                                      drop=c("^Bacteria.", "^BlobBacteria."),
                                      func.aggr="median") {

  fun <- get(func.aggr, mode="function")
  matched.feats <- unique(unlist(lapply(features, grep, getFeatureNames(x),
                                        value=TRUE)))
  drop.ext      <- c(drop, "_Bsco", "_Aggreg_")
  drop.ind      <- unique(unlist(lapply(drop.ext, grep, matched.feats)))
  if(length(drop.ind) > 0) matched.feats <- matched.feats[-drop.ind]

  if(length(matched.feats) == 0) stop("no features found.")
  molten <- extractFeatures(x, features=matched.feats)

  molten <- lapply(matched.feats, function(feat, dat) {
    res <- lapply(dat, function(type, f) {
      res <- lapply(type, function(grp, f) {
        return(grp[[f]])
      }, f)
      not.null <- which(!sapply(res, is.null))[1]
      return(res[[not.null]])
    }, feat)
    not.null <- which(!sapply(res, is.null))[1]
    return(res[[not.null]])
  }, meltData(molten))
  names(molten) <- matched.feats

  if("na.rm" %in% names(formals(fun))) {
    res <- lapply(molten, function(feat, fun) {
      return(fun(feat, na.rm=TRUE))
    }, fun)
  } else {
    res <- lapply(molten, function(feat, fun) {
      return(fun(feat))
    }, fun)
  }
  names(res) <- paste0(names(molten), paste0("_Aggreg_W_", func.aggr))
  res <- data.frame(res)

  for(i in 1:length(x$data)) {
    old <- x$data[[i]]$data.vec$Aggregate
    if(is.null(old)) {
      x$data[[i]]$data.vec$Aggregate <- res
    } else {
      x$data[[i]]$data.vec$Aggregate <- cbind(old, res)
    }
  }
  return(x)
}

#' @export
augmentAggregate.default <- function(x, ...) {
  stop("can only deal with PlateData objects.")
}