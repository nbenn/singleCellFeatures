#' Rebuild cache files for certain Data objects
#'
#' Given an existing MatData cache file and some feature files that were not
#' imported correctly, update and save it to the corresponding cache location.
#'
#' @param x        MatData object
#' @param location An optional argument for the location of the new .mat files.
#'                 If not specified, the /Errors directory is checked and if it
#'                 does not exist, the missing features are detected and
#'                 downloaded.
#' @param names    An optional argument for the name of the newly added feature.
#'                 If not specified, the filename (without extension) will be
#'                 used. It is ignored when the missing features are detected
#'                 automatically.
#' 
#' @return The updated MatData object. The updated cache file is saved to its
#'         corresponding location.
#' 
#' @examples
#' dat <- MatData(PlateLocation("J110-2C"))
#' dat <- rebuildCache(dat)
#'
#' @export
rebuildCache <- function(x, ...) {
  UseMethod("rebuildCache", x)
}

#' @export
rebuildCache.MatData <- function(x, location=NULL, names=NULL) {
  loc  <- convertToPlateLocation(x)
  if (is.null(location)) {
    location <- paste0(getLocalPath(loc), "/", "Errors")
  }
  if (dir.exists(location)) {
    message("reading from ", location)
    data <- readMatFeatureHelper(location)
    if (!is.null(names)) names(data) <- names
  } else {
    suppressWarnings(missing <- checkCompletenessFeature(x)$missing)
    if (length(missing) == 0) {
      stop("no missing features detected.")
    }
    if (length(missing) > 150) {
      # all missing features are merged into one huge or-linked regexp
      stop("consider downloading the whole dataset anew.")
    }
    message("fetching missing features:\n  ", paste(missing, collapse="\n  "))
    if ("Cells.Infection_IsInfected" %in% missing) {
      dat.infect <- dowloadFeatureHelper(loc, "dectree")
      if(is.null(dat.infect)) {
        dat.infect <- dowloadFeatureHelper(loc, "thresholded")
      }
      if(is.null(dat.infect)) {
        stop("try something other than 'dectree'/'thresholded'")
      } else {
        names(dat.infect) <- "Cells.Infection_IsInfected"
      }
      missing <- missing[-match("Cells.Infection_IsInfected", missing)]
    } else {
      dat.infect <- NULL
    }
    if (length(missing) > 0) {
      regexp <- paste0(".*/(", paste(missing, collapse="|"), ")\\.mat$")
      dat.rest <- dowloadFeatureHelper(loc, "cc", features=regexp)
    } else {
      dat.rest <- NULL
    }
    data <- c(dat.infect, dat.rest)
    dirs <- dir(getLocalPath(loc), "^HCS_ANALYSIS_CELL_", full.names=TRUE)
    unlink(dirs, recursive=TRUE)
  }
  x$data <- c(x$data, data)
  x$data <- x$data[order(names(x$data))]
  message("added features:\n  ", paste0(names(data), collapse="\n  "))
  saveToCache(x, force.write=TRUE)
  return(x)
}

#' @export
rebuildCache.default <- function(x, ...) {
  stop("can only deal with MatData objects.")
}