#' Check completeness of data
#'
#' Given a Data object (currently only MatData objects supported), determine if
#' the set of available features is complete.
#'
#' @param x Data object
#' 
#' @return A list holding two vectors of strings, one for features that are
#'         expected but missing and one for features that are present but not
#'         expected.
#' 
#' @examples
#' data  <- MatData(PlateLocation("J107-2C"))
#' checkDataCompleteness(data)
#'
#' @export
checkDataCompleteness <- function(x) {
  UseMethod("checkDataCompleteness", x)
}

#' @export
checkDataCompleteness.MatData <- function(x) {
  plate <- convertToPlateLocation(x)
  # check for completenes of single cell data
  pathogen <- getPathogen(plate)
  # load feature database
  data(featureDatabase, envir=environment())
  feat.exp <- feature.database[[pathogen]]
  feat.dat <- getFeatureNames(x)
  if(is.null(feat.exp)) {
    warning("for ", pathogen, ", currently no feature list is available. ",
            "Please add one using updateDatabaseFeatures.")
  }
  missing <- setdiff(feat.exp, feat.dat)
  superfl <- setdiff(feat.dat, feat.exp)
  if(length(missing) > 0) {
    warning("detected ", length(missing), " missing feature(s) (",
            getBarcode(plate), "):\n  ", paste(missing, collapse="\n  "))
  }
  if(length(superfl) > 0) {
    warning("detected ", length(superfl), " superfluous feature(s) (",
            getBarcode(plate), "):\n  ", paste(superfl, collapse="\n  "))
  }
  return(list(missing = missing, superfluous = superfl))
}

#' @export
checkDataCompleteness.default <- function(x, ...) {
  stop("can only deal with MatData objects.")
}