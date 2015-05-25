#' Check completeness of data
#'
#' Given a Data object, determine if the set of available features is complete.
#'
#' @param x Data object
#' 
#' @return A list holding two vectors of strings, one for features that are
#'         expected but missing and one for features that are present but not
#'         expected.
#' 
#' @examples
#' data  <- MatData(PlateLocation("J107-2C"))
#' checkCompletenessFeature(data)
#'
#' @export
checkCompletenessFeature <- function(x) {
  UseMethod("checkCompletenessFeature", x)
}

#' @export
checkCompletenessFeature.Data <- function(x) {
  ignored.feats <- c("Bacteria.SubObjectFlag", "Batch_handles",
                     "Image.ModuleError_43CreateBatchFiles")
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
  ignored <- missing[missing %in% ignored.feats]
  missing <- missing[!missing %in% ignored.feats]
  superfl <- setdiff(feat.dat, feat.exp)
  result <- TRUE
  if(length(ignored) > 0) {
    message(length(ignored), " ignored feature(s) missing (",
            getBarcode(plate), "):\n  ", paste(ignored, collapse="\n  "))
  }
  if(length(missing) > 0) {
    warning("detected ", length(missing), " missing feature(s) (",
            getBarcode(plate), "):\n  ", paste(missing, collapse="\n  "))
    result <- FALSE
  }
  if(length(superfl) > 0) {
    warning("detected ", length(superfl), " superfluous feature(s) (",
            getBarcode(plate), "):\n  ", paste(superfl, collapse="\n  "))
    result <- FALSE
  }
  return(list(missing = missing, superfluous = superfl, result=result))
}

#' @export
checkCompletenessFeature.default <- function(x, ...) {
  stop("can only deal with Data objects.")
}