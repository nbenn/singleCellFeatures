#' Enforce equal feature sets
#' 
#' Given a list of DataObject, ensure that all list enstries have the same set
#' of features.
#'
#' @param lst  A list of DataObjects.
#'
#' @return A list of DataObjects (same type as input) which is guaranteed to
#'         have the same feature set.
#'
#' @examples
#' wells <- findWells(plates=c("J110-2C", "J104-2D", "J107-2L"),
#'                    well.names=c("H2", "H6"))
#' data  <- unlist(getSingleCellData(wells), recursive=FALSE)
#' data  <- lapply(data, cleanData, "lower")
#' data  <- makeFeatureCompatible(data) 
#' 
#' @export
makeFeatureCompatible <- function(lst) {
  if(!all(sapply(lst, function(x) any(class(x) == "Data")))) {
    stop("expecting a list of Data objects")
  }
  if(!all(sapply(lst, checkConsistency))) {
    stop("Data objects have to be consistent")
  }

  features <- lapply(lst, getFeatureNames)
  intersection <- Reduce(intersect, features)
  res <- lapply(lst, extractFeatures, features=intersection)
  return(res)
}