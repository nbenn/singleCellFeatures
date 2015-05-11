#' Get directory of the package source
#' 
#' Returns the directory of the package source
#'
#' @return The package source directory as a string
#'
#' @examples
#' dir <- getPackageDir()
#' 
#' @export
getPackageDir <- function() {
  # load current settingsDatabase
  data(settingsDatabase, envir=environment())
  return(settings.database$package)
}

#' Get directory of the package source
#' 
#' Returns the directory of the package source
#'
#' @return The package source directory as a string
#'
#' @examples
#' dir <- getPackageDir()
#' 
#' @export
reloadSingleCellFeatures <- function() {
  detach("package:singleCellFeatures", unload = TRUE, character.only = TRUE)
  library("singleCellFeatures", character.only = TRUE)
}

#' Get amount of storage used by an object
#' 
#' Prints the amount of storage used by the specified object.
#'
#' @param x The object of interest
#'
#' @return The amount of storage taken up by the object as a string
#'
#' @examples
#' plate <- WellLocation("J101-2C", "B", 15)
#' size(plate)
#' 
#' @export
size <- function(x) {
  format(object.size(x), units="auto")
}