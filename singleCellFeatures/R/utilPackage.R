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

#' @export
clear <- function() {
  rm(list=ls(), envir=globalenv())
}

#' @export
size <- function(x) {
  format(object.size(x), units="auto")
}