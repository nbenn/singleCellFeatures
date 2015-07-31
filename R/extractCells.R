#' Extract cells from Data objects
#' 
#' Given a Data object, extract all infected or all uninfected cells.
#'
#' @param x        The original Data object.
#' @param infected Either NULL or a boolgen whether to return all infected (
#'                 TRUE) or all uninfected cells (FALSE).
#'
#' @return A Data object of the same type as x, but only containing a subset of
#'         all cells.
#'
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' data <- extractCells(data, infected=TRUE)
#' 
#' @export
extractCells <- function(x, ...) {
  UseMethod("extractCells", x)
}

#' @export
extractCells.PlateData <- function(x, infected=NULL) {
  x$data <- lapply(x$data, extractCells, infected)
  return(x)
}

#' @export
extractCells.WellData <- function(x, infected=NULL) {
  x$data <- lapply(x$data, extractCells, infected)
  return(x)
}

#' @export
extractCells.ImageData <- function(x, infected=NULL) {
  if(!is.null(infected)) {
    x$data.mat <- lapply(x$data.mat, function(group, infec) {
      if("Cells.Infection_IsInfected" %in% colnames(group)) {
        inf.ind <- match("Cells.Infection_IsInfected", colnames(group))
        inf.ind <- group[,inf.ind]
        if(infec) {
          return(group[inf.ind == 1,, drop=FALSE])
        } else {
          return(group[inf.ind == 0,, drop=FALSE])
        }
      } else {
        return(group)
      }
    }, infected)
  }
  return(x)
}

#' @export
extractCells.default <- function(x, ...) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}