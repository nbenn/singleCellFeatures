#' Constructor for WellLocation objects
#' 
#' Given a well specification consisting of plate barcode, row and column, an
#' S3 object of type WellLocation is created
#'
#' @param plate A plate barcode as a character string
#' @param row   Row specification (a:p, A:P, or 1:16)
#' @param col   Column specification (1:24)
#'
#' @return A WellLocation object: a list with slots plate, row, column, index
#'         (the linearized well index), space, group and experiment
#'
#' @examples
#' well <- WellLocation("J107-2C", "H", 6)
#' 
#' @export
WellLocation <- function(plate, row, col) {
  ind <- getWellIndex1D(row, col)
  row <- getWellIndex2D(ind)$wel.row
  if (!is.character(plate)) stop("plate must be character string")
  data(plateDatabase, envir=environment())
  match <- plate.database[plate.database$Barcode == plate,]
  if(nrow(match) != 1) stop("could not find the specified plate")
  result <- list(plate = plate, row=row, column=col, index=ind,
                 space=match$Space, group=match$Group,
                 experiment=match$Experiment)
  return(structure(result, class = c("WellLocation", "DataLocation")))
}

#' Constructor for PlateLocation objects
#' 
#' Given a plate specification (consisting of a plate barcode) an S3 object of
#' type PlateLocation is created
#'
#' @param barcode A plate barcode as a character string
#'
#' @return A PlateLocation object: a list with slots plate, space, group and
#'         experiment
#'
#' @examples
#' plate <- PlateLocation("J107-2C")
#' 
#' @export
PlateLocation <- function(barcode) {
  if (!is.character(barcode)) stop("barcode must be character string")
  data(plateDatabase, envir=environment())
  match <- plate.database[plate.database$Barcode == barcode,]
  if(nrow(match) != 1) stop("could not find the specified plate")
  result <- list(plate = barcode, space=match$Space, group=match$Group,
                 experiment=match$Experiment)
  return(structure(result, class = c("PlateLocation", "DataLocation")))
}