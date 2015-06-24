#' Remove some bad image data
#'
#' Given a PlateData/WellData object, remove all images where the cell count
#' does not lie within the 0.05 and 0.95 quantiles.
#'
#' @param x        The PlateData/WellData object of interest.
#' @param quants   A string for choosing between only applying lower, upper,
#'                 both one none of the quantiles for cleaning up the data.
#' @param infected A string for choosing whether to discard infected or
#'                 uninfected cells.
#' 
#' @return A cleaned up version of the input object.
#' 
#' @examples
#' data <- PlateData(PlateLocation("J101-2C"))
#' cleaned <- cleanData(data)
#' 
#' @export
cleanData <- function(x, ...) {
  UseMethod("cleanData", x)
}

#' @export
cleanData.PlateData <- function(x, ...) {
  x$data <- lapply(x$data, cleanData, ...)
  return(x)
}

#' @export
cleanData.WellData <- function(x, quants="both", infected="keepAll") {
  if(!(quants %in% c("upper", "lower", "both", "none"))) {
    stop("expecting a string in c(\"upper\", \"lower\", \"both\", \"none\") ",
         "for quants.")
  }
  if(!(infected %in% c("remUninf", "keepAll", "remInf"))) {
    stop("expecting a string in c(\"remUninf\", \"keepAll\", \"remInf\") ",
         "for infected.")
  }
  quant1 <- x$data[[1]]$counts.quantiles[1]
  quant2 <- x$data[[1]]$counts.quantiles[2]
  rem.upper <- sapply(x$data, function(img) {
    if(img$counts.cells > quant2) return(TRUE)
    else return(FALSE)
  })
  rem.lower <- sapply(x$data, function(img) {
    if(img$counts.cells < quant1) return(TRUE)
    else return(FALSE)
  })
  if(sum(rem.upper) > 0 | sum(rem.lower) > 0) {
    message("well ", getWellName(x), " (", getBarcode(x), "):")
  }
  if(sum(rem.upper) > 0) {
    discarded <- which(rem.upper)
    counts <- sapply(discarded, function(i, x) {
      return(x$data[[i]]$counts.cells)
    }, x)
    if(quants %in% c("upper", "both")) {
       message("  discarding ", length(discarded), " images (",
              paste0(discarded, collapse=", "), ") because count.cells ",
              "not in [", quant1, ", ", quant2, "] but ",
              paste0(counts, collapse=", "), ".")
      x$data <- x$data[-rem.upper]
    } else {
       message("  keeping ", length(discarded), " images (",
              paste0(discarded, collapse=", "), ") despite count.cells ",
              "not in [", quant1, ", ", quant2, "] but ",
              paste0(counts, collapse=", "), ".")
    }
  }
  if(sum(rem.lower) > 0) {
    discarded <- which(rem.lower)
    counts <- sapply(discarded, function(i, x) {
      return(x$data[[i]]$counts.cells)
    }, x)
    if(quants %in% c("lower", "both")) {
       message("  discarding ", length(discarded), " images (",
              paste0(discarded, collapse=", "), ") because count.cells ",
              "not in [", quant1, ", ", quant2, "] but ",
              paste0(counts, collapse=", "), ".")
      x$data <- x$data[-rem.lower]
    } else {
       message("  keeping ", length(discarded), " images (",
              paste0(discarded, collapse=", "), ") despite count.cells ",
              "not in [", quant1, ", ", quant2, "] but ",
              paste0(counts, collapse=", "), ".")
    }
  }

  if(infected != "keepAll") {
    x$data <- lapply(x$data, cleanData, infected)
  }
  return(x)
}

#' @export
cleanData.ImageData <- function(x, infected="keepAll") {
  if(!(infected %in% c("remUninf", "keepAll", "remInf"))) {
    stop("expecting a string in c(\"remUninf\", \"keepAll\", \"remInf\") ",
         "for infected.")
  }
  if(infected != "keepAll") {
    isInf <- as.logical(x$data.mat$Cells[,"Cells.Infection_IsInfected"])
  }

  if(infected == "remUninf") {
    x$data.mat$Cells <- x$data.mat$Cells[isInf, , drop=FALSE]
  } else if(infected == "remInf") {
    x$data.mat$Cells <- x$data.mat$Cells[!isInf, , drop=FALSE]
  }

  return(x)
}

#' @export
cleanData.default <- function(x, ...) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}