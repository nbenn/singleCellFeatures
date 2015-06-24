#' Remove some bad image data
#'
#' Given a PlateData/WellData object, remove all images where the cell count
#' does not lie within the 0.05 and 0.95 quantiles.
#'
#' @param x      The PlateData/WellData object of interest.
#' @param quants A string for choosing betwenn only applying lower, upper,
#'               both one none of the quantiles for cleaning up the data.
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
cleanData.PlateData <- function(x) {
  x$data <- lapply(x$data, cleanData)
  return(x)
}

#' @export
cleanData.WellData <- function(x, quants="both") {
  if(!(quants %in% c("upper", "lower", "both", "none"))) {
    stop("expecting a string in c(\"upper\", \"lower\", \"both\", \"none\") ",
         "for quants.")
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
  return(x)
}

#' @export
cleanData.default <- function(x, ...) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}