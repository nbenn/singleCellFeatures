#' Check if the set of wells is complete
#'
#' Given a PlateData object, check if the set of wells is complete.
#'
#' @param x The PlateData object of interest.
#' 
#' @return Logical: TRUE if complete, FALSE if incomplete.
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' checkCompletenessWell(plate)
#' 
#' @export
checkCompletenessWell <- function(...) {
  UseMethod("checkCompletenessWell")
}

#' @export
checkCompletenessWell.PlateData <- function(x) {
  formatOutput <- function(str) {
    out <- paste0("\"", str, "\"")
    colwidth <- max(nchar(out))
    out <- stri_wrap(paste(stri_pad_right(out, colwidth),
                     collapse="  ", sep=""), normalize=FALSE)
    return(paste(out, collapse="\n  "))
  }
  well.exp <- paste0(rep(LETTERS[1:16], each=24), rep(1:24, 16))
  well.dat <- unique(names(x$data))
  missing <- setdiff(well.exp, well.dat)
  superfl <- setdiff(well.dat, well.exp)
  if(length(missing) > 0) {
    message("detected ", length(missing), " missing well(s) (",
            getBarcode(x), "):\n  ", formatOutput(missing))
  }
  if(length(superfl) > 0) {
    message("detected ", length(superfl), " superfluous image(s) (",
            getBarcode(x), "):\n  ", formatOutput(superfl))
  }

  is.well <- sapply(x$data, function(img) any(class(img) == "WellData"))
  if(!all(is.well)) {
    message("not all images in well ", paste0(x$row, x$col), " are of class ",
            "ImageData:\n  ", formatOutput(names(x$data)[!is.well]))
  }
  if(length(missing) > 0 | length(superfl) > 0 | !all(is.well)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' @export
checkCompletenessWell.default <- function(x) {
  stop("can only deal with PlateData objects.\nnot with ", class(x), ".")
}