#' Check if the set of images is complete
#'
#' Given a PlateData or WellData object, check if the set of images is complete
#' over the entire object.
#'
#' @param x The object of interest.
#' 
#' @return Logical: TRUE if complete, FALSE if incomplete.
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' checkCompletenessImage(plate)
#' 
#' @export
checkCompletenessImage <- function(...) {
  UseMethod("checkCompletenessImage")
}

#' @export
checkCompletenessImage.PlateData <- function(x) {
  formatOutput <- function(str) {
    out <- paste0("\"", str, "\"")
    colwidth <- max(nchar(out))
    out <- stri_wrap(paste(stri_pad_right(out, colwidth),
                     collapse="  ", sep=""), normalize=FALSE)
    return(paste(out, collapse="\n  "))
  }
  suppressMessages(complete <- sapply(x$data, checkCompletenessImage))
  if(!all(complete)) {
    message("not all wells on plate ", getBarcode(x), " contain a complete ",
            "set of images:\n  ", formatOutput(names(x$data)[!complete]))
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' @export
checkCompletenessImage.WellData <- function(x) {
  formatOutput <- function(str) {
    out <- paste0("\"", str, "\"")
    colwidth <- max(nchar(out))
    out <- stri_wrap(paste(stri_pad_right(out, colwidth),
                     collapse="  ", sep=""), normalize=FALSE)
    return(paste(out, collapse="\n  "))
  }
  n.img <- getNoImgPerWell(x)
  img.dat <- unique(names(x$data))
  if(n.img == 9) {
    img.exp <- c("img_11", "img_12", "img_13", "img_21", "img_22", "img_23",
                 "img_31", "img_32", "img_33")
  } else if (n.img == 6) {
    img.exp <- c("img_11", "img_12", "img_21", "img_22", "img_31", "img_32")
  } else {
    stop("expecting 6 or 9 images per well, not ", n.img, ".")
  }
  missing <- setdiff(img.exp, img.dat)
  superfl <- setdiff(img.dat, img.exp)
  if(length(missing) > 0) {
    message("detected ", length(missing), " missing image(s) (",
            getBarcode(x), "):\n  ", formatOutput(missing))
  }
  if(length(superfl) > 0) {
    message("detected ", length(superfl), " superfluous image(s) (",
            getBarcode(x), "):\n  ", formatOutput(superfl))
  }
  is.image <- sapply(x$data, function(img) any(class(img) == "ImageData"))
  if(!all(is.image)) {
    message("not all images in well ", paste0(x$row, x$col), " are of class ",
            "ImageData:\n  ", formatOutput(names(x$data)[!is.image]))
  }
  if(length(missing) > 0 | length(superfl) > 0 | !all(is.image)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' @export
checkCompletenessImage.default <- function(x) {
  stop("can only deal with PlateData/WellData objects.\nnot with ",
       class(x), ".")
}