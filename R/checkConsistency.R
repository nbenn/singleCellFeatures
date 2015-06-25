#' Check consistency of Data object
#'
#' Given a PlateData or WellData object, check if all sub-objects are
#' consistent, i.e. if they contain the same features, barcodes, well names.
#'
#' @param x The object of interest.
#' 
#' @return Logical: TRUE if consistent, FALSE if inconsistent (the
#' inconsistencies are printed).
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' checkConsistency(plate)
#' 
#' @export
checkConsistency <- function(...) {
  UseMethod("checkConsistency")
}

#' @export
checkConsistency.PlateData <- function(x) {
  prototype <- buildPrototypeHelper(x$data[[1]]$data[[1]])[-2]
  if(getBarcode(x) != prototype$plate.barcode) {
    meta <- FALSE
    message("metadata and data are not consistent: ", getBarcode(x), " != ",
            prototype$plate.barcode)
  } else {
    meta <- TRUE
  }
  res <- sapply(x$data, function(well, proto) {
    compare <- buildPrototypeHelper(well$data[[1]])[-2]
    comparison <- all.equal(compare, proto)
    if(!is.logical(comparison)) {
      message("for plate ", compare$plate.barcode, ", well ",
              getWellName(well), ":\n  ", paste(comparison, collpase="\n  "))
    }
    within.well <- checkConsistency(well)
    if(within.well & is.logical(comparison)) return(TRUE)
    else return(FALSE)
  }, prototype)
  return(all(res) & meta)
}

#' @export
checkConsistency.WellData <- function(x) {
  prototype <- buildPrototypeHelper(x$data[[1]])
  meta <- TRUE
  if(getBarcode(x) != prototype$plate.barcode) {
    meta <- FALSE
    message("metadata and data are not consistent: ", getBarcode(x), " != ",
            prototype$plate.barcode)
  }
  if(getWellName(x) != prototype$well.name) {
    meta <- FALSE
    message("metadata and data are not consistent: ", getWellName(x), " != ",
            prototype$well.name, " (plate ", prototype$plate.barcode, ")")
  }

  res <- sapply(x$data, function(image, proto) {
    compare <- buildPrototypeHelper(image)
    comparison <- all.equal(compare, proto)
    if(!is.logical(comparison)) {
      message("for plate ", compare$plate.barcode, ", well ", compare$well.name,
              ", image ", image$image.index, ":\n  ",
              paste(comparison, collpase="\n  "))
      return(FALSE)
    } else {
      return(comparison)
    }
  }, prototype)
  return(all(res) & meta)
}

#' @export
checkConsistency.default <- function(x) {
  stop("can only deal with PlateData/WellData objects.\nnot with ",
       class(x), ".")
}

#' @export
buildPrototypeHelper <- function(img) {
  if(!any(class(img) == "ImageData")) {
    stop("expecting an object of class ImageData")
  }
  vec.feat <- lapply(img$data.vec, function(group) return(names(group)))
  mat.feat <- lapply(img$data.mat, function(group) return(colnames(group)))
  lst.feat <- lapply(img$data.lst, function(group) return(names(group)))
  return(list(plate.barcode = img$plate,
              well.name     = getWellName(img),
              image.total   = img$image.total,
              vec.names     = names(img$data.vec),
              vec.feat      = vec.feat,
              mat.names     = names(img$data.mat),
              mat.feat      = mat.feat,
              lst.names     = names(img$data.lst),
              lst.feat      = lst.feat))
}
