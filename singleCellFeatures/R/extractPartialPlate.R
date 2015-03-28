#' Given feature data from a whole plate, extract the data for specified rows/
#' cols/wells/images
#'
#' @param data Data for a whole plate, multiply nested list of features
#' @param rows Row(s) of the well(s), specified either by a vector consisting 
#'             of the characters A:P, a:P, or 1:16
#' @param cols Column(s) of the well(s), specified by a vector of integers in 
#'             the range of 1:24
#' @param wels Well index(es), specified by a vector of integers in the range 
#'             of 1:384
#' @param imgs Image number(s) within the well, specified by a vector of 
#'             integers in the range of 1:9
#'
#' @return A multiply nested list of features with all values for the given 
#'         plate coordinates (for structure, refer to importCompletePlate)
#'
#' @examples
#' subset.img <- extractPartialPlate(dat.plat)
#' subset.wel <- extractPartialPlate(dat.plat, cols=c(3,4), rows=c("A","B"), 
#'                                   wels=c(333, 21), imgs=NULL)
#'
#' @export

extractPartialPlate <- function(data, rows=NULL, cols=NULL, wels=NULL, 
                                imgs=5) {
  # Input validation
  # 16x24 wells: check input data, as well as vector for wells
  if (length(data$data) != 384) stop("expecting $data to be of length 384")
  if (!(is.vector(wels, mode="numeric") | is.null(wels))) {
    stop("expecting NULL or an integer vector for wels")
  }
  if (!is.null(wels)) wels <- as.integer(wels)
  if (any(ifelse(wels < 1 | wels > 384, TRUE, FALSE))) {
    stop("expecting 1 <= wels <= 384")
  }
  # 9 images per well
  if (any(sapply(data$data, function(x) length(x$data) != 9))) {
    stop("expecting $data$plate$data to be of length 9")
  }
  if (!(is.vector(imgs, mode="numeric") | is.null(imgs))) {
    stop("expecting NULL or an integer vector for imgs")
  }
  if (!is.null(imgs)) imgs <- as.integer(imgs)
  if (any(ifelse(imgs < 1 | imgs > 9, TRUE, FALSE))) {    
    stop("expecting 1 <= imgs <= 9")
  }
  # type/range of cols/rows will be checked in getLinearizedPlateIndex
  if (!(is.vector(rows) | is.null(rows))) {
    stop("expecting NULL or an integer vector for rows")
  }
  if (!(is.vector(cols) | is.null(cols))) {
    stop("expecting NULL or an integer vector for cols")
  }
  
  # for the specified rows and cols, calculate the linear well indices
  if (is.null(rows) & !is.null(cols)) {
    wells.extra <- as.vector(sapply(cols, function(i) {
      sapply(1:16, function(j, i) getLinearizedPlateIndex(j, i, NULL), i)
    }))
    wels <- c(wels, wells.extra)
  } else if (!is.null(rows) & is.null(cols)) {
    wells.extra <- as.vector(sapply(rows, function(i) {
      sapply(1:24, function(j, i) getLinearizedPlateIndex(i, j, NULL), i)
    }))
    wels <- c(wels, wells.extra)
  } else if (!is.null(rows) & !is.null(cols)) {
    wells.extra <- as.vector(sapply(rows, function(i, cols) {
      sapply(cols, function(j, i) getLinearizedPlateIndex(i, j, NULL), i)
    }, cols))
    wels <- c(wels, wells.extra)
  } else if (is.null(rows) & is.null(cols) & is.null(wels)) {
    wels <- 1:384
  }
  # drop duplicates that can occur through specifying the same well twice
  wels <- sort(unique(wels))
  
  # extract the data from input
  if (is.null(imgs)) {
    well.data <- lapply(wels, function (i, data) data$data[[i]], data)
  } else {
    imgs <- sort(unique(imgs))
    well.data <- lapply(wels, function (i, j, data) {
      images <- lapply(j, function(j, data) data$data[[j]], data$data[[i]])
      names(images) <- names(data$data[[i]]$data)[j]
      return(list(meta=data$data[[i]]$meta, data=images))
    }, imgs, data) 
  }
  names(well.data) <- names(data$data)[wels]
    
  return(list(meta=data$meta, data=well.data))
}