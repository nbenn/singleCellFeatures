#' Linearize the index of a well on a 16x24 well plate
#'
#' Given row, column and possible image coordinates on a plate, get the
#' corresponding index of a linearized representation
#'
#' @param row  Row of the well, specified either by A:P, a:P, or 1:16
#' @param col  Column number of the well, specified by 1:24
#' @param img  NULL/Image number within the well, 1:9
#'
#' @return The linearized index; if img==NULL, range 1:384, else range: 1:3456
#'
#' @examples
#' ind1 <- getWellIndex1D("B", 12, NULL)
#' ind2 <- getWellIndex1D(2, 12, 5)
#'
getWellIndex1D <- function(row, col, img=NULL) {
  # input validation
  col <- as.integer(col)
  # 24 columns per plate  
  if (col < 1 | col > 24) 
    stop("col should be in the range 1:24, but is " , col, ".")
  # 16 rows indicated by A/a:P/p
  if (is.integer(row)) {
    # do nothing
  } else if (is.numeric(row)) {
    # coerce to integer
    row <- as.integer(row)
  } else if (is.character(row)) {
    # get corresponding integer
    upper <- which(LETTERS[1:16]==row)
    lower <- which(letters[1:16]==row)
    if (length(upper) == 1) row <- upper
    else if (length(lower) == 1) row <- lower
    else stop("row should be in the range a:p or A:P, but is ", row, ".")
  }
  # check if row range is in 1:16
  if (row < 1 | row > 16)
    stop("row has to be in the range 1:16, but is ", row, ".")
  
  if (!is.null(img)) {
    # 9 images per well
    img <- as.integer(img)
    if(img < 1 | img > 9) stop("img has to be in the range 1:9")
    # the data is in row major order, in groups of 9 images
    index <- ((row-1)*24+(col-1))*9+img
  } else {
    # the data is in row major order
    index <- (row-1)*24+col
  }
  
  return(index)    
}

#' Get the 2D index of a well on a 16x24 well plate
#'
#' Given the well index and a logical indication of whether image of well level
#' resolution is desired, calculate row, column and possibly image indices
#'
#' @param index 1D index of a well, either in 1:384 or 1:3456 depending on the
#'              image setting
#' @param image Well (FALSE) or image (TRUE) level resolution
#'
#' @return The 2D index of the specified well as a list; if image==FALSE,
#'         img=NULL
#'
#' @examples
#' index <- getWellIndex1D("H", 12, 5)
#' well  <- getWellIndex2D(index, TRUE)
#'
getWellIndex2D <- function(index, image=FALSE) {
  # input validation
  index <- as.integer(index)
  image <- as.logical(image)
  if(!image) {
    if (index < 1 | index > 384)
      stop("col has to be in the range 1:384, but is ", col, ".")
    # integer division
    row <- ((index-1) %/% 24) + 1
    # modulo division
    col <- ((index-1) %% 24)  + 1
    row <- LETTERS[1:16][row]
    img <- NULL
  } else {
    if (index < 1 | index > 3456)
      stop("col has to be in the range 1:3456, but is ", col, ".")
    img <- ((index-1) %% 9) + 1
    rem <- (index-1) %/% 9
    row <- (rem %/% 24) + 1
    col <- (rem %% 24)  + 1
    row <- LETTERS[1:16][row]
  }
  return(list(row=row, col=col, img=img))    
}