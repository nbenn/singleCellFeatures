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
#' @export
getWellIndex1D <- function(row, col, img=NULL, n.img=9) {
  # input validation
  if(!n.img %in% c(6, 9)) stop("n.img has to be 6 or 9")
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
    upper <- which(LETTERS[1:16] == row)
    lower <- which(letters[1:16] == row)
    if (length(upper) == 1) row <- upper
    else if (length(lower) == 1) row <- lower
    else stop("row should be in the range a:p or A:P, but is ", row, ".")
  }
  # check if row range is in 1:16
  if (row < 1 | row > 16)
    stop("row has to be in the range 1:16, but is ", row, ".")

  if (!is.null(img)) {
    # 6/9 images per well
    img <- as.integer(img)
    if(img < 1 | img > n.img) stop("img has to be in the range 1:", n.img)
    # the data is in row major order, in groups of 6/9 images
    index <- ((row - 1) * 24 + (col - 1)) * n.img + img
  } else {
    # the data is in row major order
    index <- (row - 1) * 24 + col
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
#' @param n.img Well (NULL) or image (6/9) level resolution
#'
#' @return The 2D index of the specified well as a list; if image == NULL,
#'         img.ind, img.row and img.col = NULL
#'
#' @examples
#' index <- getWellIndex1D("H", 12, 5)
#' well  <- getWellIndex2D(index, 9)
#'
#' @export
getWellIndex2D <- function(index, n.img=NULL) {
  # input validation
  index <- as.integer(index)
  if(!is.null(n.img)) {
    if(!n.img %in% c(6, 9)) stop("n.img has to be NULL, 6 or 9")
  }
  if(is.null(n.img)) {
    if (index < 1 | index > 384)
      stop("col has to be in the range 1:384, but is ", col, ".")
    # integer division
    row <- ((index - 1) %/% 24) + 1
    # modulo division
    col <- ((index - 1) %% 24)  + 1
    row <- LETTERS[1:16][row]
    img <- NULL
    i.c <- NULL
    i.r <- NULL
  } else {
    tot.nimgs <- n.img * 384
    if (index < 1 | index > tot.nimgs)
      stop("col has to be in the range 1:", tot.nimgs, ", but is ", col, ".")
    img <- ((index - 1) %% n.img) + 1
    rem <- (index - 1) %/% n.img
    row <- (rem %/% 24) + 1
    col <- (rem %% 24)  + 1
    row <- LETTERS[1:16][row]
    img.cols <- n.img / 3
    i.r <- ((img - 1) %/% img.cols) + 1
    i.c <- ((img - 1) %% img.cols)  + 1
  }
  return(list(wel.row=row, wel.col=col, img.ind=img,
              img.row=i.r, img.col=i.c))
}