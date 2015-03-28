#' Linearize the index of a 16x24 well plate
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
#' getLinearizedPlateIndex("B", 12, NULL)
#' getLinearizedPlateIndex(2, 12, 5)

getLinearizedPlateIndex <- function(row, col, img=NULL) {
  # Input validation
  col <- as.integer(col)
  # 24 columns per plate  
  if (col < 1 | col > 24) stop("col has to be in the range 1:24")
  # 16 rows indicated by A/a:P/p
  if (!is.integer(row)) {
    if (row == "a" | row == "A") row <- 1
    else if (row == "b" | row == "B") row <- 2
    else if (row == "c" | row == "C") row <- 3
    else if (row == "d" | row == "D") row <- 4
    else if (row == "e" | row == "E") row <- 5
    else if (row == "f" | row == "F") row <- 6
    else if (row == "g" | row == "G") row <- 7
    else if (row == "h" | row == "H") row <- 8
    else if (row == "i" | row == "I") row <- 9
    else if (row == "j" | row == "J") row <- 10
    else if (row == "k" | row == "K") row <- 11
    else if (row == "l" | row == "L") row <- 12
    else if (row == "m" | row == "M") row <- 13  
    else if (row == "n" | row == "N") row <- 14 
    else if (row == "o" | row == "O") row <- 15 
    else if (row == "p" | row == "P") row <- 16
    else row <- as.integer(row)
  }
  # or 1:16
  if (row < 1 | row > 16) stop("row has to be in the range 1:16, a:p or A:P")
  
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