#' Enforce full rank design matrix
#' 
#' Given a design matrix (if the response is also included it has to be named 
#' "Response"), rank deficiency is corrected by dropping zero variance variables
#' and one of each pair of highly correlated variables to make the design matrix
#' invertible.
#'
#' @param data Matrix/data frame holding the design matrix and optionally a
#'             response vector named "Response"
#'
#' @return A subset of the original data structure (some of the variables are
#'         removed)
#'
#' @examples
#' # get gene locations
#' mtor.bruc.du <- findPlateWellsFromGene("brucella", "du-k", "MTOR")
#' scra.bruc.du <- findPlateWellsFromGene("brucella", "du-k", "SCRAMBLED")
#' # combine for faster fetching
#' targets <- rbind(mtor.bruc.du[1:6,], scra.bruc.du[c(6,18,78,90,150,162),])
#' # fetch data
#' data <- getPlateWells(targets)
#' # prepare data for glm
#' j101_2c <- prepareDataforGlm(data[["J101-2C"]]$H6$data, 
#'                              data[["J101-2C"]]$G23$data)
#' # enforce full rank
#' training <- makeRankFull(j101_2c$train)
#' # run glm
#' model <- glm("Response ~ .", binomial, training)
#' 
#' @export

makeRankFull <- function(data) {
  design.mat <- subset(data, select=-c(Response))
  orig.rank <- qr(design.mat)$rank
  if(orig.rank < length(design.mat)) {
    # rank deficiency
    cat("\nrank:", orig.rank, "full:", length(design.mat), "\n")    
    col.var <- apply(design.mat, 2, var)
    zero.var <- names(which(col.var==0))
    if(length(zero.var) > 0) {
      cat("\nremoving zero variance variables:\n")
      print(zero.var)
      data <- subset(data, select=c(names(which(col.var!=0)), "Response"))
      warning("removed ", length(zero.var), " zero variance variables.")  
    }
    design.mat <- subset(data, select=-c(Response))
    new.rank <- qr(design.mat)$rank
    if(new.rank < length(design.mat)) {
      # still rank deficient      
      cat("\nrank:", new.rank, "full:", length(design.mat), "\n")
      one.cor <- which(cor(design.mat) > 0.9999, arr.ind = TRUE)
      off.diag <- one.cor[one.cor[,1]!=one.cor[,2],]
      remo.ind <- unique(pmax(off.diag[,1], off.diag[,2]))
      if(length(remo.ind) > 0) {
        remo.nam <- names(design.mat)[remo.ind]
        keep.nam <- names(design.mat)[-remo.ind]
        cat("\nremoving highly correlated (>0.9999) variables:\n")
        print(remo.nam)        
        data <- subset(data, select=c(keep.nam, "Response"))  
        warning("removed ", length(remo.ind), " highly correlated (>0.9999) ",
                "variables.")  
      }
    }
    design.mat <- subset(data, select=-c(Response))
    new.rank <- qr(design.mat)$rank
    if(new.rank < length(design.mat)) {
      # still rank deficient      
      stop("\nrank: ", new.rank, " full: ", length(design.mat), "\n")
    } else cat("\nfull rank for design matrix\n\n")
  }
  return(data)
}