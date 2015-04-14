#' Prepare a dataset for GLM analysis
#' 
#' Given two design matrices (one corresponding to a knockdown and one to a
#' control well), create a concatenated design matrix with a column "Response"
#' encoding for original membership. This unified data matrix is then divided
#' into 90% training and 10% testing. To avoid missingness problems, all
#' variables containing NA/NaN are dropped.
#'
#' @param active  Data coming from a well where a gene knockdown occurred
#' @param control Data belonging to a control well
#'
#' @return A list with entries "test" and "train" each holding a data frame
#'         containing a design matrix and a response vector
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
#' # run glm (enforcing full rank might be advisbale)
#' model <- glm("Response ~ .", binomial, j101_2c$train)
#' 
#' @export

prepareDataforGlm <- function(active, control) {

  prepareSet <- function(x, isActive) {
    lapply(x, function(img) {
      len <- unlist(lapply(img, function(x) {
        return(length(x))
      }))
      if(length(unique(len)) != 1) {
        stop("non-equal number of ojects within image")
      }
    })
    result <- do.call("rbind", lapply(x, as.data.frame))
    if(isActive) {
      result$Response <- factor(rep("active", nrow(result)),
                                levels=c("active", "control"))
    } else {
      result$Response <- factor(rep("control", nrow(result)),
                                levels=c("active", "control"))
    }
    return(result)
  }
  
  act <- prepareSet(active, TRUE)
  ctr <- prepareSet(control, FALSE)
  data <- rbind(act,ctr)
  complete <- complete.cases(t(data))
  if(sum(!complete) > 0) {
    cat("\nfound incomplete cases. removing vars:\n")
    print(names(data[!complete]))
    data <- data[,complete.cases(t(data))]
    warning("removed ", sum(!complete), " variables containing Na/NaN.")
  }
  set.seed(7)
  test.ind <- sort(sample.int(nrow(data), nrow(data)/10))
  testing <- data[test.ind,]
  trainin <- data[-test.ind,]
  return(list(test=testing, train=trainin))
}