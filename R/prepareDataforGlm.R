#' Prepare a dataset for GLM analysis
#' 
#' Given two design matrices (one corresponding to a knockdown and one to a
#' control well), create a concatenated design matrix with a column "Response"
#' encoding for original membership. This unified data matrix is then divided
#' into 90% training and 10% testing. To avoid missingness problems, all
#' variables containing NA/NaN are dropped.
#'
#' @param active  Data coming from a well where a gene knockdown occurred.
#' @param control Data belonging to a control well.
#' @param drop    A vector of strings of column names that will be dropped.
#'
#' @return A list with entries "test" and "train" each holding a data frame
#'         containing a design matrix and a response vector.
#'
#' @examples
#' # get gene locations
#' mtor.loc <- findWells(experiments="brucella-du-k1", contents="MTOR")
#' scra.loc <- findWells(plates=sapply(mtor.loc, getBarcode),
#'                       contents="SCRAMBLED", well.names="G23")
#' # combine for faster fetching
#' data     <- getSingleCellData(list(mtor.loc[[1]], scra.loc[[1]]))
#' mtor.dat <- meltData(cleanData(data[[1]]$H6))
#' scra.dat <- meltData(cleanData(data[[1]]$G23))
#' # prepare data for glm
#' data <- prepareDataforGlm(mtor.dat$mat$Cells, scra.dat$mat$Cells)
#' data <- makeRankFull(data)
#' # run glm
#' model <- glm("Response ~ .", binomial, data$train)
#' 
#' @export

prepareDataforGlm <- function(active, control, drop=NULL) {

  if(!(is.data.frame(active) & is.data.frame(control))) {
    stop("expecting data frames as aguments active and control.")
  }

  active$Response <- factor(rep("active", nrow(active)),
                                levels=c("active", "control"))
  control$Response <- factor(rep("control", nrow(control)),
                                levels=c("active", "control"))

  data.all <- rbind(active, control)

  drop <- c(drop, "Image.Index", "Well.Index", "Well.Name", "Plate.Barcode")
  data.all <- data.all[, !names(data.all) %in% drop]

  complete <- complete.cases(t(data.all))
  if(sum(!complete) > 0) {
    message("found incomplete cases; removing vars:\n  ",
            paste(names(data.all[!complete]), collapse="\n  "))
    data.all <- data.all[, complete]
    warning("removed ", sum(!complete), " variables containing Na/NaN.")
  }
  set.seed(7)
  test.ind <- sort(sample.int(nrow(data.all), nrow(data.all) / 10))
  testing <- data.all[test.ind,]
  trainin <- data.all[-test.ind,]
  return(list(test=testing, train=trainin))
}

#' Enforce full rank design matrix
#' 
#' Given a design matrix (if the response is also included it has to be named 
#' "Response"), rank deficiency is corrected by dropping zero variance variables
#' and one of each pair of highly correlated variables to make the design matrix
#' invertible.
#'
#' @param data    Matrix/data frame holding the design matrix and optionally a
#'                response vector named "Response"
#'
#' @return A subset of the original data structure (some of the variables are
#'         removed)
#'
#' @examples
#' # get gene locations
#' mtor.loc <- findWells(experiments="brucella-du-k1", contents="MTOR")
#' scra.loc <- findWells(plates=sapply(mtor.loc, getBarcode),
#'                       contents="SCRAMBLED", well.names="G23")
#' # combine for faster fetching
#' data     <- getSingleCellData(list(mtor.loc[[1]], scra.loc[[1]]))
#' mtor.dat <- meltData(cleanData(data[[1]]$H6))
#' scra.dat <- meltData(cleanData(data[[1]]$G23))
#' # prepare data for glm
#' data <- prepareDataforGlm(mtor.dat$mat$Cells, scra.dat$mat$Cells)
#' data <- makeRankFull(data)
#' # run glm (enforcing full rank might be advisbale)
#' model <- glm("Response ~ .", binomial, data$train)
#' 
#' @export

makeRankFull <- function(data) {
  if(!is.list(data)) stop("expecting a LIST")
  if(length(data) != 2) stop("expecting a list of LENGTH 2")
  if(!all(names(data) %in% c("test", "train"))) {
    stop("expecting a list with slots \"test\" and \"train\".")
  }
  temp <- data$train
  temp <- temp[, !names(temp) %in% "Response"]
  remove <- NULL
  orig.rank <- qr(temp)$rank
  if(orig.rank < ncol(temp)) {
    # rank deficiency
    message("rank: ", orig.rank, " full: ", ncol(temp))
    col.var <- apply(temp, 2, var)
    remove <- names(which(col.var == 0))
    if(length(remove) > 0) {
      message("removing zero variance variables:\n  ",
              paste(remove, collapse="\n  "))
      temp <- temp[, !names(temp) %in% remove]
      warning("removed ", length(remove), " zero variance variables.")
    }
    new.rank <- qr(temp)$rank
    if(new.rank < ncol(temp)) {
      # still rank deficient
      message("rank: ", new.rank, " full: ", ncol(temp))
      one.cor <- which(cor(temp) > 0.9999, arr.ind = TRUE)
      off.diag <- one.cor[one.cor[,1] != one.cor[,2],]
      remo.ind <- unique(pmax(off.diag[,1], off.diag[,2]))
      if(length(remo.ind) > 0) {
        remo.nam <- names(temp)[remo.ind]
        message("removing one of highly correlated (>0.9999) variable pairs:",
                "\n  ", paste(remo.nam, collapse="\n  "))
        temp <- temp[, !names(temp) %in% remo.nam]
        warning("removed ", length(remo.ind), " variables due to highly ",
                "correlation (>0.9999) ")
        remove <- c(remove, remo.nam)
      }
    }
    new.rank <- qr(temp)$rank
    if(new.rank < ncol(temp)) {
      # still rank deficient
      stop("rank: ", new.rank, " full: ", ncol(temp))
    }
  }
  if(length(remove) > 0) {
    data$train <- data$train[, !names(data$train) %in% remove]
    data$test  <- data$test[, !names(data$test) %in% remove]
  }
  return(data)
}

#' Compare result of binary classification to truth
#'
#' Given two vectors of factors (two identical levels) of equal length,
#' calculate confusion matrix marginals
#'
#' @param estim Vector holding a two level factor with the classification
#'              result.
#' @param truth Vector holding a two level factor with the truth.
#'
#' @return A list with various key charactereistics of the resulting confusion
#'         matrix
#'
#' @examples
#' # get gene locations
#' mtor.loc <- findWells(experiments="brucella-du-k1", contents="MTOR")
#' scra.loc <- findWells(plates=sapply(mtor.loc, getBarcode),
#'                       contents="SCRAMBLED", well.names="G23")
#' # combine for faster fetching
#' data     <- getSingleCellData(list(mtor.loc[[1]], scra.loc[[1]]))
#' mtor.dat <- meltData(cleanData(data[[1]]$H6))
#' scra.dat <- meltData(cleanData(data[[1]]$G23))
#' # prepare data for glm
#' data <- prepareDataforGlm(mtor.dat$mat$Cells, scra.dat$mat$Cells)
#' data <- makeRankFull(data)
#' # run glm
#' model <- glm("Response ~ .", binomial, data$train)
#' # compare to testing data
#' predi <- as.factor(round(predict(model, newdata=data$test, type="response")))
#' levels(predi) <- c("active", "control")
#' comparison <- compareModeltoTruth(predi, data$test$Response)
#'
#' @export

compareModeltoTruth <- function(estim, truth) {
  if(length(estim) != length(truth)) {
    stop("comparing vectors of unequal length")
  }
  if(!identical(levels(truth), levels(estim))) {
    stop("levels don't match")
  }
  po <- levels(estim)[1]
  ne <- levels(estim)[2]
  tp <- sum(estim == po & truth == po)
  tn <- sum(estim == ne & truth == ne)
  fp <- sum(estim == po & truth == ne)
  fn <- sum(estim == ne & truth == po)
  tpr <- tp / (tp + fn)
  tnr <- tn / (fp + tn)
  fpr <- fp / (fp + tn)
  fnr <- fn / (fn + tp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  fdr <- fp / (fp + tp)
  fur <- fn / (fn + tn)
  f1s <- 2 * tp / (2 * tp + fp + fn)
  acc <- (tp + tn) / (tp + tn + fp + fn)
  n <- tn + tp + fn + fp
  s <- (tp + fn) / n
  p <- (tp + fp) / n
  mcc <- (tp / n - s * p) / sqrt(p * s * (1 - s) * (1 - p))
  return(list(POS=po, NEG=ne, TP=tp, TN=tn, FP=fp, FN=fn, TPR=tpr, TNR=tnr,
              FPR=fpr, FNR=fnr, PPV=ppv, NPV=npv, FDR=fdr, FOR=fur, F1S=f1s,
              ACC=acc, MCC=mcc)
  )
}