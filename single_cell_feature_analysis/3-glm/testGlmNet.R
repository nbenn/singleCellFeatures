## Load functions
library(singleCellFeatures)

## Define functions
prepareData <- function(active, scram) {
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
                                levels=c("active", "scrambled"))
    } else {
      result$Response <- factor(rep("scrambled", nrow(result)),
                                levels=c("active", "scrambled"))
    }
    return(result)
  }
  act <- prepareSet(active, TRUE)
  scr <- prepareSet(scram, FALSE)
  data <- rbind(act,scr)
  complete <- complete.cases(t(data))
  if(sum(!complete) > 0) {
    cat("\nfound incomplete cases. removing vars:\n")
    print(names(data[!complete]))
    data <- data[,complete.cases(t(data))]
  }
  set.seed(7)
  test.ind <- sort(sample.int(nrow(data), nrow(data)/10))
  testing <- data[test.ind,]
  trainin <- data[-test.ind,]
  return(list(test=testing, train=trainin))
}

compareModeltoTruth <- function(estim, truth) {
  tp <- sum(estim==truth)
  fp <- sum(estim!=truth)
  al <- length(estim==truth)
  cat("true positive rate: ", tp/al, "\nfalse positive rate:", fp/al)
}

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

glm.regular <- function(data) {
  design.mat <- makeRankFull(data$train)
  model <- glm("Response ~ .", binomial, design.mat)
  predi <- as.factor(round(predict(model, newdata=data$test, type="response")))
  levels(predi) <- c("active", "scrambled")
  compareModeltoTruth(predi, data$test$Response)
}

glm.glm2 <- function(data) {
  library(glm2)
  design.mat <- makeRankFull(data$train)
  model <- glm2("Response ~ .", binomial, design.mat)
  predi <- as.factor(round(predict(model, newdata=data$test, type="response")))
  levels(predi) <- c("active", "scrambled")
  compareModeltoTruth(predi, data$test$Response)
}

glm.safe <- function(data) {
  library(safeBinaryRegression)
  design.mat <- makeRankFull(data$train)
  #design.mat <- data$train
  tryCatch({
    model <- glm("Response ~ .", binomial, design.mat)
    predi <- as.factor(round(predict(model, newdata=data$test,
                                     type="response")))
    levels(predi) <- c("active", "scrambled")
    compareModeltoTruth(predi, data$test$Response)
  }, finally = {
      detach("package:safeBinaryRegression", unload=TRUE)
  })
}

glm.brglm <- function(data) {
  library(brglm)  
  design.mat <- makeRankFull(data$train)  
  model <- brglm("Response ~ .", data=design.mat)
  predi <- as.factor(round(predict(model, newdata=data$test, type="response")))
  levels(predi) <- c("active", "scrambled")
  compareModeltoTruth(predi, data$test$Response)
}

glm.glmnet <- function(data) {
  library(glmnet)
  x <- as.matrix(subset(data$train, select=-c(Response)))
  y <- data$train$Response
  model <- glmnet(x, y, family="binomial")
  x.test <- as.matrix(subset(data$test, select=-c(Response)))
  predi <- as.factor(round(predict(model, newx=x.test, type="response")))
  levels(predi) <- c("active", "scrambled")
  compareModeltoTruth(predi, data$test$Response)
}

glm.bestglm <- function(data) {
  library(bestglm)
  full.rank <- makeRankFull(data$train)
  X <- as.matrix(subset(full.rank, select=-c(Response)))
  y <- as.numeric(data$train$Response)-1
  Xy <- as.data.frame(cbind(X,y))
  names(Xy) <- c(paste("X", 1:ncol(X) ,sep=""), "y")
  model <- bestglm(Xy, family=binomial)
  predi <- as.factor(round(predict(model, newdata=data$test, type="response")))
  levels(predi) <- c("active", "scrambled")
  compareModeltoTruth(predi, data$test$Response)
}

# get data
mtor.bruc.du <- findPlateWellsFromGene("brucella-du-k", "MTOR")
scra.bruc.du <- findPlateWellsFromGene("brucella-du-k", "SCRAMBLED")

targets <- rbind(mtor.bruc.du[1:6,], scra.bruc.du[c(6,18,78,90,150,162),])

data <- getPlateWells(targets)

# prepare input for same plate analysis
j101_2c.mtor <- data[["J101-2C"]]$H6$data
j101_2c.scra <- data[["J101-2C"]]$G23$data

j101_2c <- prepareData(j101_2c.mtor, j101_2c.scra)

# carry out same plate analysis: complete sepration problems, good predictions
# probably mostly due to well effects
glm.regular(j101_2c)
glm.glm2(j101_2c)
glm.safe(j101_2c)
glm.brglm(j101_2c)
glm.glmnet(j101_2c)

# prepare input for same replication analysis
j10x_2c.mtor <- unlist(list(data[["J101-2C"]]$H6$data, 
                            data[["J104-2C"]]$H6$data,
                            data[["J107-2C"]]$H6$data), recursive=FALSE)
j10x_2c.scra <- unlist(list(data[["J101-2C"]]$G23$data, 
                            data[["J104-2C"]]$G23$data,
                            data[["J107-2C"]]$H6$data), recursive=FALSE)
                     
j10x_2c <- prepareData(j10x_2c.mtor, j10x_2c.scra)

# carry out same replication analysis: complete sepration problems gone,
# mediocre prediction quality
glm.regular(j10x_2c)
glm.glmnet(j10x_2c)
glm.bestglm(j10x_2c)

# prepare input for all available data
j10x_2x.mtor <- unlist(list(data[["J101-2C"]]$H6$data, 
                            data[["J104-2C"]]$H6$data,
                            data[["J107-2C"]]$H6$data,
                            data[["J101-2D"]]$H6$data, 
                            data[["J104-2D"]]$H6$data,
                            data[["J107-2D"]]$H6$data), recursive=FALSE)
j10x_2x.scra <- unlist(list(data[["J101-2C"]]$G23$data, 
                            data[["J104-2C"]]$G23$data,
                            data[["J107-2C"]]$H6$data,
                            data[["J101-2D"]]$H6$data, 
                            data[["J104-2D"]]$H6$data,
                            data[["J107-2D"]]$H6$data), recursive=FALSE)

j10x_2x <- prepareData(j10x_2x.mtor, j10x_2x.scra)

# carry out full analysis: complete sepration problems gone, bad prediction
# quality
glm.regular(j10x_2x)
glm.glmnet(j10x_2x)


