## Load functions
library(singleCellFeatures)

## Define functions
compareModeltoTruth <- function(estim, truth) {
  tp <- sum(estim==truth)
  fp <- sum(estim!=truth)
  al <- length(estim==truth)
  cat("true positive rate: ", tp/al, "\nfalse positive rate:", fp/al)
}

glm.regular <- function(data) {
  model <- glm("Response ~ .", binomial, data$train)
  predi <- as.factor(round(predict(model, newdata=data$test, type="response")))
  levels(predi) <- c("active", "control")
  compareModeltoTruth(predi, data$test$Response)
}

glm.glm2 <- function(data) {
  library(glm2)
  design.mat <- makeRankFull(data$train)
  model <- glm2("Response ~ .", binomial, design.mat)
  predi <- as.factor(round(predict(model, newdata=data$test, type="response")))
  levels(predi) <- c("active", "control")
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
    levels(predi) <- c("active", "control")
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
  levels(predi) <- c("active", "control")
  compareModeltoTruth(predi, data$test$Response)
}

glm.glmnet <- function(data) {
  library(glmnet)
  x <- as.matrix(subset(data$train, select=-c(Response)))
  y <- data$train$Response
  model <- glmnet(x, y, family="binomial")
  x.test <- as.matrix(subset(data$test, select=-c(Response)))
  predi <- as.factor(round(predict(model, newx=x.test, type="response")))
  levels(predi) <- c("active", "control")
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
  levels(predi) <- c("active", "control")
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

j101_2c <- prepareDataforGlm(j101_2c.mtor, j101_2c.scra)

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
                     
j10x_2c <- prepareDataforGlm(j10x_2c.mtor, j10x_2c.scra)

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

j10x_2x <- prepareDataforGlm(j10x_2x.mtor, j10x_2x.scra)

# carry out full analysis: complete sepration problems gone, bad prediction
# quality
glm.regular(j10x_2x)
glm.glmnet(j10x_2x)


