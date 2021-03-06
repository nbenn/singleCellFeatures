% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareDataforGlm.R
\name{compareModeltoTruth}
\alias{compareModeltoTruth}
\title{Compare result of binary classification to truth}
\usage{
compareModeltoTruth(estim, truth)
}
\arguments{
\item{estim}{Vector holding a two level factor with the classification
result.}

\item{truth}{Vector holding a two level factor with the truth.}
}
\value{
A list with various key charactereistics of the resulting confusion
        matrix
}
\description{
Given two vectors of factors (two identical levels) of equal length,
calculate confusion matrix marginals
}
\examples{
# get gene locations
mtor.loc <- findWells(experiments="brucella-du-k1", contents="MTOR")
scra.loc <- findWells(plates=sapply(mtor.loc, getBarcode),
                      contents="SCRAMBLED", well.names="G23")
# combine for faster fetching
data     <- getSingleCellData(list(mtor.loc[[1]], scra.loc[[1]]))
mtor.dat <- meltData(cleanData(data[[1]]$H6))
scra.dat <- meltData(cleanData(data[[1]]$G23))
# prepare data for glm
data <- prepareDataforGlm(mtor.dat$mat$Cells, scra.dat$mat$Cells)
data <- makeRankFull(data)
# run glm
model <- glm("Response ~ .", binomial, data$train)
# compare to testing data
predi <- as.factor(round(predict(model, newdata=data$test,
                                 type="response")))
levels(predi) <- c("active", "control")
comparison <- compareModeltoTruth(predi, data$test$Response)

}

