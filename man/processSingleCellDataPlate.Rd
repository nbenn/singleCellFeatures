% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSingleCellData.R
\name{processSingleCellDataPlate}
\alias{processSingleCellDataPlate}
\title{Read all single cell features for wells/a plate}
\usage{
processSingleCellDataPlate(locations, select = NULL, drop = NULL,
  images = NULL)
}
\arguments{
\item{locations}{A list of PlateLocation/WellLocation objects for which the
data is to be returned. They all have to lie on the same
plate.}

\item{select}{A vector of strings that are matched with features to be
kept.}

\item{drop}{A vector of strings that are matched with features to be
dropped.}

\item{images}{A vector of image indices (1:9) indicating which images to
return for each well.}
}
\value{
A list with well/plate data
}
\description{
Read all available single cell feature data for multiple wells and/or the
whole plate.
}
\examples{
plate <- PlateLocation("J101-2C")
well1 <- WellLocation("J101-2C", "A", 24)
well2 <- WellLocation("J101-2C", "B", 12)
well3 <- WellLocation("J101-2C", "P", 24)

locations <- list(well1, well2, well3, plate)

data <- processSingleCellDataPlate(locations)

}

