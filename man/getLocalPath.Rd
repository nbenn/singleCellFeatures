% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLocalPath.R
\name{getLocalPath}
\alias{getLocalPath}
\title{Get the local path of a PlateLocation/WellLocation object}
\usage{
getLocalPath(x)
}
\arguments{
\item{x}{A PlateLocation/WellLocation object}
}
\value{
The local path of the PlateLocation/WellLocation object (as a
        string)
}
\description{
Given a PlateLocation/WellLocation object, return the local path
}
\examples{
plate <- PlateLocation("J101-2C")
path <- getLocalPath(plate)

}

