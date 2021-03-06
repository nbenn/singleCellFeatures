% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getPathogen.R
\name{getPathogen}
\alias{getPathogen}
\title{Get the pathogen name of a PlateLocation/WellLocation object}
\usage{
getPathogen(x)
}
\arguments{
\item{x}{A PlateLocation/WellLocation object}
}
\value{
The pathogen name of the PlateLocation/WellLocation object (as a
        string)
}
\description{
Give a PlateLocation/WellLocation object, return the pathogen name
}
\examples{
plate <- PlateLocation("J101-2C")
pathogen <- getPathogen(plate)

}

