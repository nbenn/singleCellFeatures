% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilPackage.R
\name{size}
\alias{size}
\title{Get amount of storage used by an object}
\usage{
size(x)
}
\arguments{
\item{x}{The object of interest}
}
\value{
The amount of storage taken up by the object as a string
}
\description{
Prints the amount of storage used by the specified object.
}
\examples{
plate <- WellLocation("J101-2C", "B", 15)
size(plate)

}

