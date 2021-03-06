% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkCompletenessWell.R
\name{checkCompletenessWell}
\alias{checkCompletenessWell}
\title{Check if the set of wells is complete}
\usage{
checkCompletenessWell(...)
}
\arguments{
\item{x}{The PlateData object of interest.}
}
\value{
Logical: TRUE if complete, FALSE if incomplete.
}
\description{
Given a PlateData object, check if the set of wells is complete.
}
\examples{
plate <- PlateLocation("J101-2C")
checkCompletenessWell(plate)

}

