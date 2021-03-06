% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataMeta.R
\name{PlateMetadata}
\alias{PlateMetadata}
\title{Constructor for PlateMetadata objects}
\usage{
PlateMetadata(plate, quantiles = NULL)
}
\arguments{
\item{plate}{A PlateLocation object}

\item{quantiles}{A list containing the upper and lower 5% quantiles of cell
counts.}
}
\value{
A PlateMetadata object: a list with slots describing conditions on
        a plate.
}
\description{
Given a PlateLocation object, the corresponding metadata is compiled
}
\examples{
meta <- PlateMetadata(PlateLocation("J107-2C"))

}

