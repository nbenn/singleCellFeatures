% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCacheFilenameData.R
\name{getCacheFilenameData}
\alias{getCacheFilenameData}
\title{Get path/filename of a PlateLocation/WellLocation object single cell data
cache}
\usage{
getCacheFilenameData(x)
}
\arguments{
\item{x}{PlateLocation/WellLocation object}
}
\value{
The path/filename of the single cell data cache of
        PlateLocation/WellLocation object (as a string)
}
\description{
Upon downloading from openBIS and importing to R, all single cell data
associated with a PlateLocation/WellLocation is cached in an .rds file. This
function retrieves the filepath of this cache file.
}
\examples{
plate <- PlateLocation("J101-2C")
path <- getCacheFilenameData(plate)

}

