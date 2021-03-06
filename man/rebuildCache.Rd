% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rebuildCache.R
\name{rebuildCache}
\alias{rebuildCache}
\title{Rebuild cache files for certain Data objects}
\usage{
rebuildCache(x, ...)
}
\arguments{
\item{x}{MatData object}

\item{location}{An optional argument for the location of the new .mat files.
If not specified, the /Errors directory is checked and if it
does not exist, the missing features are detected and
downloaded.}

\item{names}{An optional argument for the name of the newly added
feature. If not specified, the filename (without extension)
will be used. It is ignored when the missing features are
detected automatically.}
}
\value{
The updated MatData object. The updated cache file is saved to its
        corresponding location.
}
\description{
Given an existing MatData cache file and some feature files that were not
imported correctly, update and save it to the corresponding cache location.
}
\examples{
dat <- MatData(PlateLocation("J110-2C"))
dat <- rebuildCache(dat)

}

