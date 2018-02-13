% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findHotspots.R
\name{findPVals}
\alias{findPVals}
\title{Calculate a vector of p-values expressing the difference between two
coverages.}
\usage{
findPVals(x, y, wins = 10000)
}
\arguments{
\item{x}{First coverage.}

\item{y}{Second coverage.}

\item{wins}{Size of the window.}
}
\description{
Calculate a vector of p-values expressing the difference between two
coverages. Works by windows.
}
\author{
Ricard Illa \email{ricard.illa@irbbarcelona.org}
}
\keyword{manip}