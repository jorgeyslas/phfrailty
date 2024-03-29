% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4.bivph.R
\name{moment,bivphasetype-method}
\alias{moment,bivphasetype-method}
\title{Moment method for bivphasetype class}
\usage{
\S4method{moment}{bivphasetype}(x, k = c(1, 1))
}
\arguments{
\item{x}{An object of class \linkS4class{bivphasetype}.}

\item{k}{A vector with the location.}
}
\value{
An real value.
}
\description{
Moment method for bivphasetype class
}
\examples{
obj <- bivphasetype(dimensions = c(3, 3))
moment(obj, c(1, 1))
}
