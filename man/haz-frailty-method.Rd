% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2.phfrailty.R
\name{haz,frailty-method}
\alias{haz,frailty-method}
\title{Hazard rate method for univariate phase type frailty models}
\usage{
\S4method{haz}{frailty}(x, y, X = numeric(0))
}
\arguments{
\item{x}{An object of class \linkS4class{frailty}.}

\item{y}{A vector of locations.}

\item{X}{A matrix of covariates.}
}
\value{
A list containing the locations and corresponding hazard rate evaluations.
}
\description{
Hazard rate method for univariate phase type frailty models
}
\examples{
obj <- frailty(phasetype(structure = "general"), bhaz = "weibull", bhaz_pars = 2)
haz(obj, c(1, 2, 3))
}
