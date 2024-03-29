% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ph_density}
\alias{ph_density}
\title{Phase-type density}
\usage{
ph_density(x, alpha, S)
}
\arguments{
\item{x}{Non-negative values.}

\item{alpha}{Vector of initial probabilities.}

\item{S}{Sub-intensity matrix.}
}
\value{
The density at \code{x}.
}
\description{
Computes the density of a phase-type distribution with parameters
\code{alpha} and \code{S} at \code{x}.
}
