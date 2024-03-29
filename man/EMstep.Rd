% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{EMstep}
\alias{EMstep}
\title{EM for phase-type distributions using Pade for matrix exponential}
\usage{
EMstep(alpha, S, obs, weight)
}
\arguments{
\item{alpha}{Initial probabilities.}

\item{S}{Sub-intensity.}

\item{obs}{The observations.}

\item{weight}{The weights for the observations.}
}
\value{
Fitted alpha and S after one iteration.
}
\description{
EM for phase-type distributions using Pade for matrix exponential
}
