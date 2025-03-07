% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cn_Methods.R
\name{get_founderWeights}
\alias{get_founderWeights}
\title{Retrieve the Founder Weights}
\usage{
get_founderWeights(procPed, subtype_weights = NULL)
}
\arguments{
\item{procPed}{The information returned by ProcPed}

\item{subtype_weights}{A vector of length 2.  The weights for individuals with subtypes A and B, respectively.  By default, \code{subtype_weights = NULL} so that no founder weights are applied}
}
\value{
A vector of founder weights
}
\description{
Retrieve the Founder Weights
}
\keyword{internal}
