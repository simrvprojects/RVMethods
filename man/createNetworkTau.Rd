% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rgain_setup_functions.R
\name{createNetworkTau}
\alias{createNetworkTau}
\title{create bayesian network from processed pedigree}
\usage{
createNetworkTau(procPed, tau, prior = c(1, 2, 1))
}
\arguments{
\item{procPed}{processed Pedigree object}

\item{prior}{prior on number of alleles for founders}
}
\value{
bayesian network from gRain package
}
\description{
create bayesian network from processed pedigree
}
\details{
Creates a bayesian network using the gRain package.
The network is built based on the information in a pedigree object
that has been processed using \code{processPedigree}; as well as tau,
the probability that a heterozygous parent transmits a cRV
}
\keyword{internal}
