% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GlobalConfigFunctions.R
\name{compute_globalWeights}
\alias{compute_globalWeights}
\title{Compute Global weights for a set of families}
\usage{
compute_globalWeights(peds, subtypes, tau_increment = 0.05,
  subtype_weights = NULL, identify_trueConfig = FALSE)
}
\arguments{
\item{peds}{A ped object containing two or more pedigrees.}

\item{subtypes}{A list of length 2. Contains character labels for the two subtypes that occur in the pedigrees contained in the \code{famStudy_obj}. The label for the more the genetically-complelling subtype must be listed first.}

\item{tau_increment}{Numeric. The width of the grid for the taus}

\item{subtype_weights}{A vector of length 2.  When assuming an informative prior for the founders; the weights the individuals with subtypes A and B, respectively.  By default, \code{subtype_weights = NULL} so that no founder weights are applied (i.e. flat founder prior).}

\item{identify_trueConfig}{Logical. Identify true sharing configuration when \code{ped} includes familial genotypes at the disease locus. By default \code{identify_trueConfigs = FALSE}.}
}
\value{
A data frame of global configurations anf weights for the set of provided families.
}
\description{
Compute Global weights for a set of families
}
\keyword{internal}
