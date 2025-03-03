% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FamAnalysisFunctions.R
\name{compute_famWeights_internal}
\alias{compute_famWeights_internal}
\title{Compute sharing weights for a pedigree given the results of test_allcombos}
\usage{
compute_famWeights_internal(famLike_results, tau_grid, ped, subtypes,
  subtype_weights = NULL, identify_trueConfig = FALSE)
}
\arguments{
\item{famLike_results}{The output returned by the \code{\link{test_allcombos}} function.}

\item{ped}{An object of class \code{ped}. A pedigree generated by \code{sim_ped} or \code{sim_RVped}, or an object created by the function \code{new.ped}, for details please see documantation and vignette for package \code{SimRVPedigree}.}

\item{subtypes}{A list of length 2. Contains character labels for the two subtypes that occur in the pedigrees contained in the \code{famStudy_obj}. The label for the more the genetically-complelling subtype must be listed first.}

\item{subtype_weights}{A vector of length 2.  When assuming an informative prior for the founders; the weights the individuals with subtypes A and B, respectively.  By default, \code{subtype_weights = NULL} so that no founder weights are applied (i.e. flat founder prior).}

\item{identify_trueConfig}{Logical. Identify true sharing configuration when \code{ped} includes familial genotypes at the disease locus. By default \code{identify_trueConfigs = FALSE}.}
}
\value{
A data frame containing the weights for each configuration in the pedigree. See Details.
}
\description{
Compute sharing weights for a pedigree given the results of test_allcombos
}
\keyword{internal}
