% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LR_by_obsSNV_Functions.R
\name{condition_globalDist_zeroConfig2}
\alias{condition_globalDist_zeroConfig2}
\title{(Option 1) Conditon the probabilty of each global configuration on the event that the variant was observed in at least one study member.}
\usage{
condition_globalDist_zeroConfig2(likeGrids_byFam, famID_index, tau_grid,
  obs_bins = NULL)
}
\arguments{
\item{likeGrids_byFam}{A list of results from test_allcombos, one for each family.}

\item{famID_index}{The family IDs, indexed in the same order as likeGrids_byFam and weights_byFam}

\item{tau_grid}{The grid of tau values}

\item{weights_byFam}{A list of results from compute_famWeights, one for each family.}

\item{seq_option}{Logical. If false ass}
}
\value{
A data frame of global configurations anf weights for the set of provided families.
}
\description{
Compute MLE's for tau_a and tau_b under global configuration framework: here configurations that are identically zero contribute to the estimation of the MLEs; also, probabilities are conditioned on the variant being observed in at least one disease-affected study member.
}
\details{
This will be an internal function
}
\keyword{internal}
