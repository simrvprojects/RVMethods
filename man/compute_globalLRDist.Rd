% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zeroConfig_Functions.R
\name{compute_globalLRDist}
\alias{compute_globalLRDist}
\title{Compute Global weights for a set of families}
\usage{
compute_globalLRDist(peds, subtypes, carrier_prob = 0.002,
  tau_increment = 0.05, subtype_weights = NULL)
}
\arguments{
\item{peds}{A ped object containing two or more pedigrees.}

\item{subtypes}{A list of length 2. Contains character labels for the two subtypes that occur in the pedigrees contained in the \code{famStudy_obj}. The label for the more the genetically-complelling subtype must be listed first.}

\item{carrier_prob}{The cumulative carrier probabilty of all crvs as a group.}

\item{tau_increment}{Numeric. The width of the grid for the taus}

\item{subtype_weights}{A vector of length 2.  When assuming an informative prior for the founders; the weights the individuals with subtypes A and B, respectively.  By default, \code{subtype_weights = NULL} so that no founder weights are applied (i.e. flat founder prior).}
}
\value{
A data frame of global configurations anf weights for the set of provided families.
}
\description{
Compute Global weights for a set of families
}
\examples{
library(RVMethods)
data(study_pedigrees)

library(SimRVPedigree)
plot(study_pedigrees[study_pedigrees$FamID == 304, ])

peds = study_pedigrees[study_pedigrees$FamID \%in\%c(58, 304), ]
subtypes = c("HL", "NHL")


myCombos = test_allcombos_includeZero(ped = study_pedigrees[study_pedigrees$FamID == 58, ],
                                      subtypes = c("HL", "NHL"),
                                      tau_grid = make_tauGrid(),
                                      carrier_prob = 0.002)

my_GDist = compute_globalLRDist(peds = study_pedigrees[study_pedigrees$FamID \%in\%c(58, 304), ],
                                subtypes = c("HL", "NHL"))

order(my_GDist$LR) == order(my_GDist$w_LR)

}
\keyword{internal}
