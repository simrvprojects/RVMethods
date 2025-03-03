% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StudyAnalysisFunctions.R
\name{tabulate_binID_by_SNV}
\alias{tabulate_binID_by_SNV}
\title{Tabulate each families configurations stored in binary format}
\usage{
tabulate_binID_by_SNV(Ahaplos, Amap, FamIDs)
}
\arguments{
\item{Ahaplos}{the haplotype matrix for the disease-affected relatives only}

\item{Amap}{the pedigree to haplotype matrix map for the disease affected relatives only}

\item{FamIDs}{a list of FamIDs for pedigrees in the study}
}
\value{
A matrix of binIDs.  Rows are families, columns are SNVs, the $(i,j)^{th}$ entry is the binID for the sharing configuration for the ith family and the jth SNV.
}
\description{
Tabulate each families configurations stored in binary format
}
\keyword{internal}
