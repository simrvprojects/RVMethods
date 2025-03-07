
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RVMethods

RVMethods provides routines to compute statistics for priortization of
RVs observed in family-based studies. This software is considered a work
in progress.

## Installation

You can install RVMethods from github with:

``` r
# install.packages("devtools")
devtools::install_github("cnieuwoudt/RVMethods")
```

## Example Usage

``` r
library(RVMethods)
#load example founder data
data("pop_fdat")

#load example pedigrees
data("study_pedigrees")

#load SimRVPedigree - to plot pedigrees. 
library(SimRVPedigree)
#pedigree number 58
plot(study_pedigrees[study_pedigrees$FamID==58, ])
#pedigree number 116
plot(study_pedigrees[study_pedigrees$FamID==116, ])

#compute distributions 
FBSdists = compute_distributions(peds = study_pedigrees[study_pedigrees$FamID %in% c(58, 116), ],
                                 subtypes = c("HL", "NHL"),
                                 carrier_prob = 0.002)
    
#View first few rows global likelihood ratio statistics over all possible configurations
head(FBSdists$GlobalDist)

#View first few rows global transmission statistics over all possible configurations
head(FBSdists$GlobalTransDist)

#View first few rows of the local likelihood ratio, RVS, and modified RVS statistics over all possible conifgurations
head(FBSdists$LocalDist)
```
