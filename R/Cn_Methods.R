#' Find all descendents of a pedigree member
#'
#' @param id Numeric. The ID of the member whose offspring we would like to identify
#' @param parents The parent matrix returned by ProcessPed
#'
#' @return a list of descendent IDs
#' @keywords internal
#'
findDescendents <- function(id, parents){


  #find the offspring of the pedigree member with id
  offspring <- descendents <- findOffspring(id, parents)

  cont = TRUE
  while (cont == TRUE) {
    offspring <- unlist(lapply(offspring, function(x){ findOffspring(x, parents) }))

    if (length(offspring) == 0){
      cont = FALSE
    } else {
      descendents <- c(descendents, offspring)
    }

  }

  return(descendents)

}


#' Find offspring of a pedigree member
#'
#' @param id Numeric. The ID of the member whose offspring we would like to identify
#' @param parents The parent matrix returned by ProcessPed
#'
#' @return a list of offspring IDs
#' @keywords internal
findOffspring <- function(id, parents){
  offspring_index <- which(apply(parents, 2, function(x){id %in% x}))

  # sex <- ped$sex[ped$ID == id]
  #
  # if (sex == 0) {
  #   offspring_ids <- ped$ID[which(ped$dadID == id)]
  # } else {
  #   offspring_ids <- ped$ID[which(ped$momID == id)]
  # }

  return(offspring_index)
}

#' Retrieve the Founder Weights
#'
#' @param procPed The information returned by ProcPed
#' @param subtype_weights A vector of length 2.  The weights for individuals with subtypes A and B, respectively.  By default, \code{subtype_weights = NULL} so that no founder weights are applied
#'
#' @return A vector of founder weights
#' @keywords internal
#'
get_founderWeights <- function(procPed, subtype_weights = NULL){

  if (is.null(subtype_weights)){
    #compute equal weights for founders
    f_weight <- rep(1/length(procPed$founders), length(procPed$founders))

  } else {
    #For each founder: find the number of descendants affected by disease-subtype A
    f_countA <- sapply(procPed$founders, function(x){
      length(intersect(findDescendents(x, procPed$parents), procPed$typeA))
    })

    #For each founder: find the number of descendants affected by disease-subtype B
    f_countB <- sapply(procPed$founders, function(x){
      length(intersect(findDescendents(x, procPed$parents), procPed$typeB))
    })

    #weight and normalize
    f_weight <- (f_countA*subtype_weights[1] + f_countB*subtype_weights[2])/sum(f_countA*subtype_weights[1] + f_countB*subtype_weights[2])

  }

  return(f_weight)
}



#' numerator of sharing probability
#' @keywords internal
#'
#' @description calculates the numerator of the sharing probability
#' outline in section 2.1 of Bureau et al.
#' @param gRain bayesian network
#' @param procPed pedigree object that has been process with processPedigree
#' @return numerator value
numerProb <- function(net, procPed)
{
  #create a list of RVstates (?) for each carrier.  The RVstates take on
  #values 1 or 2; i.e. heterozygous and homozygous carriers??
  rvInCarriers <- sapply(simplify=FALSE, FUN=function(dummy) 1:2,
                         X=as.character(procPed$carriers))

  #For the non-carriers, set each RVstate to 0.  The non-carriers are determined
  #by the set difference of the carriers and the affecteds. If all affecteds
  #are also carriers this will return an empty named list ???
  noRvInNonCarriers <- sapply(simplify=FALSE, FUN=function(dummy) 0,
                              X=as.character(setdiff(procPed$affected, procPed$carriers)))

  return(marginalProb(net, c(rvInCarriers, noRvInNonCarriers)))
}

#' Calculate the probability of an observed sharing configuration given specified transmission probabilities.
#' @param ped An object of class \code{ped}.
#' @param carriers subjects in pedigree share the variant
#' @param tau A list of two transmission probabilites, with the transmission probability for the genetically-complelling disease-subtype listed first. These represent the transmission probability from a heterozygous parent to an offspring.
#' @param subtypes  A list of two subtypes that identify the subtype ranking.  The label for the more the genetically-complelling subtype is listed first.
#' @return The transmssion-based probability of the sharing configuration for the transmission probabilities specified in the argument \code{tau}.
#' @importFrom gRain setEvidence
#' @importFrom gRain retractEvidence
#' @keywords internal
compute_sharingProb <- function(ped, subtypes, tau,
                                carriers = NULL,
                                subtype_weights = NULL){

  #process the pedigree
  procPed <- processPed(ped, carriers, subtypes)

  # #define carriers
  # if (!is.null(carriers)) {
  #   #identify indices of carriers, when provided.
  #   procPed$carriers <- which(procPed$origID %in% carriers)
  # } else {
  #   #otherwise default to disease-affected relatives
  #   procPed$carriers <- procPed$affected
  # }


  net <- createNetworkTau(procPed, tau)
  net <- gRain::setEvidence(net, as.character(procPed$founders),
                            rep('0', length(procPed$founders)))
  #NOTE: after we execute the command above, all founders are set to
  #have RV status 0.  We consider the possiblility that each founder
  #introduces the RV (one at a time) in the for loop below.

  # sum over probs, conditioning on each founder introducing variant
  numer <- 0

  #get founder weights
  fW <- get_founderWeights(procPed, subtype_weights)
  for (f in 1:length(procPed$founders)) { #TODO: use sapply here
    # condition on founder and calculate distribution
    condNet <- gRain::retractEvidence(net, as.character(procPed$founders[f]))
    condNet <- gRain::setEvidence(condNet, as.character(procPed$founders[f]), '1')

    # compute probability
    #denom <- denom + denomProb(condNet, procPed)
    numer <- numer + numerProb(condNet, procPed)*fW[f]
  }

  return(numer)
}
