#' create bayesian network from processed pedigree
#'
#' create bayesian network from processed pedigree
#'
#' Creates a bayesian network using the gRain package.
#' The network is built based on the information in a pedigree object
#' that has been processed using \code{processPedigree}; as well as tau,
#' the probability that a heterozygous parent transmits a cRV
#'
#' @param procPed processed Pedigree object
#' @param prior prior on number of alleles for founders
#' @importFrom gRain cptable
#' @importFrom gRain grain
#' @return bayesian network from gRain package
#' @keywords internal
createNetworkTau <- function(procPed, tau, prior = c(1, 2, 1)) {
  if (length(tau) == 1) tau <- rep(tau, 2)
  # process founders
  # Create a cptable for each founder.  Since founders do not have any ancestors
  # these cptables depend only on the this creates a cptable for the offspring node
  # given the parent nodes.
  founderNodes <- lapply(procPed$founders, gRain::cptable,
                         values = prior, levels = 0:2)


  # process non-founders
  # For each non-founder, this creates a cptable for the offspring node
  # given the parent nodes.

  nonFounderNodesA <- lapply(procPed$typeA[!procPed$typeA %in% procPed$founders], function(nf) {
    gRain::cptable(c(nf, procPed$parents[1,nf], procPed$parents[2,nf]),
                   values = TauProbTable(tau[1]), levels=0:2)
  })

  nonFounderNodesB <- lapply(procPed$typeB[!procPed$typeB %in% procPed$founders], function(nf) {
    gRain::cptable(c(nf, procPed$parents[1,nf], procPed$parents[2,nf]),
                   values = TauProbTable(tau[2]), levels=0:2)
  })

  # create bayesian network
  condProbTable <- gRain::compileCPT(c(founderNodes, nonFounderNodesA, nonFounderNodesB))
  return(gRain::grain(condProbTable))
}

#' Table of inheritance probabilities that depend on tau
#'
#' @param tau The probability that an offspring inherits the RV from a parent who is a carrier.
#'
#' @return  An array of dimension c(3, 3, 3).
#' @keywords internal
TauProbTable <- function(tau){
  ProbTab <- array(0, c(3,3,3))
  ProbTab[1,,] <- c(1, 1-tau, 0,
                    1-tau, (1-tau)^2, 0,
                    0, 0, 0) # 0 copies in offspring
  ProbTab[2,,] <- c(0, tau, 1,
                    tau, 2*tau*(1-tau), 1-tau,
                    1, 1-tau, 0) # 1 copy in offspring
  ProbTab[3,,] <-   c(0, 0, 0,
                      0, tau^2, tau,
                      0, tau, 1)       # 2 copies in offspring

  return(ProbTab)
}


check_tau <- function(tau){
  if(tau < 0 | tau > 1){
    stop("ERROR: tau < 0 | tau > 1 \n  The transmission probability must be between 0 and 1.")
  } else if(tau < 0.5){
    stop("WARNING: tau < 0.5 is used for protective variants; that is, variants that are enriched among unaffected relatives.")
  }

}


#' calculates the marginal probability of a set of nodes
#'
#' @description Given a bayesian network from the gRain package and a
#' named list of (nodes, states), this function returns the joint-marginal
#' probability of each node taking a value in the specified set of states.
#' @details This function calculates the probability P(A,B,C) by factoring
#' it into conditional probabilities, i.e. P(A|B,C) * P(B|C) * P(C).
#' Starting at the right side, P(C) is computed and then evidence of C
#' being true is added to the network and P(B) is computed - effectively
#' giving the probability P(B|C). This process continues from right to
#' left until the entire product has been computed.
#' @param net bayesian network from gRain package
#' @param states named list of states for each node
#' @return joint-marginal probability
#'
#' @keywords internal
#' @importFrom gRain querygrain
#' @importFrom gRain setEvidence
marginalProb <- function(net, states)
{
  prob <- 1
  for (n in names(states))
  {
    if (prob > 0) # prevents conditioning on zero prob events
    {
      # calculate probability for this node
      p <- unname(gRain::querygrain(net, n, exclude=FALSE)[[1]])
      prob <- prob * sum(p[states[[n]]+1])

      # condition on this node being in the correct states
      net <- gRain::setEvidence(net, evidence=sapply(simplify=FALSE,
                                                     X=n, FUN=function(x) as.numeric(0:2 %in% states[[n]])))
    }
  }
  return(prob)
}
