#' Compute sharing probabilities for every configuration of disease-affected relatives and initialization of tau
#'
#' Compute sharing probabilities for every configuration of disease-affected relatives and initialization of tau
#'
#' @param tau_grid A grid of tau values produced by \code{expand.grid}.  \code{test_allcombos} will compute each sharing probability at every specified value of tau in \code{tau_grid}.
#' @inheritParams compute_familyWeights
#'
#' @return  A list containing the following:
#' @return \item{\code{like_mat} }{A matrix of values that represent the likelihood given a configuration (column) and relized value of tau (row).}
#' @return \item{\code{configs} }{A matrix of sharing configurations.  The columns are named for the IDs of the disease-affected relatives in the pedigree. When \code{TRUE} the individual is a carrier of the cRV, when \code{FALSE} the individual is not a carrier.}
#' @keywords internal
#' @examples
#' library(RVMethods)
#' data(study_pedigrees)
#'
#' library(SimRVPedigree)
#' plot(study_pedigrees[study_pedigrees$FamID == 58, ])
#' test_allcombos(ped = study_pedigrees[study_pedigrees$FamID == 58, ],
#'                subtypes = c("HL", "NHL"),
#'                tau_grid = make_tauGrid())
#'
test_allcombos <- function(ped, subtypes, tau_grid,
                           subtype_weights = NULL){

  #determine the total number of affected in this pedigree
  num_affected <- sum(ped$affected, na.rm = TRUE)

  #store the IDs of the affected relatives
  aff_Ids <- ped$ID[ped$affected & ped$available]

  #determine all possible sharing configurations
  config = as.matrix(expand.grid(rep(list(c(T, F)), num_affected)))
  config = config[rowSums(config) > 0, ]
  colnames(config) = as.character(aff_Ids)

  # compute and store the likelihood results in a matrix
  # ROW: realizations of tau
  # COLUMNS: different sharing configurations
  like_mat <- lapply(1:nrow(config), function(i){
    sapply(1:nrow(tau_grid), function(x){
      compute_sharingProb(ped, subtypes, tau = as.numeric(tau_grid[x, ]),
                          carriers = aff_Ids[as.vector(config[i, ])],
                          subtype_weights)})
  })
  like_mat <- do.call(cbind, like_mat)


  #compute the binID, this is a numeric identifier for the configuration.
  #The configuration is interpreted as a base 2 number and then converted to
  #base 10.  I.e. the configuration (0, 1, 0) in base 2 --> 010
  #the binID is the base 2 converted to base 10 --> 2.
  #
  #For configuration (0, 0, 1) ---> binID = 1
  #For configuration (0, 1, 0) ---> binID = 2
  #For configuration (0, 1, 1) ---> binID = 3, etc.
  binIDs <- as.numeric(apply(config*1, 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  }))

  #The binID will come in handy to quickly refer to and find
  #configurations, especially when dealing with
  #global configurations.

  return(list(like_mat = like_mat, configs = config, binIDs = binIDs))
}


#' Compute sharing weights for a pedigree given the results of test_allcombos
#'
#' @param famLike_results The output returned by the \code{\link{test_allcombos}} function.
#' @inheritParams compute_familyWeights
#'
#' @return A data frame containing the weights for each configuration in the pedigree. See Details.
#' @keywords internal
#'
compute_famWeights_internal <- function(famLike_results, tau_grid, ped,
                               subtypes, subtype_weights = NULL,
                               identify_trueConfig = FALSE){

  # find the values of tau that maximize the likelihood for each configuration
  MLEs <- get_MLE(famLike_results$like_mat, tau_grid)

  #create a data.frame that combines all of the information
  famWeights <- as.data.frame(famLike_results$configs*1)
  famWeights$binID <- famLike_results$binIDs

  #calculate K, the number of carriers in the configuration
  famWeights$K <- rowSums(famLike_results$configs)

  #find the column IDs of the relatives with
  #the genetically-complelling subtypes
  GC_colID <- which(colnames(famWeights) %in% ped$ID[which(ped$subtype == subtypes[[1]])])
  if(length(GC_colID) > 1){
    #calculate K_sub, the number of carriers in the configuration
    #who have the genetically compelling subtype
    famWeights$K_sub <- rowSums(famWeights[, GC_colID])
  } else {
    famWeights$K_sub <- famWeights[, GC_colID]
  }

  # find the maximum value of the likelihood
  famWeights$L_max <- apply(famLike_results$like_mat, 2, function(x){
    x[which.max(x)]
  })

  # value of the likelihood under H_0
  # NOTE: this is now the first row in the likelihood matrix
  # YOU HARDCODED THE NULL FOR THE FIRST ROW
  # See function make_tauGrid() if you forget.
  famWeights$L_null <- famLike_results$like_mat[1, ]

  #compute the composite LR statistic
  famWeights$LR <- famWeights$L_max/famWeights$L_null


  # Find the MLEs.
  # NOTE: they may not be unique, HOWEVER, since you ordered
  # tau_grid from smallest to largest (see make_tauGrid), choosing
  # the values in the first row chooses the smallest values of
  # tau that maximixe the likelihood.
  # i.e. values that tend toward the null in case of ties.
  famWeights$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  famWeights$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  # For certain sharing configurations and subtype assinments, the likelihood
  # may have multiple maxima, resulting in non-unique MLEs
  # We identify these in the output.
  famWeights$unique_tau <- sapply(MLEs, function(x){nrow(x) == 1})

  #When requested, via the argument identify_trueConfig, and when the data is available
  #identify the true configuration in the supplied pedigree
  if(identify_trueConfig){
    if (all(c("DA1", "DA2") %in% colnames(ped))) {
      actual_config <- as.numeric(rowSums(ped[ped$affected & ped$available, c("DA1", "DA2")]))

      famWeights$true_config <- sapply(1:nrow(famWeights), function(x){
        all(actual_config == as.numeric(famWeights[x, c(1:ncol(famLike_results$configs))]))
      })
    } else {
      warning("Cannot identify the true configuration at the DSL.
              Variables 'DA1' and 'DA2' are missing from argument ped.")
    }
  }


  #compute the numerator of the RVsharing probablity
  #NOTE this is the first row of like_mat in the argument famLike_results
  #i.e. the likelihood when tau_A = tau_B = 0.5
  RVSres <- famLike_results$like_mat[1, ]

  #compute the RVsharing probability
  famWeights$RVsharing = as.vector(RVSres/sum(RVSres))


  #compute the likelihood-based weight
  famWeights$w_LR <- sapply(1:nrow(famWeights), function(x){
    compute_transmission_weight(all_stats = famWeights, config_index = x)
  })

  #calculate the modified RVS weight
  famWeights$w_LRbasic <- sapply(1:nrow(famWeights), function(x){
    compute_LRB_weight(all_stats = famWeights, config_index = x)
  })

  #compute the RVS weight
  famWeights$w_RVS <- sapply(1:nrow(famWeights), function(x){
    compute_RVS_weight(all_stats = famWeights, config_index = x)
  })


  #calculate the modified RVS weight
  famWeights$w_RVS2 <- sapply(1:nrow(famWeights), function(x){
    compute_modified_RVS_weight(all_stats = famWeights, config_index = x)
  })


  return(famWeights)
}


#' Compute all weights for a pedigree
#'
#' Compute all weights for a pedigree
#'
#'
#' Note: the parameter space of  \eqn{\tau_a} and \eqn{\tau_b} is given by
#' \eqn{\Omega_{con} = ((\tau_{a}, \tau_{b}): 0.5 \le \tau_{b} < \tau_{a} \le 1 )}.
#'
#' The data frame returned by \code{compute_familyWeights} contains the following:
#' Sharing configuration are identified in the first $n$ columns (i.e. for the $n$ disease-affected relatives), by the identification number of the disease-affected realtive in the pedigree. Following the first $n$ columns the remaining variables are described as follows.
#'
#' \tabular{ll}{
#' \strong{name} \tab \strong{description} \cr
#' \code{binID} \tab This is a numeric identifier for the configuration. Computed by considering the \cr
#'              \tab configuration as a number in base 2 ands then converting from base 2 to base 10. \cr
#'              \tab For example, the configuration (1, 0, 1) would have binID = 5 since 101 converted \cr
#'              \tab from base 2 to base 10 is 5. \cr
#' \code{K}  \tab The total number of disease-affected relatives who carry a variant in the configuration.  \cr
#' \code{K_sub}  \tab the number of disease-affected relatives who carry a variant and who are affected \cr
#'               \tab by the genetically-compelling subtype in the specified configuration \cr
#' \code{L_max}  \tab The maximum value of the likelihood \cr
#' \code{L_null}  \tab The likelihood under the null \cr
#' \code{LR}  \tab The likelihood ratio statistic\cr
#' \code{tau_A}  \tab \eqn{\hat{\tau}_a}, the MLE of \eqn{\tau_a}\cr
#' \code{tau_A}  \tab \eqn{\hat{\tau}_b}, the MLE of \eqn{\tau_b} \cr
#' \code{unique_tau}  \tab Does the likelihood for the configuration have a unique maximum?  \cr
#' \code{RVsharing}  \tab The RVS probability for the configuration \cr
#' \code{w_LR}  \tab The transmission-based weight for the specified configuration \cr
#' \code{w_RVS}  \tab The unmodified RVS weight for the specified configuration \cr
#' \code{w_RVS2}  \tab The modified RVS weight for the specified configuration \cr
#' }
#'
#' @inheritParams compute_studyWeights
#' @param ped An object of class \code{ped}. A pedigree generated by \code{sim_ped} or \code{sim_RVped}, or an object created by the function \code{new.ped}, for details please see documantation and vignette for package \code{SimRVPedigree}.
#' @param subtype_weights A vector of length 2.  When assuming an informative prior for the founders; the weights the individuals with subtypes A and B, respectively.  By default, \code{subtype_weights = NULL} so that no founder weights are applied (i.e. flat founder prior).
#' @param identify_trueConfig Logical. Identify true sharing configuration when \code{ped} includes familial genotypes at the disease locus. By default \code{identify_trueConfigs = FALSE}.
#'
#'
#'
#'
#' @return A data frame containing the weights for each configuration in the pedigree. See Details.
#' @keywords internal
#'
#' @examples
#' library(SimRVPedigree)
#' # for more information about specifying pedigrees
#' # execute help(new.ped) in the console
#'
#' fam133 = data.frame(FamID = rep(133, 11),
#'                     ID = seq(1, 11, by  = 1),
#'                     dadID = c(NA, NA, NA, 1, 1, 1, NA, 3, 3, 3, 6),
#'                     momID = c(NA, NA, NA, 2, 2, 2, NA, 4, 4, 4, 7),
#'                     sex = c(0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1),
#'                     affected = c(F, F, F, F, T, F, F, T, T, T, T),
#'                     DA1 =   c( 1,  0,  0, 1, 0, 1,  0, 0, 0, 0, 1),
#'                     DA2 =   c( 0,  0,  0, 0, 0, 0,  0, 1, 0, 1, 0),
#'                     Gen =   c( 1,  1,  2, 2, 2, 2, 2, 3, 3, 3, 3),
#'                     subtype = c(NA, NA, NA, NA, "NHL", NA, NA, "HL", "NHL", "HL", "HL"))
#' fam133 = new.ped(fam133)
#' plot(fam133, cex = 0.9)
#'
#' compute_familyWeights(ped = fam133,
#'                       subtypes = c("HL", "NHL"))
#'
compute_familyWeights <- function(ped, subtypes,
                                  tau_increment = 0.05,
                                  subtype_weights = NULL,
                                  identify_trueConfig = FALSE){
  #create tau grid
  #note removing constranined = false option, for now
  #will need to alter output method when non-unique tau occurs
  tau_grid <- make_tauGrid(increment_width = tau_increment, constrained = TRUE)
  #order taus to save time later.
  tau_grid <- tau_grid[order(tau_grid[, 1], tau_grid[, 2], decreasing = FALSE), ]
  row.names(tau_grid) = NULL

  famWeights <- compute_famWeights_internal(famLike_results = test_allcombos(ped, subtypes,
                                                                             tau_grid = tau_grid,
                                                                             subtype_weights),
                                            tau_grid = tau_grid, ped,
                                            subtypes, subtype_weights,
                                            identify_trueConfig)

  #remove invalid configurations
  famWeights <- famWeights[!is.nan(famWeights$LR), ]

  return(famWeights)
}
