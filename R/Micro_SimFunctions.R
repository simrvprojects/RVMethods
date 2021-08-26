#' (Option 0) Compute sharing probabilities for a single global configuration (includes the zero configuration) and initialization of tau
#'
#' Compute sharing probabilities for every global configuration (includes the zero configuration) and initialization of tau.  The probability of the zero configuration will be used to condition on the event that the RV was observed in at least one study member.  This function requires an assummed value for the carrier probability of the RV
#'
#' @param tau_grid A grid of tau values produced by \code{expand.grid}.  \code{test_allcombos} will compute each sharing probability at every specified value of tau in \code{tau_grid}.
#' @param carrier_prob The cumulative carrier probabilty of all crvs as a group.
#' @inheritParams compute_familyWeights
#'
#' @return  A list containing the following:
#' @return \item{\code{like_mat} }{A matrix of values that represent the likelihood given a configuration (column) and relized value of tau (row).}
#' @return \item{\code{configs} }{A matrix of sharing configurations.  The columns are named for the IDs of the disease-affected relatives in the pedigree. When \code{TRUE} the individual is a carrier of the cRV, when \code{FALSE} the individual is not a carrier.}
#' @keywords internal
#'
test_onecombo_withZero <- function(ped, subtypes, pedcarriers, tau_grid, carrier_prob,
                                    subtype_weights = NULL){

  #determine the total number of affected in this pedigree
  num_affected <- sum(ped$affected, na.rm = TRUE)

  #store the IDs of the affected relatives
  aff_Ids <- ped$ID[ped$affected & ped$available]

  #determine all possible sharing configurations
  config = as.matrix(expand.grid(rep(list(c(T, F)), num_affected)))
  #config = config[rowSums(config) > 0, ]
  colnames(config) = as.character(aff_Ids)

  keep_rows = apply(config, 1, function(x){
    if(length(aff_Ids[as.vector(x)]) != length(pedcarriers)){FALSE}
    else {
      all(aff_Ids[as.vector(x)] == pedcarriers)
    }
  })
  keep_rows[length(keep_rows)] = TRUE #keep the zero configuration
  config = config[keep_rows, ]

  if (!is.null(nrow(config))) {
    # compute and store the likelihood results in a matrix
    # ROW: realizations of tau
    # COLUMNS: different sharing configurations
    like_mat <- lapply(1:nrow(config), function(i){
      sapply(1:nrow(tau_grid), function(x){
        compute_sharingProb(ped, subtypes, tau = as.numeric(tau_grid[x, ]),
                            carriers = aff_Ids[as.vector(config[i, ])],
                            subtype_weights)})
    })
  } else {
    like_mat <- list(sapply(1:nrow(tau_grid), function(x){
      compute_sharingProb(ped, subtypes, tau = as.numeric(tau_grid[x, ]),
                          carriers = aff_Ids[as.vector(config)],
                          subtype_weights)}))
  }


  num_found = length(unique(ped$ID[which(is.na(ped$dadID) & is.na(ped$momID))]))
  fint_prob = num_found*carrier_prob

  #calculate probability of zero config
  #NOTE: this is P(0|introduced)*carrier_prob +  P(0|not introduced)(1 - carrier_prob),
  # where P(0|not introduced) = 1
  for(i in 1:length(like_mat)){
    if(i < length(like_mat)){
      #incorporate carrier_prob appropriately for non-zero configurations
      like_mat[[i]] <- like_mat[[i]]*fint_prob
    } else {
      #calculate probability of zero config
      #NOTE: this is P(0|introduced)*carrier_prob +  P(0|not introduced)(1 - carrier_prob),
      # P(0|not introduced) = 1
      like_mat[[i]] <- like_mat[[i]]*fint_prob + (1 - fint_prob)
    }
  }

  like_mat <- do.call("cbind", like_mat)



  #compute the binID, this is a numeric identifier for the configuration.
  #The configuration is interpreted as a base 2 number and then converted to
  #base 10.  I.e. the configuration (0, 1, 0) in base 2 --> 010
  #the binID is the base 2 converted to base 10 --> 2.
  #
  #For configuration (0, 0, 1) ---> binID = 1
  #For configuration (0, 1, 0) ---> binID = 2
  #For configuration (0, 1, 1) ---> binID = 3, etc.
  if (!is.null(nrow(config))) {
    # compute and store the likelihood results in a matrix
    # ROW: realizations of tau
    # COLUMNS: different sharing configurations
    binIDs <- as.numeric(apply(config*1, 1, function(x){base::strtoi(paste0(x, collapse = ""), base = 2)}))
  } else {
    binIDs <- as.numeric(base::strtoi(paste0(config*1, collapse = ""), base = 2))
  }

  #The binID will come in handy to quickly refer to and find
  #configurations, especially when dealing with
  #global configurations.

  return(list(like_mat = like_mat,
              configs = config, binIDs = binIDs))
}

#' Compute all distributions - used for power simulation
#'
#' used for power simulation - avoids having to re-run this function multiple times
#'
#' @inheritParams compute_familyWeights
#' @param peds A ped object containing two or more pedigrees.
#' @param carrier_prob The cumulative carrier probabilty of all crvs as a group.
#' @param assumption The assumption to use for computing the LRs (can be "global" or "semi-global")
#'
#' @return A data frame of global configurations anf weights for the set of provided families.
#' @keywords internal
#' @examples
#' library(RVMethods)
#' data(study_pedigrees)
#'
#' library(SimRVPedigree)
#' plot(study_pedigrees[study_pedigrees$FamID == 304, ])
#' plot(study_pedigrees[study_pedigrees$FamID == 58, ])
#'
#' peds = study_pedigrees[study_pedigrees$FamID %in%c(58, 304), ]
#' subtypes = c("HL", "NHL")
#'
#'
#' myCombos = test_allcombos_includeZero(ped = study_pedigrees[study_pedigrees$FamID == 58, ],
#'                                       subtypes = c("HL", "NHL"),
#'                                       tau_grid = make_tauGrid(),
#'                                       carrier_prob = 0.002)
#'
#' my_GDist = compute_globalDist(peds = study_pedigrees[study_pedigrees$FamID %in%c(58, 304), ],
#'                                 subtypes = c("HL", "NHL"))
#'
#' order(my_GDist$LR) == order(my_GDist$w_LR)
#'
compute_configWeigths <- function(peds, subtypes, carriers,
                                  carrier_prob = 0.002,
                                  tau_increment = 0.05,
                                  subtype_weights = NULL,
                                  useK = FALSE){


  #determine probability cRV is introduced and probability not introduced
  #create tau grid
  #NOTE: make_tauGrid now orders taus to save time later.
  tau_grid <- make_tauGrid(increment_width = tau_increment, constrained = TRUE)
  #determine the FamIDs for all families
  study_FamIDs <- unique(peds$FamID)
  famIndex <- c(1:length(study_FamIDs))

  if(length(famIndex) < 2){
    stop("peds only contians one pedigree. Please use compute_familyWeights instead.")
  }

  #carrier_probs = compute_carrierProbs(carrier_probs),
  print(Sys.time())
  #compute all pedigree specific statistics
  message("Compute pedigree-specific statistics.")
  pb <- txtProgressBar(min = 0, max = length(study_FamIDs), style = 3)
  fam_likeGrids <- list()
  for (i in 1:length(study_FamIDs)){
    fam_likeGrids[[i]] <- test_onecombo_withZero(ped = peds[peds$FamID == study_FamIDs[[i]], ],
                                                  subtypes, pedcarriers = carriers[[i]], tau_grid, carrier_prob,
                                                  subtype_weights)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # This section computes the estimates tauA and tauB, and computes the LR and probablity of
  # the configuration under the null.  NOTE: that when the "conditioned semi-global" option
  # is selected the probability of the configuration unser the null is the RVS probability.
  #
  # Additionally, when the "conditioned semi-global" option is selected the returned
  # dataframe will contain a dist ID variable, which is used to identify each marginal
  # distribution.  This will be a critical variable for computing p-values under this
  # assumption.

  #compute probabilities and statistics for the global distribution
  D1_global_sharing_byBinID <- condition_globalDist_zeroConfig(likeGrids_byFam = fam_likeGrids,
                                                              famID_index = famIndex,
                                                              tau_grid)
  GlobalLikelihoodRatio = D1_global_sharing_byBinID[1, which(colnames(D1_global_sharing_byBinID) %in% c("LR", "tau_A",  "tau_B"))]

  #semi-global
  print(paste0("Start Semi-Global Distribution: ", Sys.time()))
  D2_global_sharing_byBinID <- approx_globalDist(likeGrids_byFam = fam_likeGrids,
                                                 famID_index = famIndex,
                                                 tau_grid)
  GlobalTransmission = D2_global_sharing_byBinID[1, which(colnames(D2_global_sharing_byBinID) %in% c("LR", "tau_A",  "tau_B"))]

  #conditioned semi-global
  D3_global_sharing_byBinID <- conditioned_semiGlobalDist(likeGrids_byFam = fam_likeGrids,
                                                          famID_index = famIndex,
                                                          tau_grid)
  LocalLikelihoodRatio = D1_global_sharing_byBinID[1, which(colnames(D3_global_sharing_byBinID) %in% c("LR", "tau_A",  "tau_B"))]


  return(list(GlobalLikelihoodRatio = GlobalLikelihoodRatio,
              GlobalTransmission = GlobalTransmission,
              LocalLikelihoodRatio = LocalLikelihoodRatio))

}
