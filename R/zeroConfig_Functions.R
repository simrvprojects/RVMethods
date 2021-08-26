#' Find the values of tau that maximize the likelihood function
#'
#' Find the values of tau that maximize the likelihood function. NOTE: if there is not a unique MCLE, this funciton orders the results conservatively.  That is, the smallest values of tau_A and tau_B are the first values returned.
#'
#' @inheritParams test_allcombos
#' @param like_mat The matrix of likelihood values given a sharing configuation and values of tau, i.e. the first item returned by \code{\link{test_allcombos}}.
#'
#' @return  A list of taus that maximize the likelihood for each configuration
#' @keywords internal
get_MLE_zeroConfig <- function(like_mat, tau_grid){

  #condition on cRV being sen in at least one study member
  like_mat_zero = like_mat[, ncol(like_mat)]
  like_mat = like_mat[, -ncol(like_mat)]

  like_mat_noZero = do.call(rbind, lapply(1:nrow(like_mat), function(x){
    like_mat[x, ]/like_mat_zero[[x]]
  }))

  # max_like <- apply(like_mat, 2, function(x){
  #   x[which.max(x)]
  # })
  #
  # max_taus <- lapply(1:ncol(like_mat), function(x){
  #   tau_grid[like_mat[, x] == like_mat[max_like[[x]], x], ]
  # })

  max_taus <- apply(like_mat_noZero, 2, function(x){
    tau_grid[x == max(x), ]
  })

  return(max_taus)
}


#' Compute sharing probabilities for every configuration of disease-affected relatives and initialization of tau
#'
#' Compute sharing probabilities for every configuration of disease-affected relatives and initialization of tau
#'
#' @param tau_grid A grid of tau values produced by \code{expand.grid}.  \code{test_allcombos} will compute each sharing probability at every specified value of tau in \code{tau_grid}.
#' @param carrier_prob The cumulative carrier probabilty of all crvs as a group.
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
#' myCombos = test_allcombos_includeZero(ped = study_pedigrees[study_pedigrees$FamID == 58, ],
#'                                       subtypes = c("HL", "NHL"),
#'                                       tau_grid = make_tauGrid(),
#'                                       carrier_prob = 0.002)
#'
#' myCombos
#'
#' rowSums(myCombos$like_mat)
#'
test_allcombos_includeZero <- function(ped, subtypes, tau_grid, carrier_prob,
                           subtype_weights = NULL){

  #determine the total number of affected in this pedigree
  num_affected <- sum(ped$affected, na.rm = TRUE)

  #store the IDs of the affected relatives
  aff_Ids <- ped$ID[ped$affected & ped$available]

  #determine all possible sharing configurations
  config = as.matrix(expand.grid(rep(list(c(T, F)), num_affected)))
  #config = config[rowSums(config) > 0, ]
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

  #calculate probability of zero config
  #NOTE: this is P(0|introduced)*carrier_prob +  P(0|not introduced)(1 - carrier_prob),
  # P(0|not introduced) = 1
  for(i in 1:length(like_mat)){
    if(i < length(like_mat)){
      #incorporate carrier_prob appropriately for non-zero configurations
      like_mat[[i]] <- like_mat[[i]]*carrier_prob
      #like_mat[[i]] <- like_mat[[i]]*carrier_probs[[1]]
    } else {
      #calculate probability of zero config
      #NOTE: this is P(0|introduced)*carrier_prob +  P(0|not introduced)(1 - carrier_prob),
      # P(0|not introduced) = 1
      #like_mat_zero_crvProb  = like_mat_zero*carrier_prob + (1 - carrier_prob)
      like_mat[[i]] <- like_mat[[i]]*carrier_prob + (1 - carrier_prob)
      #like_mat[[i]] <- like_mat[[i]]*carrier_probs[[1]] + (1 - carrier_probs[[2]])
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
  binIDs <- as.numeric(apply(config*1, 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  }))

  #The binID will come in handy to quickly refer to and find
  #configurations, especially when dealing with
  #global configurations.




  return(list(like_mat = like_mat,
              configs = config, binIDs = binIDs))
}



#' Compute the sharing probablities under the null, i.e. when all taus are 1/2
#'
#' @param famLike_results The output returned by the \code{\link{test_allcombos}} function.
#' @inheritParams compute_familyWeights
#'
#' @return A data frame containing the weights for each configuration in the pedigree. See Details.
#' @keywords internal
#'
compute_famLR_zeroConfig <- function(famLike_results, tau_grid,
                                     ped, subtypes){

  # find the values of tau that maximize the likelihood for each configuration
  MLEs <- get_MLE(famLike_results$like_mat, tau_grid)

  #create a data.frame that combines all of the information
  famWeights <- as.data.frame(famLike_results$configs*1)
  famWeights$binID <- famLike_results$binIDs


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
  #famWeights$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  #famWeights$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  # For certain sharing configurations and subtype assinments, the likelihood
  # may have multiple maxima, resulting in non-unique MLEs
  # We identify these in the output.
  famWeights$unique_tau <- sapply(MLEs, function(x){nrow(x) == 1})


  #compute the numerator of the RVsharing probablity
  #NOTE this is the first row of like_mat in the argument famLike_results
  #i.e. the likelihood when tau_A = tau_B = 0.5
  famWeights$null_configProb <- famLike_results$like_mat[1, ]

  return(famWeights)
}


#' Compute Global weights for a set of families
#'
#' @param likeGrids_byFam  A list of results from test_allcombos, one for each family.
#' @param weights_byFam A list of results from compute_famWeights, one for each family.
#' @param tau_grid The grid of tau values
#' @param famID_index The family IDs, indexed in the same order as likeGrids_byFam and weights_byFam
#'
#' @return A data frame of global configurations anf weights for the set of provided families.
#' @keywords internal
#'
compute_globalConfigDist <- function(likeGrids_byFam, weights_byFam, famID_index, tau_grid) {
  #determine the family configuration binIDs
  fam_bin_configs = lapply(likeGrids_byFam, `[[`, 3)

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their binary representations
  global_configs = expand.grid(fam_bin_configs)

  #this list holds the family likelihoods for each configuration
  #evaluated at each value of tau
  fam_likes = lapply(likeGrids_byFam, `[[`, 1)

  #get the family config index, i.e. the column of the config
  #in fam_likes.  NOTE: these are already ordered
  fam_bin_configs_index <- lapply(fam_bin_configs, function(x){
    seq(1:length(x))
  })

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their binary representations
  global_configs_index = expand.grid(fam_bin_configs_index)

  #expand each familial likelihhod matrix accordingly
  expanded_like_by_fam = list()
  expanded_prob_by_fam = list()
  for(i in 1:length(likeGrids_byFam)){
    expanded_like_by_fam[[i]] <- do.call(cbind, lapply(global_configs_index[, i], function(x){
      fam_likes[[i]][, x]
    }))
    expanded_prob_by_fam[[i]] <- weights_byFam[[i]]$null_configProb[global_configs_index[, i]]
  }

  #EACH global configuration will become a column in the likelihood matrix grio
  #EACH value of tau_A and tau_B have already been computed for
  #each family.  When families are unrelated we can multiply the likelihoods
  #to get the likelihood of the global configuration for the specified values of
  #tau_A and tau_B
  global_likeMat <- Reduce("*", expanded_like_by_fam)


  # find the values of tau that maximize the likelihood for each configuration
  #MLEs <- get_MLE_zeroConfig(like_mat = global_likeMat, tau_grid)
  MLEs <- get_MLE(like_mat = global_likeMat, tau_grid)

  # find the maximum value of the likelihood
  L_max <- apply(global_likeMat, 2, function(x){
    x[which.max(x)]
  })

  # value of the likelihood under H_0
  L_null <- global_likeMat[1, ]

  #compute the composite LR statistic
  LR <- L_max/L_null
  configProb <- Reduce("*", expanded_prob_by_fam)


  globalWeights <- global_configs
  colnames(globalWeights) <- famID_index

  globalWeights$LR <- L_max/L_null
  globalWeights$configProb <- Reduce("*", expanded_prob_by_fam)
  # Find the MLEs.  NOTE: they may not be unique???
  globalWeights$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  globalWeights$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  return(globalWeights)
}


#' Compute Global weights for a set of families
#'
#' @inheritParams compute_familyWeights
#' @param peds A ped object containing two or more pedigrees.
#' @param carrier_prob The cumulative carrier probabilty of all crvs as a group.
#'
#' @return A data frame of global configurations anf weights for the set of provided families.
#' @keywords internal
#' @examples
#' library(RVMethods)
#' data(study_pedigrees)
#'
#' library(SimRVPedigree)
#' plot(study_pedigrees[study_pedigrees$FamID == 304, ])
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
#' my_GDist = compute_globalLRDist(peds = study_pedigrees[study_pedigrees$FamID %in%c(58, 304), ],
#'                                 subtypes = c("HL", "NHL"))
#'
#' order(my_GDist$LR) == order(my_GDist$w_LR)
#'
compute_globalLRDist <- function(peds, subtypes, carrier_prob = 0.002,
                                 tau_increment = 0.05,
                                 subtype_weights = NULL){


  #create tau grid
  #NOTE: values are ordered in make_tauGrid return.
  tau_grid <- make_tauGrid(increment_width = tau_increment, constrained = TRUE)
  #determine the FamIDs for all families
  study_FamIDs <- unique(peds$FamID)
  famIndex <- c(1:length(study_FamIDs))

  if(length(famIndex) < 2){
    stop("peds only contians one pedigree. Please use compute_familyWeights instead.")
  }

  #compute all pedigree specific statistics
  message("Compute pedigree-specific statistics... this may take a few minutes.")
  pb <- txtProgressBar(min = 0, max = length(study_FamIDs), style = 3)
  fam_likeGrids <- list()
  for (i in 1:length(study_FamIDs)){
    fam_likeGrids[[i]] <- test_allcombos_includeZero(ped = peds[peds$FamID == study_FamIDs[[i]], ],
                                                     subtypes, tau_grid, carrier_prob)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  #carrier_probs = compute_carrierProbs(carrier_probs),

  # Remove any invalid configurations from the all combos result.
  # Invalid configurations are configurations that are not possible
  # under the assumption of a single founder introducing one copy of the
  # RV.
  # NOTE: this is to save space when we conpute gloabl sharing configurations
  #
  # An invalid configuration occurs when mutliple founders must
  # introduce an SNV in order to observe a sharing configuration
  fam_likeGrids <- lapply(fam_likeGrids, function(x){
    remove_invalidConfigs(x)
  })


  #For each family compute sharing weights for SNVs observed
  #in only one family
  fam_LRdist <- lapply(1:length(study_FamIDs), function(x){
    compute_famLR_zeroConfig(famLike_results = fam_likeGrids[[x]], tau_grid,
                             ped = peds[peds$FamID == study_FamIDs[x], ],
                             subtypes)
  })


  global_sharing_byBinID <- compute_globalConfigDist(likeGrids_byFam = fam_likeGrids,
                                                     weights_byFam = fam_LRdist,
                                                     famID_index = famIndex,
                                                     tau_grid)




  #expand familial configuratiosn and re-lable so that individuals
  #are identified by family ID and individual ID
  fam_configs <- list()
  for(i in famIndex){
    fam_configs[[i]] <- 1*fam_likeGrids[[i]]$configs[match(global_sharing_byBinID[, i],
                                                           fam_likeGrids[[i]]$binIDs), ]
    colnames(fam_configs[[i]]) <- paste0(study_FamIDs[[i]], ":", colnames(fam_configs[[i]]))
  }

  fam_configs <- do.call(cbind, fam_configs)

  #combine the configurations with their statistics
  global_dist <- cbind(fam_configs,
                          global_sharing_byBinID[, -c(1:length(study_FamIDs))])

  #remove the all-zero configuration and re-calculate the config probability
  #given that it was observed in at least one study member.
  #NOTE: the all zero config is always in the very last row but adding
  #next line in case I'm wrong
  zero_config_row = which(apply(global_sharing_byBinID[ ,c(1:length(study_FamIDs))],
                                1, function(x){(all(x == 0))}))
  zero_config_prob = global_dist$configProb[zero_config_row]

  global_dist = global_dist[-zero_config_row, ]
  global_dist$configProb = global_dist$configProb/(sum(global_dist$configProb))


  #calculate the modified RVS weight
  global_dist$w_LR <- sapply(1:nrow(global_dist), function(x){
    1/sum(global_dist$configProb[global_dist$LR >= global_dist$LR[x]],
          na.rm = TRUE)
  })


  global_dist$binID <- apply(global_dist[, 1:ncol(fam_configs)], 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })

  return(global_dist)

}
