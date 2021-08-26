#' Compute Global weights for a set of families
#'
#' @inheritParams compute_familyWeights
#' @param peds A ped object containing two or more pedigrees.
#'
#' @return A data frame of global configurations anf weights for the set of provided families.
#' @keywords internal
#'
compute_globalWeights <- function(peds, subtypes,
                                  tau_increment = 0.05,
                                  subtype_weights = NULL,
                                  identify_trueConfig = FALSE){
  #likeGrids_byFam, weights_byFam, famID_index, tau_grid) {

  #create tau grid
  #NOTE: make_tauGrid now orders taus to save time later.
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
    fam_likeGrids[[i]] <- test_allcombos(ped = peds[peds$FamID == study_FamIDs[[i]], ],
                                         subtypes, tau_grid, subtype_weights)
    setTxtProgressBar(pb, i)
  }
  close(pb)


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
  fam_sharingWeights <- lapply(1:length(study_FamIDs), function(x){
    compute_famWeights_internal(famLike_results = fam_likeGrids[[x]], tau_grid,
                                ped = peds[peds$FamID == study_FamIDs[x], ],
                                subtypes, subtype_weights)
  })


  global_sharing_byBinID <- compute_global_weight_internal(likeGrids_byFam = fam_likeGrids,
                                                           weights_byFam = fam_sharingWeights,
                                                           famID_index = c(1:length(study_FamIDs)),
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
  global_weights <- cbind(fam_configs,
                          global_sharing_byBinID[, -c(1:length(study_FamIDs))])


  # #compute the likelihood-based weight
  # global_weights$w_LR <- sapply(1:nrow(global_weights), function(x){
  #   compute_transmission_weight(all_stats = global_weights, config_index = x)
  # })
  #
  # #calculate the modified RVS weight
  # global_weights$w_LRbasic <- sapply(1:nrow(global_weights), function(x){
  #   compute_LRB_weight(all_stats = global_weights, config_index = x)
  # })

  #compute the RVS weight
  global_weights$w_RVS <- sapply(1:nrow(global_weights), function(x){
    compute_RVS_weight(all_stats = global_weights, config_index = x)
  })


  #calculate the modified RVS weight
  global_weights$w_RVS2 <- sapply(1:nrow(global_weights), function(x){
    compute_modified_RVS_weight(all_stats = global_weights, config_index = x)
  })

  global_weights$binID <- apply(global_weights[, 1:ncol(fam_configs)], 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })


return(global_weights)

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
compute_global_weight_internal <- function(likeGrids_byFam, weights_byFam, famID_index, tau_grid) {
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
  expanded_RVS_by_fam = list()
  expanded_K_by_fam = list()
  expanded_K_sub_by_fam = list()
  for(i in 1:length(likeGrids_byFam)){
    expanded_like_by_fam[[i]] <- do.call(cbind, lapply(global_configs_index[, i], function(x){
      fam_likes[[i]][, x]
    }))
    expanded_RVS_by_fam[[i]] <- weights_byFam[[i]]$RVsharing[global_configs_index[, i]]
    expanded_K_by_fam[[i]] <- weights_byFam[[i]]$K[global_configs_index[, i]]
    expanded_K_sub_by_fam[[i]] <- weights_byFam[[i]]$K_sub[global_configs_index[, i]]
  }

  #EACH global configuration will become a column in the likelihood matrix grio
  #EACH value of tau_A and tau_B have already been computed for
  #each family.  When families are unrelated we can multiply the likelihoods
  #to get the likelihood of the global configuration for the specified values of
  #tau_A and tau_B
  global_likeMat <- Reduce("*", expanded_like_by_fam)


  # find the values of tau that maximize the likelihood for each configuration
  MLEs <- get_MLE(global_likeMat, tau_grid)

  # find the maximum value of the likelihood
  L_max <- apply(global_likeMat, 2, function(x){
    x[which.max(x)]
  })

  # value of the likelihood under H_0
  L_null <- global_likeMat[1, ]

  #compute the composite LR statistic
  LR <- L_max/L_null
  RVsharing <- Reduce("*", expanded_RVS_by_fam)
  K <- Reduce("+", expanded_K_by_fam)
  K_sub <- Reduce("+", expanded_K_sub_by_fam)


  globalWeights <- global_configs
  colnames(globalWeights) <- famID_index

  globalWeights$LR <- L_max/L_null
  globalWeights$RVsharing <- Reduce("*", expanded_RVS_by_fam)
  globalWeights$K <- Reduce("+", expanded_K_by_fam)
  globalWeights$K_sub <- Reduce("+", expanded_K_sub_by_fam)

  # Find the MLEs.  NOTE: they may not be unique???
  globalWeights$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  globalWeights$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  # #compute the likelihood-based weight
  # globalWeights$w_LR <- sapply(1:nrow(globalWeights), function(x){
  #   sum(RVsharing[LR >= LR[x] & K >= K[x]],
  #       na.rm = TRUE)
  # })
  # globalWeights$w_LR <- 1/globalWeights$w_LR
  #
  # #compute the RVS weight
  # globalWeights$w_RVS <- sapply(1:nrow(globalWeights), function(x){
  #   sum(RVsharing[RVsharing <= RVsharing[x] & K >= K[x]],
  #       na.rm = TRUE)
  # })
  # globalWeights$w_RVS <- 1/globalWeights$w_RVS
  #
  # #calculate the modified RVS weight
  # globalWeights$w_RVS2 <- sapply(1:nrow(globalWeights), function(x){
  #   sum(RVsharing[K >= K[x]
  #                 & K_sub >= K_sub[x]
  #                 & RVsharing <= RVsharing[x]],
  #       na.rm = TRUE)
  # })
  # globalWeights$w_RVS2 <- 1/globalWeights$w_RVS2


  return(globalWeights)
}


