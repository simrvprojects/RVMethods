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
compute_global_weight <- function(likeGrids_byFam, weights_byFam, famID_index, tau_grid) {
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

  # Find the MLEs.  NOTE: they may not be unique???
  #globalWeights$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  #globalWeights$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  globalWeights <- global_configs
  colnames(globalWeights) <- famID_index

  globalWeights$LR <- L_max/L_null
  globalWeights$RVsharing <- Reduce("*", expanded_RVS_by_fam)
  globalWeights$K <- Reduce("+", expanded_K_by_fam)
  globalWeights$K_sub <- Reduce("+", expanded_K_sub_by_fam)


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


