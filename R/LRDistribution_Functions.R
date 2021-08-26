#' (Option 1) Compute sharing probabilities for every global configuration (includes the zero configuration) and initialization of tau
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
test_allcombos_withZero <- function(ped, subtypes, tau_grid, carrier_prob,
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
  # for(i in 1:length(like_mat)){
  #   if(i < length(like_mat)){
  #     #incorporate carrier_prob appropriately for non-zero configurations
  #     like_mat[[i]] <- like_mat[[i]]*carrier_prob
  #   } else {
  #     #calculate probability of zero config
  #     #NOTE: this is P(0|introduced)*carrier_prob +  P(0|not introduced)(1 - carrier_prob),
  #     # P(0|not introduced) = 1
  #     like_mat[[i]] <- like_mat[[i]]*carrier_prob + (1 - carrier_prob)
  #   }
  # }


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


#' (Option 1) Conditon the probabilty of each global configuration on the event that the variant was observed in at least one study member.
#'
#' Compute MLE's for tau_a and tau_b under global configuration framework: here configurations that are identically zero contribute to the estimation of the MLEs; also, probabilities are conditioned on the variant being observed in at least one disease-affected study member.
#'
#' This will be an internal function
#'
#' @param likeGrids_byFam  A list of results from test_allcombos, one for each family.
#' @param weights_byFam A list of results from compute_famWeights, one for each family.
#' @param tau_grid The grid of tau values
#' @param famID_index The family IDs, indexed in the same order as likeGrids_byFam and weights_byFam
#'
#' @return A data frame of global configurations anf weights for the set of provided families.
#' @keywords internal
#'
condition_globalDist_zeroConfig <- function(likeGrids_byFam, famID_index, tau_grid) {

  #Pull the familial binID from likeGrids_byFam and store in one list
  fam_bin_configs = lapply(likeGrids_byFam, `[[`, 3)

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their binary representations
  global_configs = expand.grid(fam_bin_configs)

  #Pull the probability matricies from likeGrids_byFam and store in one list
  fam_likes = lapply(likeGrids_byFam, `[[`, 1)

  #get the family config index, i.e. the column of the config
  #in fam_likes.  NOTE: these are already ordered
  fam_bin_configs_index <- lapply(fam_bin_configs, function(x){
    seq(1:length(x))
  })

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their index number in
  # fam_bin_configs
  global_configs_index = expand.grid(fam_bin_configs_index)

  #expand each familial likelihhod matrix accordingly
  expanded_like_by_fam = list()

  # The probability of global configuration g is defined to the product of the
  # individual configurations (since families are assumend to be unrelated,
  # and thus indpendent). For example, with two configurations:
  # P(g|t) = P(c1, c2 |t) = P(c1|t)*P(c2|t)
  #
  # Based on the global configurations defined in global config index we pull
  # the approporate columns (multiple times) from the familial probability grid.
  #
  # Then, after each families probability grid has been expanded approproately
  # we perform element wise multiplication to get P(c1|t)*P(c2|t)
  #

  for(i in 1:length(likeGrids_byFam)){
    expanded_like_by_fam[[i]] <- do.call(cbind, lapply(global_configs_index[, i], function(x){
      fam_likes[[i]][, x]
    }))
  }

  # Perform element wise multiplication to get P(c1|t)*P(c2|t)
  # NOTE: columns are global configrations
  # rows are probabilities computed under different values of tau
  global_probMat <- Reduce("*", expanded_like_by_fam)


  #find the row of tau_A and tau_B, under the null.
  #NOTE: by design, this is always teh first row of tau grid and in
  #fam_likeResults, but determing row explicitly just in case I've screwed up.
  #TODO: remove this after you implement a test to check that
  #the null taus are always in the first row.
  null_taus <- which(tau_grid[, 1] == 0.5 & tau_grid[, 2] == 0.5)

  #determine the zero configuration.
  #NOTE: by design this is the row in global_configs and last column in global_probMat,
  #but implementing catch incase I'm mistaken
  #TODO: remove after implementing test to ensure this asumption is correct.
  zero_config <- which(rowSums(global_configs) == 0)

  #The likelihood matrix has the following form:
  # NOTE: to keep this simple only showing tau value, but techinically there are
  # two taus. The tau values in rows of tau_grid describe the assumed taus in
  # row of the likelihood matrix.
  #
  # Column configurations are descibed by the rows in famLike_results$configs
  #
  # | p(g1|t = 1/2)  p(g2|t = 1/2) ...  p(g = 0|t = 1/2) |
  # | p(g1|t = 3/4)  p(g2|t = 3/4) ...  p(g = 0|t = 3/4) |
  # |    ...             ...                  ...        |
  # | p(g1|t = 1)     p(g2|t = 1) ...    p(g = 0|t = 1)  |
  #
  #
  #

  # Determine the probablity that configuration is non-zero for each tau
  # i.e. 1 - P(g = 0|t)
  prob_nonzero_given_tau <- 1 - global_probMat[, zero_config]


  #remove the zero configuration from the probability grid
  nonZero_probGrid <- global_probMat[, -zero_config]

  #this is needed to re-configure when we are computing the probability
  #for unequally weighted founders under the alternative
  if(!is.matrix(nonZero_probGrid)){
    nonZero_probGrid = matrix(nonZero_probGrid, nrow = 1)
  }

  #condition on the variant being seen in at least one familiy member
  #first create a matrix with the following form
  # | (1- p(g = 0|t = 1/2))  (1 - p(g = 0|t = 1/2)) ...  (1 - p(g = 0|t = 1/2)) |
  # | (1- p(g = 0|t = 3/4))  (1 - p(g = 0|t = 3/4)) ...  (1 - p(g = 0|t = 3/4)) |
  # |     ...                   ...                          ...                |
  # | (1- p(g = 0|t = 1))    (1 - p(g = 0|t = 1))   ...  (1 - p(g = 0|t = 1))   |
  #
  # Then multiply nonZero_probGrid by reciprocal so that each entry is of the form
  # p(g1|t)*(1/p(g != 0|t))
  #
  denomProb <- matrix(rep(prob_nonzero_given_tau, ncol(nonZero_probGrid)),
                     ncol = ncol(nonZero_probGrid))

  #compute conditioned probs, note does NOT include zero configuration
  conditoned_Probs <- nonZero_probGrid*(1/denomProb)


  # find the values of tau that maximize the likelihood for each configuration
  #MLEs <- get_MLE_zeroConfig(like_mat = global_likeMat, tau_grid)
  MLEs <- get_MLE(like_mat = conditoned_Probs, tau_grid)



  #Store the results in a single matrix
  globalDist <- global_configs[-zero_config, ] #NOTE these are the true binIDs, not their index
  colnames(globalDist) <- famID_index

  # find the maximum value of the likelihood for each configuration
  globalDist$L_max <- apply(conditoned_Probs, 2, function(x){
    x[which.max(x)]
  })

  # value of the likelihood under H_0, i.e. when both taus are 1/2
  globalDist$L_null <- conditoned_Probs[null_taus, ]

  globalDist$LR <- globalDist$L_max/globalDist$L_null
  globalDist$null_configProb <- globalDist$L_null
  # Find the MLEs.  NOTE: they may not be unique???
  globalDist$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  globalDist$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  return(globalDist)
}




#' (Option 2) Determine the MLEs and LR in the semi global framework.
#'
#' This will be an internal function
#'
#' @param likeGrids_byFam  A list of results from test_allcombos, one for each family.
#' @param weights_byFam A list of results from compute_famWeights, one for each family.
#' @param tau_grid The grid of tau values
#' @param famID_index The family IDs, indexed in the same order as likeGrids_byFam and weights_byFam
#'
#' @return A data frame of global configurations anf weights for the set of provided families.
#' @keywords internal
#'
approx_globalDist <- function(likeGrids_byFam, famID_index, tau_grid) {

  #find the row of tau_A and tau_B, under the null.
  #NOTE: by design, this is always the first row of tau grid and in
  #fam_likeResults, but determing row explicitly just in case I've screwed up.
  #TODO: remove this after you implement a test to check that
  #the null taus are always in the first row.
  null_taus <- which(tau_grid[, 1] == 0.5 & tau_grid[, 2] == 0.5)

  #Pull the familial binID from likeGrids_byFam and store in one list
  fam_bin_configs = lapply(likeGrids_byFam, `[[`, 3)

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their binary representations
  global_configs = expand.grid(fam_bin_configs)

  #Pull the probability matricies from likeGrids_byFam and store in one list
  fam_likes = lapply(likeGrids_byFam, `[[`, 1)

  #this will be used to compute the null probability of the configuration
  #i.e. approximated by the global distribution null probabilites
  #
  # Store null probabilities under tau_a = 0.5 and tau_b = 0.5
  fam_likes_null = lapply(fam_likes, function(x){
    x[null_taus, ]
  })

  #If the familial configuraiton is identically zero set probability to 1
  # for every realization of tau_a and tau_b.
  #
  # i.e. if c1 = 0 then P(c1|t) = 1
  for(i in 1:length(fam_likes)){
    fam_likes[[i]][, which(fam_bin_configs[[i]] == 0)] = 1
  }

  #get the family config index, i.e. the column of the config
  #in fam_likes.  NOTE: these are already ordered
  fam_bin_configs_index <- lapply(fam_bin_configs, function(x){
    seq(1:length(x))
  })

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their index number in
  # fam_bin_configs
  global_configs_index = expand.grid(fam_bin_configs_index)

  #expand each familial likelihhod matrix accordingly
  expanded_like_by_fam = list()

  # The probability of global configuration g is defined to the product of the
  # individual configurations (since families are assumend to be unrelated,
  # and thus indpendent). For example, with two configurations:
  # P(g|t) = P(c1, c2 |t) = P(c1|t)*P(c2|t)
  #
  # Based on the global configurations defined in global config index we pull
  # the approporate columns (multiple times) from the familial probability grid.
  #
  # Then, after each families probability grid has been expanded approproately
  # we perform element wise multiplication to get P(c1|t)*P(c2|t)
  #
  # NOTE: if c1 = 0 then P(c1|t) = 1, then P(c1|t)*P(c2|t) = P(c2|t)
  # That is, configurations that are identically zero DO NOT contribute to the likelihood
  #
  for(i in 1:length(likeGrids_byFam)){
    expanded_like_by_fam[[i]] <- do.call(cbind, lapply(global_configs_index[, i], function(x){
      fam_likes[[i]][, x]
    }))
  }


  # Perform element wise multiplication to get P(c1|t)*P(c2|t)
  # NOTE: columns are global configrations
  # rows are probabilities computed under different values of tau
  global_probMat <- Reduce("*", expanded_like_by_fam)


  #determine the zero configuration.
  #NOTE: by design this is the last row in global_configs and last column in global_probMat,
  #but implementing catch incase I'm mistaken
  #TODO: remove after implementing test to ensure this asumption is correct.
  zero_config <- which(rowSums(global_configs) == 0)

  #The likelihood matrix has the following form:
  # NOTE: to keep this simple only showing tau value, but techinically there are
  # two taus. The tau values in rows of tau_grid describe the assumed taus in
  # row of the likelihood matrix.
  #
  # Column configurations are descibed by the rows in famLike_results$configs
  #
  # | p(g1|t = 1/2)  p(g2|t = 1/2) ...  p(g = 0|t = 1/2) |
  # | p(g1|t = 3/4)  p(g2|t = 3/4) ...  p(g = 0|t = 3/4) |
  # |    ...             ...                  ...        |
  # | p(g1|t = 1)     p(g2|t = 1) ...    p(g = 0|t = 1)  |
  #
  #
  #

  #remove the zero configuration from the probability grid
  nonZero_probGrid <- global_probMat[, -zero_config]

  #this is needed to re-configure when we are computing the probability
  #for unequally weighted founders under the alternative
  if(!is.matrix(nonZero_probGrid)){
    nonZero_probGrid = matrix(nonZero_probGrid, nrow = 1)
  }

  # find the values of tau that maximize the likelihood for each configuration
  #MLEs <- get_MLE_zeroConfig(like_mat = global_likeMat, tau_grid)
  MLEs <- get_MLE(like_mat = nonZero_probGrid, tau_grid)



  #Store the results in a single matrix
  globalDist <- global_configs[-zero_config, ] #NOTE these are the true binIDs, not their index
  colnames(globalDist) <- famID_index

  # find the maximum value of the likelihood for each configuration
  globalDist$L_max <- apply(nonZero_probGrid, 2, function(x){
    x[which.max(x)]
  })

  # value of the likelihood under H_0, i.e. when both taus are 1/2
  globalDist$L_null <- nonZero_probGrid[null_taus, ]

  globalDist$LR <- globalDist$L_max/globalDist$L_null

  #compute the probability under the null (NOTE: this incorporates the carrier probablitity and
  #is thus the probabilty under the null according to the global distribution)
  config_null_probs = Reduce("*", lapply(1:length(fam_likes_null), function(x){
    fam_likes_null[[x]][global_configs_index[, x]]
  }))
  #condition on variant being observed in at least one sequenced disease-affected study member
  config_null_probs = config_null_probs[-zero_config]/(1 - config_null_probs[zero_config])

  globalDist$null_configProb <- config_null_probs

  # Find the MLEs.  NOTE: they may not be unique
  globalDist$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  globalDist$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  return(globalDist)
}

#' (Option 3) Determine the MLEs and LR in the semi global framework, with conditioning imposed.
#'
#' This will be an internal function
#'
#' @param likeGrids_byFam  A list of results from test_allcombos, one for each family.
#' @param weights_byFam A list of results from compute_famWeights, one for each family.
#' @param tau_grid The grid of tau values
#' @param famID_index The family IDs, indexed in the same order as likeGrids_byFam and weights_byFam
#'
#' @return A data frame of global configurations anf weights for the set of provided families.
#' @keywords internal
#'
conditioned_semiGlobalDist <- function(likeGrids_byFam, famID_index, tau_grid) {

  #find the row of tau_A and tau_B, under the null.
  #NOTE: by design, this is always the first row of tau grid and in
  #fam_likeResults, but determing row explicitly just in case I've screwed up.
  #TODO: remove this after you implement a test to check that
  #the null taus are always in the first row.
  null_taus <- which(tau_grid[, 1] == 0.5 & tau_grid[, 2] == 0.5)

  #Pull the familial binID from likeGrids_byFam and store in one list
  fam_bin_configs = lapply(likeGrids_byFam, `[[`, 3)

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their binary representations
  global_configs = expand.grid(fam_bin_configs)

  #Pull the probability matricies from likeGrids_byFam and store in one list
  fam_likes = lapply(likeGrids_byFam, `[[`, 1)

  #store location of each familial zero configuration.
  fam_zero_config <- sapply(fam_bin_configs, function(x){
    which(x == 0)
  })

  #compute the probablility that the variant was observed in each family for each
  #value of tau_a and tau_b
  cond_mats <- lapply(1:length(fam_likes), function(x){
    matrix(1/(1 - fam_likes[[x]][, fam_zero_config[[x]]]), ncol = ncol(fam_likes[[x]]), nrow = nrow(fam_likes[[x]]))
  })

  #condition on the variant being observed in each family
  cond_fam_likes <- lapply(1:length(fam_likes), function(x){
    fam_likes[[x]]*cond_mats[[x]]
  })

  #If the familial configuraiton is identically zero set probability to 1
  # for every realization of tau_a and tau_b.
  #
  # i.e. if c1 = 0 then P(c1|t) = 1
  for(i in 1:length(cond_fam_likes)){
    cond_fam_likes[[i]][, which(fam_bin_configs[[i]] == 0)] = 1
  }

  #get the family config index, i.e. the column of the config
  #in fam_likes.  NOTE: these are already ordered
  fam_bin_configs_index <- lapply(fam_bin_configs, function(x){
    seq(1:length(x))
  })

  #determine all global configs for the supplied families
  #NOTE: the configurations are indicated by their index number in
  # fam_bin_configs
  global_configs_index = expand.grid(fam_bin_configs_index)

  #expand each familial likelihhod matrix accordingly
  expanded_like_by_fam = list()

  # The probability of global configuration g is defined to the product of the
  # individual configurations (since families are assumend to be unrelated,
  # and thus indpendent). For example, with two configurations:
  # P(g|t) = P(c1, c2 |t) = P(c1|t)*P(c2|t)
  #
  # Based on the global configurations defined in global config index we pull
  # the approporate columns (multiple times) from the familial probability grid.
  #
  # Then, after each families probability grid has been expanded approproately
  # we perform element wise multiplication to get P(c1|t)*P(c2|t)
  #
  # NOTE: if c1 = 0 then P(c1|t) = 1, then P(c1|t)*P(c2|t) = P(c2|t)
  # That is, configurations that are identically zero DO NOT contribute to the likelihood
  #
  for(i in 1:length(likeGrids_byFam)){
    expanded_like_by_fam[[i]] <- do.call(cbind, lapply(global_configs_index[, i], function(x){
      cond_fam_likes[[i]][, x]
    }))
  }


  # Perform element wise multiplication to get P(c1|t)*P(c2|t)
  # NOTE: columns are global configrations
  # rows are probabilities computed under different values of tau
  global_probMat <- Reduce("*", expanded_like_by_fam)


  #determine the zero configuration.
  #NOTE: by design this is the last row in global_configs and last column in global_probMat,
  #but implementing catch incase I'm mistaken
  #TODO: remove after implementing test to ensure this asumption is correct.
  zero_config <- which(rowSums(global_configs) == 0)

  #The likelihood matrix has the following form:
  # NOTE: to keep this simple only showing tau value, but techinically there are
  # two taus. The tau values in rows of tau_grid describe the assumed taus in
  # row of the likelihood matrix.
  #
  # Column configurations are descibed by the rows in famLike_results$configs
  #
  # | p(g1|t = 1/2)  p(g2|t = 1/2) ...  p(g = 0|t = 1/2) |
  # | p(g1|t = 3/4)  p(g2|t = 3/4) ...  p(g = 0|t = 3/4) |
  # |    ...             ...                  ...        |
  # | p(g1|t = 1)     p(g2|t = 1) ...    p(g = 0|t = 1)  |
  #
  #
  #

  #remove the zero configuration from the probability grid
  nonZero_probGrid <- global_probMat[, -zero_config]

  #this is needed to re-configure when we are computing the probability
  #for unequally weighted founders under the alternative
  if(!is.matrix(nonZero_probGrid)){
    nonZero_probGrid = matrix(nonZero_probGrid, nrow = 1)
  }

  # find the values of tau that maximize the likelihood for each configuration
  #MLEs <- get_MLE_zeroConfig(like_mat = global_likeMat, tau_grid)
  MLEs <- get_MLE(like_mat = nonZero_probGrid, tau_grid)

  #Store the results in a single matrix
  globalDist <- global_configs[-zero_config, ] #NOTE these are the true binIDs, not their index
  colnames(globalDist) <- famID_index

  # find the maximum value of the likelihood for each configuration
  globalDist$L_max <- apply(nonZero_probGrid, 2, function(x){
    x[which.max(x)]
  })

  # value of the likelihood under H_0, i.e. when both taus are 1/2
  globalDist$L_null <- nonZero_probGrid[null_taus, ]

  globalDist$LR <- globalDist$L_max/globalDist$L_null

  globalDist$null_configProb <- globalDist$L_null

  # Find the MLEs.  NOTE: they may not be unique
  globalDist$tau_A <- sapply(MLEs, function(x){x[1, 1]})
  globalDist$tau_B <- sapply(MLEs, function(x){x[1, 2]})

  #find distribution ID of each configuration - this will be used to compute p-values
  dist_inds <- do.call(cbind, lapply(1:ncol(global_configs_index),
                                     function(x){1 - global_configs_index[, x] %in% fam_zero_config[[x]]}))
  #remove zero configuration
  dist_inds <- dist_inds[-zero_config, ]

  globalDist$distID <- sapply(1:nrow(dist_inds), function(x){
    base::strtoi(paste0(dist_inds[x, ], collapse = ""), base = 2)
  })

  return(globalDist)
}


#' Compute the distribution of the familial configurations for a set of families
#'
#' NOTE: TO compute the distrubtion of the LR in this setting, we have to make assumptions regarding the carrier probability of the RV.
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
compute_globalDist <- function(peds, subtypes, carrier_prob = 0.002,
                               tau_increment = 0.05,
                               assumption = "global",
                               subtype_weights = NULL){


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

  #compute all pedigree specific statistics
  message("Compute pedigree-specific statistics... this may take a few minutes.")
  pb <- txtProgressBar(min = 0, max = length(study_FamIDs), style = 3)
  fam_likeGrids <- list()
  for (i in 1:length(study_FamIDs)){
    fam_likeGrids[[i]] <- test_allcombos_withZero(ped = peds[peds$FamID == study_FamIDs[[i]], ],
                                                  subtypes, tau_grid, carrier_prob,
                                                  subtype_weights)
    setTxtProgressBar(pb, i)
  }
  close(pb)

  #compute the familial probabilies under the null hypothesis of equally weighted founders
  if(!is.null(subtype_weights)){
    #compute all pedigree specific statistics - under the null hypothesis of equally-weighted founders
    fam_likeGrids_null <- list()
    for (i in 1:length(study_FamIDs)){
      fam_likeGrids_null[[i]] <- test_allcombos_withZero(ped = peds[peds$FamID == study_FamIDs[[i]], ],
                                                         subtypes, tau_grid = tau_grid[1, ], carrier_prob,
                                                         subtype_weights = NULL)
    }

  }


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

  if(!is.null(subtype_weights)){
    fam_likeGrids_null <- lapply(fam_likeGrids_null, function(x){
      remove_invalidConfigs(x)
    })
  }

  # This section computes the estimates tauA and tauB, and computes the LR and probablity of
  # the configuration under the null.  NOTE: that when the "conditioned semi-global" option
  # is selected the probability of the configuration unser the null is the RVS probability.
  #
  # Additionally, when the "conditioned semi-global" option is selected the returned
  # dataframe will contain a dist ID variable, which is used to identify each marginal
  # distribution.  This will be a critical variable for computing p-values under this
  # assumption.
  if(assumption == "global"){
    #compute probabilities and statistics for the global distribution
    global_sharing_byBinID <- condition_globalDist_zeroConfig(likeGrids_byFam = fam_likeGrids,
                                                              famID_index = famIndex,
                                                              tau_grid)

  } else if (assumption == "semi-global"){
    global_sharing_byBinID <- approx_globalDist(likeGrids_byFam = fam_likeGrids,
                                              famID_index = famIndex,
                                              tau_grid)


  } else if (assumption == "conditioned semi-global"){
    global_sharing_byBinID <- conditioned_semiGlobalDist(likeGrids_byFam = fam_likeGrids,
                                                         famID_index = famIndex,
                                                         tau_grid)

  }



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

  global_dist$K = apply(fam_configs, 1, sum)

  #create a list of names of individuals with the genetically compelling subtype
  GCsub_list <- apply(peds[which(peds$available & peds$subtype == subtypes[[1]]), c("FamID", "ID")],
                      1, function(x){paste0(x, collapse = ":")})

  global_dist$K_sub <- sapply(1:nrow(fam_configs), function(x){
    sum(colnames(fam_configs)[fam_configs[x, ] == 1] %in% GCsub_list)
  })


  if (assumption != "conditioned semi-global"){
  #calculate the p-value for the likelihood ratio statistic
  global_dist$LR_pvalue <- sapply(1:nrow(global_dist), function(x){
    sum(global_dist$null_configProb[global_dist$LR >= global_dist$LR[x]],
        na.rm = TRUE)
  })
  }

  #calculate RVS and modified RVS pavlues, NOTE:can only be computed when
  #assumption == "conditioned semi-global"
  if (assumption == "conditioned semi-global"){

    #calculate the p-value for the likelihood ratio statistic
    global_dist$LR_pvalue <- sapply(1:nrow(global_dist), function(x){
      sum(global_dist$null_configProb[global_dist$LR >= global_dist$LR[x]
                                      & global_dist$distID == global_dist$distID[x]],
          na.rm = TRUE)
    })

    global_dist$RVS_pvalue <- sapply(1:nrow(fam_configs), function(x){
      sum(global_dist$null_configProb[global_dist$null_configProb <= global_dist$null_configProb[x]
                                      & global_dist$K >= global_dist$K[x]
                                      & global_dist$distID == global_dist$distID[x]],
          na.rm = TRUE)
    })

    global_dist$modRVS_pvalue <- sapply(1:nrow(fam_configs), function(x){
      sum(global_dist$null_configProb[global_dist$null_configProb <= global_dist$null_configProb[x]
                                      & global_dist$K >= global_dist$K[x]
                                      & global_dist$K_sub >= global_dist$K_sub[x]
                                      & global_dist$distID == global_dist$distID[x]],
          na.rm = TRUE)
    })

  }


  global_dist$binID <- apply(global_dist[, 1:ncol(fam_configs)], 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })

  return(global_dist)

}



#' Compute prioritization statistics and probabilities
#'
#' Computes all considered prioritization statistics and probabilities
#'
#' @inheritParams compute_familyWeights
#' @param peds A ped object containing two or more pedigrees.
#' @param carrier_prob The cumulative carrier probabilty of all cRVs as a group.
#'
#' @return A list of three \code{data.frames}: \code{GlobalDist}, \code{GlobalTransDist}, and \code{LocalDist}; which include:
#' @return \item{RV status indicators}{Binary indicators of RV status for each disease-affected relative, notated as \code{FamID:ID}, where \code{FamID} is the family identification number and \code{ID} is the individual identification number.}
#' @return \item{L_max/L_null}{The likelihood under the alternative and null hypotheses}
#' @return \item{LR}{The likelihood ratio statistic}
#' @return \item{null_configProb}{The null probability of the sharing configuration.}
#' @return \item{tau_A/tau_B}{The values of tau_A and tau_B under the alternative hypothesis.}
#' @return \item{K}{The number of relatives that were observed to carry the RV}
#' @return \item{LR_pvalue}{the p-value if the LR statistic}
#' @return \item{binID}{a unique ID for each sharing configuration.}
#'
#' @export
#' @examples
#' library(RVMethods)
#' data(study_pedigrees)
#'
#' library(SimRVPedigree)
#' plot(study_pedigrees[study_pedigrees$FamID == 304, ])
#' plot(study_pedigrees[study_pedigrees$FamID == 58, ])
#'
#'
#' FBSDists = compute_distributions(ped = study_pedigrees[study_pedigrees$FamID %in% c(58, 304), ],
#'                                  subtypes = c("HL", "NHL"),
#'                                  carrier_prob = 0.002)
#'
#'
compute_distributions <- function(peds, subtypes, carrier_prob = 0.002,
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
  message("Compute pedigree-specific statistics... this may take a few minutes.")
  pb <- txtProgressBar(min = 0, max = length(study_FamIDs), style = 3)
  fam_likeGrids <- list()
  for (i in 1:length(study_FamIDs)){
    fam_likeGrids[[i]] <- test_allcombos_withZero(ped = peds[peds$FamID == study_FamIDs[[i]], ],
                                                  subtypes, tau_grid, carrier_prob,
                                                  subtype_weights)
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

  # This section computes the estimates tauA and tauB, and computes the LR and probablity of
  # the configuration under the null.  NOTE: that when the "conditioned semi-global" option
  # is selected the probability of the configuration unser the null is the RVS probability.
  #
  # Additionally, when the "conditioned semi-global" option is selected the returned
  # dataframe will contain a dist ID variable, which is used to identify each marginal
  # distribution.  This will be a critical variable for computing p-values under this
  # assumption.

  print(paste0("Generate Results: ", Sys.time()))
  #compute probabilities and statistics for the global distribution
  D1_global_sharing_byBinID <- condition_globalDist_zeroConfig(likeGrids_byFam = fam_likeGrids,
                                                            famID_index = famIndex,
                                                            tau_grid)
  #semi-global
#  print(paste0("Start Semi-Global Distribution: ", Sys.time()))
  D2_global_sharing_byBinID <- approx_globalDist(likeGrids_byFam = fam_likeGrids,
                                               famID_index = famIndex,
                                               tau_grid)

  #conditioned semi-global
#  print(paste0("Start Local Distribution: ", Sys.time()))
  D3_global_sharing_byBinID <- conditioned_semiGlobalDist(likeGrids_byFam = fam_likeGrids,
                                                       famID_index = famIndex,
                                                       tau_grid)

  #expand familial configuratiosn and re-lable so that individuals
  #are identified by family ID and individual ID
  fam_configs <- list()
  for(i in famIndex){
    fam_configs[[i]] <- 1*fam_likeGrids[[i]]$configs[match(D1_global_sharing_byBinID[, i],
                                                           fam_likeGrids[[i]]$binIDs), ]
    colnames(fam_configs[[i]]) <- paste0(study_FamIDs[[i]], ":", colnames(fam_configs[[i]]))
  }

  fam_configs <- do.call(cbind, fam_configs)


  #combine the configurations with their statistics
  global_dist <- cbind(fam_configs,
                       D1_global_sharing_byBinID[, -c(1:length(study_FamIDs))])
  global_dist$K = apply(fam_configs, 1, sum)

#  print(paste0("Global P-value: ", Sys.time()))
  if(useK){
    #compute p_value
    global_dist$LR_pvalue <- sapply(1:nrow(global_dist), function(x){
      sum(global_dist$null_configProb[global_dist$LR >= global_dist$LR[x]
                                      & global_dist$K >= global_dist$K[x]],
          na.rm = TRUE)
    })
  } else{
    global_dist$LR_pvalue <- sapply(1:nrow(global_dist), function(x){
      sum(global_dist$null_configProb[global_dist$LR >= global_dist$LR[x]],
          na.rm = TRUE)
    })
  }


  #combine the configurations with their statistics
  semiglobal_dist <- cbind(fam_configs,
                           D2_global_sharing_byBinID[, -c(1:length(study_FamIDs))])
  semiglobal_dist$K = apply(fam_configs, 1, sum)


#  print(paste0("Semi-Global P-value: ", Sys.time()))
  if(useK){
    #compute p_value
    semiglobal_dist$LR_pvalue <- sapply(1:nrow(global_dist), function(x){
      sum(semiglobal_dist$null_configProb[semiglobal_dist$LR >= semiglobal_dist$LR[x]
                                          & semiglobal_dist$K >= semiglobal_dist$K[x]],
          na.rm = TRUE)
    })
  } else {
    #compute p_value
    semiglobal_dist$LR_pvalue <- sapply(1:nrow(global_dist), function(x){
      sum(semiglobal_dist$null_configProb[semiglobal_dist$LR >= semiglobal_dist$LR[x]],
          na.rm = TRUE)
    })
  }


  #combine the configurations with their statistics
  condsemiglobal_dist <- cbind(fam_configs,
                               D3_global_sharing_byBinID[, -c(1:length(study_FamIDs))])

  condsemiglobal_dist$K = apply(fam_configs, 1, sum)

  #create a list of names of individuals with the genetically compelling subtype
  GCsub_list <- apply(peds[which(peds$available & peds$subtype == subtypes[[1]]), c("FamID", "ID")],
                      1, function(x){paste0(x, collapse = ":")})

  condsemiglobal_dist$K_sub <- sapply(1:nrow(fam_configs), function(x){
    sum(colnames(fam_configs)[fam_configs[x, ] == 1] %in% GCsub_list)
  })

#  print(paste0("Local P-value: ", Sys.time()))
  if(useK){
    #calculate the p-value for the likelihood ratio statistic
    condsemiglobal_dist$LR_pvalue <- sapply(1:nrow(condsemiglobal_dist), function(x){
      sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$LR >= condsemiglobal_dist$LR[x]
                                              & condsemiglobal_dist$K >= condsemiglobal_dist$K[x]
                                              & condsemiglobal_dist$distID == condsemiglobal_dist$distID[x]],
          na.rm = TRUE)
    })
  } else {
    condsemiglobal_dist$LR_pvalue <- sapply(1:nrow(condsemiglobal_dist), function(x){
      sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$LR >= condsemiglobal_dist$LR[x]
                                              & condsemiglobal_dist$distID == condsemiglobal_dist$distID[x]],
          na.rm = TRUE)
    })
  }


#  print(paste0("RVS-Based P-value: ", Sys.time()))
  condsemiglobal_dist$RVS_pvalue <- sapply(1:nrow(fam_configs), function(x){
    sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$null_configProb <= condsemiglobal_dist$null_configProb[x]
                                    & condsemiglobal_dist$K >= condsemiglobal_dist$K[x]
                                    & condsemiglobal_dist$distID == condsemiglobal_dist$distID[x]],
        na.rm = TRUE)
  })

  condsemiglobal_dist$modRVS_pvalue <- sapply(1:nrow(fam_configs), function(x){
    sum(condsemiglobal_dist$null_configProb[condsemiglobal_dist$null_configProb <= condsemiglobal_dist$null_configProb[x]
                                    & condsemiglobal_dist$K >= condsemiglobal_dist$K[x]
                                    & condsemiglobal_dist$K_sub >= condsemiglobal_dist$K_sub[x]
                                    & condsemiglobal_dist$distID == condsemiglobal_dist$distID[x]],
        na.rm = TRUE)
  })


  global_dist$binID <- apply(global_dist[, 1:ncol(fam_configs)], 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })

  semiglobal_dist$binID <- apply(semiglobal_dist[, 1:ncol(fam_configs)], 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })

  condsemiglobal_dist$binID <- apply(condsemiglobal_dist[, 1:ncol(fam_configs)], 1, function(x){
    base::strtoi(paste0(x, collapse = ""), base = 2)
  })


  return(list(GlobalDist = global_dist,
              GlobalTransDist = semiglobal_dist,
              LocalDist = condsemiglobal_dist))

}


