
#' Compute sharing probabilities for all configuration (includes the zero configuration) and initialization of tau
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
#'                                       tau_grid = make_tauGrid())
#'
#' myCombos
#'
#' rowSums(myCombos$like_mat)
#'
test_allcombos_noCarProb <- function(ped, subtypes, tau_grid,
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




#' Inocrporate the carrier probability into the probability of configuraiton
#'
#' Compute sharing probabilities for every global configuration (includes the zero configuration) and initialization of tau.  The probability of the zero configuration will be used to condition on the event that the RV was observed in at least one study member.  This function requires an assummed value for the carrier probability of the RV
#'
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
allcombos_addCarProb <- function(ped, like_mat, carrier_prob){

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


  return(like_mat)
}




#' Compute all distributions - used for power simulation - avoids haveing to re-run this function multiple times
#'
#' used for power simulation - avoids haveing to re-run this function multiple times
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
quick_compute_distributions <- function(peds, subtypes, carrier_prob = 0.002,
                                        tau_increment = 0.05, file_path, ped_indicies,useK = FALSE){


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

    load(paste0(file_path, ped_indicies[[i]], ".rda", sep = ""))
    fam_likeGrids[[i]] = fam_res

    fam_likeGrids[[i]][[1]] <- allcombos_addCarProb(ped = peds[peds$FamID == study_FamIDs[[i]], ],
                                                    like_mat = fam_likeGrids[[i]][[1]],
                                                    carrier_prob)
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

  print(paste0("Start Global Distribution: ", Sys.time()))
  #compute probabilities and statistics for the global distribution
  D1_global_sharing_byBinID <- condition_globalDist_zeroConfig(likeGrids_byFam = fam_likeGrids,
                                                               famID_index = famIndex,
                                                               tau_grid)
  #semi-global
  print(paste0("Start Semi-Global Distribution: ", Sys.time()))
  D2_global_sharing_byBinID <- approx_globalDist(likeGrids_byFam = fam_likeGrids,
                                                 famID_index = famIndex,
                                                 tau_grid)

  #conditioned semi-global
  print(paste0("Start Local Distribution: ", Sys.time()))
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

  print(paste0("Global P-value: ", Sys.time()))
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


  print(paste0("Semi-Global P-value: ", Sys.time()))
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

  print(paste0("Local P-value: ", Sys.time()))
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


  print(paste0("RVS-Based P-value: ", Sys.time()))
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


  return(list(globDist = global_dist,
              semiglobDist = semiglobal_dist,
              condsemiglobDist = condsemiglobal_dist))

}



