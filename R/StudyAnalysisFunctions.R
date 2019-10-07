#' Determine the families that share variants in the simulated data
#'
#' @param binID_mat The matrix of binIDs for the families in the study. The $(i,j)^{th}$ entry represents the binary representation of the configuration for the $j^{th}$ variant of the $i^{th}$ family.
#'
#' @return A list of sharing combinations, i.e. a list of FamID lists.
#' @keywords internal
find_famSharing <- function(binID_mat){

  #return all families that share each variant
  fam_combos_by_SNV <- lapply(1:ncol(binID_mat), function(x){
    which(binID_mat[, x] != 0)
  })

  #determine which SNVs are shared by multiple families
  fam_combos_is_sharing <- sapply(fam_combos_by_SNV, function(x){
    length(x) > 1
  })

  #determine all unique combinations of families that share SNVs
  fam_combos <- unique(fam_combos_by_SNV[fam_combos_is_sharing])

  return(fam_combos)
}


#' Remove invalid configurations from the items returned by test_allcombos
#'
#' @param test_allCombos_result The items returned by \code{test_allcombos}, i.e. a list containing \code{like_mat}, \code{configs}, and \code{binIDs}.
#'
#' @return A list containing \code{like_mat}, \code{configs}, and \code{binID}, with the same format as \code{test_allcombos}, but with invalid configurations removed.
#' @keywords internal
remove_invalidConfigs <- function(test_allCombos_result){
  remove_configs <- which(sapply(1:ncol(test_allCombos_result$like_mat), function(x){
    all(test_allCombos_result$like_mat[, x] == 0)
  }))

  if (length(remove_configs) >= 1){
    return(list(like_mat = test_allCombos_result$like_mat[, -remove_configs],
                configs = test_allCombos_result$configs[-remove_configs, ],
                binIDs = test_allCombos_result$binIDs[-remove_configs]))
  } else {
    return(list(like_mat = test_allCombos_result$like_mat,
                configs = test_allCombos_result$configs,
                binIDs = test_allCombos_result$binIDs))
  }
}

#' Compute weights for variants shared among disease-affected study participants
#'
#' Compute weights for variants shared among disease-affected study participants
#'
#' The parameter space of  \eqn{\tau_a} and \eqn{\tau_b} is given by
#' \eqn{\Omega_{con} = ((\tau_{a}, \tau_{b}): 0.5 \le \tau_{b} < \tau_{a} \le 1 )}.
#'
#'
#'
#' @param famStudy_obj An object of class \code{famStudy}, i.e. returned by the \code{sim_RVstudy} function (included in the \code{SimRVSequences} package).
#' @param tau_increment Numeric. The width of the grid for the taus
#' @param subtypes  A list of length 2. Contains character labels for the two subtypes that occur in the pedigrees contained in the \code{famStudy_obj}. The label for the more the genetically-complelling subtype must be listed first.
#' @param subtype_weights A vector of length 2.  The weights to applied to individuals with the subtypes listed in \code{subtypes}. By default, \code{subtype_weights = NULL} so that no weights are applied.
#' @importFrom Matrix colSums
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'
#' @return A dataframe cataloging all of the SNVs, with weights included as columns w_LR, w_RVS, and w_RVS2.  See \code{\link{compute_familyWeights}} documentation for weight definitions.
#'
#' @examples
#' library(RVMethods)
#'
#' # load example pedigrees
#' data(study_pedigrees)
#' head(study_pedigrees)
#'
#' #load example population haplotype data
#' data(pop_fdat)
#' str(pop_fdat)
#'
#' #simulate sequence data for pedigrees with the SimRVSequences package
#' library(SimRVSequences)
#' ped_seq = sim_RVstudy(ped_files = study_pedigrees,
#'                       SNV_data = pop_fdat)
#' class(ped_seq)
#' # for more information about the sim_RVstudy function and it's output,
#' # i.e. objects of class famStudy, please refer to the vigentte
#' # for the SimRVSequences R package
#'
#'
#' # compute weights for SNVs shared among the
#' # disease-affected relatives in the study
#' SNV_weights = compute_studyWeights(famStudy_obj = ped_seq,
#'                                    subtypes = c("HL", "NHL"))
#'
#' head(SNV_weights)
#'
compute_studyWeights <- function(famStudy_obj, subtypes,
                                 tau_increment = 0.05,
                                 subtype_weights = NULL){

  #create tau grid
  #NOTE: make_tauGrid now orders taus to save time later.
  tau_grid <- make_tauGrid(increment_width = tau_increment, constrained = TRUE)
  #determine the FamIDs for all families in the study
  study_FamIDs <- unique(famStudy_obj$ped_files$FamID)


  message("Processing pedigree and SNV data")
  # since we only compute weights for variants shared by diseased relatives
  # we subset the data to for these individuals.
  # NOTE: currently, the data includes SNVs which may only be carried by
  # unaffected relatives

  #remove any SNVs that are introduced by more than one founder
  #or that are not carried by disease affected relatives
  famStudy_obj <- preProcessSNVs_studyDat(famStudy_obj)

  # The item named "haplo_map" in the famStudy object contains a variable
  # named affected, which identified the rows of data belonging to a disease affected relative.
  Affected_ped_haplos <- famStudy_obj$ped_haplos[famStudy_obj$haplo_map$affected, ]
  Affected_haplo_map <- famStudy_obj$haplo_map[famStudy_obj$haplo_map$affected, ]
  SNV_map <- famStudy_obj$SNV_map



  # tabulate binIDs for the study
  binMat <- tabulate_binID_by_SNV(Ahaplos = Affected_ped_haplos,
                                  Amap = Affected_haplo_map,
                                  FamIDs = study_FamIDs)

  # Occasionally, a SNV will recieve a binIDs of zero when the SNV is seen in an
  # affected relative contained in the study.  This occurs when the relative
  # carries multiple copies of the SNV, which is not possible under our model
  # as this requires the SNV to be introduced by multiple founders.
  # Check to ensure that this has not occured in binMat, if
  # it has, we need to remove these SNVs from the data.
  any_invalid_binCols <- apply(binMat, 2, function(x){all(x == 0)})
  if (any(any_invalid_binCols)){
    #NOTE: this entire section can likely be removed.
    #Since adding the preProcessSNVs_studyDat function
    #there should not be any configuraitons to remove
    #HOWEVER, keeping this section in for now, to check that
    #preProcessSNVs_studyDat is performing as expected.
    warning("PreProcess missed an SNV *please tell CN if you see this warning*")
    remove_SNVs <- which(any_invalid_binCols)

    #remove from haplotype matrix
    Affected_ped_haplos <- Affected_ped_haplos[, -c(remove_SNVs)]
    #remove from SNV_map
    SNV_map <- SNV_map[-c(remove_SNVs), ]
    #remove from binMat
    binMat <- binMat[, -c(remove_SNVs)]
  }
  #rename colIDs in SNV_map
  SNV_map$colID <- c(1:nrow(SNV_map))

  #TODO: for the love of all that is good and holy,
  #please get rid of the variable colID in the SNVmap
  #This needs to be done in the SimRVSequences package.


  message("Compute pedigree-specific statistics... this may take a few minutes.")
  pb <- txtProgressBar(min = 0, max = length(study_FamIDs), style = 3)
  fam_likeGrids <- list()
  for (i in 1:length(study_FamIDs)){
    fam_likeGrids[[i]] <- test_allcombos(ped = famStudy_obj$ped_files[famStudy_obj$ped_files$FamID == study_FamIDs[[i]], ],
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
                       ped = famStudy_obj$ped_files[famStudy_obj$ped_files$FamID == study_FamIDs[x], ],
                       subtypes, subtype_weights)
  })

  #determine all combinations of families that share SNVs
  #NOTE: this list is a list of FamIDs
  #not their position in fam_likeGrids
  global_famCombos_index <- find_famSharing(binMat)
  if(length(global_famCombos_index) > 0){
    message("Computing statistics for global configurations")
    global_Weights = list()

    for(i in 1:length(global_famCombos_index)){
      global_Weights[[i]] <- compute_global_weight(likeGrids_byFam = fam_likeGrids[global_famCombos_index[[i]]],
                                                   weights_byFam = fam_sharingWeights[global_famCombos_index[[i]]],
                                                   famID_index = global_famCombos_index[[i]],
                                                   tau_grid)
    }
  }

  message("Tabulating weights")
  #initialize weight variables in SNV_map
  SNV_map$LR = NA
  SNV_map$w_LR = NA
  SNV_map$w_RVS = NA
  SNV_map$w_RVS2 = NA

  for(i in 1:nrow(SNV_map)){
   #determine which familes share the RV
   loop_fams_index = which(binMat[, i] != 0)
   loop_bins = binMat[loop_fams_index, i]

   if(length(loop_fams_index) == 1){
     #THIS OPTION IS USED WHEN ONLY A SINGLE FAMILY CARRIES THE VARIANT
     #check to see if this binID is a valid configuration, if not
     #all weights will be set to zero
     if(loop_bins %in% fam_sharingWeights[[loop_fams_index]]$binID){
       #If only a single family segregates the RV, we can simply pull their weights
       #from the results that were computed for the family. These are stored in the
       #fam_sharingWeights list
       SNV_map[i, c("LR", "w_LR", "w_RVS", "w_RVS2")] <- fam_sharingWeights[[loop_fams_index]][fam_sharingWeights[[loop_fams_index]]$binID == loop_bins, c("LR", "w_LR", "w_RVS", "w_RVS2")]
     } else {
       #check to see if this binID is a valid configuration, if not
       #weights will be set to zero
       SNV_map[i, c("LR", "w_LR", "w_RVS", "w_RVS2")] <- c(0, 0, 0, 0)
     }
   } else if (length(loop_fams_index) > 1){
     #THIS OPTION IS USED FOR GLOBAL CONFIGURATIONS SPANNING MULTIPLE FAMILIES


     #find the index of the global configuration in global_famCombos_index
     global_set_index <- which(unlist(lapply(global_famCombos_index, identical, loop_fams_index)))
     if(length(global_set_index) > 1){
        stop("global index broken **tell CN if you see this message**")
     } else if (length(global_set_index) == 1) {
       Config_line <- which(apply(global_Weights[[global_set_index]][, 1:length(loop_bins)], 1,
                                  function(x){all(x == loop_bins)}))
       if(length(Config_line) != 0){

         #store the LR value
         SNV_map$LR[i] <- global_Weights[[global_set_index]][Config_line, "LR"]


         #compute transmission based-weight
         SNV_map$w_LR[i] <- compute_transmission_weight(all_stats = global_Weights[[global_set_index]],
                                                        config_index = Config_line)

         #compute RVS weight
         SNV_map$w_RVS[i] <- compute_RVS_weight(all_stats = global_Weights[[global_set_index]],
                                               config_index = Config_line)

         #compute modified RVS weight
         SNV_map$w_RVS2[i] <- compute_modified_RVS_weight(all_stats = global_Weights[[global_set_index]],
                                                         config_index = Config_line)

       } else {
         SNV_map[i, c("LR", "w_LR", "w_RVS", "w_RVS2")] <- c(0, 0, 0, 0)
       }
     }

   }
  }

  #identify familial crvs as well as which family simualted the crv
  loop_cRVs <- as.character(unique(famStudy_obj$haplo_map$FamCRV))
  SNV_map$FamCRV <- ifelse(SNV_map$marker %in% loop_cRVs, TRUE, FALSE)
  SNV_map$FamIDs <- sapply(1:nrow(SNV_map), function(x){
    ifelse(SNV_map$FamCRV[[x]] == TRUE,
           paste0(unique(famStudy_obj$haplo_map$FamID[which(famStudy_obj$haplo_map$FamCRV == SNV_map$marker[[x]])]),
                  collapse = ", "),
           NA)
  })


  return(SNV_map)

  # return(list(SNV_map = SNV_map,
  #             affected_haplos = Affected_ped_haplos,
  #             affected_haplo_map = Affected_haplo_map,
  #             ped_files = famStudy_obj$ped_files))
}




#' Tabulate each families configurations stored in binary format
#'
#' @param Ahaplos the haplotype matrix for the disease-affected relatives only
#' @param Amap the pedigree to haplotype matrix map for the disease affected relatives only
#' @param FamIDs a list of FamIDs for pedigrees in the study
#'
#' @return A matrix of binIDs.  Rows are families, columns are SNVs, the $(i,j)^{th}$ entry is the binID for the sharing configuration for the ith family and the jth SNV.
#' @keywords internal
#'
tabulate_binID_by_SNV <- function(Ahaplos, Amap, FamIDs){
  #want to populate a matrix with familial binIDs (i.e. configs)
  #for each SNV in the simulated data
  #ROWS = families
  #Columns = SNVs
  #i,jth entry is the binID for family i and SNV j


  bin_mat <- matrix(NA, nrow = length(FamIDs), ncol = ncol(Ahaplos))
  for (i in 1:length(FamIDs)) {
    loop_dat <- Ahaplos[Amap$FamID == FamIDs[[i]], ]
    bin_vec <- apply(loop_dat, 2, function(x){
      compute_binID(x)
    })

    bin_mat[i, ] <- bin_vec

  }

  return(bin_mat)
}

#' Compute the binID for the observed familial genotypes
#'
#' Compute the binID for the observed familial genotypes
#'
#' Note that under our model, cRVs may only be introduced by a single founder. If an affected relative has 2 copies of the variant then it must have been introduced by multiple founders.  When this occurs, we set the  the binID to 0 since these variants are not applicable under our model.
#'
#' @param affected_genotypes_at_marker The genotype for the disease affected relatives at a marker.
#'
#' @return The binID, which is the base 2 representation of the sharing configuration.
#' @keywords internal
#'
compute_binID <- function(affected_genotypes_at_marker){
  binID <- base::strtoi(paste0(rowSums(matrix(affected_genotypes_at_marker,
                                              nrow = length(affected_genotypes_at_marker)/2,
                                              byrow = TRUE)),
                               collapse = ""), base = 2)

  #NA binIDs are not possible under our model, for example
  #if a relative has two copied of a variant the configuration
  #will result in an NA binID.  We set these to zero, which
  #will result in a zero weight.
  ifelse(!is.na(binID), binID, 0)
}
