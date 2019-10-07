#' Process the SNVs in the famStudy object
#'
#' This willl become an internal function. This is a processing step whereby we filter SNVs for specific exclusion criteria.
#'
#' In this pre-processing step we filter any variants that are introduced by multiple founders in the same pedigree or that are introduced by the same founder mulitple times as this violates our model assumptions.  Additionally, we remove any SNVs that are not observed in at least one disease-affected relative in the study.
#'
#' @param famStudy_obj An object of class \code{famStudy}.
#'
#' @return An object of class \code{famStudy}. The SNVs in the returned \code{famStudy} object have been filtered based on exclusion criteria.
#' @keywords internal
#'
preProcessSNVs_studyDat <- function(famStudy_obj){

  FamIDs <- unique(famStudy_obj$haplo_map$FamID)

  #determine the number of affected ped family
  nAff_perFam <- sapply(1:length(FamIDs), function(x){
    sum(famStudy_obj$ped_files$affected[famStudy_obj$ped_files$FamID == FamIDs[[x]]])
  })

  #remove SNVs that are not carried by at least one DISEASE-AFFECTED relative
  #in the study

  Aff_dat <- famStudy_obj$ped_haplos[famStudy_obj$haplo_map$affected, ]
  remove_SNVs <- colSums(Aff_dat) == 0

  if (length(which(remove_SNVs)) > 0){
    #remove SNVs from ped_haplos
    famStudy_obj$ped_haplos <- famStudy_obj$ped_haplos[, !remove_SNVs]

    #remove SNVs from SNV_map
    famStudy_obj$SNV_map <- famStudy_obj$SNV_map[!remove_SNVs, ]
    famStudy_obj$SNV_map$colID <- 1:nrow(famStudy_obj$SNV_map)
    row.names(famStudy_obj$SNV_map) = NULL

  }


  #remove SNVs introduced by multiple founders
  #store famID and individual IDs of all founders
  founderIDs <- famStudy_obj$ped_files[is.na(famStudy_obj$ped_files$dadID),
                                       c("FamID", "ID")]

  founderIDs_by_fam <- lapply(1:length(FamIDs), function(x){
    famStudy_obj$ped_files$ID[is.na(famStudy_obj$ped_files$dadID) &
                             famStudy_obj$ped_files$FamID == FamIDs[[x]]]
  })

  # founder_rows <- match(paste(founderIDs$FamID, founderIDs$ID),
  #                       paste(famStudy_obj$haplo_map$FamID, famStudy_obj$haplo_map$ID))
  #
  # #this is a sparse matrix of founder haplotypes
  # #ROWs are people, each person has two rows
  # #COLUMNs are SNVs
  # founder_haps <- famStudy_obj$ped_haplos[founder_rows, ]


  remove_SNVs <- unlist(lapply(1:length(FamIDs), function(x){
    which(colSums(famStudy_obj$ped_haplos[famStudy_obj$haplo_map$FamID == FamIDs[[x]] &
                                            famStudy_obj$haplo_map$ID %in% founderIDs_by_fam[[x]], ]) > 1)
  }))

#
#
#   #identify any SNVs that are introduced by the same founder more than once
#   remove_SNVs <- apply(founder_haps, 2, function(x){
#     any(rowSums(matrix(x, ncol = 2, byrow = TRUE)) == 2)
#   })

  if (length(remove_SNVs) > 0){
    #remove SNVs from ped_haplos
    famStudy_obj$ped_haplos <- famStudy_obj$ped_haplos[, -c(remove_SNVs)]

    #remove SNVs from SNV_map
    famStudy_obj$SNV_map <- famStudy_obj$SNV_map[-c(remove_SNVs), ]
    famStudy_obj$SNV_map$colID <- 1:nrow(famStudy_obj$SNV_map)
    row.names(famStudy_obj$SNV_map) = NULL
  }





  # #determine the FamIDs of families who share each retained SNV
  # prodSpace_size <- sapply(1:ncol(Aff_share_df), function(x){
  #   prod(2^nAff_perFam[Aff_share_df[, x]])
  # })


  #return(list(maxCRV = prodSpace_size[which(famStudy_obj$SNV_map$is_CRV)],
  #            maxSNV = max(prodSpace_size)))

  return(famStudy_obj)

}
