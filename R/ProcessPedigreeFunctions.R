#' Assign generation number based on oldest founder
#'
#' @param x an object of class ped
#'
#' @return a list of generation numbers for pedigree members, in the order listed in \code{x}.
#' @importFrom kinship2 kindepth
#' @keywords internal
assign_gen <- function(x){
  Gen <- NA
  mates <- cbind(x$dadID, x$momID)
  mates <- unique(mates)
  #remove rows with only NAs, these represent the founders
  mates <- mates[!is.na(mates[, 1]), ]

  #get kindepth and set Gen to kindepth when kindepth is non-zero
  kd <- kinship2::kindepth(ped2pedigree(x))
  Gen[kd != 0] <- kd[kd != 0]

  # NOTE: when the pedigree reduces to only 2 generations,
  # e.g. as in a family-trio, mates will NOT be a matrix.
  if(class(mates) == "matrix"){
    for(i in 1:nrow(mates)){
      mate_gens <-  kd[x$ID %in% mates[i, ]]
      Gen[x$ID %in% mates[i, ]] <- max(mate_gens)
    }
  } else {
    Gen[x$ID %in% mates] <- 0
  }

  return(Gen + 1)
}

#' Reduce pedigree to contain only disease-affecteds, obligate carriers, and founders.
#'
#' Reduce pedigree to contain only disease-affecteds, obligate carriers, and founders.
#'
#' @param ped_file An object of class \code{ped}, inherits from a data.frame. The pedigree to reduce to disease-affected relatives.
#'
#' @return \code{retA_ped} A pedigree containing only affected members, obligate carriers, and founders.
#' @keywords internal
affected_onlyPed = function(ped_file){

  #assign generation number if not included in ped_file
  if(!"Gen" %in% colnames(ped_file)){
    ped_file$Gen <- assign_gen(ped_file)
  }

  #create new ped file with affecteds only
  retA_ped <- ped_file[ped_file$affected & ped_file$available, ]

  if (nrow(retA_ped) == 0) {
    stop(paste0("Disease-affected relatives are not present in pedigree with FamID ",
                   sep = "", ped_file$FamID[1]))
  } else {
    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(retA_ped$dadID,
                               retA_ped$ID[which(retA_ped$sex == 0)])
      readd_dad <- retA_ped$dadID[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(retA_ped$momID,
                               retA_ped$ID[which(retA_ped$sex == 1)])
      readd_mom <- retA_ped$momID[miss_mom]
      readd_mom <- unique(readd_mom[!is.na(readd_mom)])

      #check to see if we need to readd anyone
      if (length(c(readd_dad, readd_mom)) == 0) {
        d <- 1
      } else {
        #Now pull the rows containing the required parents
        # from the original ped_file
        readd <- ped_file[which(ped_file$ID %in% c(readd_dad, readd_mom)), ]

        #combine with affected ped file
        retA_ped <- rbind(retA_ped, readd)
      }
    }
  }

 return(retA_ped)
}

#' Assign tau to obligate carriers
#'
#' @param id the id (index) of the obligate carrier
#' @param parents the parent matrix
#' @param typeA the list of indicies of the type A's
#' @param typeB the list of indicies of the type B's
#'
#' @return character, either "A" or "B
#' @keywords internal
assign_obgcarrier_tau <- function(id, parents, typeA, typeB){
  #find descendents
  descendents <- findDescendents(id, parents)

  #count number of typeA descendents
  numA <- length(intersect(descendents, typeA))

  #count number of typeB descendents
  numB <- length(intersect(descendents, typeB))

  if (numA == 0){
    #if
    obg_tau <- "B"
  } else if (numB == 0){
    obg_tau <- "A"
  } else if (numA == 0 & numB == 0){
    stop("error: this person is not an obligate carrier")
    #remove this when you are sure this will never happen
  } else if (numA/numB >= 1) {
    obg_tau <- "A"
  } else if (numA/numB < 1) {
    obg_tau <- "B"
  }

  return(obg_tau)
}


#' Process Pedigree
#'
#' @inheritParams compute_sharingProb
#' @importFrom SimRVPedigree ped2pedigree
#'
#' @return A list of information for the pedigree
#' @keywords internal
processPed <- function(ped, carriers, subtypes){
  #if Gen is included but some generation numbers are missing,
  #remove this variable entirely
  if("Gen" %in% colnames(ped) & any(is.na(ped$Gen))){
    ped <- ped[, -which(colnames(ped) == "Gen")]
  }

  #re-label the subtypes generically, the genetically-compelling subtype is
  #labeled "A" and the other subtype is labelled "B".
  ped$subtype <- ifelse(ped$subtype == subtypes[1], "A",
                        ifelse(ped$subtype == subtypes[2], "B", NA))

  #reduce pedigree to contain only disease-affected relatives,
  #their parents, and common ancestors
  ped <- affected_onlyPed(ped)
  ped <- ped[order(ped$ID), ]

  # #store index of type A and type B affecteds
  typeA <- which(ped$subtype == "A")
  typeB <- which(ped$subtype == "B")

  #convert to a pedigree object
  pdgree <- ped2pedigree(ped)


  # relabel subject IDs and store basic pedigree information
  origID <- pdgree$id
  #relabelling IDs sequentially for easy reference
  pdgree$id <- 1:length(pdgree$id)
  parents <- sapply(pdgree$id, function(i) c(pdgree$findex[i], pdgree$mindex[i]))
  founders <- which(parents[1, ] == 0)
  nonfounders <- which(parents[1, ] != 0)
  finalDescendants <- which(!(pdgree$id %in% c(parents[1,], parents[2,])))

  # NOTE: pedigrees simulated by SimRVPedigree and
  # converted by ped2pedigree will have additional info
  # in this matrix. We want to retain only affection statuses
  if (inherits(pdgree$affected,"matrix")) {
    pdgree$affected <- pdgree$affected[, 1]
  }

  #store the IDs of the disease-affected relatives
  affected <- which(pdgree$affected == 1) #affected provided

  if (!is.null(carriers)) {
    #identify indices of carriers, when provided.
    carriers <- which(origID %in% carriers)
  } else {
    #otherwise default to disease-affected relatives
    carriers <- affected
  }


  # check pedigree is valid
  #if (sum(affected %in% founders) > 0)
  #    stop('some founders are affected')
  if (length(affected) < 2)
    stop('need at least 2 affected subjects')
  if (sum(!carriers %in% affected) > 0)
    stop('carriers must be a subset of affected')


  #find obligate carriers, i.e. non-founders who are unaffected
  obgCarriers <- setdiff(x =  pdgree$id, y = c(founders, affected))

  #assign taus to the obligate carriers
  obg_taus <- unlist(lapply(obgCarriers, function(x){assign_obgcarrier_tau(x, parents, typeA, typeB)}))

  #append the type A obligate carriers to the typeA list,
  # and the type B obligate carriers to the typeB list
  typeA <- c(typeA, obgCarriers[which(obg_taus == "A")])
  typeB <- c(typeB, obgCarriers[which(obg_taus == "B")])

  # save info in list
  return(list('origID' = origID, 'ped' = pdgree, 'parents' = parents,
              'founders' = founders, 'affected' = affected,
              'size' = length(pdgree$id), 'carriers' = carriers,
              'id' = pdgree$id, 'finalDescendants' = finalDescendants,
              'typeA' = typeA, 'typeB' = typeB))
}
