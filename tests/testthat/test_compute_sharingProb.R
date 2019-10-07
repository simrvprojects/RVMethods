library(testthat)
library(RVMethods)
library(SimRVPedigree)
context("compute_sharingProb")


# create some test families
test_fam1 <- new.ped(data.frame(FamID = rep(133, 11),
                                ID = seq(1, 11, by  = 1),
                                dadID = c(NA, NA, NA, 1, 1, 1, NA, 3, 3, 3, 6),
                                momID = c(NA, NA, NA, 2, 2, 2, NA, 4, 4, 4, 7),
                                sex = c(0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1),
                                affected = c(F, F, F, F, T, F, F, T, T, T, T),
                                subtype = c(NA, NA, NA, NA, "B", NA, NA, "A", "B", "A", "A")))
#number affected
nt1 <- sum(test_fam1$affected, na.rm = TRUE)
#store the IDs of the affected relatives
aID1 <- test_fam1$ID[test_fam1$affected]
#determine all possible sharing configurations
config1 = as.matrix(expand.grid(rep(list(c(T, F)), nt1)))
config1 = config1[rowSums(config1) > 0, ]
colnames(config1) = as.character(aID1)


test_fam2 <- new.ped(data.frame(FamID = rep(133, 11),
                                ID = seq(1, 11, by  = 1),
                                dadID = c(NA, NA, NA, 1, 1, 1, NA, 3, 3, 3, 6),
                                momID = c(NA, NA, NA, 2, 2, 2, NA, 4, 4, 4, 7),
                                sex = c(0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1),
                                affected = c(F, F, F, F, T, F, F, F, F, T, T),
                                subtype = c(NA, NA, NA, NA, "B", NA, NA, NA, NA, "A", "A")))
#number affected
nt2 <- sum(test_fam2$affected, na.rm = TRUE)
#store the IDs of the affected relatives
aID2 <- test_fam2$ID[test_fam2$affected]
#determine all possible sharing configurations
config2 = as.matrix(expand.grid(rep(list(c(T, F)), nt2)))
config2 = config2[rowSums(config2) > 0, ]
colnames(config2) = as.character(aID2)


test_fam3 <- new.ped(data.frame(FamID = rep(133, 11),
                                ID = seq(1, 11, by  = 1),
                                dadID = c(NA, NA, NA, 1, 1, 1, NA, 3, 3, 3, 6),
                                momID = c(NA, NA, NA, 2, 2, 2, NA, 4, 4, 4, 7),
                                sex = c(0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1),
                                affected = c(F, F, F, T, T, F, F, F, F, T, T),
                                subtype = c(NA, NA, NA, "A", "B", NA, NA, NA, NA, "A", "A")))
#number affected
nt3 <- sum(test_fam3$affected, na.rm = TRUE)
#store the IDs of the affected relatives
aID3 <- test_fam3$ID[test_fam3$affected]
#determine all possible sharing configurations
config3 = as.matrix(expand.grid(rep(list(c(T, F)), nt3)))
config3 = config3[rowSums(config3) > 0, ]
colnames(config3) = as.character(aID3)


test_that("probability is always non-negative", {

  expect_true(all(apply(config1, 1, function(x){
    compute_sharingProb(ped = test_fam1, subtypes = c("A", "B"), carriers = aID1[x], tau = c(0.5, 0.5))
  }) >= 0))

  expect_true(all(apply(config2, 1, function(x){
    compute_sharingProb(ped = test_fam2, subtypes = c("A", "B"), carriers = aID2[x], tau = c(0.5, 0.5))
  }) >= 0))

  expect_true(all(apply(config3, 1, function(x){
    compute_sharingProb(ped = test_fam3, subtypes = c("A", "B"), carriers = aID3[x], tau = c(0.5, 0.5))
  }) >= 0))

})
