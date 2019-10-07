library(testthat)
library(RVMethods)
context("compute_familyWeights")
#ADD a test for randomly simulated pedigree


data("study_pedigrees")
sampID = sample(unique(study_pedigrees$FamID), size = 1)
famWeights <- compute_familyWeights(study_pedigrees[study_pedigrees$FamID == sampID, ],
                                    subtypes = c("HL", "NHL"))
head(famWeights)

test_that("RVSharing probabilty output meets expectations", {

  expect_equal(sum(famWeights$RVsharing), 1)

  expect_true(all(famWeights$RVsharing > 0))

  expect_true(all(famWeights$RVsharing <= 1))

})

test_that("LR statisic behaves as expected", {

  expect_true(all(famWeights$LR >= 1))

})

test_that("Likelihood max/null behaves as expected", {

  expect_true(all(famWeights$L_max  >= 0) & all(famWeights$L_max <= 1))
  expect_true(all(famWeights$L_null >= 0) & all(famWeights$L_null <= 1))

})

test_that("tau ests behave as expected", {

  expect_true(all(famWeights$tau_A >= 0.5) & all(famWeights$tau_A <= 1))
  expect_true(all(famWeights$tau_B >= 0.5) & all(famWeights$tau_B <= 1))

  expect_true(all(famWeights$tau_A >= famWeights$tau_B))

})

