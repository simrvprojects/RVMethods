library(testthat)
library(RVMethods)
context("compute_binID")

test_that("when RVstatuses are 0/1 always returns a non-zero binID", {

  for (l in 1:30){
    my_n = sample(x = 2:8, size = 1)
    geno_vec <- c()
    for (i in 1:my_n){
      #sample genotypes for each individual
      geno_vec <- c(geno_vec,
                    c(sample(x = c(0, 1), size = 2, replace = FALSE)))
    }

    expect_true(compute_binID(geno_vec) > 0)
  }
})

test_that("when a person carries multiple copies returns a zero binID", {

  for (l in 1:30){
    my_n = sample(x = 2:8, size = 1)
    my_multi_carrier = sample(x = 1:my_n, size = 1)
    geno_vec <- c()
    for (i in 1:my_n){
      if(i != my_multi_carrier){
        #sample genotypes for each individual
        geno_vec <- c(geno_vec,
                      c(sample(x = c(0, 1), size = 2, replace = FALSE)))
      } else {
        geno_vec <- c(geno_vec, c(1, 1))
      }
    }
    geno_vec

    expect_true(compute_binID(geno_vec) == 0)
  }
})
