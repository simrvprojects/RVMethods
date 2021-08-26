#' Create grid of tau values
#'
#' @param increment_width The distance between consecutive tau values
#' @inheritParams compute_familyWeights
#'
#' @return tau_grid a matrix of taus
#' @keywords internal
make_tauGrid <- function(increment_width = 0.05, constrained = TRUE){

  #create grid of tau values
  tau_grid <- expand.grid(seq(0.5, 1, by = increment_width), seq(0.5, 1, by = increment_width))

  #if we are considering the constrined space, reduce
  #taus to combinations for which tau_A > tau_B
  if (constrained) {
    tau_grid <- tau_grid[which(tau_grid[, 1] > tau_grid[, 2]), ]
    tau_grid[(nrow(tau_grid) + 1), ] <- c(0.5, 0.5)
  }

  tau_grid <- tau_grid[order(tau_grid[, 1], tau_grid[, 2], decreasing = FALSE), ]

  row.names(tau_grid) = NULL
  colnames(tau_grid) = c("tau_A", "tau_B")

  return(tau_grid)
}


#' Find the values of tau that maximize the likelihood function
#'
#' Find the values of tau that maximize the likelihood function. NOTE: if there is not a unique MCLE, this funciton orders the results conservatively.  That is, the smallest values of tau_A and tau_B are the first values returned.
#'
#' @inheritParams test_allcombos
#' @param like_mat The matrix of likelihood values given a sharing configuation and values of tau, i.e. the first item returned by \code{\link{test_allcombos}}.
#'
#' @return  A list of taus that maximize the likelihood for each configuration
#' @keywords internal
get_MLE <- function(like_mat, tau_grid){

  # max_like <- apply(like_mat, 2, function(x){
  #   x[which.max(x)]
  # })
  #
  # max_taus <- lapply(1:ncol(like_mat), function(x){
  #   tau_grid[like_mat[, x] == like_mat[max_like[[x]], x], ]
  # })

  max_taus <- apply(like_mat, 2, function(x){
    tau_grid[x == max(x), ]
  })

  return(max_taus)
}
