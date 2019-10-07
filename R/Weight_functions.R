#' Compute the transmission based weight
#'
#' @param all_stats The data frame of stats for all configurations
#' @param config_index The index of the observed confiugration
#'
#' @return  the transmission based weight
#' @keywords internal
compute_transmission_weight <- function(all_stats, config_index){
  #find a better fix. Implementing this right now
  #to fix floating point errors
  all_stats$LR <- round(all_stats$LR, digits = 8)

  1/sum(all_stats$RVsharing[all_stats$LR >= all_stats$LR[config_index]
                            & all_stats$K >= all_stats$K[config_index]],
        na.rm = TRUE)
}


#' Compute the modified RVS weight
#'
#' @inheritParams compute_transmission_weight
#'
#' @return the modified RVS weight
#' @keywords internal
compute_modified_RVS_weight <- function(all_stats, config_index){
  1/sum(all_stats$RVsharing[all_stats$RVsharing <= all_stats$RVsharing[config_index]
                            & all_stats$K >= all_stats$K[config_index]
                            & all_stats$K_sub >= all_stats$K_sub[config_index]],
        na.rm = TRUE)
}

#' Compute the RVS weight
#'
#' @inheritParams compute_transmission_weight
#'
#' @return the RVS weight
#' @keywords internal
compute_RVS_weight <- function(all_stats, config_index){
  1/sum(all_stats$RVsharing[all_stats$RVsharing <= all_stats$RVsharing[config_index]
                            & all_stats$K >= all_stats$K[config_index]],
        na.rm = TRUE)
}
