#' Fit a Bayesian analysis of networks
#'
#' @param fit Object with a particular class that will dispatch to the respective package functions
#' @param ... Additional arguments to be passed onto the respective fitting functions


bgm_fit <- function(fit, ...){

  UseMethod("bgm_fit", fit)

}
