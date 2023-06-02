#' Fit a Bayesian analysis of networks
#'
#' @param x Object with a particular class that will dispatch to the respective package functions
#' @param ... Additional arguments to be passed onto the respective fitting functions
#' 
set_prior_defaults <- function(x, ...) {
  UseMethod("set_prior_defaults", x)
}