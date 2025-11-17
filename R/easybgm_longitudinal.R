#' @title Fit a Bayesian analysis of longitudinal networks
#'
#' @description Easy estimation of a Bayesian analysis of networks to obtain conditional (in)dependence relations between variables in a network.
#'
#' @name easybgm_longitudinal
#'
#' @param data A data frame or matrix containing the time series data of a
#'   single subject. The data should be in long format,  where the columns
#'   represent the variables and the rows represent the time points. See the
#'   example data \code{\link{ts_data}} for the correct format.
#' @param package The R-package that should be used for fitting the network model; supports currently only tsnet. 
#' @param type The type of the data to be fit.
#' @param iter number of iterations for the sampler.
#' @param centrality Logical. Should the centrality measures be extracted (default = FALSE)? 
#' @param progress Logical. Should a progress bar be shown (default = TRUE)?
#' @param ... Additional arguments that are handed to the fitting functions of the packages, e.g., informed prior specifications.
#'
#'
#' @return The returned object of \code{easybgm_longitudinal} contains several elements:
#'
#' \itemize{
#'
#' \item \code{parameters} A p x p matrix containing partial associations.
#'
#' \item \code{inc_probs} A p x p matrix containing the posterior inclusion probabilities.
#'
#' \item \code{BF} A p x p matrix containing the posterior inclusion Bayes factors.
#'
#'  \item \code{structure} Adjacency matrix of the median probability model (i.e., edges with a posterior probability larger 0.5).
#' }
#'
#'
#'
#' @return For all packages, when setting `save = TRUE` and `centrality = TRUE`, the function will return the following objects respectively:
#'
#' \itemize{
#'
#' \item \code{samples_posterior} A k x iter matrix containing the posterior samples for each parameter (i.e., k = (p/(p-1))/2) at each iteration (i.e., iter) of the sampler.
#'
#' \item \code{centrality} A p x iter matrix containing the centrality of a node at each iteration of the sampler.
#' }
#'
#' @details
#'
#' Users may oftentimes wish to deviate from the default, usually uninformative, prior specifications of the
#' packages to informed priors. This can be done by simply adding additional arguments to the \code{easybgm} function.
#' Depending on the package that is running the underlying network estimation, researcher can specify different prior
#' arguments. We give an overview of the prior arguments per package below.
#'
#'
#' \strong{tsnet}:
#'
#' \itemize{
#'
#' \item \code{xx} xx
#'
#' \item \code{xx} xx
#'
#' }
#'
#' We would always encourage researcher to conduct prior robustness checks.
#'
#' @export
#'
#' @importFrom tsnet stan_gvar stan_fit_convert get_centrality
#'
#' 
#'




easybgm_longitudinal <- function(data, 
                                 type = "continuous", 
                                 package = NULL, 
                                 iter = 1e4,
                                 centrality = FALSE, 
                                 progress = TRUE, 
                                 ...){
  
  
  # Set default values for fitting if package is unspecified
  if(is.null(package)){
    package <- "package_tsnet"
  } else {
    if(package == "tsnet") package <- "package_tsnet"
  }
  
  
  fit <- list()
  class(fit) <- c(package, "easybgm_longitudinal", "easybgm")
  
  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, iter = iter,
                    progress = progress, centrality = centrality, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop(paste("Error meassage: ", e$message, "Please consult the original message for more information.") )
    })
  
  # Extract the results
  res <- bgm_extract(fit, type = type,
                     data = data, centrality = centrality,
                     ...)
  
  # Output results
  class(res) <- c(package, "easybgm_longitudinal", "easybgm")
  return(res)
}
