#' @title Fit a Bayesian analysis of networks
#'
#' @description Easy estimation of a Bayesian analysis of networks to obtain conditional (in)dependence relations between variables in a network.
#'
#' @name easybgm
#'
#' @param data An n x p matrix or dataframe containing the variables for n independent observations on p variables.
#' @param type What is the data type? Options: continuous, mixed, ordinal, binary
#' @param package The R-package that should be used for fitting the network model; supports BGGM, BDgraph, and bgms. Optional argument;
#'     default values are specified depending on the datatype.
#' @param not.cont If data-type is mixed, a vector of length p, specifying the not-continuous
#'     variables (1 = not continuous, 0 = continuous).
#' @param iter number of iterations for the sampler.
#' @param save Logical. Should the posterior samples be obtained (default = FALSE)?
#' @param centrality Logical. Should the centrality measures be extracted (default = FALSE)? Note, that it will significantly increase the computation time.
#' @param progress Logical. Should a progress bar be shown (default = TRUE)?
#' @param edge.prior Default is 0.5. Single value or p x p matrix with one value per edge.
#' @param ... Additional arguments that are handed to the fitting functions of the packages, e.g., informed prior specifications.
#'
#'
#' @return The returned object of \code{easybgm} contains a lot of information that
#'         is used for visualizing and interpreting the results of a Bayesian analysis of network.
#'         The output includes:
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
#' @return In addition, for `BDgraph` and `bgms`, the function returns:
#'
#' \itemize{
#'
#' \item \code{structure_probabilities} A vector containing the posterior probabilities of all visited structures, between 0 and 1.
#'
#' \item \code{graph_weights} A vector containing the number of times a particular structure was visited.
#'
#' \item \code{sample_graphs} A vector containing the indexes of a particular structure.
#' }
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
#' @export
#'
#' @importFrom bgms bgm
#' @importFrom BDgraph bdgraph bdgraph.mpl plinks
#' @importFrom BGGM explore select
#'
#' @examples
#' \donttest{
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- Wenchuan
#'
#' # Fitting the Wenchuan PTSD data
#'
#' fit <- easybgm(data, type = "ordinal",
#'                 iter = 1000 # note for demonstrative purposes, normally 1e5 is recommended
#'                 )
#'
#' # To extract the posterior parameter distribution
#' # and centrality measures
#'
#' fit <- easybgm(data, type = "ordinal",
#'                 iter = 1000, save = TRUE,
#'                 centrality = TRUE)
#' }



easybgm <- function(data, type, package = NULL, not.cont = NULL, iter = 1e4,
                    save = FALSE, centrality = FALSE, progress = TRUE,
                    edge.prior = 0.5, ...){


  if(type == "mixed" & is.null(not.cont)){
    stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
         call. = FALSE)
  }

  # Set default values for fitting if package is unspecified
  if(is.null(package)){
    if(type == "continuous") package <- "package_bdgraph"
    if(type == "mixed") package <- "package_bdgraph"
    if(type == "ordinal") package <- "package_bgms"
    if(type == "binary") package <- "package_bgms"
  } else {
    if(package == "BDgraph") package <- "package_bdgraph"
    if(package == "BGGM") package <- "package_bggm"
    if(package == "bgms") package <- "package_bgms"
  }

  if(type =="continuous" & package == "package_bdgraph" & any(is.na(data))){
    warning("The data contains missing values which cannot be handled as continuous data by BDgraph.
            Note that we switched the type to \"mixed\", which estimates a GCGM and can impute missing data.")
    type <- "mixed"
    not.cont <- rep(0, ncol(data))
  }

  if((package == "package_bgms") & (type %in% c("continuous", "mixed"))){
    warning("bgms can only fit ordinal or binary datatypes. For continuous or mixed data,
           choose either the BDgraph or BGGM package. By default we have changed the package to BDgraph",
            call. = FALSE)
    package <- "package_bdgraph"

  }


  fit <- list()
  class(fit) <- c(package, "easybgm")

  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, not.cont = not.cont, iter = iter,
                    save = save, centrality = centrality, progress = progress, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop(paste("Error meassage: ", e$message, "Please consult the original message for more information.") )
    })

  # Extract the results
  res <- bgm_extract(fit, model = fit$model,
                     edge.prior = edge.prior,
                     save = save, not.cont = not.cont,
                     data = data, centrality = centrality, ...)

  # Output results
  class(res) <- c(package, "easybgm")
  return(res)
}
