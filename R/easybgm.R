#' @title Fit a Bayesian analysis of networks
#'
#' @description Easy estimation of a Bayesian analysis of networks to obtain conditional (in)dependence relations between variables in a network.
#'
#' @name easybgm
#'
#' @param data An n x p matrix or dataframe containing the variables for n independent observations on p variables.
#' @param type What is the data type? Options: continuous, mixed, ordinal, binary, or blume-capel.
#' @param package The R-package that should be used for fitting the network model; supports BGGM, BDgraph, and bgms. Optional argument;
#'     default values are specified depending on the datatype.
#' @param not_cont If data-type is mixed, a vector of length p, specifying the not-continuous
#'     variables (1 = not continuous, 0 = continuous).
#' @param iter number of iterations for the sampler.
#' @param save Logical. Should the posterior samples be obtained (default = FALSE)?
#' @param centrality Logical. Should the centrality measures be extracted (default = FALSE)? Note, that it will significantly increase the computation time.
#' @param progress Logical. Should a progress bar be shown (default = TRUE)?
#' @param posterior_method Determines how the posterior samples of the edge weight parameters are obtained for models fit with BDgraph. The argument can be either MAP for the maximum-a-posteriori or model-averaged. If MAP, samples are obtained for the edge weights only for the most likely structure. If model-averaged, samples are obtained for all plausible structures weighted by their posterior probability. Default is model-averaged.
#' @param ... Additional arguments that are handed to the fitting functions of the packages, e.g., informed prior specifications.
#'
#' @return The returned object of \code{easybgm} contains several elements:
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
#'
#' For the \code{bgms} package, when \code{edge_prior = "Stochastic-Block"},
#' the function will also return an object \code{sbm} which contains:
#'   \itemize{
#'
#'    \item \code{posterior_num_blocks} A data frame with the estimated posterior
#'    probability of the possible number of clusters.
#'
#'    \item \code{posterior_mean_allocations} The posterior mean of the cluster assignments of the nodes.
#'
#'    \item \code{posterior_mode_allocations} The posterior mode of the cluster assignments of the nodes.
#'
#'    \item \code{posterior_mean_coclustering_matrix} A p x p matrix containing the estimated
#'       pairwise proportions of cluster occurrence of every variable. This matrix
#'       can be plotted to visually inspect the estimated number of clusters
#'       and visually inspect nodes that tend to switch clusters.
#'
#'   }
#' }
#'
#' If using version 0.1.6.1 or higher of the \code{bgms} package, the function also returns the
#' the Gelman-Rubin convergence statistic for each edge weight parameter. As well as the
#'  95% Monte Carlo confidence interval for the inclusion Bayes factor.
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
#' \strong{bgms}:
#'
#' \itemize{
#'
#' \item \code{interaction_scale} the scale of the Cauchy distribution that is used as a
#' prior for the pairwise interaction parameters. The default is 2.5.
#'
#' \item \code{edge_prior} prior on the graph structure, which can be either "Bernoulli", "Beta-Bernoulli" or "Stochastic Block". The default is "Bernoulli".
#'
#' \item \code{inclusion_probability} prior edge inclusion probability for the "Bernoulli" distribution. The default is 0.5.
#'
#' \item \code{beta_bernoulli_alpha} and \code{beta_bernoulli_beta} the parameters of the "Beta-Bernoulli" or "Stochastic Block" priors. The default is 1 for both.
#'
#' \item \code{beta_bernoulli_alpha_between} and \code{beta_bernoulli_beta_between} the parameters of the "Stochastic Block" prior for edges between blocks.
#' This is currently only available in a developer version of bgms and will be available in version 0.1.6.2 or higher.
#'
#' \item \code{dirichlet_alpha} The shape of the Dirichlet prior on the node-to-block
#' allocation parameters for the Stochastic Block prior on the graph structure.
#'
#' \item \code{threshold_alpha} and \code{threshold_beta} the parameters of the beta-prime distribution for the threshold parameters. The defaults are both set to 1.
#'
#' \item \code{variable_type} What kind of variables are there in \code{x}? Can be a
#' single character string specifying the variable type of all \code{p}
#' variables at once or a vector of character strings of length \code{p}
#' specifying the type for each variable in \code{x} separately. Currently, bgm
#' supports ``ordinal'' and ``blume-capel''. Binary variables are automatically
#' treated as ``ordinal’’. Defaults to \code{variable_type = "ordinal"}.
#'
#' }
#'
#' \strong{BDgraph}:
#'
#' \itemize{
#'
#' \item \code{df.prior} prior on the parameters (i.e., inverse covariance matrix), degrees of freedom of the prior G-Wishart distribution. The default is set to 3.
#'
#' \item \code{g.prior} prior probability of edge inclusion. This can be either a scalar, if it is the same for all edges, or a matrix, if it should be different among the edges. The default is set to 0.5.
#'
#' }
#' \strong{BGGM}:
#'
#' \itemize{
#'
#' \item \code{prior_sd} the standard deviation of the prior distribution of the interaction parameters, approximately the scale of a beta distribution. The default is 0.25.

#' }
#'
#' We would always encourage researcher to conduct prior robustness checks.
#'
#' @export
#'
#' @import bgms
#' @importFrom BDgraph bdgraph bdgraph.mpl plinks
#' @importFrom BGGM explore select
#' @importFrom utils packageVersion
#'
#' @examples
#'
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#'
#' # Fitting the Wenchuan PTSD data
#'
#' fit <- easybgm(data, type = "continuous",
#'                 iter = 100 # for demonstration only (> 5e4 recommended)
#'                 )
#'
#' summary(fit)
#'
#' \donttest{
#' # To extract the posterior parameter distribution
#' # and centrality measures
#'
#' fit <- easybgm(data, type = "continuous",
#'                 iter = 100, 
#'                 centrality = TRUE, save = TRUE)
#' }



easybgm <- function(data, type, package = NULL, not_cont = NULL, iter = 1e4, save = FALSE,
                    centrality = FALSE, progress = TRUE, posterior_method = "model-averaged",
                    ...){


  if(type == "mixed" & is.null(not_cont)){
    stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
         call. = FALSE)
  }

  dots <- list(...)
  has_reference <- "reference_category" %in% names(dots)
  has_baseline  <- "baseline_category" %in% names(dots)

  # Example: If type == "blume-capel", then at least one must be present
  if (type == "blume-capel" && !(has_reference || has_baseline)) {
    stop("For the Blume-Capel model, a reference category needs to be specified.
         If type is 'blume-capel' it specifies the reference category in the Blume-Capel model.
         Should be an integer within the range of integer scores observed for the
         'blume-capel' variable. Can be a single number specifying the reference
         category for all Blume-Capel variables at once, or a vector of length
         p where the i-th element contains the reference category for
         variable i if it is Blume-Capel, and bgm ignores its elements for
         other variable types. The value of the reference category is also recoded
         when bgm recodes the corresponding observations. Only required if there is at
         least one variable of type ``blume-capel''.
         For bgms version smaller than 0.1.6, use the reference_category argument.
         For all package versions including and older than 0.1.6., the baseline_category argument.",
         call. = FALSE)
  }


  # Set default values for fitting if package is unspecified
  if(is.null(package)){
    if(type == "continuous") package <- "package_bggm"
    if(type == "mixed") package <- "package_bggm"
    if(type == "ordinal") package <- "package_bgms"
    if(type == "binary") package <- "package_bgms"
    if(type == "blume-capel") package <- "package_bgms"
  } else {
    if(package == "BDgraph") package <- "package_bdgraph"
    if(package == "BGGM") package <- "package_bggm"
    if(package == "bgms") package <- "package_bgms"
    if(type == "binary") package <- "package_bgms"
  }

  if(type =="continuous" & package == "package_bdgraph" & any(is.na(data))){
    warning("The data contains missing values which cannot be handled as continuous data by BDgraph.
            Note that we switched the type to \"mixed\", which estimates a GCGM and can impute missing data.")
    type <- "mixed"
    not_cont <- rep(0, ncol(data))
  }

  if((package == "package_bgms") & (type %in% c("continuous", "mixed"))){
    warning("bgms can only fit ordinal, binary, or blume-capel datatypes. For continuous or mixed data,
           choose either the BDgraph or BGGM package. By default we have changed the package to BDgraph",
            call. = FALSE)
    package <- "package_bdgraph"

  }


  fit <- list()
  class(fit) <- c(package, "easybgm")

  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, not_cont = not_cont, iter = iter,
                    save = save, centrality = centrality, progress = progress, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop(paste("Error meassage: ", e$message, "Please consult the original message for more information.") )
    })

  # Extract the results
  res <- bgm_extract(fit, type = type,
                     save = save, not_cont = not_cont,
                     data = data, centrality = centrality,
                     posterior_method = posterior_method,
                     iter = iter,
                     ...)

  # Output results
  class(res) <- c(package, "easybgm")
  return(res)
}
