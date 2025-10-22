#' @title Compare networks across groups using Bayesian inference
#'
#' @description Easy comparison of networks using Bayesian inference to extract differences across conditional (in)dependence across groups.
#' 
#' @name easybgm_compare
#'
#' @param data A list with two n x p matrices or dataframes containing the variables for n independent observations on p variables for two groups. Note that the variables need to be the same in the two different dataframes.
#' @param type What is the data type? Options: continuous, mixed, ordinal, binary, or blume-capel.
#' @param package The R-package that should be used for fitting the network model; supports BGGM and bgms. Optional argument;
#'     default values are specified depending on the datatype.
#' @param not_cont If data-type is mixed, a vector of length p, specifying the not-continuous
#'     variables (1 = not continuous, 0 = continuous).
#' @param iter number of iterations for the sampler.
#' @param save Logical. Should the posterior samples be obtained (default = TRUE)?
#' @param progress Logical. Should a progress bar be shown (default = TRUE)?
#' @param ... Additional arguments that are handed to the fitting functions of the packages, e.g., informed prior specifications.
#'
#'
#' @return The returned object of \code{easybgm} contains several elements:
#'
#' \itemize{
#'
#' \item \code{parameters} A p x p matrix containing difference across partial associations.
#'
#' \item \code{inc_probs} A p x p matrix containing the posterior inclusion probabilities of subgroup differences.
#'
#' \item \code{BF} A p x p matrix containing the posterior inclusion Bayes factors of subgroup differences.
#'
#'  \item \code{structure} Adjacency matrix of the median probability model (i.e., edges with a posterior probability larger 0.5).
#' }
#'
#'
#' @return In addition, for `bgms`, the function returns:
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
#' @return For both packages, when setting `save = TRUE`, the function will also return the following object:
#'
#' \itemize{
#'
#' \item \code{samples_posterior} A k x iter matrix containing the posterior samples of parameter differences (i.e., k = (p/(p-1))/2) at each iteration (i.e., iter) of the sampler.

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
#' \item \code{pairwise_difference_prior} prior on the graph structure, which can be either "Bernoulli" or "Beta-Bernoulli". The default is "Bernoulli".
#'
#' \item \code{pairwise_difference_probability} prior edge inclusion probability for the "Bernoulli" distribution. an be a single probability or a matrix of p rows and p columns specifying the probability of a difference for each edge pair. The default is 0.5.
#'
#' \item \code{pairwise_beta_bernoulli_alpha} and \code{pairwise_beta_bernoulli_beta} the parameters of the "Beta-Bernoulli" priors. The default is 1 for both.
#' 
#' \item \code{threshold_alpha} and \code{threshold_beta} the parameters of the beta-prime distribution for the threshold parameters. The defaults are both set to 1.
#'
#' }
#'
#' \strong{BGGM}:
#'
#' \itemize{
#'
#' \item \code{prior_sd} the standard deviation of the prior distribution of the interaction parameters, approximately the scale of a beta distribution. The default is 0.25.

#' }
#'
#' We always encourage researcher to conduct prior robustness checks.
#'
#'
#' @export
#'
#' @importFrom bgms bgm
#' @importFrom BGGM explore select
#' @importFrom utils packageVersion
#'





easybgm_compare <- function(data, type, package = NULL, not_cont = NULL, iter = 1e4,
                    save = TRUE, progress = TRUE,
                    ...){

  if(class(data) != "list"){
    stop("Please provide your two datasets in a list containing only the two datasets.",
         call. = FALSE)
  }
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
    if(type == "continuous") package <- "package_bggm_compare"
    if(type == "mixed") package <- "package_bggm_compare"
    if(type == "ordinal") package <- "package_bgms_compare"
    if(type == "binary") package <- "package_bgms_compare"
    if(type == "blume-capel") package <- "package_bgms_compare"
  } else {
    if(package == "BGGM") package <- "package_bggm_compare"
    if(package == "bgms") package <- "package_bgms_compare"
    if(type == "binary") package <- "package_bgms_compare"
  }


  if((package == "package_bgms_compare") & (type %in% c("continuous", "mixed"))){
    warning("bgms can only fit ordinal or binary datatypes. For continuous or mixed data,
           choose the BGGM package. By default we have changed the package to BGGM",
            call. = FALSE)
    package <- "package_bggm_compare"

  }


  fit <- list()
  class(fit) <- c(package, "easybgm")

  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, not_cont = not_cont, iter = iter,
                    save = save, progress = progress, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop(paste("Error meassage: ", e$message, "Please consult the original message for more information.") )
    })

  # Extract the results
  res <- bgm_extract(fit, type = type,
                     save = save, not_cont = not_cont,
                     data = data, 
                     ...)

  # Output results
  class(res) <- c(package, "easybgm")
  return(res)
}
