#' @title Compare networks across groups using Bayesian inference
#'
#' @description Easy comparison of networks using Bayesian inference to extract
#'   differences in conditional (in)dependence relations across groups.
#'
#' @name easybgm_compare
#'
#' @param data The data can be provided in two formats:
#'
#'   \strong{1. A list of two dataframes} (two-group comparison): Each list element
#'   is an n x p matrix or dataframe for one group. The variables (columns) must
#'   be the same across both dataframes. This format is supported by both
#'   \code{bgms} and \code{BGGM}.
#'
#'   \strong{2. A single matrix or dataframe} (multi-group comparison): An n x p
#'   matrix containing responses from all groups combined. Requires the
#'   \code{group_indicator} argument to specify which rows belong to which
#'   group. This format supports two or more groups and is only available with
#'   the \code{bgms} package.
#'
#' @param type Specifies the data type. Can be used in two ways:
#'
#'   \strong{1. A single string}:
#'   \itemize{
#'     \item \code{"continuous"}: For continuous data. Default package: BGGM.
#'     \item \code{"ordinal"}: For ordinal (Likert-type) data. Default package:
#'       bgms.
#'     \item \code{"binary"}: For binary (0/1) data. Default package: bgms.
#'     \item \code{"blume-capel"}: For Blume-Capel ordinal data. Requires
#'       \code{baseline_category}. Default package: bgms.
#'     \item \code{"mixed"}: For mixed data, requires \code{not_cont}. Default 
#'        package: BGGM.
#'   }
#'
#'   \strong{2. A character vector of length p} (per-variable specification):
#'   Each element specifies the type of the corresponding column. Valid values
#'   are \code{"ordinal"}, \code{"blume-capel"}, and \code{"binary"}; note that
#'   \code{"continuous"} is \emph{not} supported for group comparison with
#'   \code{bgms}. For example: \code{type = c("ordinal", "ordinal", "blume-capel")}.
#'
#'   Per-variable vectors are fitted by \code{bgms} and require \code{bgms}
#'   version 0.2.0.0 or later.
#'
#' @param package The R-package used for fitting the comparison model, either
#'   \code{"bgms"} or \code{"BGGM"}. If not specified, \code{bgms} is used for
#'   ordinal, binary, and blume-capel data, and \code{BGGM} for continuous and
#'   mixed data.
#'
#' @param not_cont A binary vector of length p, required when
#'   \code{type = "mixed"}. Each element indicates whether the corresponding
#'   variable is not continuous (\code{1} = not continuous/ordinal,
#'   \code{0} = continuous).
#'
#' @param group_indicator An integer vector of length n specifying group
#'   membership for each row in \code{data}. Required when \code{data} is a
#'   single matrix/dataframe (multi-group comparison). Supports two or more
#'   groups (e.g., \code{rep(c(1, 2, 3), each = 50)} for three groups of 50
#'   observations each). Only available with the \code{bgms} package.
#'
#' @param iter Number of iterations for the sampler. The default is 1e4.
#'   The recommended number of iterations depends on the data and model
#'   complexity. Check convergence diagnostics in the output.
#'
#' @param save Logical. Should the posterior samples be obtained
#'   (default = \code{TRUE})? If \code{TRUE}, the output includes a
#'   \code{samples_posterior} matrix with posterior samples for each difference
#'   parameter.
#'
#' @param progress Logical. Should a progress bar be shown
#'   (default = \code{TRUE})?
#'
#' @param ... Additional arguments passed to the fitting functions of the
#'   underlying packages (e.g., prior specifications). Consult the documentation
#'   of \code{bgms} and \code{BGGM} for the specific options available.
#'
#' @return An object of class \code{easybgm_compare} with the following
#'   elements:
#'
#' \strong{Always returned:}
#' \itemize{
#'   \item \code{parameters}: A p x p matrix of posterior mean differences in
#'     partial associations across groups.
#'   \item \code{inc_probs}: A p x p matrix of posterior inclusion probabilities
#'     for group differences (i.e., the probability that an edge differs between
#'     groups).
#'   \item \code{inc_BF}: A p x p matrix of inclusion Bayes factors for group
#'     differences.
#'   \item \code{structure}: A p x p adjacency matrix of the median probability
#'     model for differences (edges with posterior inclusion probability > 0.5).
#'   \item \code{model}: A string indicating the data type used.
#' }
#'
#' \strong{Returned for bgms only:}
#' \itemize{
#'   \item \code{structure_probabilities}: Posterior probabilities of all
#'     visited graph structures.
#'   \item \code{graph_weights}: Number of times each structure was visited.
#'   \item \code{sample_graph}: Identifiers for each visited structure.
#'   \item \code{convergence_parameter}: The Gelman-Rubin (R-hat) convergence
#'     statistic for each difference parameter. Values close to 1 indicate
#'     good convergence.
#' }
#'
#' \strong{Returned when save = TRUE:}
#' \itemize{
#'   \item \code{samples_posterior}: A k x iter matrix of posterior samples for
#'     each difference parameter (k = p*(p-1)/2 edges).
#' }
#'
#' @details
#'
#' \strong{Data types and package support for group comparison}
#'
#' \tabular{lcc}{
#'   \strong{Data type}  \tab \strong{bgms} \tab \strong{BGGM} \cr
#'   ordinal              \tab Yes (default)  \tab No            \cr
#'   binary               \tab Yes (default)  \tab Yes            \cr
#'   blume-capel          \tab Yes (default)  \tab No            \cr
#'   continuous            \tab No             \tab Yes (default) \cr
#'   mixed                 \tab No             \tab Yes (default) \cr
#' }
#'
#'
#'
#' \strong{Prior specification}
#'
#' Users may wish to adjust priors via the \code{...} argument. We summarize
#' the bgms options here and refer to \code{\link[bgms]{bgmCompare}} and
#' \code{\link[BGGM]{explore}} for full details.
#'
#' \emph{bgms} (>= 0.2.0.0)
#' \itemize{
#'   \item \code{interaction_prior}: Prior on the baseline pairwise
#'     interactions. Use \code{\link[bgms]{normal_prior}(scale)} (default
#'     \code{normal_prior(scale = 1)}), \code{\link[bgms]{cauchy_prior}(scale)},
#'     or \code{\link[bgms]{beta_prime_prior}(alpha, beta)}.
#'     For example, a cauchy prior with scale 1 would be specified with adding the 
#'     argument \code{threshold_prior = cauchy_prior(1)} to the easybgm call.
#'   \item \code{threshold_prior}: Prior on threshold parameters. Default
#'     \code{beta_prime_prior(0.5, 0.5)}, for example, specified by adding 
#'     \code{threshold_prior = beta_prime_prior(0.5, 0.5)} to the easybgm call.
#'   \item \code{difference_prior}: Indicator prior on group differences.
#'     Use \code{bernoulli_prior(inclusion_probability)} (default
#'     \code{bernoulli_prior(0.5)}) or
#'     \code{beta_bernoulli_prior(alpha, beta)}. For example, a beta bernoulli 
#'     prior with alpha 1 and beta 3 can be specified by adding 
#'     \code{difference_prior = beta_bernoulli_prior(1, 3)} to the easybgm call.
#'   \item \code{difference_family}: The family of the prior on the magnitude of
#'     the pairwise differences, either \code{"Normal"} (the default) or
#'     \code{"Cauchy"}. Versions of \code{bgms} before 0.2.0.0 had no such
#'     argument and always used a Cauchy, so comparison Bayes factors obtained
#'     with those versions correspond to \code{difference_family = "Cauchy"}.
#'   \item \code{difference_scale}: Scale of the prior on the magnitude of
#'     pairwise differences, on the family set by \code{difference_family}.
#'     Default 1, for example, specified by adding the argument 
#'     \code{difference_scale = 1} to the easybgm call. 
#'   \item \code{difference_selection}: Logical, whether to perform Bayesian
#'     selection on group differences. Default \code{TRUE}.
#' }
#'
#' For backwards compatibility of \emph{bgms} (< 0.2.0.0), the previous
#' prior specifications are still accepted and translated into the relevant constructs. 
#' Check the previous bgms version for its prior arguments. 
#'
#' We always encourage researchers to conduct prior sensitivity checks.
#'
#' @export
#' @import bgms
#' @importFrom BGGM explore select
#' @importFrom utils packageVersion
#'
#' @examples
#'
#' \dontrun{
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(ADHD)
#'
#' # --- Two-group comparison (list input) ---
#' group1 <- data[1:10, 1:3]
#' group2 <- data[11:20, 1:3]
#'
#' fit <- easybgm_compare(list(group1, group2),
#'                 type = "binary", save = TRUE,
#'                 iter = 100   # for demonstration only
#'                 )
#' summary(fit)
#'
#' # --- Multi-group comparison (single dataframe + group_indicator) ---
#' fit_multi <- easybgm_compare(data[1:80, 1:5],
#'                 group_indicator = rep(c(1, 2, 3, 4), each = 20),
#'                 type = "binary", save = TRUE,
#'                 iter = 100   # for demonstration only
#'                 )
#' summary(fit_multi)
#' }


easybgm_compare <- function(data, 
                            type, 
                            package = NULL, 
                            not_cont = NULL, 
                            group_indicator = NULL, 
                            iter = 1e3,
                            save = TRUE, 
                            progress = TRUE,
                            ...){
  
  if(!is.list(data) && is.null(group_indicator)){
    stop("Your data can't be read. There are two options of providing your data: 
         1) Provide two datasets in a list containing only the two datasets, or 
         for ordinal data with the bgms package > 0.1.6. 2) provide the data as 
         a matrix or data.frame together with specifying the 'group_indicator' 
         argument, which then also allows for multi-group comparison.",
         call. = FALSE)
  }
  
  # Per-variable 'type' vectors require bgms >= 0.2.0.0. Resolve the fitting
  # package first, then validate 'type' against it.
  bgms_supports_vector_type <- packageVersion("bgms") >= "0.2.0.0"
  is_vector_type <- length(type) > 1

  if(is_vector_type && !bgms_supports_vector_type){
    stop("Specifying 'type' as a vector with one entry per variable requires ",
         "bgms version 0.2.0.0 or later (installed: ", packageVersion("bgms"), "). ",
         "Please either update bgms, or pass a single type for all variables.",
         call. = FALSE)
  }

  # --- Resolve the fitting package ---------------------------------------
  if(is.null(package)){
    package <- if(is_vector_type){
      "package_bgms_compare"
    } else if(type %in% c("continuous", "mixed")){
      "package_bggm_compare"
    } else {
      "package_bgms_compare"
    }
  } else {
    package <- switch(package,
                      "BGGM" = "package_bggm_compare",
                      "bgms" = "package_bgms_compare",
                      stop("Unrecognized 'package': '", package, "'. ",
                           "Valid options are 'bgms' and 'BGGM'.",
                           call. = FALSE))

    # BGGM only fits 'continuous' and 'mixed' data. Anything else can only be
    # fitted by bgms, so it overrides an explicit package choice.
    if(package == "package_bggm_compare"){
      override_reason <- if(is_vector_type){
        "a per-variable 'type' vector"
      } else if(!type %in% c("continuous", "mixed")){
        paste0("type = '", type, "'")
      } else NULL

      if(!is.null(override_reason)){
        warning("BGGM can only fit 'continuous' or 'mixed' data, so it cannot ",
                "fit ", override_reason, "; the 'package' argument was ",
                "overridden and bgms will be used instead.",
                call. = FALSE)
        package <- "package_bgms_compare"
      }
    }
    # bgms can not compare continuous or mixed data. change to BGGM instead.
    if(package == "package_bgms_compare"){
      override_reason <- if(is_vector_type){
        "a per-variable 'type' vector"
      } else if(!type %in% c("binary", "ordinal", "blume-capel")){
        paste0("type = '", type, "'")
      } else NULL
      
      if(!is.null(override_reason)){
        warning("bgms can only fit 'binary', 'ordinal' or 'blume-capel' data, 
                so it cannot ", "fit ", override_reason, "; 
                the 'package' argument was ",
                "overridden and BGGM will be used instead.",
                call. = FALSE)
        package <- "package_bggm_compare"
      }
    }
  }

  # --- Validate 'type' against the resolved package ------------------------
  if(package == "package_bgms_compare"){
    valid_types <- c("ordinal", "blume-capel", "binary")
    invalid <- type[!type %in% valid_types]
    if(length(invalid) > 0) {
      warning("The following variable type(s) are not recognized for group comparison: ",
              paste0("'", unique(invalid), "'", collapse = ", "), ". ",
              "Valid types are: ", paste(valid_types, collapse = ", "), ". ",
              "Note: 'continuous' is not supported for group comparison via bgms.",
              call. = FALSE)
      stop("Invalid variable types detected. See the warning message for more details.",
           call. = FALSE)
    }
    ncols <- if(is.list(data) && !is.data.frame(data)) ncol(data[[1]]) else ncol(data)
    if(is_vector_type && length(type) != ncols) {
      stop("When 'type' is a vector, its length (", length(type), ") must equal ",
           "the number of columns in 'data' (", ncols, ").",
           call. = FALSE)
    }
    type[type == "binary"] <- "ordinal"
  } else {
    # BGGM takes a single 'continuous' or 'mixed' type for all variables.
    if(is_vector_type || !type %in% c("continuous", "mixed")){
      stop("BGGM can only fit 'continuous' or 'mixed' data. ",
           "Please use package = 'bgms' for other data types.",
           call. = FALSE)
    }
    if(type == "mixed" && is.null(not_cont)){
      stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
           call. = FALSE)
    }
  }


  dots <- list(...)
  has_reference <- "reference_category" %in% names(dots)
  has_baseline  <- "baseline_category" %in% names(dots)

  if (any(type == "blume-capel") && !(has_reference || has_baseline)) {
    stop("For the Blume-Capel model, a reference category needs to be specified.
         Use the baseline_category argument to specify the reference category.",
         call. = FALSE)
  }


  if((is.data.frame(data) || is.matrix(data)) && package == "package_bggm_compare"){
    stop("Your data can't be read. For continuous data fit with BGGM, 
         you can only provide two datasets in a list.",
         call. = FALSE)
  }
  
  fit <- list()
  class(fit) <- c(package, "easybgm")
  
  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, not_cont = not_cont, 
                    group_indicator = group_indicator, 
                    iter = iter,
                    save = save, progress = progress, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop(paste("Error meassage: ", e$message, "Please consult the original message for more information.") )
    })
  
  # Extract the results
  res <- bgm_extract(fit, type = type,
                     save = save, not_cont = not_cont,
                     group_indicator = group_indicator, 
                     data = data, 
                     ...)
  
  # Output results
  class(res) <- c(package, "easybgm_compare", "easybgm")
  return(res)
}
