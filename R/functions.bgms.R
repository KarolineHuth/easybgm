# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not_cont, centrality, progress, ...){

  if(!save && centrality){
    save <- TRUE
  }


  bgms_fit <- do.call(
    bgm, c(list(x = data, iter = iter, save = save,
                display_progress = progress,
                ...))
  )


  fit$model <- type
  fit$packagefit <- bgms_fit
  if(is.null(colnames(data))){
    fit$var_names <- paste0("V", 1:ncol(data))
  } else {
    fit$var_names <- colnames(data)
  }
  class(fit) <- c("package_bgms", "easybgm")
  return(fit)
}

# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_extract.package_bgms <- function(fit, type, save,
                                     not_cont, data, centrality, ...){
  if(any(class(fit) != "bgms")){
    varnames <- fit$var_names
    fit <- fit$packagefit
    class(fit) <- "bgms"
  } else if (any(class(fit) == "bgms")){
    varnames <- fit$arguments$data_columnnames
    if(is.null(varnames)){
      varnames <- paste0("V", 1:fit$arguments$no_variables)
    }
  }

  args <- bgms::extract_arguments(fit)

  if (args$edge_prior[1] == "Bernoulli") {
    edge.prior <- args$inclusion_probability
  } else { # if BB or SBM
    edge.prior <- calculate_edge_prior(alpha = args$beta_bernoulli_alpha,
                                       beta = args$beta_bernoulli_beta)
    # otherwise it saves the wrong values (could be done more elegantly)
    args$inclusion_probability <- edge.prior
  }

  bgms_res <- list()

  if(args$save){
    p <- args$no_variables
    if(packageVersion("bgms") < "0.1.4"){
      pars <- extract_pairwise_interactions(fit)
    } else {
      pars <- extract_pairwise_interactions(fit)}
    bgms_res$parameters <- vector2matrix(colMeans(pars), p = p)
    if(packageVersion("bgms") < "0.1.4"){
      bgms_res$thresholds <- bgms::extract_pairwise_thresholds(fit)
    } else {
      bgms_res$thresholds <- extract_category_thresholds(fit)}
    colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters),
                                 nrow = nrow(bgms_res$parameters))

    if(args$edge_selection){
      bgms_res$inc_probs <- bgms::extract_posterior_inclusion_probabilities(fit)
      bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1 - edge.prior))
      bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)
      if(packageVersion("bgms") < "0.1.4"){
        gammas <- bgms::extract_edge_indicators(fit)
      } else {
        gammas <- extract_indicators(fit)}
      structures <- apply(gammas, 1, paste0, collapse="")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[,2]/nrow(gammas)
      bgms_res$graph_weights <- table_structures[,2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
      #EDGE SELECTION + SAVE
      if (args$edge_prior[1] == "Stochastic-Block") {
        extract_sbm <- summarySBM(fit)
        bgms_res$cluster_probabilities <- extract_sbm$components
        bgms_res$cluster_allocations <- extract_sbm$allocations
      }
    }
  } else {
    if(packageVersion("bgms") < "0.1.4"){
      bgms_res$parameters <- extract_pairwise_interactions(fit)
    } else {
      bgms_res$parameters <- extract_pairwise_interactions(fit)}
    if(packageVersion("bgms") < "0.1.4"){
      bgms_res$thresholds <- bgms::extract_pairwise_thresholds(fit)
    } else {
      bgms_res$thresholds <- extract_category_thresholds(fit)}
    colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters),
                                 nrow = nrow(bgms_res$parameters))
    if(args$edge_selection){
      bgms_res$inc_probs <- bgms::extract_posterior_inclusion_probabilities(fit)
      bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1-edge.prior))
      bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)
      if (args$edge_prior[1] == "Stochastic-Block") {
        bgms_res$cluster_probabilities <- fit$components
        bgms_res$cluster_allocations <- fit$allocations
      }
    }

  }
  if(args$save){
    if(packageVersion("bgms") < "0.1.4"){
      bgms_res$samples_posterior <- extract_pairwise_interactions(fit)
    } else {
      bgms_res$samples_posterior <- extract_pairwise_interactions(fit)}

    if(centrality){
      bgms_res$centrality <- centrality(bgms_res)
    }
  }

  if(args$edge_selection){
    # Adapt column names of output
    colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
    colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters)
  }
  bgms_res$model <- type
  bgms_res$fit_arguments <- args
  output <- bgms_res
  class(output) <- c("package_bgms", "easybgm")
  return(output)

}

# --------------------------------------------------------------------------------------------------
# 3. Function for calulating Clustering Bayes factors for Stochastic Block Model
# --------------------------------------------------------------------------------------------------
#' Calculate Clustering Bayes Factors for when using the Stochastic Block Model
#' as an edge prior
#'
#' This function calculates Bayes factors to evaluate evidence inf favor of
#' clustering for models fitted with the \code{bgms} packae with
#' \code{edge_prior = "Stochastic-Block"}. It supports two types of calculations:
#' simple Bayes factors between two hypothesized number of clusters (`b1` and `b2`),
#' and Bayes factor of the hypothesis of clustering against the null hypothesis of
#' no clustering.
#'
#' @param fit A fitted object of class \code{easybgm} or \code{bgms} containing
#' the clustering results.
#' @param type A character string specifying the type of Bayes factor to calculate.
#'   Options are `"simple"` or `"complement"`. Defaults to `"simple"`.
#' @param b1 Indicates the number of clusters according to the first simple hypothesis,
#'  required for `type = "simple"`.
#' @param b2 Indicates the number of clusters according to the first simple hypothesis,
#'  required for `type = "simple"`.

#' @return A numeric value representing the Bayes factor.
#'
#' @export
clusterBayesfactor <- function(fit,
                                type = "complement",
                                b1 = NULL,
                                b2 = NULL) {

  # check if the type argument is valid
  if (!type %in% c("simple", "complement")) {
    stop("The type argument must be either 'simple' or 'complement'.")
  }

  # Check the class of fit (if it is a bgms object, rename components)
  if (inherits(fit, "bgms")) {
    names(fit)[names(fit) == "components"] <- "cluster_probabilities"
    names(fit)[names(fit) == "allocations"] <- "cluster_allocations"
    names(fit)[names(fit) == "arguments"] <- "fit_arguments"
  }

  lambda <- fit$fit_arguments$lambda

  if (type == "simple") {
    if (is.null(b1) || is.null(b2)) {
      stop("For the simple type, both b1 and b2, indicating the number of clusters to be tested, must be provided.")
    }
    # Calculate prior odds in favor of b1 against b2
    prO <- (lambda^(b1 - b2) * factorial(b2)) / factorial(b1)

    # Calculate the posterior odds in favor of b1 against b2
    poO <-  unname(fit$cluster_probabilities[b1, 2]) / unname(fit$cluster_probabilities[b2, 2])

    bayesFactor <- poO / prO

  } else if (type == "complement") {
    # In favor of the complement
    prO <- (exp(lambda) - 1 - lambda) / lambda
    poO <- sum(fit$cluster_probabilities[-1, 2]) / unname(fit$cluster_probabilities[1, 2])
    bayesFactor <- poO / prO
  }

  return(round(bayesFactor, 1))
}


