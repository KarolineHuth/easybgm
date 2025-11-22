# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not_cont, centrality, progress, ...){

  if(!save && centrality){
    save <- TRUE
  }

  if(type == "binary") {
    type <- "ordinal"
  }

  if(packageVersion("bgms") > "0.1.4.2"){
    bgms_fit <- do.call(
      bgm, c(list(x = data, iter = iter,
                  variable_type = type,
                  display_progress = progress,
                  ...))
    )}
  if(packageVersion("bgms") < "0.1.6"){
    bgms_fit <- do.call(
      bgm, c(list(x = data, iter = iter, save = T,
                  variable_type = type,
                  display_progress = progress,
                  ...))
    )}


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
bgm_extract.package_bgms <- function(fit, type, save, iter, 
                                     not_cont, data, centrality, ...){
  if (packageVersion("bgms") < "0.1.4") {
    stop("easybgm now requires bgms version 0.1.4 or higher.")
  }
  # --- Ensure proper bgms object and variable names ---
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
    class(fit) <- "bgms"
  } else {
    varnames <- fit$arguments$data_columnnames
    if (is.null(varnames)) {
      varnames <- paste0("V", 1:fit$arguments$no_variables)}
  }
  if(packageVersion("bgms") > "0.1.4.2"){
    class(fit) <- "bgms"
  }

  # --- Extract model arguments and edge priors ---
  args <- bgms::extract_arguments(fit)
  args$save <- save
  # extract SBM information
  bgms_res <- list()
  if (args$edge_prior[1] == "Stochastic-Block" && packageVersion("bgms") > "0.1.6") {
    bgms_res$sbm <- extract_sbm(fit)
  }
  # extract the prior inclusion probabilities
  if (args$edge_prior[1] == "Bernoulli") {

    # Single global inclusion probability
    edge.prior <- args$inclusion_probability

  } else if (args$edge_prior[1] == "Beta-Bernoulli") {

    # Single global Beta–Bernoulli inclusion probability
    edge.prior <- args$beta_bernoulli_alpha /
      (args$beta_bernoulli_alpha + args$beta_bernoulli_beta)

    args$inclusion_probability <- edge.prior

  } else if (args$edge_prior[1] == "Stochastic-Block") {

    # Cluster assignments: a length-p vector like c(1,1,1,2,2)
    cl <-bgms_res$sbm$posterior_mean_allocations
    p  <- length(cl)

    # --- WITHIN-cluster inclusion probability ---
    pi_within <- args$beta_bernoulli_alpha /
      (args$beta_bernoulli_alpha + args$beta_bernoulli_beta)

    # --- BETWEEN-cluster inclusion probability ---
    # If no separate "between" hyperparameters given, fall back to within-cluster value
    if (is.null(args$beta_bernoulli_alpha_between) ||
        is.null(args$beta_bernoulli_beta_between)) {

      pi_between <- pi_within

    } else {

      pi_between <- args$beta_bernoulli_alpha_between /
        (args$beta_bernoulli_alpha_between + args$beta_bernoulli_beta_between)

    }

    # Matrix indicating whether node pairs share the same cluster
    same_cluster <- outer(cl, cl, "==")

    # Edge-specific prior inclusion probabilities (p × p matrix)
    prior_mat <- ifelse(same_cluster, pi_within, pi_between)

    # Vectorize consistent with bgms edge ordering (lower triangle)
    edge.prior <- prior_mat[lower.tri(prior_mat)]

    args$inclusion_probability <- edge.prior

  } else {

    stop("Unknown edge prior type in args$edge_prior[1].")
  }
  # --- Main extraction ---
  if (args$save) {
    p <- args$no_variables
    pars <- extract_pairwise_interactions(fit)
    bgms_res$parameters <- vector2matrix(colMeans(pars), p = p)
    bgms_res$samples_posterior <- extract_pairwise_interactions(fit)
    bgms_res$thresholds <- extract_category_thresholds(fit)
    colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = p, nrow = p)
    if (args$edge_selection) {
      bgms_res$inc_probs <- extract_posterior_inclusion_probabilities(fit)
      if (args$edge_prior[1] == "Bernoulli") {
        bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
          (edge.prior / (1 - edge.prior))
      } else if (args$edge_prior[1] == "Beta-Bernoulli") {
        bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
          (args$beta_bernoulli_alpha / args$beta_bernoulli_beta)
      } else if (args$edge_prior[1] == "Stochastic-Block") {
        # when there is only one set of hyperparameters for the beta bernoulli
        if (is.null(args$beta_bernoulli_alpha_between) | packageVersion("bgms") < "0.1.6") {
          bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
            (args$beta_bernoulli_alpha / args$beta_bernoulli_beta)
        } else {
          # when there are within and between hyperparameters
          cl <- bgms_res$sbm$posterior_mean_allocations
          p <- length(cl)
          # Matrix saying whether two nodes are in the same cluster
          same_cluster <- outer(cl, cl, "==")
          # Prior odds: within-cluster vs between-cluster
          prior_within  <- args$beta_bernoulli_alpha        / args$beta_bernoulli_beta
          prior_between <- args$beta_bernoulli_alpha_between / args$beta_bernoulli_beta_between
          prior_odds <- ifelse(same_cluster, prior_within, prior_between)
          # Finally compute BFs edgewise
          post_odds <- bgms_res$inc_probs / (1 - bgms_res$inc_probs)
          bgms_res$inc_BF <- post_odds / prior_odds
        }
      } else {
        stop("Unknown edge prior type.")
      }
      # this is both for the BB and SBM (however, when the SBM has
      # within and between beta hyperparameters, this needs to be adjusted
      bgms_res$structure <- 1 * (bgms_res$inc_probs > 0.5)
      gammas <- extract_indicators(fit)
      structures <- apply(gammas, 1, paste0, collapse = "")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[, 2] / nrow(gammas)
      bgms_res$graph_weights <- table_structures[, 2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
    }
  } else {
    p <- args$no_variables
    pars <- extract_pairwise_interactions(fit)
    bgms_res$parameters <- vector2matrix(colMeans(pars), p = p)
    bgms_res$thresholds <- extract_category_thresholds(fit)
    colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters),
                                 nrow = nrow(bgms_res$parameters))
    if (args$edge_selection) {
      bgms_res$inc_probs <- extract_posterior_inclusion_probabilities(fit)
      if (args$edge_prior[1] == "Bernoulli") {
        bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
          (edge.prior / (1 - edge.prior))
      } else if (args$edge_prior[1] == "Beta-Bernoulli") {
        bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
          (args$beta_bernoulli_alpha / args$beta_bernoulli_beta)
      } else if (args$edge_prior[1] == "Stochastic-Block") {
        # when there is only one set of hyperparameters for the beta bernoulli
        if (is.null(args$beta_bernoulli_alpha_between)) {
          bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
            (args$beta_bernoulli_alpha / args$beta_bernoulli_beta)
        } else {
          # when there are within and between hyperparameters
          cl <- bgms_res$sbm$posterior_mean_allocations
          p <- length(cl)
          # Matrix saying whether two nodes are in the same cluster
          same_cluster <- outer(cl, cl, "==")
          # Prior odds: within-cluster vs between-cluster
          prior_within  <- args$beta_bernoulli_alpha        / args$beta_bernoulli_beta
          prior_between <- args$beta_bernoulli_alpha_between / args$beta_bernoulli_beta_between
          prior_odds <- ifelse(same_cluster, prior_within, prior_between)
          # Finally compute BFs edgewise
          post_odds <- bgms_res$inc_probs / (1 - bgms_res$inc_probs)
          bgms_res$inc_BF <- post_odds / prior_odds
        }
      } else {
        stop("Unknown edge prior type.")
      }
      gammas <- extract_indicators(fit)
      structures <- apply(gammas, 1, paste0, collapse = "")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[, 2] / nrow(gammas)
      bgms_res$graph_weights <- table_structures[, 2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
      # finalize output
      colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
      colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters)
    }
  }
  # --- Optionally compute centrality ---
  if (centrality) {
    bgms_res$centrality <- centrality(bgms_res)
  }
  # --- For newer version compute convergence ---
  if (packageVersion("bgms") > "0.1.4.2") {
    # extract the Rhat
    bgms_res$convergence_parameter <-  fit$posterior_summary_pairwise$Rhat
    # calculate MC uncertainty
    bgms_res$MCSE_BF <-BF_MCSE(gamma_mat = extract_indicators(fit),
                               BF_vec = bgms_res$inc_BF[lower.tri(bgms_res$inc_BF)],
                               ess = fit$posterior_summary_indicator$n_eff,
                               return = "ci",
                               smooth_bf = FALSE)
  }
  # --- Finalize output ---
  bgms_res$model <- type
  bgms_res$fit_arguments <- args
  output <- bgms_res
  class(output) <- c("package_bgms", "easybgm")
  return(output)
}
