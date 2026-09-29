# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not_cont, centrality, progress, 
                                 baseline_category, 
                                 ...){
  
  if(packageVersion("bgms") >= "0.2.0.0"){
    # Store original easybgm type before mapping
    original_type <- type
    
    # Map easybgm type to bgms variable_type
    if(length(type) > 1) {
      # Vector type: per-variable specification, map binary -> ordinal, pass to bgms
      variable_type <- type
      variable_type[variable_type == "binary"] <- "ordinal"
      # Store a simplified label for display
      original_type <- if(length(unique(type)) == 1) unique(type) else "mixed"
    } else if(type == "binary") {
      variable_type <- "ordinal"
    } else if(type == "mixed") {
      variable_type <- ifelse(not_cont == 1, "ordinal", "continuous")
    } else {
      variable_type <- type
    }
    
    bgm_args <- translate_bgm_prior_args(list(...))
    
    bgms_fit <- do.call(
      bgm, c(list(x = data, iter = iter,
                  variable_type = variable_type,
                  display_progress = progress, 
                  baseline_category = baseline_category), 
             bgm_args
      )
    )
  } else {
    if(type == "binary") {
      type <- "ordinal"
    }
    
    original_type <- type
    bgms_fit <- do.call(
      bgm, c(list(x = data, iter = iter,
                  variable_type = type,
                  display_progress = progress,
                  baseline_category = baseline_category,
                  ...))
    )
  }
  
  
  fit$model <- original_type
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
  
  if(packageVersion("bgms") >= "0.2.0.0"){
    # --- Ensure proper bgms object and variable names ---
    # Determine display-friendly model label
    if(length(type) > 1) {
      model_label <- if(length(unique(type)) == 1) unique(type) else "mixed"
    } else {
      model_label <- type
    }
  } else {
    model_label <- type
  }
  
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
  } else {
    args_tmp <- extract_arguments(fit)
    varnames <- args_tmp$data_columnnames
    if (is.null(varnames)) {
      varnames <- paste0("V", 1:args_tmp$num_variables)
    }
  }
  
  # Do not overwrite the class unconditionally: bgms (>= 0.2.0.0) fit objects
  # are S7 (class c("bgms", "S7_object")) and rely on the S7_object class for
  # `$`/`[[` access via S7::prop(). Stripping it to plain "bgms" breaks those
  # accessors. The object already inherits "bgms" in every path, so only set
  # the class in the (defensive) case where it somehow does not.
  if(!inherits(fit, "bgms")) class(fit) <- "bgms"
  
  # --- Extract model arguments and edge priors ---
  args <- bgms::extract_arguments(fit)
  args$save <- save
  # bgms names pairwise quantities after its own column names, which differ
  # from varnames for data without column names ("Variable 1" vs "V1")
  pair_varnames <- if (is.null(args$data_columnnames)) varnames else args$data_columnnames
  # extract SBM information
  bgms_res <- list()
  dots <- list(...)
  if("edge_selection" %in% names(dots)){
    if(dots$edge_selection == FALSE){
      args$edge_prior <- "None"
      args$edge_selection <- FALSE
    }
  }
  if (args$edge_prior[1] == "Stochastic-Block") {
    bgms_res$sbm <- extract_sbm(fit)
  }
  # extract the prior inclusion probabilities
  if (args$edge_prior[1] == "Bernoulli") {
    
    # Single global inclusion probability
    edge.prior <- args$inclusion_probability
    
    # needed for prior sensitivity check plot
    bgms_res$edge.prior <- edge.prior
    
  } else if (args$edge_prior[1] == "Beta-Bernoulli") {
    
    # Single global Beta–Bernoulli inclusion probability
    edge.prior <- args$beta_bernoulli_alpha /
      (args$beta_bernoulli_alpha + args$beta_bernoulli_beta)
    
    args$inclusion_probability <- edge.prior
    
    # needed for prior sensitivity check plot
    bgms_res$edge.prior <- edge.prior
    
  } else if (args$edge_prior[1] == "Stochastic-Block") {
    
    # Cluster assignments: a length-p vector like c(1,1,1,2,2)
    cl <- bgms_res$sbm$posterior_mean_allocations
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

    # The matrix above conditions on the estimated posterior partition. bgms
    # marginalizes over the prior partition instead, so its probabilities are
    # used where available; the matrix above is kept for older bgms.
    if (packageVersion("bgms") >= "0.2.0.0") {
      prior_mat <- extract_prior_inclusion_probabilities(fit)
    }
    
    # Vectorize consistent with bgms edge ordering (lower triangle)
    edge.prior <- prior_mat[lower.tri(prior_mat)]
    
    args$inclusion_probability <- edge.prior
    
    # needed for prior sensitivity check plot
    bgms_res$edge.prior <- edge.prior
    
  } else if(args$edge_prior[1] == "None"){
    # Does nothing but we just need to include it to avoid it running into the following else statement
  } else {
    
    stop("Unknown edge prior type in args$edge_prior[1].")
  }
  
  # --- Main extraction ---
  if (args$save) {
    if(packageVersion("bgms") >= "0.2.0.0"){
      p <- args$num_variables
    } else {
      p <- args$no_variables
    }
    pars <- extract_pairwise_interactions(fit)
    bgms_res$parameters <- vector2matrix_named(colMeans(pars), pair_varnames)
    bgms_res$samples_posterior <- extract_pairwise_interactions(fit)
    if(packageVersion("bgms") >= "0.2.0.0"){
      bgms_res$thresholds <- extract_main_effects(fit)
    } else {
      bgms_res$thresholds <- extract_category_thresholds(fit)
    }
    rownames(bgms_res$parameters) <- colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = p, nrow = p)
    if (args$edge_selection) {
      bgms_res$inc_probs <- extract_posterior_inclusion_probabilities(fit)
      if (packageVersion("bgms") >= "0.2.0.0") {
        # bgms marginalizes the prior inclusion odds over the prior partition of
        # the stochastic block prior (not the estimated posterior partition used
        # below), and runs a prior-only chain where no closed form exists, so
        # its Bayes factor is correct for every edge prior. The arithmetic below
        # is kept for older bgms, which lacks this extractor.
        bgms_res$inc_BF <- extract_inclusion_bf(fit)
        diag(bgms_res$inc_BF) <- 0
      } else if (args$edge_prior[1] == "Bernoulli") {
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
    if(packageVersion("bgms") >= "0.2.0.0"){
      p <- args$num_variables
    } else {
      p <- args$no_variables
    }
    pars <- extract_pairwise_interactions(fit)
    bgms_res$parameters <- vector2matrix_named(colMeans(pars), pair_varnames)
    if(packageVersion("bgms") >= "0.2.0.0"){
      bgms_res$thresholds <- extract_main_effects(fit)
    } else {
      bgms_res$thresholds <- extract_category_thresholds(fit)
    }
    rownames(bgms_res$parameters) <- colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters),
                                 nrow = nrow(bgms_res$parameters))
    if (args$edge_selection) {
      bgms_res$inc_probs <- extract_posterior_inclusion_probabilities(fit)
      if (packageVersion("bgms") >= "0.2.0.0") {
        # bgms marginalizes the prior inclusion odds over the prior partition of
        # the stochastic block prior (not the estimated posterior partition used
        # below), and runs a prior-only chain where no closed form exists, so
        # its Bayes factor is correct for every edge prior. The arithmetic below
        # is kept for older bgms, which lacks this extractor.
        bgms_res$inc_BF <- extract_inclusion_bf(fit)
        diag(bgms_res$inc_BF) <- 0
      } else if (args$edge_prior[1] == "Bernoulli") {
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
      # median probability model; inc_probs has a zero diagonal, so this does too
      bgms_res$structure <- 1 * (bgms_res$inc_probs > 0.5)
      gammas <- extract_indicators(fit)
      structures <- apply(gammas, 1, paste0, collapse = "")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[, 2] / nrow(gammas)
      bgms_res$graph_weights <- table_structures[, 2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
      # finalize output
      rownames(bgms_res$inc_probs) <- colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
      rownames(bgms_res$inc_BF) <- colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters)
    }
  }
  # --- Optionally compute centrality ---
  if (centrality) {
    # draws are placed on their edges by name: bgms stores the pairs of mixed
    # models grouped by variable type rather than in lower-triangle order
    bgms_res$centrality <- centrality(bgms_res, bycolumn = FALSE,
                                      varnames = pair_varnames)
  }
  
  # --- Compute convergence diagnostics ---
  if (args$edge_selection == TRUE) {
    # extract the Rhat
    bgms_res$convergence_parameter <- extract_rhat(fit)$pairwise
    # Calculate the MC uncertainty of the inclusion BF.
    #
    # inc_probs above is bgms's Rao-Blackwellized inclusion probability, so the
    # Monte Carlo error has to come from that same estimator. bgms reports it in
    # posterior_summary_indicator$mcse, already on the RB scale; it is placed
    # by name into the lower-triangle order that inc_BF[lower.tri()] uses below
    # (bgms lists the pairs of mixed models by variable type). Passing the raw
    # indicator average to BF_MCSE() instead (as easybgm did before) combined a
    # binomial variance belonging to the raw average with an effective sample
    # size computed for the smoother RB chain, which inflated the interval by up
    # to ~1.4x. BF_MCSE() is kept for the BDgraph path and for older bgms.
    bf_vec <- bgms_res$inc_BF[lower.tri(bgms_res$inc_BF)]
    ind_summary <- tryCatch(fit$posterior_summary_indicator,
                            error = function(e) NULL)
    mcse_rb <- ind_summary$mcse

    if (!is.null(mcse_rb) && length(mcse_rb) == length(bf_vec)) {
      # bgms lists the pairs of mixed models grouped by variable type, so the
      # MCSEs are placed into the lower-triangle order of bf_vec by name
      pair_index <- which(lower.tri(bgms_res$inc_BF), arr.ind = TRUE)
      if (!is.null(rownames(ind_summary))) {
        mcse_mat <- vector2matrix_named(stats::setNames(mcse_rb, rownames(ind_summary)),
                                        pair_varnames)
        mcse_rb <- mcse_mat[lower.tri(mcse_mat)]
      }
      p_rb <- bgms_res$inc_probs[lower.tri(bgms_res$inc_probs)]
      # delta method onto the log-BF scale: se(logBF) = se(p) / (p * (1 - p))
      se_log <- mcse_rb / (p_rb * (1 - p_rb))
      se_log[!is.finite(se_log)] <- NA_real_
      # bgms reports NA where the RB draws are constant to double precision,
      # which is routine on decisive edges; those cells stay NA here.
      z <- stats::qnorm(1 - (1 - 0.95) / 2)
      valid <- is.finite(bf_vec) & bf_vec > 0 & !is.na(se_log)
      bgms_res$MCSE_BF <- data.frame(
        lower = ifelse(valid, exp(log(bf_vec) - z * se_log), NA_real_),
        upper = ifelse(valid, exp(log(bf_vec) + z * se_log), NA_real_)
      )
      if (!is.null(rownames(ind_summary))) {
        rownames(bgms_res$MCSE_BF) <- paste(pair_varnames[pair_index[, "col"]],
                                            pair_varnames[pair_index[, "row"]], sep = "-")
      }
    } else {
      # bgms < 0.2.0.0 does not report an RB Monte Carlo error
      bgms_res$MCSE_BF <- BF_MCSE(gamma_mat = extract_indicators(fit),
                                  BF_vec = bf_vec,
                                  ess = extract_ess(fit)$indicator,
                                  return = "ci",
                                  smooth_bf = FALSE)
    }
  }
  
  if(packageVersion("bgms") >= "0.2.0.0"){
    # --- Extract interpretable parameter scales ---
    # These return NULL when the model type doesn't support them.
    # tryCatch guards against edge cases (e.g., mixed models with only
    # one variable of a given type, where the sub-matrix is a scalar).
    bgms_res$partial_correlations <- tryCatch(
      extract_partial_correlations(fit), error = function(e) NULL)
    bgms_res$precision_matrix <- tryCatch(
      extract_precision(fit), error = function(e) NULL)
    bgms_res$log_odds <- tryCatch(
      extract_log_odds(fit), error = function(e) NULL)

    # --- Blume-Capel main effects ---
    # The linear and quadratic effects of a Blume-Capel variable are of
    # substantive interest rather than nuisance parameters, so they get their
    # own table. bgms names them "cat (1)" and "cat (2)" in the main-effect
    # matrix, the same headers it uses for genuine category thresholds, so
    # relabel what is stored in $thresholds as well.
    bc_names <- bc_variable_names(args, varnames)
    if(length(bc_names) > 0){
      bgms_res$thresholds <- label_bc_thresholds(bgms_res$thresholds, bc_names)
      bc <- tryCatch(extract_blume_capel(fit, varnames, args, save = save),
                     error = function(e) NULL)
      bgms_res$blume_capel_parameters <- bc$table
      if(save) bgms_res$samples_blume_capel <- bc$samples
    }
  }
  
  # --- Finalize output ---
  bgms_res$model <- model_label
  bgms_res$fit_arguments <- args
  # the underlying bgms fit, so bgms extractors can be called without refitting
  bgms_res$packagefit <- fit
  output <- bgms_res
  class(output) <- c("package_bgms", "easybgm")
  return(output)
}
