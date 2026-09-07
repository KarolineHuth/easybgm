# 1. Turns vector into matrix
vector2matrix <- function(vec, p, diag = FALSE, bycolumn = FALSE) {
  m <- matrix(0, p, p)
  
  if(!bycolumn){
    m[lower.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
  } else {
    m[upper.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
  }
  return(m)
}

# 2. Transform precision into partial correlations for interpretation
pr2pc <- function(K) {
  D.Prec = diag(diag(K)^(-.5))
  R <- diag(2,dim(K)[1])-D.Prec%*%K%*%D.Prec
  colnames(R) <- colnames(K)
  rownames(R) <- rownames(K)
  return(R)
}

# 3. BDgraph stores graphs as byte strings for efficiency
string2graph <- function(Gchar, p) {
  Gvec <- rep(0, p*(p-1)/2)
  edges <- which(unlist(strsplit(as.character(Gchar), "")) == 1)
  Gvec[edges] = 1
  G <- matrix(0, p, p)
  G[upper.tri(G)] <- Gvec
  G <- G + t(G)
  diag(G) <- 0
  return(G)
}

# 4. BDgraph extract posterior distribution for estimates
extract_posterior <- function(fit, data, iter = 10000,method = c("ggm", "gcgm"), posterior_method = c("maximum-posterior", "model-averaged"), not_cont){
  m <- length(fit$all_graphs)
  n <- nrow(data)
  p <- ncol(data)
  
  
  if(method == "gcgm") {
    S <- BDgraph::get_S_n_p(data, method = method, n = n, not.cont = not_cont)$S
  } else {
    S <- t(data) %*% data
  }
  
  if(posterior_method == "MAP"){
    Rs = matrix(0, nrow = iter, ncol = (p*(p-1))/2)
    index <- which.max(fit$graph_weights)
    graph_ix <- fit$sample_graphs[index]
    G <- string2graph(graph_ix, p)
    K <- BDgraph::rgwish(n=iter, adj=G, b=3+n, D=diag(p) + S)
    Rs <- list2matrix(K, p, convert = T)
  }
  if(posterior_method == "model-averaged"){
    
    n_structures <- length(fit$graph_weights)
    sum_weights <- sum(fit$graph_weights)
    structure_weights <- fit$graph_weights/sum_weights
    
    for (i in 1:n_structures) {
      graph_ix <- fit$sample_graphs[i]
      n_samples <- round(structure_weights[i]*iter)
      G <- string2graph(graph_ix, p)
      K <- BDgraph::rgwish(n=n_samples, adj=G, b=3+n, D=diag(p) + S)
      samples_ix <- list2matrix(K, p, convert = T)
      if(i ==1){
        Rs = as.matrix(samples_ix)
      } else{
        Rs <- rbind(Rs, samples_ix)
      }
      
    }
  }
  
  return(list(Rs))
}

# 5. Samples from the G-wishart distribution
gwish_samples <- function(G, S, nsamples=1000) {
  n <- nrow(S)
  p <- ncol(S)
  #Rs <- array(0, dim=c(nsamples, p, p))
  Rs = matrix(0, nrow = nsamples, ncol = (p*(p-1))/2)
  
  for (i in 1:nsamples) {
    K <- BDgraph::rgwish(n=1, adj=G, b=3+n, D=diag(p) + S)*(G + diag(p))
    Rs[i,] <- as.vector(pr2pc(K)[upper.tri(pr2pc(K))])
    #Rs[i,,] <- .pr2pc(K)
  }
  return(Rs)
}


# 6. Centrality of weighted graphs

# Strength centrality only ## FASTER CODE
# `bycolumn` says how the posterior sample vectors are laid out and is passed on
# to vector2matrix(). FALSE (the default) fills the lower triangle column-wise,
# the order bgms stores pairwise interactions in ((1,2), (1,3), (1,4), (2,3),
# ...). TRUE fills the upper triangle, which is the order BGGM uses ((1,2),
# (1,3), (2,3), (1,4), ...). The caller has to state this: the two orders
# differ, and getting it wrong silently permutes the per-node strengths.
centrality <- function(res, bycolumn = FALSE){
  Nsamples <- nrow(res$samples_posterior)
  p <- nrow(res$parameters)
  strength_samples <- matrix(0, nrow = Nsamples, ncol = p)
  for(i in 1:Nsamples){
    strength_samples[i, ] <- rowSums(abs(vector2matrix(res$samples_posterior[i,],
                                                       p, bycolumn = bycolumn)))
  }
  return(strength_samples)
}

## NOTE NOT USED
# Strength, betweenness and closeness centrality ## SLOWER CODE
# centrality_all <- function(res){
#   Nsamples <- nrow(res$samples_posterior)
#   p <- as.numeric(nrow(res$parameters))
#   samples <- res$samples_posterior
#   for(i in 1:Nsamples){
#     
#     #Strength
#     strength_samples <- rowSums(abs(vector2matrix(samples[i, ], p, bycolumn = TRUE)))
#     #EI
#     influence_samples <- rowSums(vector2matrix(samples[i, ], p, bycolumn = TRUE))
#     
#     DistMat <- 1/(ifelse(abs(vector2matrix(samples[i, ], p, bycolumn = TRUE))==0,0,abs(vector2matrix(samples[i, ], p, bycolumn = T))))
#     igraphObject <- igraph::graph.adjacency(DistMat, weighted = TRUE, mode = "undirected")
#     # Closeness
#     closeness_samples <- igraph::closeness(igraphObject)
#     # Betweenness
#     betweenness_samples <- igraph::estimate_betweenness(igraphObject,cutoff = 1/1e-10)
#     
#     if(i > 1){
#       centrality_samples <- rbind(centrality_samples, cbind(c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
#                                                             rbind(strength_samples, closeness_samples, betweenness_samples, influence_samples)))
#     } else {
#       centrality_samples <- cbind(c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
#                                   rbind(strength_samples, closeness_samples, betweenness_samples, influence_samples))
#     }
#   }
#   colnames(centrality_samples) <- c("Centrality", colnames(res$parameters))
#   centrality_samples <- as.data.frame(centrality_samples)
#   centrality_samples[, 2:(p+1)] <- sapply(centrality_samples[, 2:(p+1)], as.numeric)
#   return(centrality_samples)
# }
# 

# 7. turn list into matrix
list2matrix <- function(obj, p, convert = FALSE) {
  nlist <- length(obj)/(p*p)
  m <- obj[, , 1]
  nest <- sum(lower.tri(m))
  res <- matrix(0, nrow = nlist, ncol = nest)
  for(i in 1:nlist){
    m <- obj[, , i]
    if(convert == T){
      m <- pr2pc(m)
    }
    res[i, ] <- as.vector(m[lower.tri(m)])
  }
  return(res)
}

# 8. Set defaults of a function
set_defaults <- function(args, ...) {
  dots <- list(...)
  def_args <- setdiff(names(args), names(dots))
  dots[def_args] <- args[def_args]
  
  return(dots)
}

# 9. Check for empty ...

dots_check <- function(...){
  if(...length() > 0){
    warning("Arguments specified with ... are unused. ")
  }
}

# 10. Translate legacy flat prior args into bgms (>= 0.2.0.0) prior objects.
# Why: bgms now expects prior-class objects (cauchy_prior(), beta_prime_prior(),
# bernoulli_prior(), etc.). Passing the old flat args still works but triggers
# lifecycle::deprecate_warn() on every call. Folding them here keeps the
# easybgm user-facing API unchanged while suppressing those warnings.
translate_bgm_prior_args <- function(dots) {
  out <- dots

  # interaction_prior: from pairwise_scale, or from the older interaction_scale
  # that easybgm 0.4.0 documented. bgms maps both to a Cauchy in its own
  # deprecation shim; intercepting them here keeps that mapping but avoids
  # emitting bgms's deprecation warning for an argument easybgm told users to use.
  if(!"interaction_prior" %in% names(out)) {
    scale_arg <- out$pairwise_scale %||% out$interaction_scale
    if(!is.null(scale_arg)) {
      out$interaction_prior <- bgms::cauchy_prior(scale = scale_arg)
    }
  }
  out$pairwise_scale <- NULL
  out$interaction_scale <- NULL

  # threshold_prior: from main_alpha / main_beta (or legacy threshold_alpha / threshold_beta)
  has_main <- any(c("main_alpha", "main_beta",
                    "threshold_alpha", "threshold_beta") %in% names(out))
  if(!"threshold_prior" %in% names(out) && has_main) {
    a <- out$main_alpha %||% out$threshold_alpha %||% 0.5
    b <- out$main_beta  %||% out$threshold_beta  %||% 0.5
    out$threshold_prior <- bgms::beta_prime_prior(alpha = a, beta = b)
  }
  out$main_alpha <- out$main_beta <- NULL
  out$threshold_alpha <- out$threshold_beta <- NULL

  # edge_prior: from character string + flat hyperparameters
  if("edge_prior" %in% names(out) && is.character(out$edge_prior)) {
    ep_str <- out$edge_prior
    ip <- out$inclusion_probability %||% 0.5
    bba <- out$beta_bernoulli_alpha %||% 1
    bbb <- out$beta_bernoulli_beta %||% 1
    bbab <- out$beta_bernoulli_alpha_between %||% 1
    bbbb <- out$beta_bernoulli_beta_between %||% 1
    da  <- out$dirichlet_alpha %||% 1
    lam <- out$lambda %||% 1

    out$edge_prior <- switch(ep_str,
      "Bernoulli"        = bgms::bernoulli_prior(inclusion_probability = ip),
      "Beta-Bernoulli"   = bgms::beta_bernoulli_prior(alpha = bba, beta = bbb),
      "Stochastic-Block" = bgms::sbm_prior(alpha = bba, beta = bbb,
                                           alpha_between = bbab, beta_between = bbbb,
                                           dirichlet_alpha = da, lambda = lam),
      stop("Unknown edge_prior '", ep_str, "'.")
    )
  } else if(!"edge_prior" %in% names(out) &&
            any(c("inclusion_probability", "beta_bernoulli_alpha",
                  "beta_bernoulli_beta") %in% names(out))) {
    # User supplied flat hyperparameters without naming the family: assume Bernoulli
    # if inclusion_probability is given, else Beta-Bernoulli.
    if("inclusion_probability" %in% names(out)) {
      out$edge_prior <- bgms::bernoulli_prior(
        inclusion_probability = out$inclusion_probability
      )
    } else {
      out$edge_prior <- bgms::beta_bernoulli_prior(
        alpha = out$beta_bernoulli_alpha %||% 1,
        beta  = out$beta_bernoulli_beta  %||% 1
      )
    }
  }
  out$inclusion_probability <- NULL
  out$beta_bernoulli_alpha <- out$beta_bernoulli_beta <- NULL
  out$beta_bernoulli_alpha_between <- out$beta_bernoulli_beta_between <- NULL
  out$dirichlet_alpha <- out$lambda <- NULL

  out
}

# Same for bgmCompare: difference_prior is the indicator prior on group differences.
translate_bgmcompare_prior_args <- function(dots) {
  out <- dots

  # interaction_prior (baseline) and threshold_prior: same translation as bgm()
  if(!"interaction_prior" %in% names(out)) {
    scale_arg <- out$pairwise_scale %||% out$interaction_scale
    if(!is.null(scale_arg)) {
      out$interaction_prior <- bgms::cauchy_prior(scale = scale_arg)
    }
  }
  out$pairwise_scale <- NULL
  out$interaction_scale <- NULL

  has_main <- any(c("main_alpha", "main_beta",
                    "threshold_alpha", "threshold_beta") %in% names(out))
  if(!"threshold_prior" %in% names(out) && has_main) {
    a <- out$main_alpha %||% out$threshold_alpha %||% 0.5
    b <- out$main_beta  %||% out$threshold_beta  %||% 0.5
    out$threshold_prior <- bgms::beta_prime_prior(alpha = a, beta = b)
  }
  out$main_alpha <- out$main_beta <- NULL
  out$threshold_alpha <- out$threshold_beta <- NULL

  # difference_prior: from a character string + flat hyperparameters
  if("difference_prior" %in% names(out) && is.character(out$difference_prior)) {
    dp_str <- out$difference_prior
    dprob <- out$difference_probability %||% 0.5
    bba <- out$beta_bernoulli_alpha %||% 1
    bbb <- out$beta_bernoulli_beta %||% 1

    out$difference_prior <- switch(dp_str,
      "Bernoulli"      = bgms::bernoulli_prior(inclusion_probability = dprob),
      "Beta-Bernoulli" = bgms::beta_bernoulli_prior(alpha = bba, beta = bbb),
      stop("Unknown difference_prior '", dp_str, "'.")
    )
  } else if(!"difference_prior" %in% names(out) &&
            any(c("difference_probability", "beta_bernoulli_alpha",
                  "beta_bernoulli_beta") %in% names(out))) {
    if("difference_probability" %in% names(out)) {
      out$difference_prior <- bgms::bernoulli_prior(
        inclusion_probability = out$difference_probability
      )
    } else {
      out$difference_prior <- bgms::beta_bernoulli_prior(
        alpha = out$beta_bernoulli_alpha %||% 1,
        beta  = out$beta_bernoulli_beta  %||% 1
      )
    }
  }
  out$difference_probability <- NULL
  out$beta_bernoulli_alpha <- out$beta_bernoulli_beta <- NULL

  out
}

# Null-coalescing helper used by the translators above.
`%||%` <- function(a, b) if(is.null(a)) b else a


# Given alpha and beta parameters computes the (log)probability of the beta Bernoulli distribution
# for the bgms package returns probability for a specific complexity
# @args alpha, beta parameters, c the complexity for which the probability has to b computed
# p the number of nodes
beta_bernoulli_prob <- function(c, alpha, beta, p) {
  
  k <- p * (p - 1) / 2
  
  log_nom <- lbeta(alpha + c, beta + k - c)
  log_denom <- lbeta(alpha, beta)
  
  log_choose <- lchoose(k, c)
  
  
  log_prob <- log_choose + log_nom - log_denom
  
  return(log_prob)
}

# Replacement for fits produced under bgms 0.2.0.0's shifted-Poisson SBM.
# Source this file to use the replacement in the current R session, or put the
# function AND its two internal helpers into the easybgm package source.
# This file does not modify installed packages.

# P(B=b | T=t), using the unbounded shifted-Poisson prior. Rows 1:p are
# returned, but each column is normalized over the full support b >= t.
.sbm_count_kernels <- function(p, lambda, alpha, tol = 1e-12) {
  logsumexp <- function(x) {
    m <- max(x)
    if (m == -Inf) return(-Inf)
    m + log(sum(exp(x - m)))
  }
  top <- max(p + 10, ceiling(lambda + 10 * sqrt(lambda) + 10))
  for (attempt in seq_len(20)) {
    if (!is.finite(top) || top > 1e6)
      stop("Cannot evaluate the count-prior tail accurately within the numerical limit.")
    b <- seq_len(top)
    base <- stats::dpois(b - 1, lambda, log = TRUE) +
      lgamma(alpha * b) - lgamma(alpha * b + p)
    log_full <- matrix(-Inf, p, p)
    truncated <- matrix(0, p, p)
    tail_bound <- numeric(p)
    for (t in seq_len(p)) {
      keep <- b >= t
      log_weight <- rep(-Inf, top)
      log_weight[keep] <- base[keep] + lgamma(b[keep] + 1) -
        lgamma(b[keep] - t + 1)
      log_Z <- logsumexp(log_weight)
      log_Z_p <- logsumexp(log_weight[seq_len(p)])
      log_full[, t] <- log_weight[seq_len(p)] - log_Z
      truncated[, t] <- exp(log_weight[seq_len(p)] - log_Z_p)
      
      # Successive unnormalized weights beyond b=top have ratios bounded by
      # u = lambda/top * (top+1)/(top+1-t). The remaining rising-factorial
      # ratio is <= 1. This bound decreases with b, so the tail is geometric.
      u <- lambda / top * (top + 1) / (top + 1 - t)
      tail_bound[t] <- if (u >= 1) Inf else
        exp(log_weight[top] + log(u) - log1p(-u) - log_Z)
    }
    if (all(tail_bound <= tol))
      return(list(log_full = log_full, truncated = truncated))
    top <- 2 * top
  }
  stop("Count-prior normalizers did not reach the requested numerical accuracy.")
}

# The old summary is r = Q_truncated %*% P(T | data). Q is lower triangular,
# because t occupied blocks require at least t available components.
.sbm_occupied_probabilities <- function(reported, Q, allocations = NULL) {
  p <- length(reported)
  if (!is.null(allocations)) {
    if (is.matrix(allocations)) allocations <- list(allocations)
    if (!is.list(allocations) || !length(allocations))
      stop("Saved cluster allocations have an unsupported structure.")
    counts <- unlist(lapply(allocations, function(z) {
      if (!is.matrix(z) || ncol(z) != p || !nrow(z) || anyNA(z))
        stop("Saved cluster allocations do not match the SBM summary.")
      apply(z, 1, function(row) length(unique(row)))
    }), use.names = FALSE)
    weights <- tabulate(counts, nbins = p) / length(counts)
  } else {
    if (rcond(Q) < 1e-12)
      stop("The count summary is too ill-conditioned to invert; supply the raw bgms fit with saved allocations.")
    weights <- as.numeric(forwardsolve(Q, reported))
    if (any(!is.finite(weights)))
      stop("Could not recover the occupied-count posterior from this summary.")
    slack <- max(-min(c(weights, 0)), abs(sum(weights) - 1))
    if (slack > 1e-4)
      stop("This summary is not consistent with bgms 0.2.0.0's shifted-Poisson ",
           "count convention (discrepancy ", format(slack, digits = 2),
           "). If it came from another bgms version, recompute it from that fit.")
    if (slack > 1e-8)
      warning("The summary appears rounded (discrepancy ",
              format(slack, digits = 2), "); projecting onto the simplex.")
    weights <- pmax(weights, 0)
    weights <- weights / sum(weights)
  }
  if (max(abs(as.numeric(Q %*% weights) - reported)) > 1e-4)
    stop("Saved allocations/summary disagree with the assumed count-prior convention.")
  weights
}

#' Test whether a network splits into clusters
#'
#' For a network fitted with the Stochastic Block Model (SBM) edge prior, this
#' gives a Bayes factor for the number of clusters. With `type = "complement"`
#' it weighs more than one cluster against exactly one; with `type = "point"` it
#' weighs `b1` clusters against `b2`. Values above 1 favour the first of the
#' two; take the reciprocal to read the evidence the other way round.
#'
#' The count is the number of clusters the model has available, which can be
#' larger than the number that actually hold variables, since a cluster may come
#' out empty. For the memberships themselves, see [bgms::extract_sbm()].
#' Evidence for clustering concerns the network's edge structure and is not by
#' itself evidence of multidimensionality.
#'
#' @param fit A fit of class `easybgm` or `bgms`, fitted with
#'   `edge_prior = bgms::sbm_prior()`. Raw draws are used when the fit carries
#'   them, otherwise the stored summary, which must be at full precision.
#' @param type Either `"complement"`, the default, or `"point"`.
#' @param b1,b2 Whole numbers between 1 and the number of variables, required
#'   when `type = "point"`.
#' @return A single unrounded Bayes factor. `NA` with a warning if neither point
#'   hypothesis appears in the posterior. A result of 0 or `Inf` means the
#'   sampler never visited one of the two, so run more iterations rather than
#'   reading it as decisive.
#' @export

clusterBayesfactor <- function(fit, type = "complement", b1 = NULL, b2 = NULL) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("point", "complement"))
    stop("The type argument must be either 'point' or 'complement'.")
  
  allocations <- NULL
  if (inherits(fit, "easybgm")) {
    args <- fit$fit_arguments
    num_blocks <- fit$sbm$posterior_num_blocks
  } else if (inherits(fit, "bgms")) {
    args <- bgms::extract_arguments(fit)
    num_blocks <- bgms::extract_sbm(fit)$posterior_num_blocks
    allocations <- tryCatch(fit$raw_samples$allocations, error = function(e) NULL)
  } else {
    stop("fit must be a fitted object of class 'easybgm' or 'bgms'.")
  }
  if (is.null(num_blocks))
    stop("The fit has no SBM count summary. Fit the model with an SBM edge prior.")
  
  lambda <- args$lambda
  alpha <- args$dirichlet_alpha
  positive_scalar <- function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  if (!positive_scalar(lambda) || !positive_scalar(alpha))
    stop("The fit must record positive 'lambda' and 'dirichlet_alpha' values.")
  
  if (!(is.data.frame(num_blocks) || is.matrix(num_blocks)) ||
      ncol(num_blocks) != 1L || nrow(num_blocks) < 1L)
    stop("Expected a single-column posterior_num_blocks matrix or data frame.")
  post <- num_blocks[, 1]
  if (!is.numeric(post) || any(!is.finite(post)) || any(post < 0) ||
      abs(sum(post) - 1) > 1e-8)
    stop("The count summary must contain finite probabilities summing to one.")
  p <- length(post)
  if (type == "point") {
    valid_count <- function(b)
      is.numeric(b) && length(b) == 1L && is.finite(b) &&
      b >= 1 && b <= p && b == floor(b)
    if (!valid_count(b1) || !valid_count(b2))
      stop("For type='point', b1 and b2 must be integers between 1 and ", p, ".")
    if (b1 == b2) return(1)
  }
  
  kernels <- .sbm_count_kernels(p, lambda, alpha)
  weights <- .sbm_occupied_probabilities(post, kernels$truncated, allocations)
  log_probability <- function(b) {
    x <- log(weights) + kernels$log_full[b, ]
    m <- max(x)
    if (m == -Inf) return(-Inf)
    m + log(sum(exp(x - m)))
  }
  if (type == "point") {
    log_post_b1 <- log_probability(b1)
    log_post_b2 <- log_probability(b2)
    if (log_post_b1 == -Inf && log_post_b2 == -Inf) {
      warning("Neither point hypothesis has mass in the saved posterior; this Bayes factor cannot be estimated.")
      return(NA_real_)
    }
    log_prior_odds <- stats::dpois(b1 - 1, lambda, log = TRUE) -
      stats::dpois(b2 - 1, lambda, log = TRUE)
    log_bf <- log_post_b1 - log_post_b2 - log_prior_odds
  } else {
    log_p1 <- min(log_probability(1), 0)
    log_post_odds <- log(-expm1(log_p1)) - log_p1
    # log(expm1(lambda)), evaluated without overflowing when lambda is large.
    log_prior_odds <- lambda + log(-expm1(-lambda))
    log_bf <- log_post_odds - log_prior_odds
    if (!is.finite(log_bf))
      warning("The saved posterior gives (near) no draws to B = 1, so this Bayes ",
              "factor is unbounded. Report the posterior probabilities over the ",
              "number of blocks instead.")
  }
  unname(exp(log_bf))
}


# function for calculating the MC uncertainty for the inclusion BF

BF_MCSE <- function(gamma_mat,
                    BF_vec,
                    ess = NULL,
                    smooth_bf = FALSE,
                    return = c("mcse_log", "mcse_bf", "ci"),
                    conf_level = 0.95) {
  
  # Match argument
  return <- match.arg(return)
  
  T_samples <- nrow(gamma_mat)
  
  # Compute raw posterior inclusion probabilities
  p_raw <- apply(gamma_mat, 2, mean)
  
  # Determine p_hat for BF calculation
  if (smooth_bf) {
    pseudo <- c(0.5, 0.5)
    a <- pseudo[1]
    b <- if (length(pseudo) > 1) pseudo[2] else pseudo[1]
    p_for_bf <- (colSums(gamma_mat) + a) / (T_samples + a + b)
  } else {
    p_for_bf <- p_raw
  }
  
  # Compute Bayes factors
  #post_odds <- p_for_bf / (1 - p_for_bf)
  #BF_vec <- post_odds / prior_odds
  #BF_vec <- BF_vec[lower.tri(BF_vec)]
  
  # Compute effective sample size
  if (is.null(ess)) {
    ess_vec <- apply(gamma_mat, 2, function(x) {
      ess <- as.numeric(coda::effectiveSize(x))
      # Handle invalid ESS values
      if (is.na(ess) || !is.finite(ess) || ess <= 0) {
        return(1) # as a safe option (should never happen)
      }
      return(ess)
    })
  } else {
    ess_vec <- ess
  }
  
  if (smooth_bf) {
    a <- pseudo[1]
    b <- if (length(pseudo) > 1) pseudo[2] else pseudo[1]
    p_for_variance <- (colSums(gamma_mat) + a) / (T_samples + a + b)
  } else {
    p_for_variance <- p_raw
  }
  
  # Variance of p_hat
  var_p <- p_for_variance * (1 - p_for_variance) / ess_vec
  
  # Delta method to log BF scale
  denom <- p_for_variance^2 * (1 - p_for_variance)^2
  var_logBF <- var_p / denom
  se_logBF <- sqrt(var_logBF)
  
  # Set MCSE to NA where calculations are invalid
  se_logBF[!is.finite(se_logBF)] <- NA_real_
  
  # Return based on argument
  if (return == "mcse_log") {
    # Return MCSE on log scale
    return(se_logBF)
    
  } else if (return == "mcse_bf") {
    # Return MCSE on BF scale
    se_BF <- BF_vec * se_logBF
    se_BF[!is.finite(se_BF)] <- NA_real_
    
    return(se_BF)
    
  } else if (return == "ci") {
    # Return confidence intervals as data frame
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    
    # Compute CI on log scale, then exponentiate
    log_BF <- log(BF_vec)
    ci_lower <- exp(log_BF - z * se_logBF)
    ci_upper <- exp(log_BF + z * se_logBF)
    
    # Set CI to NA where BF or MCSE is invalid
    ci_lower[!is.finite(BF_vec) | is.na(se_logBF)] <- NA_real_
    ci_upper[!is.finite(BF_vec) | is.na(se_logBF)] <- NA_real_
    
    # Return as data frame with two columns
    ci_df <- data.frame(
      lower = ci_lower,
      upper = ci_upper
    )
    return(ci_df)
  }
}

# --------------------------------------------------------------------------------------------------
# Blume-Capel main effects (bgms >= 0.2.0.0)
# --------------------------------------------------------------------------------------------------
# A Blume-Capel variable has two main-effect parameters rather than one
# threshold per category: the threshold for category x of a variable with
# baseline b is mu(x) = linear * x + quadratic * (x - b)^2. bgms labels the
# columns of its main-effect matrix "cat (k)" for every discrete variable, so
# for Blume-Capel variables those headers name the wrong quantity. These
# helpers relabel what easybgm stores and pull the pair into its own table.

# Which variables were fitted as Blume-Capel? bgms returns 'variable_type' with
# one entry per variable when the user supplied a vector, and a single string
# when the same type was applied to all of them.
bc_variable_names <- function(args, varnames){
  variable_type <- rep_len(args$variable_type, length(varnames))
  varnames[variable_type == "blume-capel"]
}

# bgms recodes discrete scores to start at 0 and shifts the baseline category
# with them, so args$baseline_category is on the recoded scale. The offset is
# args$blume_capel_shift, which is NA for variables that are not Blume-Capel.
# Adding the two recovers the baseline on the scale the user supplied.
#
# Mind the indexing: bgms indexes 'baseline_category' and 'blume_capel_shift'
# over the discrete variables only, whereas 'variable_type' has one entry per
# variable. The two therefore differ in length as soon as the model holds a
# continuous variable, and lining them up positionally silently misreports the
# baseline. The returned vector is named, so callers can index it by variable.
bc_baseline_categories <- function(args, varnames){
  variable_type <- rep_len(args$variable_type, length(varnames))
  is_discrete <- variable_type != "continuous"
  discrete_names <- varnames[is_discrete]
  n <- length(discrete_names)

  align <- function(x, default){
    if(is.null(x)) return(rep(default, n))
    if(length(x) == n) return(x)
    # Tolerate a per-variable vector as well, in case bgms ever reports one.
    if(length(x) == length(varnames)) return(x[is_discrete])
    rep_len(x, n)
  }

  baseline <- align(args$baseline_category, 0)
  shift <- align(args$blume_capel_shift, 0)
  shift[is.na(shift)] <- 0

  stats::setNames(baseline + shift, discrete_names)
}

# Rename the Blume-Capel columns of the stored main-effect matrix. When
# Blume-Capel and ordinal variables share one matrix the column headers cannot
# describe both, so record the per-row meaning as an attribute instead.
label_bc_thresholds <- function(thresholds, bc_names){
  relabel <- function(m){
    if(is.null(m) || !is.matrix(m) || is.null(rownames(m))) return(m)
    is_bc <- rownames(m) %in% bc_names
    if(!any(is_bc)) return(m)
    if(all(is_bc) && ncol(m) == 2){
      colnames(m) <- c("linear", "quadratic")
    } else {
      attr(m, "variable_type") <- ifelse(is_bc, "blume-capel", "ordinal")
    }
    m
  }

  if(is.list(thresholds) && !is.null(thresholds$discrete)){
    thresholds$discrete <- relabel(thresholds$discrete)
    return(thresholds)
  }
  relabel(thresholds)
}

# Build the Blume-Capel parameter table, and optionally the posterior samples,
# from a fitted bgms object. Returns NULL when the fit holds no Blume-Capel
# variables or does not carry the posterior summary these fields are built from.
extract_blume_capel <- function(fit, varnames, args, save = FALSE){
  bc_names <- bc_variable_names(args, varnames)
  if(length(bc_names) == 0) return(NULL)

  summ <- fit$posterior_summary_main
  if(is.null(summ) || is.null(rownames(summ))) return(NULL)

  # Match on the full parameter label rather than on a suffix, so a variable
  # that happens to be named "x (linear)" cannot be mistaken for one.
  labels <- as.vector(rbind(paste0(bc_names, " (linear)"),
                            paste0(bc_names, " (quadratic)")))
  idx <- match(labels, rownames(summ))
  if(anyNA(idx)) return(NULL)

  baseline <- bc_baseline_categories(args, varnames)

  # Credible intervals come from the draws themselves; posterior_summary_main
  # reports only mean, sd, mcse, n_eff and Rhat.
  lower <- upper <- rep(NA_real_, length(labels))
  samples <- NULL
  raw <- tryCatch(fit$raw_samples, error = function(e) NULL)
  par_names <- raw$parameter_names$main
  if(!is.null(raw$main) && !is.null(par_names)){
    cols <- match(labels, par_names)
    if(!anyNA(cols)){
      draws <- do.call(rbind, raw$main)
      if(ncol(draws) == length(par_names)){
        draws <- draws[, cols, drop = FALSE]
        colnames(draws) <- labels
        quantiles <- apply(draws, 2, stats::quantile,
                           probs = c(0.025, 0.975), names = FALSE)
        lower <- quantiles[1, ]
        upper <- quantiles[2, ]
        if(save) samples <- draws
      }
    }
  }

  table <- data.frame(
    variable = rep(bc_names, each = 2),
    effect = rep(c("linear", "quadratic"), times = length(bc_names)),
    baseline = rep(unname(baseline[bc_names]), each = 2),
    estimate = summ$mean[idx],
    sd = summ$sd[idx],
    lower = lower,
    upper = upper,
    convergence = summ$Rhat[idx],
    row.names = NULL
  )
  colnames(table) <- c("Variable", "Effect", "Baseline Category", "Estimate",
                       "Posterior SD", "Lower 2.5%", "Upper 97.5%",
                       "Convergence")

  list(table = table, samples = samples)
}
