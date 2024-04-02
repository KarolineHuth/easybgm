#' @name sparse.or.dense
#' @title Test for sparse against dense topologies
#'
#' @description Used to test for a given data set if the corresponding network is sparse or dense.
#'
#' @param x data set 
#' @param type variable type, currently only binary and ordinal are supported
#' @param ... additional arguments of the bgms function
#'
#' @return List containing results of the analysis, the (log) Bayes Factor of sparse against dense, the 
#'  mean relative number of edges under both hypotheses and the uniform hypotheses, and the number of hypotheses.
#' @export
sparse.or.dense <- function(x, type, ...) {

  if ((type != "binary") & (type != "ordinal")) {
    stop("Wrong input provided. This function is only compatible with the package bgms yet,
         and therefore only for binary and ordinal variables.")
  }

  p <- ncol(x)
  k <- p * (p - 1) / 2

  default_args1 <- list(
    edge_prior = "Beta-Bernoulli",
    alpha_s = 1,
    beta_s = k^1.1 / 4,
    alpha_d = k^1.1 / 4,
    beta_d = 1,
    iter = 1e4
  )
  
  dots <- list(...)
  
  
  args <- set_defaults(default_args1, ... )
  
  if ((args$alpha_s != args$beta_d) | (args$alpha_d != args$beta_s)){
    warning("The hypotheses are not symmetric")
  }
  args_extra <- args[names(args) %in% 
                       c("alpha_s", "alpha_d", "beta_s", "beta_d") == FALSE]
  
  args_sparse <- c(list(x = x, save = TRUE, beta_bernoulli_alpha = args$alpha_s, 
                       beta_bernoulli_beta = args$beta_s), args_extra)
  
  args_dense <- c(list(x = x, save = TRUE, beta_bernoulli_alpha = args$alpha_d,
                       beta_bernoulli_beta = args$beta_d), args_extra)
  
  args_uniform <- c(list(x = x, save = TRUE, beta_bernoulli_alpha = 1,
                         beta_bernoulli_beta = 1), args_extra)
  
  print("Running the model under the sparse hypothesis")
  res_sparse <- do.call(bgm, args_sparse)
  
  print("Running the model under the dense hypothesis")
  res_dense <- do.call(bgm, args_dense)
  
  print("Running the model under the uniform hypothesis")
  res_uniform <- do.call(bgm, args_uniform)
  
  gamma_sums_sparse <- rowSums(res_sparse$gamma)
  gamma_sums_dense <- rowSums(res_dense$gamma)
  gamma_sums_uniform <- rowSums(res_uniform$gamma)
  
  tab_gamma_sparse <- tabulate(gamma_sums_sparse + 1, k + 1)
  tab_gamma_dense <- tabulate(gamma_sums_dense + 1, k + 1)
  tab_gamma_uniform <- tabulate(gamma_sums_uniform + 1, k + 1)
  
  tab_gamma_sparse <- list(tab = tab_gamma_sparse, 
                           alpha = args$alpha_s, beta = args$beta_s)
  tab_gamma_dense <- list(tab = tab_gamma_dense, 
                          alpha = args$alpha_d, beta = args$beta_d)
  tab_gamma_uniform <- list(tab = tab_gamma_uniform, 
                            alpha = 1, beta = 1)
  
  gamma_list <- list(tab_gamma_sparse, tab_gamma_uniform, tab_gamma_dense)
  overlap_check <- is_overlap(gamma_list)
  #compute necessary bridge hypotheses
  iter = 1
  while (is.list(overlap_check)) {
    alpha_in <- overlap_check$alpha
    beta_in <- overlap_check$beta
    before_index <- overlap_check$before_pos
    
    new_args <- args_uniform
    new_args$beta_bernoulli_alpha <- alpha_in
    new_args$beta_bernoulli_beta <- beta_in
    
    gamma_new <- rowSums(do.call(bgm, new_args)$gamma)
    
    tab_gamma <- tabulate(gamma_new + 1, k + 1)
    new_el <- list(tab = tab_gamma, alpha = alpha_in, beta = beta_in)
    new_el <- list(new_el)

    gamma_list <- append(gamma_list, new_el, after = before_index)
    overlap_check <- is_overlap(gamma_list)
    
    if (iter == 5) {
      break
    }
    iter = iter + 1
  }
  
  BF <- compute_bayes_factor(gamma_list, k)
  mean_complexity_uniform <- mean(gamma_sums_uniform)
  mean_complexity_sparse <- mean(gamma_sums_sparse)
  mean_complexity_dense <- mean(gamma_sums_dense)
  
  return(list(log.BF = BF, BF = exp(BF), 
              relative.complexity.sparse = mean_complexity_sparse / k, 
              relative.complxity.dense = mean_complexity_dense / k,
              relative.complexity.uniform = mean_complexity_uniform / k,
              no.hypotheses = length(gamma_list)))
}

# The function tests if we have overlap in the posterior distributions,
# and with that if we need a(nother) bridge hypothesis.
is_overlap <- function(ordered_list) {
  
  for (i in 1: (length(ordered_list) - 1)) {
    #check all pairs of hypotheses
    this_el <- ordered_list[[i]]
    next_el <- ordered_list[[i + 1]]
    
    overlap <- which(this_el$tab != 0 & next_el$tab != 0)
    
    if (length(overlap) == 0) {
      before_position <- i
      if (this_el$alpha == next_el$alpha) {
        alpha <- this_el$alpha
        beta <- this_el$beta / 2
      }
      else if (this_el$beta == next_el$beta) {
        alpha <- next_el$alpha / 2
        beta <- this_el$beta
      }
      else {
        alpha <- this_el$alpha
        beta <- this_el$beta / 2
      }
      return(list(before_pos = before_position,
                  alpha = alpha, beta = beta))
    }
  }
  return(1)
}

# Given a list of the results for all needed hypotheses, the function
# computes the log BF of sparse against dense.
# @args ordered list of the results, with the outer two the hypotheses of 
# interest, and in between the bridge hypotheses. k the number of potential edges.
compute_bayes_factor <- function(ordered_list, k) {
  bf <- 0
  c <- 0: k
  for (i in 1: (length(ordered_list) - 1)) {
    el1 <- ordered_list[[i]]
    el2 <- ordered_list[[i + 1]]
    
    alpha1 <- el1$alpha
    alpha2 <- el2$alpha
    beta1 <- el1$beta
    beta2 <- el2$beta
    tab1 <- el1$tab
    tab2 <- el2$tab
    
    log_prior1 <- lchoose(k, c) - lbeta(alpha1, beta1) + lfactorial(alpha1 + c - 1) +
      lfactorial(beta1 + k - c - 1) - lfactorial(alpha1 + beta1 + k - 1)
    log_prior2 <- lchoose(k, c) - lbeta(alpha2, beta2) + lfactorial(alpha2 + c - 1) +
      lfactorial(beta2 + k - c - 1) - lfactorial(alpha2 + beta2 + k - 1)
    
    prob1 <- tab1 / sum(tab1)
    prob2 <- tab2 / sum(tab2)
    
    log_prob1 <- log(prob1)
    log_prob2 <- log(prob2)
    
    odds1 <- log_prior1 - log_prob1
    odds2 <- log_prob2 - log_prior2

    log_bf <- odds1 + odds2
    log_bf[is.infinite(log_bf)] <- NA
    log_bf <- mean(log_bf, na.rm = TRUE)
    
    bf <- bf + log_bf
        
  }
  return(bf)
}