#' @name sparse_or_dense
#' @title Test for sparse against dense topologies
#'
#' @description The function tests if a network is sparse (i.e., few edges in between nodes) or 
#' dense (i.e., a lot of edges between nodes). It estimates
#' the network model of a given data set under the hypothesis
#' that it is sparse and that it is dense, and computes th Bayes factor.
#'
#' @param x An n x p matrix or dataframe containing the variables for n independent observations on p variables.
#' @param type What is the data type? Options: currently only binary and ordinal are supported
#' @param ... additional arguments of the bgms function
#'
#' @return List containing results of the analysis:
#' 
#' \itemize{
#'
#' \item \code{log.BF} The log of the Bayes factor of the test of a sparse against a dense network structure. A value larger 0 indicates evidence for a sparse structure. Contrarily, a value smaller 0 indicates evidence for a dense structure. 
#'
#' \item \code{BF} The Bayes factor of the test of a sparse against a dense network structure. A value larger 1 indicates evidence for a sparse structure. Contrarily, a value between 0 and 1 indicates evidence for a dense structure. 
#'
#' \item \code{relative.complexity.sparse} The relative complexity under a sparse prior hypothesis, which is the 
#' proportion of estimated included edges relative to the total possible edges
#' under the different hypotheses.
#' 
#' \item \code{relative.complexity.dense} The relative complexity under a dense prior hypothesis.
#' 
#' \item \code{relative.complexity.uniform} The relative complexity under a uniform prior hypothesis.
#' 
#'  \item \code{no.hypotheses} The number of hypotheses that are computed in the analysis. At least the sparse,
#' dense and uniform are used, but sometimes additional bridge hypotheses have to be 
#' computed to be able to determine the Bayes factor. 
#' }
#'
#' @examples
#' \donttest{
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#'
#' # Fitting the Wenchuan PTSD data
#'
#' fit <- sparse_or_dense(data, type = "ordinal",
#'                 iter = 1000 # for demonstration only (> 5e4 recommended)
#'                 )
#' }  
#'              
#' @export
sparse_or_dense <- function(x, type, ...) {

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
  
  if(packageVersion("bgms") < "0.1.4"){
    gammas_sparse <- bgms::extract_edge_indicators(res_sparse)
    gammas_dense <- bgms::extract_edge_indicators(res_dense)
    gammas_uniform <- bgms::extract_edge_indicators(res_uniform)
  } else {
    gammas_sparse <- extract_indicators(res_sparse)
    gammas_dense <- extract_indicators(res_dense)
    gammas_uniform <- extract_indicators(res_uniform)
    }
  
  gamma_sums_sparse <- rowSums(gammas_sparse)
  gamma_sums_dense <- rowSums(gammas_dense)
  gamma_sums_uniform <- rowSums(gammas_uniform)
  
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
              relative.complexity.dense = mean_complexity_dense / k,
              relative.complexity.uniform = mean_complexity_uniform / k,
              no.hypotheses = length(gamma_list)))
}
