# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bggm_compare <- function(fit, type, data, iter, save,
                                 not_cont, progress, ...){
  
  prior_defaults <- list(
    prior_sd = .25
  )
  
  args <- set_defaults(prior_defaults, ...)
  # Fit the model
  bggm_fit <- do.call(
    BGGM::ggm_compare_explore, c(list(data[[1]], data[[2]], type = type, mixed_type = not_cont,
                          iter = iter, progress = progress, seed = NULL), args)
  )
  fit$model <- type
  fit$packagefit <- bggm_fit
  fit$fit_arguments  <- args
  if(is.null(colnames(data[[1]]))){
    fit$var_names <- paste0("V", 1:ncol(data[[1]]))
  } else {
    fit$var_names <- colnames(data[[1]])
  }
  class(fit) <- c("package_bggm_compare", "easybgm")
  return(fit)
}



# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_extract.package_bggm_compare <- function(fit, type, save,
                                     not_cont, data, ...){
  bggm_res <- list()
  varnames <- fit$var_names
  fit <- fit$packagefit
  bggm_res$fit_arguments  <- fit$args
  
  bggm_res$n <- fit$info$dat_info$n
  bggm_res$p <- fit$info$dat_info$p[1]
  bggm_res$parameters <- fit$pcor_diff
  diag(bggm_res$parameters) <- 0
  colnames(bggm_res$parameters) <- varnames
  bggm_res$inc_BF <- 1/fit$BF_01
  diag(bggm_res$inc_BF) <- 0
  bggm_res$inc_probs <- bggm_res$inc_BF/(bggm_res$inc_BF + 1)
  diag(bggm_res$inc_probs) <- 0
  bggm_res$structure <- (1*bggm_res$inc_probs > 0.5)
  
  colnames(bggm_res$inc_probs) <- colnames(bggm_res$parameters)
  colnames(bggm_res$inc_BF) <- colnames(bggm_res$parameters)
  
  if(save){
    samples_group1 <- fit$samp[[1]]$post_samp$pcors
    samples_group2 <- fit$samp[[2]]$post_samp$pcors
    samples_dif <- matrix(0, ncol = bggm_res$p*(bggm_res$p -1)/2, nrow = fit$iter)
    for(i in 1:fit$iter){
      sample1 <- samples_group1[, , i]
      sample1_vec <- as.vector(sample1[upper.tri(sample1)])
      sample2 <- samples_group2[, , i]
      sample2_vec <- as.vector(sample2[upper.tri(sample2)])
      samples_dif[i, ] <- sample1_vec - sample2_vec
    }
    bggm_res$samples_posterior <- samples_dif
  }
  
  bggm_res$parameters_g1 <- fit$samp[[1]]$pcor_mat
  bggm_res$parameters_g2 <- fit$samp[[2]]$pcor_mat
  
  if(fit$type == "continuous"){
    bggm_res$model <- "ggm"
  } else {
    bggm_res$model <- "gcgm"
  }
  output <- bggm_res
  class(output) <- c("package_bggm", "easybgm")
  return(output)
}
