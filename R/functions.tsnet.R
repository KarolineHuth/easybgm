# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_tsnet <- function(fit, type, data, iter,
                                  progress, centrality, ...){
  
  # Fit the model
  tsnet_fit <- do.call(
    tsnet::stan_gvar, c(list(data,
                             iter_sampling = iter, ...))
  )
  
  fit$model <- type
  fit$packagefit <- tsnet_fit
  if(is.null(colnames(data))){
    fit$var_names <- paste0("V", 1:ncol(data))
  } else {
    fit$var_names <- colnames(data)
  }
  class(fit) <- c("package_tsnet", "easybgm_longitudinal", "easybgm")
  return(fit)
}



# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_extract.package_tsnet <- function(fit, type,
                                      centrality,
                                      data, ...){
  tsnet_res <- list()
  varnames <- fit$var_names
  fit <- fit$packagefit
  tsnet_res$fit_arguments  <- fit$arguments
  tsnet_res$n <- fit$arguments$n_t
  tsnet_res$p <- fit$arguments$p
  post_samps <- stan_fit_convert(fit, return_params = c("beta", "pcor"))
  
  #-------------- 
  # Summary 
  tsnet_res$parameters_temporal <- post_samps$beta_mu
  colnames(tsnet_res$parameters_temporal) <- rownames(tsnet_res$parameters_temporal) <- varnames
  tsnet_res$parameters_contemporaneous <- post_samps$pcor_mu
  colnames(tsnet_res$parameters_contemporaneous) <- rownames(tsnet_res$parameters_contemporaneous) <- varnames
  
  # Raw samples
  tsnet_res$posterior_samples_temporal <- list2matrix(post_samps$fit$beta, p = tsnet_res$p, part = "all")
  tsnet_res$posterior_samples_contemporaneous <- list2matrix(post_samps$fit$pcor, p = tsnet_res$p)
  
  #--------------
  # BF
  beta_post <- post_samps$fit$beta      
  pcor_post <- post_samps$fit$pcor    
  
  n_samples <- dim(beta_post)[3]
  n_row <- dim(beta_post)[1]
  n_col <- dim(beta_post)[2]
  prior <- rnorm(n_samples)
  
  # Beta estimates
  tsnet_res$bf_temporal <- matrix(NA, nrow = n_row, ncol = n_col)
  for(i in 1:n_row){
    for(j in 1:n_col){
      post_samples <- beta_post[i, j, ]
      tsnet_res$bf_temporal[i, j] <- as.numeric(bayestestR::bayesfactor_pointnull(
        posterior = post_samples,
        prior = prior,
        null = 0
      ))
    }
  }
  
  # Partial correlation
  tsnet_res$bf_contemporaneous <- matrix(0, nrow = n_row, ncol = n_col)
  upper_idx <- which(upper.tri(tsnet_res$bf_contemporaneous), arr.ind = TRUE)
  for(k in 1:n_row){
      i <- upper_idx[k, 1]
      j <- upper_idx[k, 2]
      post_samples <- pcor_post[i, j, ]
      tsnet_res$bf_contemporaneous[i, j] <- as.numeric(bayestestR::bayesfactor_pointnull(
        posterior = post_samples,
        prior = prior,
        null = 0
      ))
  }
  tsnet_res$bf_contemporaneous <- tsnet_res$bf_contemporaneous + t(tsnet_res$bf_contemporaneous)
  
  # ------------
  # Centrality 
  if (centrality){
    cent_samples <- get_centrality(fit)
    tsnet_res$centrality_samples <- cent_samples
    cent_samples[["instrength"]] <- colMeans(cent_samples[["instrength"]])
    cent_samples[["outstrength"]] <- colMeans(cent_samples[["outstrength"]])
    cent_samples[["strength"]] <- colMeans(cent_samples[["strength"]])
    cent_samples[["density_beta"]] <- mean(cent_samples[["density_beta"]])
    cent_samples[["density_pcor"]] <- mean(cent_samples[["density_pcor"]])
    tsnet_res$centrality <- cent_samples
  }
  
  tsnet_res$model <- "GVAR"
  
  output <- tsnet_res
  class(output) <- c("package_tsnet", "easybgm_longitudinal", "easybgm")
  return(output)
}
