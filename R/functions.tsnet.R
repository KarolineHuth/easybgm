# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_tsnet <- function(fit, type, data, iter,
                                  progress, ...){
  
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
  
  # Summary 
  tsnet_res$parameters_temporal <- post_samps$beta_mu
  colnames(tsnet_res$parameters_temporal) <- rownames(tsnet_res$parameters_temporal) <- varnames
  tsnet_res$parameters_contemporaneous <- post_samps$pcor_mu
  colnames(tsnet_res$parameters_contemporaneous) <- rownames(tsnet_res$parameters_contemporaneous) <- varnames
  
  # Raw samples
  tsnet_res$posterior_samples_temporal <- list2matrix(post_samps$fit$beta, p = tsnet_res$p, part = "all")
  tsnet_res$posterior_samples_contemporaneous <- list2matrix(post_samps$fit$pcor, p = tsnet_res$p)
  
  
  tsnet_res$model <- "GVAR"
  
  output <- tsnet_res
  class(output) <- c("package_tsnet", "easybgm_longitudinal", "easybgm")
  return(output)
}
