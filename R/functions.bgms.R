# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------

bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not.cont, centrality, progress, ...){

  if(save == FALSE & centrality == TRUE){
    save <- TRUE
  }
  
  prior_defaults <- list(
    interaction_prior = "UnitInfo",
    cauchy_scale = 2.5,
    threshold_alpha = 1,
    threshold_beta = 1,
    edge_prior = "Bernoulli"
  )
  prior_defaults <- set_defaults(prior_defaults, ...)
  if (prior_defaults$edge_prior == "Bernoulli") {
    extra_args <- list(
      inclusion_probability = .5
    )
  }
  else {
    extra_args <- list(
      beta_bernoulli_alpha = 1,
      beta_bernoulli_beta = 1
    )
  }
  prior_defaults <- append(prior_defaults, extra_args)


  args <- set_defaults(prior_defaults, ...)
  bgms_fit <- do.call(
    bgm, c(list(x = data, iter = iter, save = save, display_progress = progress),
           args)
  )

  fit$model <- type
  fit$packagefit <- bgms_fit
  class(fit) <- c("package_bgms", "easybgm")
  return(fit)
}




# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
bgm_extract.package_bgms <- function(fit, model, save,
                                     not.cont, data, centrality, ...){
  if(centrality == TRUE) save <- TRUE

  fit <- fit$packagefit
  
  defaults <- list(
    edge_prior = "Bernoulli"
  )
  
  args <- set_defaults(defaults, ...)
  if (args$edge_prior == "Bernoulli") {
    extra_defaults <- list(
      inclusion_probability = .5
    )
  }
  
  else {
    extra_defaults <- list(
      beta_bernoulli_alpha = 1,
      beta_bernoulli_beta = 1
    )
  }
  
  defaults <- append(defaults, extra_defaults)
  args <- set_defaults(defaults, ...)
  
  if (args$edge_prior == "Bernoulli") {
    edge.prior <- args$inclusion_probability
  }
  else {
    edge.prior <- beta(args$beta_bernoulli_alpha, args$beta_bernoulli_beta)
  }
  bgms_res <- list()
  p <- unlist(strsplit(colnames(fit$interactions)[ncol(fit$interactions)], ", "))[2]
  p <- as.numeric(unlist(strsplit(p, ")"))[1])
  bgms_res$parameters <- vector2matrix(colMeans(fit$interactions), p = p)
  colnames(bgms_res$parameters) <- rownames(bgms_res$parameters) <- colnames(data)
  bgms_res$inc_probs <- vector2matrix(colMeans(fit$gamma), p = p)
  bgms_res$BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1-edge.prior))
  bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)

  #Obtain structure information
  structures <- apply(fit$gamma, 1, paste0, collapse="")
  table_structures <- as.data.frame(table(structures))
  bgms_res$structure_probabilities <- table_structures[,2]/nrow(fit$gamma)
  bgms_res$graph_weights <- table_structures[,2]
  bgms_res$sample_graph <- as.character(table_structures[, 1])
  if(save == TRUE){
    bgms_res$samples_posterior <- fit$interactions
    if(centrality == TRUE){
      #bgms_res$centrality_strength <- centrality_strength(bgms_res)
      bgms_res$centrality <- centrality(bgms_res)
    }
  }
  bgms_res$model <- model
  output <- bgms_res
  return(output)
}
