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
      #Obtain structure information
      if(packageVersion("bgms") < "0.1.4"){
        gammas <- bgms::extract_edge_indicators(fit)
      } else {
        gammas <- extract_indicators(fit)}
      structures <- apply(gammas, 1, paste0, collapse="")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[,2]/nrow(gammas)
      bgms_res$graph_weights <- table_structures[,2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
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