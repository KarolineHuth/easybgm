# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms_compare <- function(fit, type, data, iter, save,
                                         not_cont, progress, ...){
  if(packageVersion("bgms") < "0.1.4"){
    stop("The package requires at least 'bgms' version of 0.1.4.",
         call. = FALSE)
  }
  if(!save){
    save <- TRUE
  }
  
  if(type == "binary") {
    type <- "ordinal"
  }
  if(packageVersion("bgms") < "0.1.6.0"){
    bgms_fit <- do.call(
      bgmCompare, c(list(data[[1]], data[[2]], iter = iter, save = T, 
                         variable_type = type, 
                         display_progress = progress, 
                         ...))
    )
  } 
  if(packageVersion("bgms") > "0.1.4.2"){
    bgms_fit <- do.call(
      bgmCompare, c(list(data[[1]], data[[2]], iter = iter, 
                         variable_type = type, 
                         display_progress = progress, 
                         ...))
    )
  }
  
  
  fit$model <- type
  fit$packagefit <- bgms_fit
  if(is.null(colnames(data[[1]]))){
    fit$var_names <- paste0("V", 1:ncol(data[[1]]))
  } else {
    fit$var_names <- colnames(data[[1]])
  }
  class(fit) <- c("package_bgms_compare", "easybgm")
  return(fit)
}




# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_extract.package_bgms_compare <- function(fit, type, save,
                                             not_cont, data, centrality, ...){
  if(any(class(fit) == "easybgm")){
    varnames <- fit$var_names
    fit <- fit$packagefit
  } else if (any(class(fit) == "bgms")){
    varnames <- extract_arguments(fit)$data_columnnames
    if(is.null(varnames)){
      varnames <- paste0("V", 1:extract_arguments(fit)$no_variables)
    }
  }
  
  # To ensure reverse-compatability
  if(packageVersion("bgms") < "0.1.6.0"){
    args <- fit$arguments
    
    if (args$pairwise_difference_prior[1] == "Bernoulli") {
      edge.prior <- args$inclusion_probability_difference
    } else { # if BB or SBM
      edge.prior <- calculate_edge_prior(alpha = args$pairwise_beta_bernoulli_alpha,
                                         beta = args$pairwise_beta_bernoulli_beta)
      
      # otherwise it saves the wrong values (could be done more elegantly)
      args$inclusion_probability_difference <- edge.prior
    }
    
    bgms_res <- list()
    
    
    p <- args$no_variables
    pars <- fit$pairwise_difference
    bgms_res$parameters <- vector2matrix(colMeans(pars), p = p)
    colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters), 
                                 nrow = nrow(bgms_res$parameters))
    bgms_res$inc_probs <- vector2matrix(colMeans(fit$pairwise_difference_indicator), p = p) 
    bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1 - edge.prior))
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)
    #Obtain structure information
    
    bgms_res$parameters_g1 <- vector2matrix(colMeans(fit$interactions), p = p) + bgms_res$parameters/2
    bgms_res$parameters_g2 <- vector2matrix(colMeans(fit$interactions), p = p) - bgms_res$parameters/2
    
    if(save){
      structures <- apply(fit$pairwise_difference_indicator, 1, paste0, collapse="")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[,2]/nrow(fit$pairwise_difference_indicator)
      bgms_res$graph_weights <- table_structures[,2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
      bgms_res$samples_posterior <- fit$pairwise_difference
    } 
  }
  if(packageVersion("bgms") > "0.1.4.2"){
    
    class(fit) <- c("bgmCompare", "bgms")
    
    args <- extract_arguments(fit)
    args$save <- TRUE
    dots <- list(...)
    if (args$difference_prior[1] == "Bernoulli") {
      if("difference_probability" %in% dots){
        edge.prior <- difference_prior
        args$inclusion_probability_difference <- edge.prior 
      } else {
        edge.prior <- 0.5
        args$inclusion_probability_difference <- edge.prior 
      }
    } else { # if BB or SBM
      edge.prior <- calculate_edge_prior(alpha = args$difference_selection_alpha,
                                         beta = args$difference_selection_beta)
      
      # otherwise it saves the wrong values (could be done more elegantly)
      args$inclusion_probability_difference <- edge.prior
    }
    
    bgms_res <- list()
    
    
    p <- args$num_variables
    pars <- fit$posterior_summary_pairwise_differences$mean
    bgms_res$parameters <- vector2matrix(pars, p = p)
    colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters), 
                                 nrow = nrow(bgms_res$parameters))
    indicators <- extract_indicators(fit)
    bgms_res$inc_probs <- vector2matrix(colMeans(indicators[, grep("\\(pairwise\\)", colnames(indicators))]), p = p) 
    bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1 - edge.prior))
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)
    
    #Obtain structure information
    bgms_res$parameters_g1 <- extract_group_params(fit_bgms_compare)$pairwise_effects_groups[, 1]
    bgms_res$parameters_g2 <- extract_group_params(fit_bgms_compare)$pairwise_effects_groups[, 2]
    
    if(save){
      structures <- apply(extract_indicators(fit)$pairwise_difference_indicator, 1, paste0, collapse="")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[,2]/nrow(extract_indicators(fit)$pairwise_difference_indicator)
      bgms_res$graph_weights <- table_structures[,2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
      bgms_res$samples_posterior <- extract_pairwise_difference(fit)
    } 
  }
  
  # Adapt column names of output
  colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
  colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters) 
  
  bgms_res$model <- type
  bgms_res$fit_arguments <- args
  bgms_res$edge.prior <- edge.prior[1, 1] # otherwise it stores a whole matrix 
  #bgms_res$n <- c(args$no_cases_gr1, args$no_cases_gr2)
  
  output <- bgms_res
  class(output) <- c("package_bgms_compare", "easybgm")
  return(output)
}