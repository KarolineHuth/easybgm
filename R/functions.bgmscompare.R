# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms_compare <- function(fit, type, data, group_indicator, iter, save,
                                         not_cont, progress, ...){
  if(packageVersion("bgms") < "0.1.4"){
    stop("The package requires at least 'bgms' version of 0.1.4.",
         call. = FALSE)
  }
  
  if(type == "binary") {
    type <- "ordinal"
  }
  if(is.list(data)){
    group_indicator <- NULL
  }
  if(packageVersion("bgms") < "0.1.6.0"){
    bgms_fit <- do.call(
      bgmCompare, c(list(data[[1]], data[[2]], iter = iter, save = TRUE, 
                         variable_type = type, 
                         display_progress = progress, 
                         ...))
    )
    
    fit$model <- type
    fit$packagefit <- bgms_fit
    if(is.null(colnames(data[[1]]))){
      fit$var_names <- paste0("V", 1:ncol(data[[1]]))
    } else {
      fit$var_names <- colnames(data[[1]])
    }
  } 
  if(packageVersion("bgms") > "0.1.4.2"){
    if(is.null(group_indicator)){
      bgms_fit <- do.call(
        bgmCompare, c(list(x = data[[1]], y = data[[2]], iter = iter, 
                           variable_type = type, 
                           display_progress = progress, 
                           ...))
      )
      
      fit$model <- type
      fit$packagefit <- bgms_fit
      if(is.null(colnames(data[[1]]))){
        fit$var_names <- paste0("V", 1:ncol(data[[1]]))
      } else {
        fit$var_names <- colnames(data[[1]])
      }
    } 
    if(!is.null(group_indicator)){
      bgms_fit <- do.call(
        bgmCompare, c(list(x = data, group_indicator = group_indicator, iter = iter, 
                           variable_type = type, 
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
    }
    
  }
  
  

  class(fit) <- c("package_bgms_compare", "easybgm")
  return(fit)
}




# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_extract.package_bgms_compare <- function(fit, type, save, group_indicator,
                                             not_cont, data, centrality, ...){
  if(any(class(fit) == "easybgm")){
    varnames <- fit$var_names
    fit <- fit$packagefit
  } else if (any(class(fit) == "bgmCompare")){
    varnames <- extract_arguments(fit)$data_columnnames
    if(is.null(varnames)){
      varnames <- paste0("V", 1:extract_arguments(fit)$no_variables)
    }
  }
  
  ######--------------------
  ## Old package two group estimation
  ######--------------------
  
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
    
    bgms_res$parameters_g1 <- vector2matrix(colMeans(fit$interactions), p = p) + bgms_res$parameters/2
    bgms_res$parameters_g2 <- vector2matrix(colMeans(fit$interactions), p = p) - bgms_res$parameters/2
    
    structures <- apply(fit$pairwise_difference_indicator, 1, paste0, collapse="")
    table_structures <- as.data.frame(table(structures))
    bgms_res$structure_probabilities <- table_structures[,2]/nrow(fit$pairwise_difference_indicator)
    bgms_res$graph_weights <- table_structures[,2]
    bgms_res$sample_graph <- as.character(table_structures[, 1])
    bgms_res$samples_posterior <- fit$pairwise_difference
  }
  
  ######--------------------
  ## Newer package 
  ######--------------------
  if(packageVersion("bgms") > "0.1.4.2"){
    ######--------------------
    ## Two group estimation
    ######--------------------
    if(is.null(group_indicator)){
      class(fit) <- c("bgmCompare")
      
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
      bgms_res$parameters_g1 <- vector2matrix(extract_group_params(fit)$pairwise_effects_groups[, 1], p = p)
      bgms_res$parameters_g2 <- vector2matrix(extract_group_params(fit)$pairwise_effects_groups[, 2], p = p)
      
      structures <- apply(extract_indicators(fit), 1, paste0, collapse="")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[,2]/nrow(extract_indicators(fit))
      bgms_res$graph_weights <- table_structures[,2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
      bgms_res$samples_posterior <- extract_pairwise_interactions(fit)
      
      bgms_res$convergence_parameter <-  fit$posterior_summary_pairwise_differences$Rhat
    }
    
    ######--------------------
    ## Multi-group estimation
    ######--------------------
    if(!is.null(group_indicator)){
      class(fit) <- c("bgmCompare")
      
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
      #Compute average group difference
      diffs <- fit$posterior_summary_pairwise_differences
      edges <- unique(sub(" .*", "", diffs$parameter))
      average_difference <- numeric(length(edges))
      for (i in seq_along(edges)) {
        edge_name <- edges[i]
        idx <- edges == edge_name
        vals <- diffs$mean[idx]
        average_difference[i] <- mean(vals, na.rm = TRUE)
      }
      bgms_res$parameters <- vector2matrix(average_difference, p = p)
      colnames(bgms_res$parameters) <- varnames
      bgms_res$overall_estimate <- vector2matrix(fit$posterior_summary_pairwise_baseline$mean, p = p) #overall group estimate
      bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters), 
                                   nrow = nrow(bgms_res$parameters))
      indicators <- extract_indicators(fit)
      bgms_res$inc_probs <- vector2matrix(colMeans(indicators[, grep("\\(pairwise\\)", colnames(indicators))]), p = p) 
      bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1 - edge.prior))
      bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)
      
      structures <- apply(extract_indicators(fit), 1, paste0, collapse="")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[,2]/nrow(extract_indicators(fit))
      bgms_res$graph_weights <- table_structures[,2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
      bgms_res$samples_posterior <- extract_pairwise_interactions(fit)
      bgms_res$convergence_parameter <-  fit$posterior_summary_pairwise_baseline$Rhat
      bgms_res$group_estimates <- extract_group_params(fit)$pairwise_effects_groups
      

    }
  }

  
  # Adapt column names of output
  colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
  colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters) 
  
  bgms_res$model <- type
  bgms_res$fit_arguments <- args
  bgms_res$edge.prior <- edge.prior # otherwise it stores a whole matrix 
  #bgms_res$n <- c(args$no_cases_gr1, args$no_cases_gr2)
  
  output <- bgms_res
  class(output) <- c("package_bgms_compare", "easybgm_compare", "easybgm")
  return(output)
}