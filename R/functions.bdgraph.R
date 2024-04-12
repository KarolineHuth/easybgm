# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------

bgm_fit.package_bdgraph <- function(fit, type, data, iter, save,
                                    not_cont, centrality, progress, ...){
  
  if(type == "binary"){
    prior_defaults <- list(
      g.prior = .5
    )
    
    args <- set_defaults(prior_defaults, ...)
  } else {
    prior_defaults <- list(
      g.prior = .5,
      df.prior = 3,
      algorithm = "rjmcmc"
    )
    
    args <- set_defaults(prior_defaults, ...)
  }
  
  
  if(type == "continuous"){
    
    bdgraph_fit <- do.call(BDgraph::bdgraph,
                           c(list(data = data, method ="ggm", iter = iter,
                                  save = TRUE), args))
    
    
    fit$model <- "ggm"
  }
  if(type %in% c("mixed", "ordinal")){
    if(type == "ordinal") not_cont <- rep(1, ncol(data))
    # fitting the model
    
    bdgraph_fit <- do.call(BDgraph::bdgraph,
                           c(list(data = data, method = "gcgm",
                                  iter = iter, save = TRUE, not.cont = not_cont
                           ), args))
    
    fit$model <- "gcgm"
    
  }
  if(type == "binary"){
    bdgraph_fit <- do.call(bdgraph.mpl,
                           c(list(data = data, method = "dgm-binary",
                                  iter = iter, save = TRUE), args))
    
    
    fit$model <-  "dgm-binary"
  }
  fit$packagefit <- bdgraph_fit
  if(is.null(colnames(data))){
    fit$var_names <- paste0("V", 1:ncol(data))
  } else {

    fit$var_names <- colnames(data)
  }
  class(fit) <- c("package_bdgraph", "easybgm")
  return(fit)
}


# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------

bgm_extract.package_bdgraph <- function(fit, type, save,
                                        not_cont, data, centrality, 
                                        posterior_method, ...){
  model <- fit$model
  if(is.null(model)){
    stop("Please specify the type of model estimated with BDgraph (e.g., ggm, gcgm, dgm-binary).",
         call. = FALSE)
  }
  varnames <- fit$var_names
  fit <- fit$packagefit
  
  defaults <- list(
    g.prior = .5
  )
  
  args <- set_defaults(defaults, ...)
  edge.prior = args$g.prior
  bdgraph_res <- list()
  if(model %in% "ggm"){
    #Bayesian model-averaged estimates
    bdgraph_res$parameters <- qgraph::wi2net(fit$K_hat)
    diag(bdgraph_res$parameters) <- 0
    colnames(bdgraph_res$parameters) <- varnames
    bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
    bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
    bdgraph_res$inc_BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior /(1-edge.prior))
    bdgraph_res$structure <- 1*(bdgraph_res$inc_probs > 0.5)
    bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
    bdgraph_res$graph_weights <- fit$graph_weights
    bdgraph_res$sample_graphs <- fit$sample_graphs
    # bdgraph_res$package <- "bdgraph"
    bdgraph_res$model <- "ggm"

    if(centrality == TRUE){
      save <- TRUE
    }
    
    if(save == TRUE){
      if(posterior_method == "MAP"){
        
        warning("Posterior samples of the BDgraph package are obtained after the fit and might show small deviations from the model fit.")
      }
          # Extract posterior samples
      data<-as.matrix(data)

      bdgraph_res$samples_posterior <- extract_posterior(fit, data=data, method = model, not_cont, posterior_method = posterior_method)[[1]]
      
      if(centrality == TRUE){
        # Centrality indices
        #bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
        bdgraph_res$centrality <- centrality(bdgraph_res)
      }
      
    }
    
    output <- bdgraph_res
  }
  
  if(model %in% c("gcgm")){
    #Bayesian model-averaged estimates
    bdgraph_res$parameters <- qgraph::wi2net(fit$K_hat)
    diag(bdgraph_res$parameters) <- 0
    colnames(bdgraph_res$parameters) <- varnames
    bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
    bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
    bdgraph_res$inc_BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior/(1-edge.prior))
    bdgraph_res$structure <- 1*(bdgraph_res$inc_probs > 0.5)
    bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
    bdgraph_res$graph_weights <- fit$graph_weights
    bdgraph_res$sample_graphs <- fit$sample_graphs
    # bdgraph_res$package <- "bdgraph"
    bdgraph_res$model <- "gcgm"

    if(centrality){
      save <- TRUE
    }
    if(save){
      if(posterior_method == "MAP"){
        
        warning("Posterior samples of the BDgraph package are obtained after the fit. 
                The samples were obtained from the structure with the 
                highest posterior probability, while BDgraph itself outputs model-averaged results. 
                To obtain model-averaged posterior samples set 'posterior_method = model-averaged'. 
                Note that the computation time will be considerable increased.")
      }
      
      if(is.null(not_cont)){
        stop("Specify a vector indicating variables that are continuous with the not_cont argument (1 indicates not continuous)",
             call. = FALSE)
      }
      
      data <- as.matrix(data)
      
      # Extract posterior samples

      bdgraph_res$samples_posterior <- extract_posterior(fit, data, method = model, not_cont = not_cont, posterior_method = posterior_method)[[1]]
      
      if(centrality){
        # Centrality indices
        # bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
        bdgraph_res$centrality <- centrality(bdgraph_res)
      }
    }
    colnames(bdgraph_res$inc_probs) <- colnames(bdgraph_res$parameters)
    colnames(bdgraph_res$inc_BF) <- colnames(bdgraph_res$parameters)
    
    output <- bdgraph_res
  }
  if(model %in% c("dgm-binary")){
    #Bayesian model-averaged estimates
    bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
    bdgraph_res$inc_BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior/(1-edge.prior))
    colnames(bdgraph_res$inc_probs) <- rownames(bdgraph_res$inc_probs) <- varnames
    colnames(bdgraph_res$inc_BF) <- rownames(bdgraph_res$inc_BF) <- varnames
    bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
    bdgraph_res$graph_weights <- fit$graph_weights
    bdgraph_res$sample_graphs <- fit$sample_graphs
    bdgraph_res$model <- "dgm-binary"

    if(save == TRUE){
      warning("Posterior samples cannot be obtained for \"dgm-binary\". Solely the aggregated results are extracted.",
              call. = FALSE)
    }
    output <- bdgraph_res
  }
  return(output)
}

