# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------

bgm_fit.package_bdgraph <- function(fit, type, data, iter, save,
                            not.cont, centrality, progress, ...){
  if(type == "continuous"){
    bdgraph_fit <- BDgraph::bdgraph(data = data,               #(M) n*p matrix of responses
                                    method ="ggm",           #(M) type of data
                                    iter = iter,           #(O) no. iterations sampler
                                    save = TRUE,           #(O) Should samples be stored
                                    ...)

    fit$model <- "ggm"
  }
  if(type %in% c("mixed", "ordinal")){
    if(type == "ordinal") not.cont <- rep(1, ncol(data))
    # fitting the model
    bdgraph_fit <- BDgraph::bdgraph(data = data,         #(M) n*p matrix of responses
                                    method="gcgm",       #(M) type of data
                                    iter=iter,           #(O) no. iterations sampler
                                    save= TRUE,          #(O) Should samples be stored
                                    not.cont = not.cont, #(O) Specifies not continuous variables if method is gcgm
                                    ...)
    fit$model <- "gcgm"

  }
  if(type == "binary"){
    bdgraph_fit <- bdgraph.mpl(data = data,
                               method = "dgm-binary",
                               iter = iter,
                               save = TRUE,
                               ...)
    fit$model <-  "dgm-binary"
  }
  fit$packagefit <- bdgraph_fit

  class(fit) <- c("package_bdgraph", "easybgm")
  return(fit)
}


# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------

bgm_extract.package_bdgraph <- function(fit, model, edge.prior, save,
                                not.cont, data, centrality, ...){
  if(is.null(model)){
    stop("Please specify the type of model estimated with BDgraph (e.g., ggm, gcgm, dgm-binary).",
         call. = FALSE)
  }
  fit <- fit$packagefit
  bdgraph_res <- list()
  if(model %in% "ggm"){
    #Bayesian model-averaged estimates
    bdgraph_res$parameters <- pr2pc(fit$K_hat)
    diag(bdgraph_res$parameters) <- 0
    bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
    bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
    bdgraph_res$BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior /(1-edge.prior))
    bdgraph_res$structure <- 1*(bdgraph_res$inc_probs > 0.5)
    bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
    bdgraph_res$graph_weights <- fit$graph_weights
    bdgraph_res$sample_graph <- fit$sample_graphs
    # bdgraph_res$package <- "bdgraph"
    bdgraph_res$model <- "ggm"

    if(centrality == TRUE){
      save <- TRUE
    }

    if(save == TRUE){
      # Extract posterior samples
      data<-as.matrix(data)
      bdgraph_res$samples_posterior <- extract_posterior(fit, data=data, method = model, not.cont)[[1]]

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
    bdgraph_res$parameters <- pr2pc(fit$K_hat)
    diag(bdgraph_res$parameters) <- 0
    bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
    bdgraph_res$inc_probs  <- bdgraph_res$inc_probs + t(bdgraph_res$inc_probs)
    bdgraph_res$BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior/(1-edge.prior))
    bdgraph_res$structure <- 1*(bdgraph_res$inc_probs > 0.5)
    bdgraph_res$structure_probabilities <- fit$graph_weights/sum(fit$graph_weights)
    bdgraph_res$graph_weights <- fit$graph_weights
    bdgraph_res$sample_graph <- fit$sample_graphs
    # bdgraph_res$package <- "bdgraph"
    bdgraph_res$model <- "gcgm"

    if(centrality == TRUE){
      save <- TRUE
    }
    if(save == TRUE){
      #warning("Posterior samples cannot be extracted for mixed models with BDgraph at the moment. Results are provided without the posterior samples.")

      if(is.null(not.cont)){
        stop("Specify a vector indicating variables are continuos with the not.cont argument (1 indicates not continuous)",
             call. = FALSE)
      }

      data <- as.matrix(data)

      # Extract posterior samples
      bdgraph_res$samples_posterior <- extract_posterior(fit, data, method = model, not.cont = not.cont)[[1]]

      if(centrality == TRUE){
      # Centrality indices
      # bdgraph_res$centrality_strength <- centrality_strength(bdgraph_res)
      bdgraph_res$centrality <- centrality(bdgraph_res)
      }
    }
    output <- bdgraph_res
  }
  if(model %in% c("dgm-binary")){
    #Bayesian model-averaged estimates
    bdgraph_res$inc_probs <- as.matrix(BDgraph::plinks(fit))
    bdgraph_res$BF <- (bdgraph_res$inc_probs / (1 - bdgraph_res$inc_probs))/(edge.prior/(1-edge.prior))
    # bdgraph_res$package <- "bdgraph"
    bdgraph_res$model <- "dgm-binary"

    if(save == TRUE){
      warning("Posterior samples cannot be obtained for \"dgm-binary\". Solely the aggregated results are extracted.",
           call. = FALSE)
    }
    output <- bdgraph_res
  }
  return(output)
}

