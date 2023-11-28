# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------

bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not_cont, centrality, progress, ...){

  if(!save && centrality){
    save <- TRUE
  }

  if(packageVersion("bgms") > "0.1.0"){
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
      bgm, c(list(x = data, iter = iter, save = TRUE, display_progress = progress),
             args)
    )
  } else {
    prior_defaults <- list(
      interaction_prior = "UnitInfo",
      cauchy_scale = 2.5,
      threshold_alpha = 1,
      threshold_beta = 1
    )
    args <- set_defaults(prior_defaults, ...)
    bgms_fit <- do.call(
      bgm, c(list(x = data, iter = iter, save = TRUE, display_progress = progress),
             args)
    )
  }

  fit$model <- type
  fit$packagefit <- bgms_fit
  if(is.null(colnames(data))){
    fit$var_names <- paste0("V", 1:ncol(data))
  } else {
    fit$var_names <-colnames(data)
  }
  class(fit) <- c("package_bgms", "easybgm")
  return(fit)
}




# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
bgm_extract.package_bgms <- function(fit, type, save,
                                     not_cont, data, centrality, ...){
  if(centrality) save <- TRUE

  if(any(class(fit) != "bgms")){
  fit <- fit$packagefit
  save <- TRUE
  }

  if(packageVersion("bgms") > "0.1.0"){
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
  } else {
    edge.prior <- 0.5
  }
  bgms_res <- list()
  if(save){
    p <- length(fit$colnames)
    bgms_res$parameters <- vector2matrix(colMeans(fit$interactions), p = p)
    bgms_res$thresholds <- as.matrix(colMeans(fit$thresholds))

    if(!is.null(data)){
    colnames(bgms_res$parameters) <- rownames(bgms_res$parameters) <- colnames(data)
    } else {
      colnames(bgms_res$parameters) <- rownames(bgms_res$parameters) <- fit$colnames
    }
    bgms_res$inc_probs <- vector2matrix(colMeans(fit$gamma), p = p)
    bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1-edge.prior))
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)

    #Obtain structure information
    structures <- apply(fit$gamma, 1, paste0, collapse="")
    table_structures <- as.data.frame(table(structures))
    bgms_res$structure_probabilities <- table_structures[,2]/nrow(fit$gamma)
    bgms_res$graph_weights <- table_structures[,2]
    bgms_res$sample_graph <- as.character(table_structures[, 1])
  } else {
    bgms_res$parameters <- fit$interactions
    bgms_res$thresholds <- fit$thresholds
    bgms_res$inc_probs <- fit$gamma
    bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1-edge.prior))
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)
  }
  if(save){
    bgms_res$samples_posterior <- fit$interactions
    if(centrality){
      #bgms_res$centrality_strength <- centrality_strength(bgms_res)
      bgms_res$centrality <- centrality(bgms_res)
    }
  }
  # Adapt column names of output
  colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
  colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters)

  bgms_res$model <- type
  output <- bgms_res
  class(output) <- c("package_bgms", "easybgm")
  return(output)
}
