# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------

bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not.cont, centrality, progress, ...){

  if(save == FALSE & centrality == TRUE){
    warning("The centrality measures can only be obtained if the posterior samples are saved. Note that we automatically set
            set save = TRUE.")
    save <- TRUE
  }
  
  prior_defaults <- list(
    interaction_prior = "UnitInfo",
    cauchy_scale = 2.5,
    threshold_alpha = 1,
    threshold_beta = 1
  )
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
bgm_extract.package_bgms <- function(fit, model, edge.prior, save,
                                     not.cont, data, centrality, ...){
  fit <- fit$packagefit
  bgms_res <- list()
  if(save == FALSE){

    bgms_res$parameters <- fit$interactions
    colnames(bgms_res$parameters) <- rownames(bgms_res$parameters) <- colnames(data)
    bgms_res$inc_probs <- fit$gamma
    bgms_res$BF <- fit$gamma/(1-fit$gamma)
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)
  }
  if(save == TRUE){
    p <- unlist(strsplit(colnames(fit$interactions)[ncol(fit$interactions)], ", "))[2]
    p <- as.numeric(unlist(strsplit(p, ")"))[1])
    bgms_res$parameters <- vector2matrix(colMeans(fit$interactions), p = p)
    colnames(bgms_res$parameters) <- rownames(bgms_res$parameters) <- colnames(data)
    bgms_res$inc_probs <- vector2matrix(colMeans(fit$gamma), p = p)
    bgms_res$BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1-edge.prior))
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)

    #Obtain structure information
    bgms_res$posterior_complexity <- table(rowSums(fit$gamma))/nrow(fit$gamma)
    structures <- apply(fit$gamma, 1, paste0, collapse="")
    table_structures <- as.data.frame(table(structures))
    bgms_res$structure_probabilities <- table_structures[,2]/nrow(fit$gamma)
    bgms_res$graph_weights <- table_structures[,2]
    bgms_res$sample_graphs <- as.character(table_structures[, 1])
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
