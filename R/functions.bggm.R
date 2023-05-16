# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------

bgm_fit.package_bggm <- function(fit, type, data, iter, save,
                                 not.cont, centrality, progress, ...){

  # Fit the model
  bggm_fit <- BGGM::explore(data,                        #(M) n*p matrix of responses
                            type = type,                 #(O) type of data
                            mixed_type = not.cont,       #(O) which data should be treated as ranks
                            iter = iter,             #(O) no. iterations sampler
                            progress = progress,            #(O) Should a progress bar be plotted?
                            seed = NULL,                     #(O) Integer for random seed
                            ...)
  fit$model <- type
  fit$packagefit <- bggm_fit
  class(fit) <- c("package_bggm", "easybgm")
  return(fit)
}



# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------

bgm_extract.package_bggm <- function(fit, model, edge.prior, save,
                                     not.cont, data, centrality, ...){
  fit <- fit$packagefit

  out_select <- BGGM::select(fit)
  bggm_res <- list()
  bggm_res$parameters <- out_select$pcor_mat
  colnames(bggm_res$parameters) <- colnames(fit$Y)
  bggm_res$BF <- out_select$BF_10
  bggm_res$inc_probs <- out_select$BF_10/(out_select$BF_10 + 1)
  bggm_res$structure <- out_select$Adj_10

  if(centrality == TRUE){
    save <- TRUE
  }
  if(save == TRUE){
    p <- ncol(bggm_res$parameters)
    samples <- matrix(0, ncol = p*(p-1)/2, nrow = fit$iter)
    for(i in 1:fit$iter){
      sample <- fit$post_samp$pcors[, , i]
      samples[i, ] <- as.vector(sample[upper.tri(sample)])
    }
    bggm_res$samples_posterior <- samples

    if(centrality == TRUE){
      # bggm_res$centrality_strength <- centrality_strength(bggm_res)
      bggm_res$centrality <- centrality(bggm_res)
    }
  }

  if(fit$type == "continuous"){
    bggm_res$model <- "ggm"
  } else {
    bggm_res$model <- "gcgm"
  }

  output <- bggm_res
  print("hi")
  return(output)
}

