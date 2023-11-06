# # --------------------------------------------------------------------------------------------------
# # 1. Fitting function
# # --------------------------------------------------------------------------------------------------
# 
# bgm_fit.package_bggm <- function(fit, type, data, iter, save,
#                                  not_cont, centrality, progress, ...){
# 
#   prior_defaults <- list(
#     prior_sd = .25
#   )
# 
#   args <- set_defaults(prior_defaults, ...)
#   # Fit the model
#   bggm_fit <- do.call(
#     BGGM::explore, c(list(Y = data, type = type, mixed_type = not_cont,
#                           iter = iter, progress = progress, seed = NULL), args)
#   )
#   fit$model <- type
#   fit$packagefit <- bggm_fit
#   class(fit) <- c("package_bggm", "easybgm")
#   return(fit)
# }
# 
# 
# 
# # --------------------------------------------------------------------------------------------------
# # 2. Extracting results function
# # --------------------------------------------------------------------------------------------------
# 
# bgm_extract.package_bggm <- function(fit, type, save,
#                                      not_cont, data, centrality, ...){
#   fit <- fit$packagefit
# 
#   out_select <- BGGM::select(fit)
#   bggm_res <- list()
#   bggm_res$parameters <- out_select$pcor_mat
#   colnames(bggm_res$parameters) <- colnames(fit$Y)
#   bggm_res$inc_BF <- out_select$BF_10
#   bggm_res$inc_probs <- out_select$BF_10/(out_select$BF_10 + 1)
#   bggm_res$structure <- out_select$Adj_10
# 
#   colnames(bggm_res$inc_probs) <- colnames(bggm_res$parameters)
#   colnames(bggm_res$inc_BF) <- colnames(bggm_res$parameters)
# 
#   if(centrality){
#     save <- TRUE
#   }
#   if(save){
#     p <- ncol(bggm_res$parameters)
#     samples <- matrix(0, ncol = p*(p-1)/2, nrow = fit$iter)
#     for(i in 1:fit$iter){
#       sample <- fit$post_samp$pcors[, , i]
#       samples[i, ] <- as.vector(sample[upper.tri(sample)])
#     }
#     bggm_res$samples_posterior <- samples
# 
#     if(centrality){
#       # bggm_res$centrality_strength <- centrality_strength(bggm_res)
#       bggm_res$centrality <- centrality(bggm_res)
#     }
#   }
# 
#   if(fit$type == "continuous"){
#     bggm_res$model <- "ggm"
#   } else {
#     bggm_res$model <- "gcgm"
#   }
# 
#   output <- bggm_res
#   return(output)
# }
