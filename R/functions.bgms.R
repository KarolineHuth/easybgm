# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not_cont, centrality, progress, ...){

  if(!save && centrality){
    save <- TRUE
  }

  if(type == "binary") {
    type <- "ordinal"
  }

  if(packageVersion("bgms") > "0.1.4.2"){
    bgms_fit <- do.call(
      bgm, c(list(x = data, iter = iter,
                  variable_type = type,
                  display_progress = progress,
                  ...))
    )}
  if(packageVersion("bgms") < "0.1.6"){
    bgms_fit <- do.call(
      bgm, c(list(x = data, iter = iter, save = T,
                  variable_type = type,
                  display_progress = progress,
                  ...))
    )}


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

  if (packageVersion("bgms") < "0.1.4") {
    stop("easybgm now requires bgms version 0.1.4 or higher.")
  }

  # --- Ensure proper bgms object and variable names ---
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
    class(fit) <- "bgms"
  } else {
    varnames <- fit$arguments$data_columnnames
    if (is.null(varnames)) {
      varnames <- paste0("V", 1:fit$arguments$no_variables)}
  }

  if(packageVersion("bgms") > "0.1.4.2"){
    class(fit) <- "bgms"
  }

  # --- Extract model arguments and edge priors ---
  args <- bgms::extract_arguments(fit)
  args$save <- save

  if (args$edge_prior[1] == "Bernoulli") {
    edge.prior <- args$inclusion_probability
  } else {
    edge.prior <- calculate_edge_prior(
      alpha = args$beta_bernoulli_alpha,
      beta = args$beta_bernoulli_beta
    )
    args$inclusion_probability <- edge.prior
  }

  bgms_res <- list()

  # --- Main extraction ---
  if (args$save) {
    p <- args$no_variables
    pars <- extract_pairwise_interactions(fit)
    bgms_res$parameters <- vector2matrix(colMeans(pars), p = p)
    bgms_res$samples_posterior <- extract_pairwise_interactions(fit)
    bgms_res$thresholds <- extract_category_thresholds(fit)
    colnames(bgms_res$parameters) <- varnames

    bgms_res$structure <- matrix(1, ncol = p, nrow = p)

    if (args$edge_selection) {
      bgms_res$inc_probs <- extract_posterior_inclusion_probabilities(fit)
      bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
        (edge.prior / (1 - edge.prior))
      bgms_res$structure <- 1 * (bgms_res$inc_probs > 0.5)

      gammas <- extract_indicators(fit)
      structures <- apply(gammas, 1, paste0, collapse = "")
      table_structures <- as.data.frame(table(structures))

      bgms_res$structure_probabilities <- table_structures[, 2] / nrow(gammas)
      bgms_res$graph_weights <- table_structures[, 2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
    }
  } else {
    p <- args$no_variables
    pars <- extract_pairwise_interactions(fit)
    bgms_res$parameters <- vector2matrix(colMeans(pars), p = p)
    bgms_res$thresholds <- extract_category_thresholds(fit)
    colnames(bgms_res$parameters) <- varnames

    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters),
                                 nrow = nrow(bgms_res$parameters))

    if (args$edge_selection) {
      bgms_res$inc_probs <- extract_posterior_inclusion_probabilities(fit)
      bgms_res$inc_BF <- (bgms_res$inc_probs / (1 - bgms_res$inc_probs)) /
        (edge.prior / (1 - edge.prior))
      bgms_res$structure <- 1 * (bgms_res$inc_probs > 0.5)
    }
  }

  if (args$edge_prior[1] == "Stochastic-Block" && packageVersion("bgms") > "0.1.6") {
    bgms_res$sbm <- extract_sbm(fit)
  }

  # --- Optionally compute centrality ---
  if (centrality) {
    bgms_res$centrality <- centrality(bgms_res)
  }

  # --- For newer version compute convergence ---
  if (packageVersion("bgms") > "0.1.4.2") {
    bgms_res$convergence_parameter <-  fit$posterior_summary_pairwise$Rhat
  }

  # --- Finalize output ---
  colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
  colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters)

  bgms_res$model <- type
  bgms_res$fit_arguments <- args
  output <- bgms_res

  class(output) <- c("package_bgms", "easybgm")
  return(output)

}
