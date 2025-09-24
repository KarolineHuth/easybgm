# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms <- function(fit, type, data, iter, save,
                                 not_cont, centrality, progress, ...){

  if(!save && centrality){    # the save argument is not needed for bgms
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
                                     not_cont, data, centrality, ...) {

  # Normalize 'fit' and variable names
  if (!inherits(fit, "bgms")) {
    varnames <- fit$var_names
    fit <- fit$packagefit
    class(fit) <- "bgms"
  } else {
    varnames <- fit$arguments$data_columnnames
    if (is.null(varnames)) {
      varnames <- paste0("V", seq_len(fit$arguments$no_variables))
    }
  }

  # Arguments from bgms -
  args <- bgms::extract_arguments(fit)
  p <- args$no_variables

  #  Edge prior
  if (identical(args$edge_prior[1], "Bernoulli")) {
    edge.prior <- args$inclusion_probability
  } else {
    edge.prior <- calculate_edge_prior(
      alpha = args$beta_bernoulli_alpha,
      beta  = args$beta_bernoulli_beta
    )
    args$inclusion_probability <- edge.prior
  }

  # Save/draws behavior across versions
  bgms_ver <- tryCatch(utils::packageVersion("bgms"), error = function(e) "0.0.0")
  has_draws <- if (bgms_ver >= "0.1.6") TRUE else {
    if ("save" %in% names(args)) isTRUE(args$save) else isTRUE(save)
  }

  # Build result
  bgms_res <- list()

  # parameters
  pars <- get_interactions(fit, bgms_ver)
  bgms_res$parameters <- to_pp_matrix(pars, p = p)
  colnames(bgms_res$parameters) <- varnames
  rownames(bgms_res$parameters) <- varnames

  # thresholds
  bgms_res$thresholds <- get_thresholds(fit, bgms_ver)

  # default structure
  bgms_res$structure <- matrix(1L, nrow = p, ncol = p,
                               dimnames = list(varnames, varnames))

  # edge selection
  if (isTRUE(args$edge_selection)) {
    inc_probs <- bgms::extract_posterior_inclusion_probabilities(fit)
    inc_probs <- to_pp_matrix(inc_probs, p = p)
    colnames(inc_probs) <- rownames(inc_probs) <- varnames
    bgms_res$inc_probs <- inc_probs

    bgms_res$inc_BF <- (inc_probs / (1 - inc_probs)) /
      (edge.prior / (1 - edge.prior))
    colnames(bgms_res$inc_BF) <- rownames(bgms_res$inc_BF) <- varnames

    bgms_res$structure <- 1L * (inc_probs > 0.5)

    if (isTRUE(has_draws)) {
      gammas <- get_indicators(fit, bgms_ver)
      structures <- apply(gammas, 1L, paste0, collapse = "")
      table_structures <- as.data.frame(table(structures))
      bgms_res$structure_probabilities <- table_structures[, 2] / nrow(gammas)
      bgms_res$graph_weights <- table_structures[, 2]
      bgms_res$sample_graph <- as.character(table_structures[, 1])
    }
  }

  # posterior samples (if draws available)
  if (isTRUE(has_draws)) {
    bgms_res$samples_posterior <- get_interactions(fit, bgms_ver)
    if (isTRUE(centrality)) {
      bgms_res$centrality <- centrality(bgms_res)
    }
  }

  # bookkeeping
  bgms_res$model <- type
  bgms_res$fit_arguments <- args
  class(bgms_res) <- c("package_bgms", "easybgm")
  bgms_res
}
