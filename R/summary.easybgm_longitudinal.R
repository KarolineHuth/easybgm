#' @name summary.easybgm_longitudinal
#' @title  Summary method for \code{easybgm_longitudinal} objects
#'
#' @description Used to create a object of easybgm results and in turn print it
#'
#' @param object easybgm_longitudinal object
#' @param evidence_thresh Bayes Factor which will be considered sufficient evidence for in-/exclusion, default is 10.
#' @param ... unused argument
#'
#' @return Creates and prints the output of a Bayesian cross-sectional network analysis. The summary output has four parts. The first part lists the package used, the number of variables, and the data type. The second part is a matrix of edge-specific information. Each edge is listed in a row. This row contains the posterior parameter estimate, the posterior inclusion probability, the inclusion Bayes factor, and the categorization of the edge. The category encodes whether an edge is included, excluded, or inconclusive based on the inclusion Bayes factor. Users can set the threshold for the Bayes factor classification with the evidence threshold. By default, the threshold is set to $10$. The third part of the summary provides aggregated edge information. It lists the number of included, excluded, and inconclusive edges in the network, as well as the number of possible edges. This gives the user a quick overview of the robustness and density of the network. The higher the number of conclusive edges (i.e., classified as either included or excluded), the more robust the network. Conversely, if the network has a high percentage of inconclusive edges, the network is not robust. Researchers should refrain from making strong inferential conclusions. The final output section is a description of the structure uncertainty. It shows the number of structures visited, the number of possible structures, and the highest posterior structure probability. This last section can only be obtained for networks fitted with 'BDgraph' and 'bgms'. 
#'
#' @export

summary.easybgm_longitudinal <- function(object, evidence_thresh = 10, ...) {
  ## -----------------------------
  ## 0. Check arguments
  ## -----------------------------
  
  dots_check(...)
  
  ## -----------------------------
  ## 1. Determine number of nodes
  ## -----------------------------
  
  p <- object$p
  
  ## -----------------------------
  ## 2. Create data frame with edge-specific results
  ## -----------------------------
  
  ## ---- 2a. Parameters temporal 
  samples <- object$posterior_samples_temporal
  hdi_intervals <- as.data.frame(apply(samples, MARGIN = 2, FUN = hdi))
  posterior_medians <- apply(samples, MARGIN = 2, FUN = median)
  
  names <- colnames(object$parameters_temporal)
  names_bycol <- matrix(rep(names, each = object$p), ncol = object$p)
  names_byrow <- matrix(rep(names, each = object$p), ncol = object$p, byrow = T)
  names_comb <- matrix(paste0(names_byrow, "->", names_bycol), ncol = object$p)
  index <- as.vector(names_comb)
  
  posterior_temporal <- cbind(index, 
                              data.frame(round(posterior_medians, 3), row.names = NULL),
                              data.frame(t(round(hdi_intervals, 3)), row.names = NULL))
  colnames(posterior_temporal) <- c("Relation", "Parameter", "HDI lower", "HDI upper")
  
  ## ---- 2b. Parameters contemporaneous
  samples <- object$posterior_samples_contemporaneous
  hdi_intervals <- as.data.frame(apply(samples, MARGIN = 2, FUN = hdi))
  posterior_medians <- apply(samples, MARGIN = 2, FUN = median)
  
  names <- colnames(object$parameters_temporal)
  names_bycol <- matrix(rep(names, each = object$p), ncol = object$p)
  names_byrow <- matrix(rep(names, each = object$p), ncol = object$p, byrow = T)
  names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = object$p)
  index <- names_comb[upper.tri(names_comb)]
  
  posterior_contemporaneous <- cbind(index, 
                                     data.frame(round(posterior_medians, 3), row.names = NULL),
                                     data.frame(t(round(hdi_intervals, 3)), row.names = NULL) 
                                     )
  colnames(posterior_contemporaneous) <- c("Relation", "Parameter", "HDI lower", "HDI upper")
  
  ## ---- 2c. Classify edges ---- (i.e., similar, different, inconclusive)
  #### NOT YET IMPLEMENTED
  # category <- character(length(BF))
  # category[(BF < evidence_thresh) & (BF > 1/evidence_thresh)] <- "inconclusive"
  # category[BF > evidence_thresh] <- "different"
  # category[BF < 1/evidence_thresh] <- "similar"
  
  
  ## -----------------------------
  ## 3. Create summary output list
  ## -----------------------------
  out <- list()
  out$parameters_temporal <- posterior_temporal
  out$parameters_contemporaneous <- posterior_contemporaneous
  out$package <- strsplit(class(object)[1], "_")[[1]][2]
  out$model <- object$model
  out$n_nodes <- p
  out$n_possible_edges <- p*(p-1)/2
  
  ## ---- 3a. Aggregate edge counts ----
  #### NOT YET IMPLEMENTED
  # out$n_inclu_edges <- sum(BF > evidence_thresh)
  # out$n_incon_edges <- sum((BF < evidence_thresh) & (BF > 1/evidence_thresh))
  # out$n_exclu_edges <- sum(BF < 1/evidence_thresh)
  
  
  ## -----------------------------
  ## 4. Save call and BF threshold info
  ## -----------------------------
  out$fit_object <- object
  out$evidence_thresh <- evidence_thresh
  
  ## -----------------------------
  ## 5. Return summary object
  ## -----------------------------
  class(out) <- class(object)
  return(out)
  print(out)
}

#' @name print.easybgm_longitudinal
#' @title  Print method for \code{easybgm_longitudinal} objects
#'
#' @description Used to print easybgm results. The nicest overview is created by first feeding it to
#' `summary()`
#'
#' @param x easybgm_longitudinal object
#' @param ... unused argument
#' 
#' @return Prints the output of a Bayesian cross-sectional network comparison fitted with 'easybgm'
#' 
#' @export
#'

print.easybgm_longitudinal <- function(x, ...){
  
  dots_check(...)
  
  if(is.null(x$n_possible_edges)){
    #NextMethod("print")
    print(summary.easybgm_longitudinal(x))
  } else {
    cat("\n BAYESIAN LONGITUDINAL NETWORK ESTIMATION",
        "\n Model type:", x$model,
        "\n Number of nodes:", x$n_nodes,
        "\n Fitting Package:", x$package,
        "\n---",
        "\n TEMPORAL EDGE OVERVIEW",
        "\n")
    print(x$parameters_temporal, quote = FALSE, right = TRUE, row.names=F)
    cat("\n ---",
        "\n CONTEMPORANEOUS EDGE OVERVIEW",
        "\n")
    print(x$parameters_contemporaneous, quote = FALSE, right = TRUE, row.names=F)
  }
}

