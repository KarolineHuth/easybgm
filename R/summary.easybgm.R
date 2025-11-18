#' @name summary.easybgm
#' @title  Summary method for \code{easybgm} objects
#'
#' @description Used to create a object of easybgm results and in turn print it
#'
#' @param object easybgm object
#' @param evidence_thresh Bayes Factor which will be considered sufficient evidence for in-/exclusion, default is 10.
#' @param ... unused argument
#'
#' @return Creates and prints the output of a Bayesian cross-sectional network analysis. The summary output has four parts. The first part lists the package used, the number of variables, and the data type. The second part is a matrix of edge-specific information. Each edge is listed in a row. This row contains the posterior parameter estimate, the posterior inclusion probability, the inclusion Bayes factor, and the categorization of the edge. The category encodes whether an edge is included, excluded, or inconclusive based on the inclusion Bayes factor. Users can set the threshold for the Bayes factor classification with the evidence threshold. By default, the threshold is set to $10$. The third part of the summary provides aggregated edge information. It lists the number of included, excluded, and inconclusive edges in the network, as well as the number of possible edges. This gives the user a quick overview of the robustness and density of the network. The higher the number of conclusive edges (i.e., classified as either included or excluded), the more robust the network. Conversely, if the network has a high percentage of inconclusive edges, the network is not robust. Researchers should refrain from making strong inferential conclusions. The final output section is a description of the structure uncertainty. It shows the number of structures visited, the number of possible structures, and the highest posterior structure probability. This last section can only be obtained for networks fitted with 'BDgraph' and 'bgms'.
#'
#' @export

summary.easybgm <- function(object, evidence_thresh = 10, ...) {
  ## -----------------------------
  ## 0. Check arguments
  ## -----------------------------

  dots_check(...)

  ## -----------------------------
  ## 1. Determine number of nodes
  ## -----------------------------
  if(is.null(object$inc_probs)){
    p <- ncol(object$parameters)
  } else {
    p <- ncol(object$inc_probs)
  }

  ## -----------------------------
  ## 2. Create data frame with edge-specific results
  ## -----------------------------

  if(object$model %in% c("dgm-binary")){
    ## ---- 2a. Prepare relation names ----
    names <- colnames(object$inc_probs)
    names_bycol <- matrix(rep(names, each = p), ncol = p)
    names_byrow <- matrix(rep(names, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    mat_names <- names_comb[lower.tri(names_comb)]

    ## ---- 2b. Extract inclusion probabilities and Bayes Factors ----
    inc_probs  <- round(object$inc_probs, 3)[lower.tri(object$inc_probs)]
    BF <- round(object$inc_BF, 3)[lower.tri(object$inc_BF)]

    ## ---- 2c. Classify edges based on Bayes Factor ----(i.e., included, excluded, inconclusive)
    category <- character(length(BF))
    category[(BF < evidence_thresh) & (BF > 1/evidence_thresh)] <- "inconclusive"
    category[BF > evidence_thresh] <- "included"
    category[BF < 1/evidence_thresh] <- "excluded"

    ## ---- 2d. Create results data frame ----
    results <-
      data.frame(
        relation = mat_names,
        inc_probs =  inc_probs,
        BF = BF,

        category = category,
        row.names = NULL
      )
    colnames(results) <- c(
      "Relation",
      "Posterior Incl. Prob.",
      "Inclusion BF",
      "Category")
  } else if(is.null(object$inc_probs)){
    ## ---- 2e. Case: Only parameter estimates (no inclusion probs) ----
    names <- colnames(object$parameters)
    names_bycol <- matrix(rep(names, each = p), ncol = p)
    names_byrow <- matrix(rep(names, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    mat_names <- names_comb[lower.tri(names_comb)]

    parameter_values <- round(object$parameters, 3)[lower.tri(object$parameters)]


    results <-
      data.frame(
        relation = mat_names,
        parameter_values = parameter_values,
        row.names = NULL
      )
    colnames(results) <- c(
      "Relation",
      "Parameter")
  } else {
    ## ---- 2f. General case: parameters + inclusion probs + Bayes Factors
    names <- colnames(object$parameters)
    names_bycol <- matrix(rep(names, each = p), ncol = p)
    names_byrow <- matrix(rep(names, each = p), ncol = p, byrow = T)
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
    mat_names <- names_comb[lower.tri(names_comb)]

    ## ---- 2g. Extract and round relevant values ----
    parameter_values <- round(object$parameters, 3)[lower.tri(object$parameters)]
    inc_probs  <- round(object$inc_probs, 3)[lower.tri(object$inc_probs)]
    BF <- round(object$inc_BF, 3)[lower.tri(object$inc_BF)]

    ## ---- 2h. Classify edges ---- (i.e., included, excluded, inconclusive)
    category <- character(length(BF))
    category[(BF < evidence_thresh) & (BF > 1/evidence_thresh)] <- "inconclusive"
    category[BF > evidence_thresh] <- "included"
    category[BF < 1/evidence_thresh] <- "excluded"

    ## ---- 2i. Create results data frame ----
    ## ----  Create results data frame with convergence (newer bgms)----
    if("package_bgms" %in% class(object) && packageVersion("bgms") > "0.1.4.2"){
      results <-
        data.frame(
          relation = mat_names,
          parameter_values = parameter_values,
          inc_probs =  inc_probs,
          BF = BF,
          category = category,
          convergence = round(object$convergence_parameter, 3)
        )
      colnames(results) <- c(
        "Relation",
        "Estimate",
        "Posterior Incl. Prob.",
        "Inclusion BF",
        "Category",
        "Convergence")
    } else {
      ## ----  Create results data frame without convergence----
      results <-
        data.frame(
          relation = mat_names,
          parameter_values = parameter_values,
          inc_probs =  inc_probs,
          BF = BF,
          category = category
        )
      colnames(results) <- c(
        "Relation",
        "Estimate",
        "Posterior Incl. Prob.",
        "Inclusion BF",
        "Category")
    }
  }

  ## -----------------------------
  ## 3. Create summary output list
  ## -----------------------------
  out <- list()
  out$parameters <- results
  out$package <- strsplit(class(object)[1], "_")[[1]][2]
  out$model <- object$model
  out$n_nodes <- p
  out$n_possible_edges <- p*(p-1)/2

  ## ---- 3a. Aggregate edge counts ----
  if(!is.null(object$inc_probs)){
    out$n_inclu_edges <- sum(BF > evidence_thresh)
    out$n_incon_edges <- sum((BF < evidence_thresh) & (BF > 1/evidence_thresh))
    out$n_exclu_edges <- sum(BF < 1/evidence_thresh)
  }

  ## ---- 3b. Structure uncertainty (only for BDgraph/bgms) ----
  if(!is.null(object$structure_probabilities)){
    out$possible_struc <- 2^(p*(p-1)/2)
    out$n_structures <- length(object$sample_graph)
    out$max_structure_prob <- max(object$structure_probabilities)
  }

  ## ---- 3c. Clustering information (Stochastic-Block prior) ----
  if(!is.null(object$fit_arguments)) {
    if(!is.null(object$fit_arguments$edge_prior) && object$fit_arguments$edge_prior == "Stochastic-Block") {
      out$posterior_num_blocks <- cbind(seq_len(nrow(object$sbm$posterior_num_blocks)), round(object$sbm$posterior_num_blocks, 3))
      out$posterior_node_allocations <- data.frame(object$fit_arguments$data_columnnames, object$sbm$posterior_mean_allocations,  object$sbm$posterior_mode_allocations)
      colnames(out$posterior_num_blocks) <- c(
        "No. of Clusters",
        "Posterior Prob.")
      colnames(out$posterior_node_allocations) <- c("Node", "Posterior Mean", "Posterior Mode")
      out$BF <- clusterBayesfactor(object, type = "complement")
    }
  }

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

#' @name print.easybgm
#' @title  Print method for \code{easybgm} objects
#'
#' @description Used to print easybgm results. The nicest overview is created by first feeding it to
#' `summary()`
#'
#' @param x easybgm object
#' @param ... unused argument
#'
#' @return Prints the output of a Bayesian cross-sectional network model fitted with 'easybgm'
#'
#' @export
#'

print.easybgm <- function(x, ...){

  dots_check(...)

  if(is.null(x$n_possible_edges)){
    #NextMethod("print")
    print(summary.easybgm(x))
  } else if(any(class(x) == "package_bggm")){
    cat("\n BAYESIAN ANALYSIS OF NETWORKS",
        "\n Model type:", x$model,
        "\n Number of nodes:", x$n_nodes,
        "\n Fitting Package:", x$package,
        "\n---",
        "\n EDGE SPECIFIC OVERVIEW",
        "\n")
    print(x$parameters, quote = FALSE, right = TRUE, row.names=F)
    cat("\n Bayes factors larger than", x$evidence_thresh, "were considered sufficient evidence for the classification.",
        "\n Bayes factors were obtained via single-model comparison.",
        "\n ---",
        "\n AGGREGATED EDGE OVERVIEW",
        "\n Number of edges with sufficient evidence for inclusion:", x$n_inclu_edges,
        "\n Number of edges with insufficient evidence:", x$n_incon_edges,
        "\n Number of edges with sufficient evidence for exclusion:", x$n_exclu_edges,
        "\n Number of possible edges:", x$n_possible_edges,
        "\n")
  } else if (ncol(x$parameters) < 3) {
    cat("\n BAYESIAN ANALYSIS OF NETWORKS",
        "\n Model type:", x$model,
        "\n Number of nodes:", x$n_nodes,
        "\n Fitting Package:", x$package,
        "\n---",
        "\n EDGE SPECIFIC OVERVIEW",
        "\n")
    print(x$parameters, quote = FALSE, right = TRUE, row.names=F)
  } else if(is.null(x$n_structures)){
    cat("\n BAYESIAN ANALYSIS OF NETWORKS",
        "\n Model type:", x$model,
        "\n Number of nodes:", x$n_nodes,
        "\n Fitting Package:", x$package,
        "\n---",
        "\n EDGE SPECIFIC OVERVIEW",
        "\n")
    print(x$parameters, quote = FALSE, right = TRUE, row.names=F)
    cat("\n Bayes Factors larger than", x$evidence_thresh, "were considered sufficient evidence for the classification",
        "\n Bayes factors were obtained using Bayesian model-averaging.",
        "\n ---",
        "\n AGGREGATED EDGE OVERVIEW",
        "\n Number of edges with sufficient evidence for inclusion:", x$n_inclu_edges,
        "\n Number of edges with insufficient evidence:", x$n_incon_edges,
        "\n Number of edges with sufficient evidence for exclusion:", x$n_exclu_edges,
        "\n Number of possible edges:", x$n_possible_edges,
        "\n")

    if(!is.null(x$posterior_num_blocks)){
      cat("\n ---",
          "\n CLUSTERING OVERVIEW",
          "\n")
      print(x$posterior_num_blocks, quote = FALSE, right = TRUE, row.names=F)
      cat("\n Estimated cluster assignments of the nodes:",
          "\n")
      print(x$posterior_node_allocations, quote = FALSE, right = TRUE, row.names=F)
      cat("\n Bayes factor in favor of more than one cluster:", x$BF,
          "\n If you wish to test hypotheses about specific number of clusters,",
          "\n please use the clusterBayesfactor function.",
          "\n")
    }
  } else {
    cat("\n BAYESIAN ANALYSIS OF NETWORKS",
        "\n Model type:", x$model,
        "\n Number of nodes:", x$n_nodes,
        "\n Fitting Package:", x$package,
        "\n---",
        "\n EDGE SPECIFIC OVERVIEW",
        "\n")
    print(x$parameters, quote = FALSE, right = TRUE, row.names=F)
    cat("\n Bayes Factors larger than", x$evidence_thresh, "were considered sufficient evidence for the classification",
        "\n Bayes factors were obtained using Bayesian model-averaging.",
        "\n ")
    if("package_bgms" %in% class(x) && packageVersion("bgms") > "0.1.4.2"){
      cat("\n Convergence indicates the R-hat (Gelman–Rubin) statistic measuring how well MCMC chains have converged to",
          "\n the same target distribution. Values greater than about 1.01–1.05 are considered concerning, indicating",
          "\n potential lack of convergence for the estimates of the pairwise interactions. ",
          "\n ---")
    }
    cat("\n AGGREGATED EDGE OVERVIEW",
        "\n Number of edges with sufficient evidence for inclusion:", x$n_inclu_edges,
        "\n Number of edges with insufficient evidence:", x$n_incon_edges,
        "\n Number of edges with sufficient evidence for exclusion:", x$n_exclu_edges,
        "\n Number of possible edges:", x$n_possible_edges,
        "\n")

    if(!is.null(x$posterior_num_blocks)){
      cat("\n ---",
          "\n CLUSTERING OVERVIEW",
          "\n")
      print(x$posterior_num_blocks, quote = FALSE, right = TRUE, row.names=F)
      cat("\n Estimated cluster assignments of the nodes:",
          "\n")
      print(x$posterior_node_allocations, quote = FALSE, right = TRUE, row.names=F)
      cat("\n Bayes factor in favor of more than one cluster:", x$BF,
          "\n If you wish to test hypotheses about specific number of clusters,",
          "\n please use the clusterBayesfactor function.",
          "\n")
    }

    cat("\n ---",
        "\n STRUCTURE OVERVIEW",
        "\n Number of visited structures:", x$n_structures,
        "\n Number of possible structures:", x$possible_struc,
        "\n Posterior probability of most likely structure:", x$max_structure_prob,
        "\n---")
  }
}
