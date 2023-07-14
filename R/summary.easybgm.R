#' Summary function for easybgm objects
#'
#' @param x easybgm object
#' @param evidence_thresh Bayes Factor which will be considered sufficient evidence for in-/exclusion, default is 10.
#'
#' @export

summary.easybgm <- function(x, evidence_thresh = 10) {


  # nodes
  p <- ncol(x$parameters)

  # names for each relation
  names <- colnames(x$parameters)
  names_bycol <- matrix(rep(names, each = p), ncol = p)
  names_byrow <- matrix(rep(names, each = p), ncol = p, byrow = T)
  names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = p)
  mat_names <- names_comb[upper.tri(names_comb)]

  # create data frame with parameter results
  if(x$model %in% c("dgm-binary")){
    inc_probs  <- round(x$inc_probs, 3)[upper.tri(x$inc_probs)]
    BF <- round(x$BF, 3)[upper.tri(x$BF)]

    results <-
      data.frame(
        relation = mat_names,
        inc_probs =  inc_probs,
        BF = BF, row.names = NULL
      )
    colnames(results) <- c(
      "Relation",
      "Posterior Incl. Prob.",
      "Inclusion BF",
      "Category")
  } else {
    parameter_values <- round(x$parameters, 3)[upper.tri(x$parameters)]
    inc_probs  <- round(x$inc_probs, 3)[upper.tri(x$inc_probs)]
    BF <- round(x$BF, 3)[upper.tri(x$BF)]
    #create the category of the edge (i.e., included, excluded, inconclusive)
    category <- character(length(BF))
    category[abs(BF) < evidence_thresh] <- "inconclusive"
    category[BF > evidence_thresh] <- "included"
    category[BF < -evidence_thresh] <- "excluded"


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
  # create list with output
  out <- list()
  out$parameters <- results
  out$n_structures <- length(x$sample_graphs)
  out$max_structure_prob <- max(x$structure_probabilities)
  out$package <- strsplit(class(x)[1], "_")[[1]][2]
  out$model <- x$model
  out$n_nodes <- p
  out$n_possible_edges <- p*(p-1)/2
  out$n_inclu_edges <- sum(BF > evidence_thresh)
  out$n_incon_edges <- sum(abs(BF) < evidence_thresh)
  out$n_exclu_edges <- sum(BF < -evidence_thresh)
  out$possible_struc <- 2^(p*(p-1)/2)
  out$fit_object <- x
  out$evidence_thresh <- evidence_thresh

  # return object
  class(out) <- class(x)
  return(out)
  print(out)
}

#' Print function for easybgm objects
#'
#' @param x easybgm object
#'
#' @export
#'

print.easybgm <- function(x){
  if(is.null(x$n_inclu_edges)){
    print(summary.easybgm(x))
  } else {
  cat("\n BAYESIAN ANALYSIS OF NETWORKS",
      "\n Model type:", x$model,
      "\n Number of nodes:", x$n_nodes,
      "\n Fitting Package:", x$package,
      "\n---",
      "\n EDGE SPECIFIC OVERVIEW",
      "\n")
  print(x$parameters, quote = FALSE, right = TRUE, row.names=F)
  cat("\n Bayes Factors larger", x$evidence_thresh, "were considered sufficient evidence for the categorization.",
      "\n ---",
      "\n AGGREGATED EDGE OVERVIEW",
      "\n Number of included edges:", x$n_inclu_edges,
      "\n Number of inconclusive edges:", x$n_incon_edges,
      "\n Number of excluded edges:", x$n_exclu_edges,
      "\n Number of possible edges:", x$n_possible_edges,
      "\n",
      "\n ---",
      "\n STRUCTURE OVERVIEW",
      "\n Number of visited structures:", x$n_structures,
      "\n Number of possible structures:", x$possible_struc,
      "\n Posterior probability of most likely structure:", x$max_structure_prob,
      "\n---")
  }
}

