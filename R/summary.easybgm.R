summary.easybgm <- function(x) {


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
      "Inclusion BF")
  } else {
    parameter_values <- round(x$parameters, 3)[upper.tri(x$parameters)]
    inc_probs  <- round(x$inc_probs, 3)[upper.tri(x$inc_probs)]
    BF <- round(x$BF, 3)[upper.tri(x$BF)]

    results <-
      data.frame(
        relation = mat_names,
        parameter_values = parameter_values,
        inc_probs =  inc_probs,
        BF = BF
      )
    colnames(results) <- c(
      "Relation",
      "Estimate",
      "Posterior Incl. Prob.",
      "Inclusion BF")
  }
  # create list with
  out <- list()
  out$parameters <- results
  out$n_structures <- length(x$sample_graphs)
  out$max_structure_prob <- max(x$structure_probabilities)
  out$package <- strsplit(class(x)[1], "_")[[1]][2]
  out$model <- x$model
  out$n_nodes <- p
  out$n_possible_edges <- p*(p-1)/2
  out$n_inc_edges <- sum(inc_probs > 0.5)
  out$possible_struc <- 2^(p*(p-1)/2)
  out$fit_object <- x

  # return object
  class(out) <- class(x)
  return(out)
  print(out)
}

print.easybgm <- function(x){
  cat("\n BAYESIAN ANALYSIS OF NETWORKS",
      "\n Model type:", x$model,
      "\n Number of nodes:", x$n_nodes,
      "\n Fitting Package:", x$package,
      "\n---",
      "\n")
  print(x$parameters, quote = FALSE, right = TRUE, row.names=F)
  cat("---\n",
      "\n Number of included edges:", x$n_inc_edges,
      "\n Number of possible edges:", x$n_possible_edges,
      "\n",
      "\n Number of visited structures:", x$n_structures,
      "\n Number of possible structures:", x$possible_struc,
      "\n Posterior probability of most likely structure:", x$max_structure_prob,
      "\n---")
}

