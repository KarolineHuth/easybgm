# 1. Turns vector into matrix
vector2matrix <- function(vec, p, diag = FALSE, bycolumn = FALSE) {
  m <- matrix(0, p, p)

  if(!bycolumn){
    m[lower.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
  } else {
    m[upper.tri(m, diag = diag)] <- vec
    m <- t(m)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
  }
  return(m)
}

# 2. Transform precision into partial correlations for interpretation
pr2pc <- function(K) {
  D.Prec = diag(diag(K)^(-.5))
  R <- diag(2,dim(K)[1])-D.Prec%*%K%*%D.Prec
  colnames(R) <- colnames(K)
  rownames(R) <- rownames(K)
  return(R)
}

# 3. BDgraph stores graphs as byte strings for efficiency
string2graph <- function(Gchar, p) {
  Gvec <- rep(0, p*(p-1)/2)
  edges <- which(unlist(strsplit(as.character(Gchar), "")) == 1)
  Gvec[edges] = 1
  G <- matrix(0, p, p)
  G[upper.tri(G)] <- Gvec
  G <- G + t(G)
  diag(G) <- 0
  return(G)
}

# 4. BDgraph extract posterior distribution for estimates
extract_posterior <- function(fit, data, method = c("ggm", "gcgm"), posterior_method = c("maximum-posterior", "model-averaged"), not_cont){
  m <- length(fit$all_graphs)
  n <- nrow(data)
  p <- ncol(data)


  if(method == "gcgm") {
    S <- BDgraph::get_S_n_p(data, method = method, n = n, not.cont = not_cont)$S
  } else {
    S <- t(data) %*% data
  }
  
  if(posterior_method == "MAP"){
    Rs = matrix(0, nrow = 10000, ncol = (p*(p-1))/2)
    index <- which.max(fit$graph_weights)
    graph_ix <- fit$sample_graphs[index]
    G <- string2graph(graph_ix, p)
    K <- BDgraph::rgwish(n=10000, adj=G, b=3+n, D=diag(p) + S)
    Rs <- list2matrix(K, p, convert = T)
  }
  if(posterior_method == "model-averaged"){

    n_structures <- length(fit$graph_weights)
    sum_weights <- sum(fit$graph_weights)
    structure_weights <- fit$graph_weights/sum_weights
    Rs = matrix(0, nrow = 1, ncol = (p*(p-1))/2)
    for (i in 1:n_structures) {
      graph_ix <- fit$sample_graphs[i]
      n_samples <- round(structure_weights[i]*100000)
      G <- string2graph(graph_ix, p)
      K <- BDgraph::rgwish(n=n_samples, adj=G, b=3+n, D=diag(p) + S)
      samples_ix <- list2matrix(K, p, convert = T)
      Rs <- rbind(Rs, samples_ix)

    }
  }

  return(list(Rs))
}

# 5. Samples from the G-wishart distribution
gwish_samples <- function(G, S, nsamples=1000) {
  n <- nrow(S)
  p <- ncol(S)
  #Rs <- array(0, dim=c(nsamples, p, p))
  Rs = matrix(0, nrow = nsamples, ncol = (p*(p-1))/2)

  for (i in 1:nsamples) {
    K <- BDgraph::rgwish(n=1, adj=G, b=3+n, D=diag(p) + S)*(G + diag(p))
    Rs[i,] <- as.vector(pr2pc(K)[upper.tri(pr2pc(K))])
    #Rs[i,,] <- .pr2pc(K)
  }
  return(Rs)
}


# 6. Centrality of weighted graphs

# Strength centrality only ## FASTER CODE
centrality <- function(res){
  Nsamples <- nrow(res$samples_posterior)
  p <- nrow(res$parameters)
  strength_samples <- matrix(0, nrow = Nsamples, ncol = p)
  for(i in 1:Nsamples){
    strength_samples[i, ] <- rowSums(abs(vector2matrix(res$samples_posterior[i,], p, bycolumn = T)))
  }
  return(strength_samples)
}

# Strength, betweenness and closeness centrality ## SLOWER CODE
centrality_all <- function(res){
  Nsamples <- nrow(res$samples_posterior)
  p <- as.numeric(nrow(res$parameters))
  samples <- res$samples_posterior
  for(i in 1:Nsamples){

    #Strength
    strength_samples <- rowSums(abs(vector2matrix(samples[i, ], p, bycolumn = TRUE)))
    #EI
    influence_samples <- rowSums(vector2matrix(samples[i, ], p, bycolumn = TRUE))

    DistMat <- 1/(ifelse(abs(vector2matrix(samples[i, ], p, bycolumn = TRUE))==0,0,abs(vector2matrix(samples[i, ], p, bycolumn = T))))
    igraphObject <- igraph::graph.adjacency(DistMat, weighted = TRUE, mode = "undirected")
    # Closeness
    closeness_samples <- igraph::closeness(igraphObject)
    # Betweenness
    betweenness_samples <- igraph::estimate_betweenness(igraphObject,cutoff = 1/1e-10)

    if(i > 1){
      centrality_samples <- rbind(centrality_samples, cbind(c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
                                                            rbind(strength_samples, closeness_samples, betweenness_samples, influence_samples)))
    } else {
      centrality_samples <- cbind(c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
                                  rbind(strength_samples, closeness_samples, betweenness_samples, influence_samples))
    }
  }
  colnames(centrality_samples) <- c("Centrality", colnames(res$parameters))
  centrality_samples <- as.data.frame(centrality_samples)
  centrality_samples[, 2:(p+1)] <- sapply(centrality_samples[, 2:(p+1)], as.numeric)
  return(centrality_samples)
}


# 7. turn list into matrix
list2matrix <- function(obj, p, convert = FALSE) {
  nlist <- length(obj)/(p*p)
  m <- obj[, , 1]
  nest <- sum(lower.tri(m))
  res <- matrix(0, nrow = nlist, ncol = nest)
  for(i in 1:nlist){
    m <- obj[, , i]
    if(convert == T){
      m <- pr2pc(m)
    }
    res[i, ] <- as.vector(m[lower.tri(m)])
  }
  return(res)
}

# 8. Set defaults of a function
set_defaults <- function(args, ...) {
  dots <- list(...)
  def_args <- setdiff(names(args), names(dots))
  dots[def_args] <- args[def_args]

  return(dots)
}

# 9. Check for empty ...

dots_check <- function(...){
  if(...length() > 0){
    warning("Arguments specified with ... are unused. ")
  }
}


# Sparse vs Dense Test
# The function tests if we have overlap in the posterior distributions,
# and with that if we need a(nother) bridge hypothesis.
is_overlap <- function(ordered_list) {
  
  for (i in 1: (length(ordered_list) - 1)) {
    #check all pairs of hypotheses
    this_el <- ordered_list[[i]]
    next_el <- ordered_list[[i + 1]]
    
    overlap <- which(this_el$tab != 0 & next_el$tab != 0)
    
    if (length(overlap) == 0) {
      before_position <- i
      if (this_el$alpha == next_el$alpha) {
        alpha <- this_el$alpha
        beta <- this_el$beta / 2
      }
      else if (this_el$beta == next_el$beta) {
        alpha <- next_el$alpha / 2
        beta <- this_el$beta
      }
      else {
        alpha <- this_el$alpha
        beta <- this_el$beta / 2
      }
      return(list(before_pos = before_position,
                  alpha = alpha, beta = beta))
    }
  }
  return(1)
}

# Given a list of the results for all needed hypotheses, the function
# computes the log BF of sparse against dense.
# @args ordered list of the results, with the outer two the hypotheses of 
# interest, and in between the bridge hypotheses. k the number of potential edges.
compute_bayes_factor <- function(ordered_list, k) {
  bf <- 0
  c <- 0: k
  for (i in 1: (length(ordered_list) - 1)) {
    el1 <- ordered_list[[i]]
    el2 <- ordered_list[[i + 1]]
    
    alpha1 <- el1$alpha
    alpha2 <- el2$alpha
    beta1 <- el1$beta
    beta2 <- el2$beta
    tab1 <- el1$tab
    tab2 <- el2$tab
    
    log_prior1 <- lchoose(k, c) - lbeta(alpha1, beta1) + lfactorial(alpha1 + c - 1) +
      lfactorial(beta1 + k - c - 1) - lfactorial(alpha1 + beta1 + k - 1)
    log_prior2 <- lchoose(k, c) - lbeta(alpha2, beta2) + lfactorial(alpha2 + c - 1) +
      lfactorial(beta2 + k - c - 1) - lfactorial(alpha2 + beta2 + k - 1)
    
    prob1 <- tab1 / sum(tab1)
    prob2 <- tab2 / sum(tab2)
    
    log_prob1 <- log(prob1)
    log_prob2 <- log(prob2)
    
    odds1 <- log_prior1 - log_prob1
    odds2 <- log_prob2 - log_prior2
    
    log_bf <- odds1 + odds2
    log_bf[is.infinite(log_bf)] <- NA
    log_bf <- mean(log_bf, na.rm = TRUE)
    
    bf <- bf + log_bf
    
  }
  return(bf)
}

# Given alpha and beta parameters computes the (log)probability of the beta Bernoulli distribution
# for the bgms package returns probability for a specific complexity
# @args alpha, beta parameters, c the complexity for which the probability has to b computed
# p the number of nodes
beta_bernoulli_prob <- function(c, alpha, beta, p) {
  
  k <- p * (p - 1) / 2

  log_nom <- lbeta(alpha + c, beta + k - c)
  log_denom <- lbeta(alpha, beta)

  log_choose <- lchoose(k, c)
  

  log_prob <- log_choose + log_nom - log_denom
  
  return(log_prob)
}

# Calculates the probability of an edge being present in the beta-bernoulli (BB) 
# or the stochastic block prior on the network structure of the bgms package
# Returns the prior inclusion probability for individual edges
# which simply calculates the expected value of the BB distribution 
# @args alpha, beta arguments of beta bernoulli prior, 


calculate_edge_prior <- function(alpha, beta) {
   alpha / (alpha + beta)
}


# include extractor functions to support bgms 0.1.3 version

extract_arguments <- function(bgms_object) {
  if(!inherits(bgms_object, what = "bgms"))
    stop(paste0("Expected an object with class bgms and not one with class ",
                class(bgms_object)))
  
  if(is.null(bgms_object$arguments)) {
    stop(paste0("Extractor functions have been defined for bgms versions 0.1.3 and up but not \n",
                "for older versions. The current fit object predates version 0.1.3."))
  } else {
    return(bgms_object$arguments)
  }
}

extract_pairwise_interactions <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  
  return(bgms_object$interactions)
}

extract_category_thresholds <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  
  return(bgms_object$thresholds)
}

extract_indicators <- function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  if(arguments$edge_selection & arguments$save) {
    if(arguments$version < "0.1.4") {
      edge_indicators = bgms_object$gamma
    } else {
      edge_indicators = bgms_object$indicator
    }
    return(edge_indicators)
  } else {
    stop(paste0("To access the sampled edge indicators the bgms package needs to be run using \n",
                "edge_selection = TRUE and save = TRUE."))
  }
}
