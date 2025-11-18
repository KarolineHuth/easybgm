#' @export
plot_structure_probabilities.easybgm_longitudinal <- function(output, as_BF = FALSE, ...) {
  stop("This plot is currently not supported for easybgm_longitudinal objects.")
}

# ---------------------------------------------------------------------------------------------------------------
#' @export
plot_complexity_probabilities.easybgm_longitudinal <- function(output, ...) {
  stop("This plot is currently not supported for easybgm_longitudinal objects.")
}

# ---------------------------------------------------------------------------------------------------------------
#' @export
plot_edgeevidence.easybgm_longitudinal <- function(output, 
                                                   evidence_thresh = 10, 
                                                   pars = "temporal", 
                                                   ...) {
  
  default_args <- list(
    edge.color = c("#36648b","#86a2b9", "#bfbfbf", "#f9d183", "#eeb004"),
    colnames = colnames(output$parameters),
    layout = qgraph::averageLayout(as.matrix(output$parameters_temporal*(output$bf_temporal>1))),
    legend = TRUE,
    vsize = 10,
    nodeNames = colnames(output$parameters),
    edge.width = 2,
    label.cex = 1,
    legend.cex = .6
  )
  
  args <- set_defaults(default_args, ...)
  print(args$edge.color)
  if(pars == "temporal"){
    graph <- output$bf_temporal
    
    # assign a color to each edge (inclusion - blue, exclusion - red, no conclusion - grey)
    graph_color <- graph
    # 1. Most evidence for inclusion
    graph_color[graph > evidence_thresh] <- args$edge.color[1]
    # 2. Moderate inclusion (BF > 3 but ≤ evidence_thresh)
    graph_color[graph > 3 & graph <= evidence_thresh] <- args$edge.color[2]
    # 3. Inconclusive (BF between 1/3 and 3)
    graph_color[graph >= 1/3 & graph <= 3] <- args$edge.color[3]
    # 4. Moderate exclusion (BF < 1/3 but > 1/evidence_thresh)
    graph_color[graph < 1/3 & graph > 1/evidence_thresh] <- args$edge.color[4]
    # 5. Strong evidence for exclusion (BF ≤ 1/evidence_thresh)
    graph_color[graph <= 1/evidence_thresh] <- args$edge.color[5]
    print(graph_color)
  } 
  if(pars == "contemporaneous"){
    graph <- output$bf_contemporaneous
    diag(graph) <- 1
    
    # assign a color to each edge (inclusion - blue, exclusion - red, no conclusion - grey)
    graph_color <- graph
    # 1. Most evidence for inclusion
    graph_color[graph > evidence_thresh] <- args$edge.color[1]
    # 2. Moderate inclusion (BF > 3 but ≤ evidence_thresh)
    graph_color[graph > 3 & graph <= evidence_thresh] <- args$edge.color[2]
    # 3. Inconclusive (BF between 1/3 and 3)
    graph_color[graph >= 1/3 & graph <= 3] <- args$edge.color[3]
    # 4. Moderate exclusion (BF < 1/3 but > 1/evidence_thresh)
    graph_color[graph < 1/3 & graph > 1/evidence_thresh] <- args$edge.color[4]
    # 5. Strong evidence for exclusion (BF ≤ 1/evidence_thresh)
    graph_color[graph <= 1/evidence_thresh] <- args$edge.color[5]
  }
  
  graph[graph > -.001] <- 1
  colnames(graph) <- colnames(output$fn_args$data)
  qgraph_plot <- qgraph::qgraph(graph,
                                edge.color = graph_color,
                                layout = args$layout, 
                                vsize = args$vsize,
                                nodeNames = args$nodeNames,
                                edge.width = args$edge.width,
                                label.cex = args$label.cex,
                                legend.cex = args$legend.cex,
                                ...
  )
}

# ---------------------------------------------------------------------------------------------------------------


#' @export
plot_network.easybgm_longitudinal <- function(output, 
                                              exc_prob = 0.5, 
                                              evidence_thresh = 10, 
                                              dashed = FALSE, 
                                              pars = "temporal", 
                                              ...) {
  
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  
  if(is.null(output$inc_probs) & dashed == TRUE){
    dashed <- FALSE
    warning("The model was fitted without edge selection and no inclusion probabilities were obtained. Therefore, edges cannot be dashed according to their PIP.",
            call. = FALSE)
  }
  
  
  if(pars == "temporal"){
    graph <- output$parameters_temporal
  } else {
    graph <- output$parameters_contemporaneous
  }
  
  default_args <- list(
    layout = qgraph::averageLayout(as.matrix(graph)),
    evidence_thresh = 10,
    theme = "TeamFortress",
    vsize = 10,
    nodeNames = colnames(graph),
    legend = T,
    label.cex = 1.2,
    legend.cex = .6, 
    edge.labels = FALSE
  )
  args <- set_defaults(default_args, ...)
  
  # Currently not supported
  # Exclude edges with an inclusion probability lower than exc_prob
  # inc_probs_m <- output$inc_probs
  # graph[inc_probs_m < exc_prob] <- 0
  # diag(graph) <- 1
  
  # Plot
  if(dashed){
    # currently not supported
    # graph_dashed <- ifelse(output$inc_BF < args$evidence_thresh, 2, 1)
    # 
    # 
    # qgraph_plot <- qgraph::qgraph(graph, layout = args$layout, 
    #                               lty = graph_dashed,
    #                               theme = args$theme, vsize = args$vsize,
    #                               nodeNames = args$nodeNames,
    #                               legend = args$legend,
    #                               label.cex = args$label.cex,
    #                               legend.cex = args$legend.cex, 
    #                               edge.labels = args$edge.labels, ...)
  } else {
    qgraph_plot <- qgraph::qgraph(graph, theme = args$theme, 
                                  layout = args$layout, vsize = args$vsize,
                                  nodeNames = args$nodeNames,
                                  legend = args$legend,
                                  label.cex = args$label.cex,
                                  legend.cex = args$legend.cex, 
                                  edge.labels = args$edge.labels, ...)
    
  }
  return(invisible(qgraph_plot))
}

# -------------------------------------------------
#' @export
plot_structure.easybgm_longitudinal <- function(output, ...) {
  stop("This plot is currently not supported for easybgm_longitudinal objects.")
}

# ---------------------------------------------------------------------------------------------------------------
#' @export
plot_parameterHDI.easybgm_longitudinal <- function(output, 
                                                   pars = "temporal", ...) {
  
  def_args <- list(
    theme_ = theme_bw(),
    geom_pointrange = geom_pointrange(position = position_dodge(width = c(0.3)),
                                      size = .3),
    xlab = "",
    ylab = "Highest Density Interval of Parameter",
    geom_hline = geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.3),
    theme = theme(
      axis.text = element_text(size=8),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 1.1),
      axis.ticks.length = unit(.2, "cm"),
      axis.ticks = element_line(linewidth = .8),
      axis.title.x = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 18, face = "bold")
    )
  )
  
  args <- set_defaults(def_args, ...)
  
  if(pars == "temporal"){
    samples <- output$posterior_samples_temporal
  } else {
    samples <- output$posterior_samples_contemporaneous
  }
  
  hdi_intervals <- as.data.frame(apply(samples, MARGIN = 2, FUN = hdi))
  posterior_medians <- apply(samples, MARGIN = 2, FUN = median)
  
  names <- colnames(output$parameters_temporal)
  names_bycol <- matrix(rep(names, each = output$p), ncol = output$p)
  names_byrow <- matrix(rep(names, each = output$p), ncol = output$p, byrow = T)
  if(pars == "temporal"){
    names_comb <- matrix(paste0(names_byrow, "->", names_bycol), ncol = output$p)
    index <- as.vector(names_comb)
  } else {
    names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = output$p)
    index <- names_comb[upper.tri(names_comb)]
  }
  
  
  posterior <- cbind(data.frame(posterior_medians, row.names = NULL),
                     data.frame(t(hdi_intervals), row.names = NULL), 
                     index)
  colnames(posterior) <- c("posterior_medians", "lower", "upper", "names")
  posterior <- posterior[order(posterior$posterior_medians, decreasing = F),]
  posterior$names <- factor(posterior$names, levels = posterior$names)
  
  
  ggplot2::ggplot(data = posterior, aes(x = .data$names, y = .data$posterior_medians, ymin = .data$lower,
                                        ymax = .data$upper, ...)) +
    args$geom_pointrange +
    args$theme_ +
    coord_flip() +
    ylab(args$ylab)+
    xlab(args$xlab) +
    args$geom_hline +
    args$theme
  
}


# ---------------------------------------------------------------------------------------------------------------
#' @export
plot_centrality.easybgm_longitudinal <- function(output, 
                                                 group_names = group_names, 
                                                 metric = "strength", 
                                                 ...){
  
  if(is.null(output$centrality)){
    stop("Centrality results are required. When estimating the model, set \"centrality = TRUE\".")
  }
  
  default_args <- list(
    theme_ = theme_minimal(),
    ylab = "Centrality",
    xlab = "Nodes",
    geom_errorbar = geom_errorbar(aes(y=.data$mean, ymin = .data$lower, ymax = .data$upper)
                                  , linewidth =.5, width=.4),
    theme = theme(
      axis.text = element_text(size=16),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 1.1),
      axis.ticks.length = unit(.2, "cm"),
      axis.ticks = element_line(linewidth = .8),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 18, face = "bold"),
      panel.grid.major = element_blank()
    )
  )

  args <- set_defaults(default_args, ...)
  cent_samples <- output$centrality_samples[[metric]]
  p <- ncol(output$parameters_cont)
  rownames(cent_samples) <- NULL
  # Creating summary statistics

  centrality_means <- colMeans(cent_samples)
  centrality_hdi <- apply(cent_samples, MARGIN = 2, FUN = hdi, allowSplit = F)
  centrality_summary <- data.frame(node = colnames(output$parameters_cont),
                                   mean = centrality_means,
                                   lower = centrality_hdi[1, ],
                                   upper = centrality_hdi[2, ])

  centrality_summary |>
    dplyr::arrange(mean) |>
    dplyr::mutate(node = factor(.data$node, levels = .data$node)) |>
    ggplot2::ggplot(aes(x = .data$node, y=.data$mean, ...))+
    args$theme_ +
    geom_point()+
    args$geom_errorbar+
    coord_flip() +
    ylab(args$ylab) +
    xlab(args$xlab) +
    args$theme
  

  
}

# ---------------------------------------------------------------------------------------------------------------
# Centrality plot for two or more groups
#' @export
plot_centrality.list <- function(output, group_names = NULL, ...){
  
  stop("This plot is currently not supported for easybgm_longitudinal objects.")
  
  # Convert to easybgm output if provided objects are bgms
  
  # Check for bgms package version 
  if(any(class(output[[1]]) == "bgms")) {
    if(packageVersion("bgms") < "0.1.3"){
      stop("Your version of the package bgms is not supported anymore. Please update.")
    }
    
    res <- list()
    for(i in 1:length(output)) {
      fit_args <- bgms::extract_arguments(output[[i]])
      
      if(!fit_args$save){
        stop("Samples of the posterior distribution required but not required for at least one fit. When estimating the model with bgm, set \"save = TRUE\".")
      }
      
      fit_args <- bgms::extract_arguments(output[[i]])
      
      res[[i]] <- bgm_extract.package_bgms(fit = output[[i]], save = fit_args$save, centrality = TRUE,
                                           type = NULL, not_cont = NULL, data = NULL,
                                           edge_prior = fit_args$edge_prior,
                                           inclusion_probability  = fit_args$inclusion_probability,
                                           beta_bernoulli_alpha = fit_args$beta_bernoulli_alpha,
                                           beta_bernoulli_beta = fit_args$beta_bernoulli_beta)
      
    }
    output <- res
  }
  
  
  default_args <- list(
    theme_ = theme_minimal(),
    ylab = "Strength Centrality",
    xlab = "Nodes",
    geom_errorbar = geom_errorbar(aes(y=.data$mean, ymin = .data$lower, ymax = .data$upper)
                                  , linewidth =.5, width=.4, position = position_dodge(width = 0.6)),
    theme = theme(
      axis.text = element_text(size=16),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 1.1),
      axis.ticks.length = unit(.2, "cm"),
      axis.ticks = element_line(linewidth = .8),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 18, face = "bold"),
      panel.grid.major = element_blank(), 
      text = element_text(size = 14)
    )
  )
  
  args <- set_defaults(default_args, ...)
  
  
  for(i in 1:length(output)){
    
    if(is.null(output[[i]]$centrality)){
      stop("At least one provided fit does not have the centrality values provided. When estimating the model, set \"centrality = TRUE\".")
    }
    
    cent_samples <- output[[i]]$centrality
    p <- ncol(output[[i]]$parameters)
    rownames(cent_samples) <- NULL
    
    if(is.null(group_names)){
      group_i <- paste0("Group ", i)
    } else {
      group_i <- group_names[i]
    }
    
    # Creating summary statistics
    centrality_means <- colMeans(cent_samples)
    centrality_hdi <- apply(cent_samples, MARGIN = 2, FUN = hdi, allowSplit = F)
    centrality_summary_i <- data.frame(node = colnames(output[[i]]$parameters),
                                       group = rep(group_i, nrow(output[[i]]$parameters)),
                                       mean = centrality_means,
                                       lower = centrality_hdi[1, ],
                                       upper = centrality_hdi[2, ])
    if(i == 1){
      centrality_summary <- centrality_summary_i
    } else {
      centrality_summary <- rbind(centrality_summary, centrality_summary_i)
    }
  }
  
  ordering <- centrality_summary %>%
    filter(group == group_i) %>%
    arrange(mean) %>%
    pull(node)
  
  centrality_summary |>
    mutate(node = factor(.data$node, levels = ordering)) |>
    ggplot2::ggplot(aes(x = .data$node, y=.data$mean, color =.data$group, ...))+
    args$theme_ +
    geom_point(position = position_dodge(width = 0.6))+
    args$geom_errorbar +
    coord_flip() +
    ylab(args$ylab) +
    xlab(args$xlab) +
    labs(color = "Group")+
    args$theme
}

# -------------------------------------------------------------------------------
#' @export
plot_prior_sensitivity.list <- function(output,
                                        evidence_thres = 10, ...) {
  
  stop("This plot is currently not supported for easybgm_longitudinal objects.")
  
  # if(any(class(output[[1]]) == "easybgm_compare")){
  #   stop("The centrality plot cannot be obtained for Bayesian network comparison fits.")
  # }
  # 
  # # Convert to easybgm output if provided objects are bgms
  # 
  # # Check for bgms package version 
  # if(any(class(output[[1]]) == "bgms")) {
  #   if(packageVersion("bgms") < "0.1.3"){
  #     stop("Your version of the package bgms is not supported anymore. Please update.")
  #   }
  #   
  #   res <- list()
  #   for(i in 1:length(output)) {
  #     fit_args <- bgms::extract_arguments(output[[i]])
  #     
  #     res[[i]] <- bgm_extract.package_bgms(fit = output[[i]], save = fit_args$save, centrality = TRUE,
  #                                          type = NULL, not_cont = NULL, data = NULL,
  #                                          edge_prior = fit_args$edge_prior,
  #                                          inclusion_probability  = fit_args$inclusion_probability,
  #                                          beta_bernoulli_alpha = fit_args$beta_bernoulli_alpha,
  #                                          beta_bernoulli_beta = fit_args$beta_bernoulli_beta)
  #     
  #   }
  #   output <- res
  # }
  # 
  # default_args <- list(
  #   theme_ = theme_minimal(),
  #   ylab = ylab("Relative no. edges"),
  #   xlab = xlab("Prior edge probability"),
  #   xlim = xlim(0,1),
  #   ylim = ylim(0,1),
  #   theme = theme(
  #     axis.text = element_text(size=16),
  #     panel.border = element_blank(),
  #     axis.line = element_line(colour = "black", linewidth = 1.1),
  #     axis.ticks.length = unit(.2, "cm"),
  #     axis.ticks = element_line(linewidth = .8),
  #     axis.title.x = element_text(size = 18, face = "bold"),
  #     axis.title.y = element_text(size = 18, face = "bold"),
  #     plot.title = element_text(size = 18, face = "bold"),
  #     panel.grid.major = element_blank(),
  #     legend.text = element_text(size = 12),
  #     
  #   ),
  #   colors = c("#36648b", "#eeb004", "#bfbfbf"),
  #   size = 1
  # )
  # 
  # args <- set_defaults(default_args, ...)
  # no_priors <- length(output)
  # edge_priors <- rep(NA, no_priors)
  # 
  # for (i in 1:no_priors) {
  #   if (!(any(class(output[[i]]) == "easybgm") | any(class(output[[i]]) == "package_bgms"))) {
  #     stop(paste0("Wrong input provided on index ", i, ". The function requires
  #                 as input a list of output values of the easybgm or the bgms function."))
  #   }
  # }
  # 
  # incl_edges = excl_edges = inconcl_edges <- rep(NA, no_priors)
  # 
  # for (i in 1:no_priors) {
  #   res <- output[[i]]
  #   
  #   if (is.null(res$edge.prior)) {
  #     stop("Wrong input provided. At least one of the output of the easybgm
  #          function does not include a specification of the edge prior. Please note that this 
  #          plot cannot be obtained with the package BGGM.")
  #   }
  #   edge_priors[i] <- res$edge.prior
  #   
  #   
  #   incl_bf <- res$inc_BF
  #   incl_bf <- incl_bf[lower.tri(incl_bf)]
  #   
  #   k <- length(incl_bf)
  #   
  #   incl_edges[i] <- length(which(incl_bf > evidence_thres)) / k
  #   excl_edges[i] <- length(which(incl_bf < (1 / evidence_thres))) / k
  #   
  # }
  # 
  # inconcl_edges <- 1 - incl_edges - excl_edges
  # 
  # data <- data.frame(cbind(edge_priors, incl_edges, excl_edges, inconcl_edges))
  # 
  # ggplot2::ggplot(data, aes(x = edge_priors, ...)) +
  #   geom_line(aes(y = incl_edges, color = "included"), size = args$size) + 
  #   geom_point(aes(y = incl_edges, color = "included"), size = args$size+ 0.5) +
  #   geom_line(aes(y = excl_edges, color = "excluded"), size = args$size) +
  #   geom_point(aes(y = excl_edges, color = "excluded"), size = args$size + 0.5) +
  #   geom_line(aes(y = inconcl_edges, color = "inconclusive"), size = args$size) +
  #   geom_point(aes(y = inconcl_edges, color = "inconclusive"), size = args$size + 0.5) +
  #   args$theme_ + args$xlab + args$ylab + scale_color_manual(values = args$colors, name = "")  +
  #   args$xlim + args$ylim
  # 
  # 
}

