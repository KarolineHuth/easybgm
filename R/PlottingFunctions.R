

#' Plot Posterior Structure Probabilities
#'
#' @param output Output object from the bgm_extract function
#' @param as.BF if TRUE plots the y-axis as Bayes factor instead of posterior structure probability
#'
#' @export
#' @import ggplot2
#' @import dplyr
#'

plot_posteriorstructure <- function(output, as.BF = FALSE, ...) {
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  if (any(class(output) == "package_bggm")) {
    stop("The plot cannot be obtained for models estimated with BGGM. Suggestion: Change the package to BDgraph.",
         call. = FALSE)
  }
  
  default_args <- list(
    xlab = "Structures",
    ylab = ifelse(as.BF == TRUE, expression(log(BF[1][s])), "Posterior Structure Probability"),
    theme = theme_classic(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 1.1),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14), 
    axis.ticks.length = unit(.25, "cm")
  )
  
  args <- set_defaults(default_args, ...)
  sorted_structure_prob <- as.data.frame(sort(output$structure_probabilities, decreasing=T))
  colnames(sorted_structure_prob) <- "posterior_prob"
  if(as.BF){
    
    BF1s <- sorted_structure_prob$posterior_prob[1] / sorted_structure_prob$posterior_prob # BF best structure vs. others
    data <- data.frame(structures = 1:length(BF1s), BayesFactor = BF1s)
    ggplot2::ggplot(data, aes(x = structures, y = BayesFactor)) +
      geom_point(size = 4, shape = 1, ...) +
      scale_y_continuous(trans = "log10") +
      args$theme +
      labs(x = args$xlab,
           y = args$ylab) +
      theme(panel.grid.major = args$panel.grid.major,
            panel.grid.minor = args$panel.grid.minor, axis.line = args$axis.line,
            axis.text = args$axis.text, axis.title = args$axis.title,
            axis.ticks.length = args$axis.ticks.length)
  } else {
    data <- data.frame(structures = 1:nrow(sorted_structure_prob), Probs = sorted_structure_prob)
    ggplot2::ggplot(data, aes(x = structures, y = posterior_prob, ...)) +
      geom_point(size = 4, shape = 1) +
      args$theme +
      labs(x = args$xlab,
           y = args$ylab) +
      theme(panel.grid.major = args$panel.grid.major,
            panel.grid.minor = args$panel.grid.minor, axis.line = args$axis.line,
            axis.text = args$axis.text, axis.title = args$axis.text,
            axis.ticks.length = args$axis.ticks.length)
  }
}

# ---------------------------------------------------------------------------------------------------------------


#' Plot posterior structure complexity
#'
#' @param output Output object from the bgm_extract function
#'
#' @export
#' @import ggplot2
#'

plot_posteriorcomplexity <- function(output, ...) {
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  if (any(class(output) == "package_bggm")) {
    stop("The plot cannot be obtained for models estimated with BGGM. Suggestion: Change the package to BDgraph.",
         call. = FALSE)
  }
  
  default_args <- list(
    xlab = "Complexity",
    ylab = "Posterior Complexity Probability",
    theme = theme_minimal(),
    legend.position = c(.85, .25),
    axis.text.size = 20,
    legend.background = element_rect(fill = NULL),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1.1),
    axis.ticks = element_line(size=.8),
    legend.text = element_text(size = 14),
    
    axis.ticks.length = unit(.2, "cm"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    panel.grid.major = element_blank()
  )
  
  args <- set_defaults(default_args, ...)
  complexity <- c()
  for(i in 1:length(output$sample_graph)){
    complexity[i] <- sum(as.numeric(unlist(strsplit(output$sample_graph[i], ""))))
  }

  data_complexity <- tibble(complexity, weights = output$graph_weights)  %>%
    group_by(complexity) %>%
    summarise(complexity_weight = sum(weights)) %>%
    mutate(complexity_weight = complexity_weight/sum(complexity_weight))

  ggplot(data_complexity, aes(x = complexity, y = complexity_weight, ...)) +
    geom_point() +
    ylab(args$ylab) +
    xlab(args$xlab)  +
    args$theme +
    theme(legend.position = args$legend.position, axis.text=element_text(size=args$axis.text.size),
          legend.background =  args$legend.background, panel.border = args$panel.border,
          axis.line = args$axis.line, axis.ticks.length=args$axis.ticks.length,
          axis.ticks = args$axis.ticks, legend.text = args$legend.text,
          axis.title.x = args$axis.title.x,
          axis.title.y = args$axis.title.y,
          panel.grid.major = args$panel.grid.major
    )
}

# ---------------------------------------------------------------------------------------------------------------

#' Edge evidence plot
#'
#' @param output Output object from the bgm_extract function
#' @param evidence_thresh BF which will be considered sufficient evidence for in-/exclusion
#' @param split if TRUE, plot is split in included and excluded edges
#' @param show specifies which edges should be shown, indicated by "all", "included", "inconclusive", "excluded"
#' @param ... Additional `qgraph` arguments
#'
#' @export
#' @import qgraph
#'
plot_edgeevidence <- function(output, evidence_thresh = 10, split = F, show = "all", ...) {
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  
  default_args <- list(
    colors = c("#36648b", "#990000", "#bfbfbf"),
    colnames = colnames(output$parameters),
    layout_avg = qgraph::averageLayout(output$parameters*output$structure),
    theme = "TeamFortress",
    legend = T,
    vsize = 10,
    nodeNames = colnames(output$parameters),
    labels = colnames(output$parameters),
    edge.width = 3,
    label.cex = 1,
    legend.cex = 1
    
  )
  args <- set_defaults(default_args, ...)
  graph <- output$BF
  diag(graph) <- 1
  
  # assign a color to each edge (inclusion - blue, exclusion - red, no conclusion - grey)
  graph_color <- graph
  graph_color <-  ifelse(graph < evidence_thresh & graph > 1/evidence_thresh, 
                         graph_color <- args$colors[3], graph_color <- args$colors[1])
  graph_color[graph < (1/evidence_thresh)] <- args$colors[2]

  if (show == "all") {
    if (split == F) {
      graph[output$inc_probs <= 1] <- 1
      diag(graph) <- 1
      colnames(graph) <- args$colnames
      qgraph::qgraph(graph,
                     edge.color = graph_color,
                     layout = args$layout_avg,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
                     labels = args$labels,
                     edge.width = args$edge.width,
                     label.cex = args$label.cex,
                     legend.cex = args$legend.cex,
                     ...
      )
    }
    if (split == T) {
      graph_inc <- graph_exc <- graph
      # plot included graph
      graph_inc[output$inc_probs >= .5] <- 1
      graph_inc[output$inc_probs < .5] <- 0
      diag(graph_inc) <- 1
      colnames(graph_inc) <- colnames(output$parameters)
      qgraph::qgraph(graph_inc,
                     edge_color = graph_color,
                     layout = args$layout_avg,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
                     labels = args$labels,
                     edge.width = args$edge.width,
                     label.cex = args$label.cex,
                     legend.cex = args$legend.cex, # specifies the color of the edges
                     ...
      )
      # Plot excluded graph
      graph_exc[output$inc_probs >= .5] <- 0
      graph_exc[output$inc_probs < .5] <- 1
      diag(graph_exc) <- 1
      colnames(graph_exc) <- colnames(output$parameters)
      qgraph::qgraph(graph_exc,
                     edge.color = graph_color,
                     # specifies the color of the edges
                     layout = args$layout_avg,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
                     labels = args$labels,
                     edge.width = args$edge.width,
                     label.cex = args$label.cex,
                     legend.cex = args$legend.cex,
                     ...
      )
    }
  }
    if(show != "all"){
      graph_show <- matrix(0, ncol = ncol(graph), nrow = nrow(graph))
      if("included" %in% show){
        graph_show[output$BF > evidence_thresh] <- 1
      }
      if("excluded" %in% show){
        graph_show[output$BF < (1/evidence_thresh)] <- 1
      }
      if("inconclusive" %in% show){
        graph_show[(output$BF > (1/evidence_thresh)) & (output$BF < evidence_thresh)] <- 1
      }
      diag(graph_show) <- 1
      colnames(graph_show) <- colnames(output$parameters)
      qgraph::qgraph(graph_show,
                     edge.color = graph_color,
                     layout = args$layout_avg,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
                     labels = args$labels,
                     label.cex = args$label.cex,
                     legend.cex = args$legend.cex,# specifies the color of the edges
                     ...
      )
    }

}

# ---------------------------------------------------------------------------------------------------------------

#' Network plot
#'
#' @param output Output object from the bgm_extract function
#' @param exc_prob threshold for excluding edges; all edges with a lower inclusion probability will not be shown
#' @param dashed binary parameter indicating whether edges with inconclusive evidence should be dashed
#' @param ... Additional `qgraph` arguments

#'
#' @export
#' @import qgraph

plot_network <- function(output, exc_prob = .5, dashed = F, ...) {
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  if(output$model == "dgm-binary"){
    stop("Plot cannot be obtained for 'dgm-binary' models. Use the package bgms instead to obtain parameter estimates.",
         call. = FALSE)
  }
  graph <- output$parameters
  default_args <- list(
    evidence_thres = 10
  )
  args <- set_defaults(default_args, ...)
  # Exclude edges with a inclusion probability lower exc_prob
  inc_probs_m <- output$inc_probs
  graph[inc_probs_m < exc_prob] <- 0
  diag(graph) <- 1

  # Plot
  if(dashed == T){
    graph_dashed <- ifelse(output$BF < evidence_thres, "dashed", "solid")
    qgraph::qgraph(graph, lty = graph_dashed, ...)
  } else {
    qgraph::qgraph(graph, ...)
  }

}

# -------------------------------------------------

#' Structure plot
#'
#' @param output Output object from the bgm_extract function
#' @param ... Additional `qgraph` arguments
#'
#' @export
#'
#' @import qgraph
#'

plot_structure <- function(output, ...) {

  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  graph <- output$structure

  # Plot
  qgraph::qgraph(graph, ...)

}

# ---------------------------------------------------------------------------------------------------------------

#' Plot of interaction parameters and their 95% highest density intervals
#'
#' @param output Output object from the bgm_extract function
#'
#' @export
#' @import ggplot2 HDInterval
#'
plot_parameterHDI <- function(output, ...) {
  
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }

  if(is.null(output$samples_posterior)){
    stop("Samples of the posterior distribution required. When estimating the model, set \"save = TRUE\".")
  }
  
  def_args <- list(
    theme_ = theme_bw(),
    geom_pointrange = geom_pointrange(position = position_dodge(width = c(0.3)),
                                      size = .3),
    xlab = "",
    ylab = "Highest Density Interval of Parameter",
    geom_hline = geom_hline(yintercept = 0, linetype = "dashed", size = 1.3),
    theme = theme(
    axis.text = element_text(size=8),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1.1),
    axis.ticks.length = unit(.2, "cm"),
    axis.ticks = element_line(size = .8),
    axis.title.x = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold")
    )
  )
  
  args <- set_defaults(def_args, ...)
  hdi_intervals <- as.data.frame(apply(output$samples_posterior, MARGIN = 2, FUN = hdi))
  posterior_medians <- apply(output$samples_posterior, MARGIN = 2, FUN = median)

  names <- colnames(output$parameters)
  names_bycol <- matrix(rep(names, each = ncol(output$parameters)), ncol = ncol(output$parameters))
  names_byrow <- matrix(rep(names, each = ncol(output$parameters)), ncol = ncol(output$parameters), byrow = T)
  names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = ncol(output$parameters))
  index <- names_comb[upper.tri(names_comb)]

  posterior <- cbind(data.frame(posterior_medians, row.names = NULL),
                     data.frame(t(hdi_intervals), row.names = NULL), index)
  colnames(posterior) <- c("posterior_medians", "lower", "upper", "names")
  posterior <- posterior[order(posterior$posterior_medians, decreasing = F),]
  posterior$names <- factor(posterior$names, levels = posterior$names)


  ggplot2::ggplot(data = posterior, aes(x = names, y = posterior_medians, ymin = lower,
                                        ymax = upper, ...)) +
    args$geom_pointrange +
    args$theme_ +
    coord_flip() +
    ylab(args$ylab)+
    xlab(args$xlab) +
    args$geom_hline +
    args$theme
}


# ---------------------------------------------------------------------------------------------------------------
# Centrality plot

#' Plot centrality measures and 95% highest density interval
#'
#' @param output Output object from the bgm_extract function
#' @param measure Centrality measures that should be plotted. Users can choose "all" or a
#'        subselection of the list: "Strength", "Closeness", "Betweenness", or "ExpectedInfluence"
#'
#' @export
#' @import tibble
#' @import tidyr
#'

plot_centrality <- function(output, measure = "Strength", ... ){

  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }

  if(is.null(output$centrality)){
    stop("Centrality results are required. When estimating the model, set \"centrality = TRUE\".")
  }
  
  default_args <- list(
    ylab = "Value",
    xlab = "Nodes",
    geom_errorbar = geom_errorbar(aes(y=value, ymin = lower, ymax = upper)
                                  , size =.5, width=.4)
    
  )
  
  args <- set_defaults(default_args, ...)
  cent_samples <- output$centrality
  p <- ncol(output$parameters)
  rownames(cent_samples) <- NULL
  # Creating summary statistics

  centrality_means <- cent_samples %>%
    as_tibble() %>%
    group_by(Centrality) %>%
    group_modify(~ as.data.frame(colMeans(.x)))
  centrality_means <- cbind(centrality_means, rep(colnames(output$parameters), 4))
  colnames(centrality_means)[2:3] <- c("value", "node")
  centrality_means <- centrality_means[order(centrality_means$Centrality, centrality_means$node), ]
  centrality_hdi <- cent_samples %>%
    as_tibble() %>%
    group_by(Centrality) %>%
    group_modify(~ as.data.frame(hdi(.x, allowSplit = F)))
  firstnode <- colnames(output$parameters)[1]
  lastnode <- colnames(output$parameters)[nrow(output$parameters)]
  centrality_hdi <- centrality_hdi %>%
    gather(node, value, firstnode:lastnode) %>%
    tibble::add_column(interval = rep(c("lower", "upper"), p*4)) %>%
    tidyr::spread(interval, value)

  centrality_summary <- merge(centrality_hdi, centrality_means, all = T)

  measure_options <- c("all", "Betweenness", "Closeness", "ExpectedInfluence", "Strength")
  if(any((measure %in% measure_options) == FALSE)) {
    stop("This centrality measure cannot be plotted. Please choose one or several of the following measures: Betweenness, Closeness, ExpectedInfluence, Strength.")
  }
  if(any(measure == "all")){
    measure <- c("Betweenness", "Closeness", "ExpectedInfluence", "Strength")
  }
  centrality_summary %>%
    filter(Centrality %in% measure) %>%
    ggplot(aes(x = node, y = value, group= Centrality, ...))+
    geom_line()+
    geom_point()+
    args$geom_errorbar+
    facet_wrap(.~ Centrality, ncol = 4, scales = "free_x") +
    coord_flip() +
    ylab(args$ylab) +
    xlab(args$xlab)
}

