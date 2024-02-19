#' @export
plot_structure_probabilities.easybgm <- function(output, as_BF = FALSE, ...) {
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  if (any(class(output) == "package_bggm")) {
    stop("The plot cannot be obtained for models estimated with BGGM. Suggestion: Change the package to BDgraph.",
         call. = FALSE)
  }

  default_args <- list(
    xlab = "Structures",
    ylab = ifelse(as_BF == TRUE, expression(log(BF[1][s])), "Posterior Structure Probability"),
    theme = theme_minimal(),
    axis.text.size = 16,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.1),
    axis.ticks = element_line(linewidth = .8),
    axis.ticks.length = unit(.2, "cm"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    panel.grid.major = element_blank()
  )

  args <- set_defaults(default_args, ...)
  sorted_structure_prob <- as.data.frame(sort(output$structure_probabilities, decreasing=T))
  colnames(sorted_structure_prob) <- "posterior_prob"
  if(as_BF){

    BF1s <- sorted_structure_prob$posterior_prob[1] / sorted_structure_prob$posterior_prob # BF best structure vs. others
    data <- data.frame(structures = 1:length(BF1s), BayesFactor = BF1s)
    ggplot2::ggplot(data, aes(x = .data$structures, y = .data$BayesFactor, ...)) +
      geom_point(size = 3) +
      scale_y_continuous(trans = "log10") +
      args$theme +
      labs(x = args$xlab,
           y = args$ylab) +
      theme(legend.position = args$legend.position, axis.text=element_text(size=args$axis.text.size),
            legend.background =  args$legend.background, panel.border = args$panel.border,
            axis.line = args$axis.line, axis.ticks.length=args$axis.ticks.length,
            axis.ticks = args$axis.ticks, legend.text = args$legend.text,
            axis.title.x = args$axis.title.x,
            axis.title.y = args$axis.title.y,
            panel.grid.major = args$panel.grid.major)
  } else {
    data <- data.frame(structures = 1:nrow(sorted_structure_prob), Probs = sorted_structure_prob)
    ggplot2::ggplot(data, aes(x = .data$structures, y = .data$posterior_prob, ...)) +
      geom_point(size = 3) +
      args$theme +
      labs(x = args$xlab,
           y = args$ylab) +
      theme(legend.position = args$legend.position, axis.text=element_text(size=args$axis.text.size),
            legend.background =  args$legend.background, panel.border = args$panel.border,
            axis.line = args$axis.line, axis.ticks.length=args$axis.ticks.length,
            axis.ticks = args$axis.ticks, legend.text = args$legend.text,
            axis.title.x = args$axis.title.x,
            axis.title.y = args$axis.title.y,
            panel.grid.major = args$panel.grid.major)
  }
}

# ---------------------------------------------------------------------------------------------------------------
#' @export
plot_complexity_probabilities.easybgm <- function(output, ...) {
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
    axis.text.size = 16,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.1),
    axis.ticks = element_line(linewidth = .8),
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

  data_complexity <- data.frame(complexity = complexity, weights = output$graph_weights) |>
    dplyr::group_by(complexity) |>
    dplyr::summarise(complexity_weight = sum(.data$weights)) |>
    dplyr::mutate(complexity_weight = .data$complexity_weight/sum(.data$complexity_weight))

  ggplot(data_complexity, aes(x = .data$complexity, y = .data$complexity_weight, ...)) +
    geom_point(size = 3) +
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
#' @export
plot_edgeevidence.easybgm <- function(output, evidence_thresh = 10, split = FALSE, show = "all", ...) {
if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }

  default_args <- list(
    colors = c("#36648b", "#990000", "#bfbfbf"),
    colnames = colnames(output$parameters),
    layout = qgraph::averageLayout(output$parameters*output$structure),
    theme = "TeamFortress",
    legend = TRUE,
    vsize = 10,
    nodeNames = colnames(output$parameters),
    edge.width = 3,
    label.cex = 1,
    legend.cex = .6
  )
  args <- set_defaults(default_args, ...)
  graph <- output$inc_BF
  diag(graph) <- 1

  # assign a color to each edge (inclusion - blue, exclusion - red, no conclusion - grey)
  graph_color <- graph
  graph_color <-  ifelse(graph < evidence_thresh & graph > 1/evidence_thresh,
                         graph_color <- args$colors[3], graph_color <- args$colors[1])
  graph_color[graph < (1/evidence_thresh)] <- args$colors[2]

  if (show == "all") {
    if (!split) {
      graph[output$inc_probs <= 1] <- 1
      diag(graph) <- 1
      colnames(graph) <- args$colnames
      qgraph_plot <- qgraph::qgraph(graph,
                     edge.color = graph_color,
                     layout = args$layout,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
                     edge.width = args$edge.width,
                     label.cex = args$label.cex,
                     legend.cex = args$legend.cex,
                     ...
      )
    }

    if (split) {

      graph_inc <- graph_exc <- graph
      # plot included graph
      graph_inc[output$inc_probs >= .5] <- 1
      graph_inc[output$inc_probs < .5] <- 0
      diag(graph_inc) <- 1
      colnames(graph_inc) <- colnames(output$parameters)
      qgraph_plot1 <- qgraph::qgraph(graph_inc,
                     edge.color = graph_color,
                     layout = args$layout,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
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
      qgraph_plot2 <- qgraph::qgraph(graph_exc,
                     edge.color = graph_color,
                     # specifies the color of the edges
                     layout = args$layout,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
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
        graph_show[output$inc_BF > evidence_thresh] <- 1
      }
      if("excluded" %in% show){
        graph_show[output$inc_BF < (1/evidence_thresh)] <- 1
      }
      if("inconclusive" %in% show){
        graph_show[(output$inc_BF > (1/evidence_thresh)) & (output$BF < evidence_thresh)] <- 1
      }
      diag(graph_show) <- 1
      colnames(graph_show) <- colnames(output$parameters)
      qgraph_plot <- qgraph::qgraph(graph_show,
                     edge.color = graph_color,
                     layout = args$layout,# specifies the color of the edges
                     theme = args$theme,
                     vsize = args$vsize,
                     nodeNames = args$nodeNames,
                     legend = args$legend,
                     label.cex = args$label.cex,
                     legend.cex = args$legend.cex,# specifies the color of the edges
                     ...
      )
    }
  if (split == TRUE) {
    return(invisible(list(qgraph_plot1, qgraph_plot2)))
  } else {
    return(invisible(qgraph_plot))
  }
}

# ---------------------------------------------------------------------------------------------------------------


#' @export
plot_network.easybgm <- function(output, exc_prob = 0.5, evidence_thresh = 10,  dashed = TRUE, ...) {

  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  if(output$model == "dgm-binary"){
    stop("Plot cannot be obtained for 'dgm-binary' models. Use the package bgms instead to obtain parameter estimates.",
         call. = FALSE)
  }
  graph <- output$parameters
  default_args <- list(
    layout = qgraph::averageLayout(output$parameters*output$structure),
    evidence_thresh = 10,
    theme = "TeamFortress",
    vsize = 10,
    nodeNames = colnames(output$parameters),
    legend = T,
    label.cex = 1.2,
    legend.cex = .6, 
    edge.labels = TRUE
  )
  args <- set_defaults(default_args, ...)

  # Exclude edges with an inclusion probability lower than exc_prob
  inc_probs_m <- output$inc_probs
  graph[inc_probs_m < exc_prob] <- 0
  diag(graph) <- 1

  # Plot
  if(dashed){
    graph_dashed <- ifelse(output$inc_BF < args$evidence_thresh, 2, 1)

    qgraph_plot <- qgraph::qgraph(graph, layout = args$layout, lty = graph_dashed,
                   theme = args$theme, vsize = args$vsize,
                   nodeNames = args$nodeNames,
                   legend = args$legend,
                   label.cex = args$label.cex,
                   legend.cex = args$legend.cex, edge.labels = args$edge.labels, ...)
  } else {
    qgraph_plot <- qgraph::qgraph(graph, theme = args$theme, layout = args$layout, vsize = args$vsize,
                   nodeNames = args$nodeNames,
                   legend = args$legend,
                   label.cex = args$label.cex,
                   legend.cex = args$legend.cex, edge.labels = args$edge.labels, ...)

  }
  return(invisible(qgraph_plot))
}

# -------------------------------------------------
#' @export
plot_structure.easybgm <- function(output, ...) {
  default_args <- list(
    layout = qgraph::averageLayout(output$parameters*output$structure),
    theme = "TeamFortress",
    vsize = 10,
    nodeNames = colnames(output$parameters),
    legend = T,
    label.cex = 1.2,
    legend.cex = .6
  )
  args <- set_defaults(default_args, ...)
  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }
  graph <- output$structure
  colnames(graph) <- colnames(output$parameters)
  # Plot
  qgraph_plot <- qgraph::qgraph(graph, layout = args$layout,
                 theme = args$theme, vsize = args$vsize,
                 nodeNames = args$nodeNames,
                 legend = args$legend,
                 label.cex = args$label.cex,
                 legend.cex = args$legend.cex, ...)
  
  return(invisible(qgraph_plot))
}

# ---------------------------------------------------------------------------------------------------------------
#' @export
plot_parameterHDI.easybgm <- function(output, ...) {

  if(is.null(output$samples_posterior)){
    stop("Samples of the posterior distribution required. When estimating the model, set \"save = TRUE\".")
  }

  if(any(class(output) == "package_bdgraph")){
    warning("Posterior samples of the interaction parameters are obtained after the estimation. Especially for ordinal and skewed data, these estimates might be sub-optimal. Please interpret with caution.")
  }

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
plot_centrality.easybgm <- function(output, ...){

  if(!any(class(output) == "easybgm")){
    stop("Wrong input provided. The function requires as input the output of the easybgm function.")
  }

  if(is.null(output$centrality)){
    stop("Centrality results are required. When estimating the model, set \"centrality = TRUE\".")
  }

  default_args <- list(
    theme_ = theme_minimal(),
    ylab = "Strength Centrality",
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
  cent_samples <- output$centrality
  p <- ncol(output$parameters)
  rownames(cent_samples) <- NULL
  # Creating summary statistics

  centrality_means <- colMeans(cent_samples)
  centrality_hdi <- apply(cent_samples, MARGIN = 2, FUN = hdi, allowSplit = F)
  centrality_summary <- data.frame(node = colnames(output$parameters),
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

