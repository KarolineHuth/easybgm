

#' Plot Posterior Structure Probabilities
#'
#' @param output Output object from the bgm_extract function
#' @param as.BF if TRUE plots the y-axis as Bayes factor instead of posterior structure probability
#'
#' @export
#' @import ggplot2
#'

plot_structure_probability <- function(output, as.BF = FALSE) {
  if (output$package == "rbinnet") {
    stop("The plot cannot be obtained for rbinnet.",
         call. = FALSE)
  }
  if (output$package == "BGGM") {
    stop("The plot cannot be obtained for BGGM.",
         call. = FALSE)
  }

  sorted_structure_prob <- as.data.frame(sort(output$structure_probabilities, decreasing=T))
  colnames(sorted_structure_prob) <- "posterior_prob"
  if(as.BF){
    BF1s <- sorted_structure_prob$posterior_prob[1] / sorted_structure_prob$posterior_prob # BF best structure vs. others
    data <- data.frame(structures = 1:length(BF1s), BayesFactor = BF1s)
    ggplot2::ggplot(data, aes(x = structures, y = BayesFactor)) +
      geom_point(size = 4, shape = 1) +
      scale_y_continuous(trans = "log10") +
      theme_classic()+
      labs(x = "Structures",
           y = expression(log(BF[1][s])))+
      geom_hline(yintercept = 1/10, linetype = "dashed", size = 1.5)  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.1),
            axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            axis.ticks.length = unit(.25, "cm"))
  } else {
    data <- data.frame(structures = 1:nrow(sorted_structure_prob), Probs = sorted_structure_prob)
    ggplot2::ggplot(data, aes(x = structures, y = posterior_prob)) +
      geom_point(size = 4, shape = 1) +
      theme_classic()+
      labs(x = "Structures",
           y = "Posterior Structure Probability")+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", linewidth = 1.1),
            axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            axis.ticks.length = unit(.25, "cm"))
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

plot_posteriorcomplexity <- function(output) {
  if (output$package == "rbinnet") {
    stop("The plot cannot be obtained for rbinnet.",
         call. = FALSE)
  }
  if (output$package == "BGGM") {
    stop("The plot cannot be obtained for BGGM.",
         call. = FALSE)
  }
  complexity <- c()
  for(i in 1:length(output$sample_graph)){
    complexity[i] <- sum(as.numeric(unlist(strsplit(res$sample_graph[i], ""))))
  }

  data_complexity <- tibble(complexity, weights = res$graph_weights)  %>%
    group_by(complexity) %>%
    summarise(complexity_weight = sum(weights)) %>%
    mutate(complexity_weight = complexity_weight/sum(complexity_weight))

  ggplot(data_complexity, aes(x = complexity, y = complexity_weight))+
    geom_point() +
    ylab("Posterior Complexity Probability") +
    xlab("Complexity")  +
    theme_minimal()+
    theme(legend.position = c(.85, 0.25), axis.text=element_text(size=20),
          legend.background = element_rect(fill = NULL), panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.2, "cm"),
          axis.ticks = element_line(size= .8), legend.text = element_text(size=14),
          axis.title.x = element_text(size=18,face="bold"),
          axis.title.y = element_text(size=18,face="bold"),
          text=element_text(  family="Times New Roman"),
          panel.grid.major = element_blank()
    )
}

# ---------------------------------------------------------------------------------------------------------------

#' Edge evidence plot
#'
#' @param output Output object from the bgm_extract function
#' @param evidence_thresh BF which will be considered sufficient evidence for in-/exclusion
#' @param layout Layout of the network; qgraph argument
#' @param edge.width Layout of the network; qgraph argument
#' @param split if TRUE, plot is split in included and excluded edges
#' @param ... Additional `qgraph` arguments
#'
#' @export
#' @import qgraph
#'
plot_edgeevidence <- function(output, evidence_thresh = 10, split = F, ...) {

  if(output$model == "dgm-binary"){
    stop("Plot cannot be obtained for 'dgm-binary' models. Use the package rbinnet instead to obtain parameter estimates for the Ising model.",
         call. = FALSE)
  }
  graph <- output$BF
  diag(graph) <- 1

  # assign a color to each edge (inclusion - blue, exclusion - red, no conclusion - grey)
  graph_color <- graph
  graph_color <-  ifelse(graph < evidence_thresh & graph > 1/evidence_thresh, graph_color <- "#bfbfbf", graph_color <- "#36648b")
  graph_color[graph < (1/evidence_thresh)] <- "#990000"

  if(split == F){
    graph[output$inc_probs <= 1] <- 1
    diag(graph) <- 1
    colnames(graph) <- colnames(output$sigma)
    qgraph::qgraph(graph,
                   edge.color = graph_color, # specifies the color of the edges
                   ...
    )
  }
  if(split==T){
    par(mfrow=c(1, 2))
    graph_inc <- graph_exc <- graph
    # plot included graph
    graph_inc[output$inc_probs >= .5] <- 1
    graph_inc[output$inc_probs < .5] <- 0
    diag(graph_inc) <- 1
    colnames(graph_inc) <- colnames(output$sigma)
    qgraph::qgraph(graph_inc,
                   edge.color = graph_color, # specifies the color of the edges
                   ...
    )
    # Plot excluded graph
    graph_exc[output$inc_probs >= .5] <- 0
    graph_exc[output$inc_probs < .5] <- 1
    diag(graph_exc) <- 1
    colnames(graph_exc) <- colnames(output$sigma)
    qgraph::qgraph(graph_exc,
                   edge.color = graph_color, # specifies the color of the edges
                   ...
    )
    par(mfrow=c(1, 1))
  }

}

# ---------------------------------------------------------------------------------------------------------------

#' Network plot
#'
#' @param output Output object from the bgm_extract function
#' @param exc_prob threshold for excluding edges; all edges with a lower inclusion probability will not be shown
#' @param ... Additional `qgraph` arguments
#'
#' @export
#' @import qgraph

plot_network <- function(output, exc_prob = .5, ...) {
  if(output$model == "dgm-binary"){
    stop("Plot cannot be obtained for 'dgm-binary' models. Use the package rbinnet instead to obtain parameter estimates.",
         call. = FALSE)
  }

  graph <- output$sigma

  # Exclude edges with a inclusion probability lower exc_prob
  inc_probs_m <- output$inc_probs
  graph[inc_probs_m < exc_prob] <- 0
  diag(graph) <- 1

  # Plot
  qgraph::qgraph(graph, ...)

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
plot_parameterHDI <- function(output) {

  package <- output$package
  if(output$package == "BDgraph" & output$model == "gcgm"){
    stop("Plot cannot be obtained for GCGMs.")
  }
  if(is.null(output$samples_posterior)){
    stop("Samples of the posterior distribution required. When extracting the results, set \"posterior_samples = TRUE\".")
  }

  hdi_intervals <- as.data.frame(apply(output$samples_posterior, MARGIN = 2, FUN = hdi))
  posterior_medians <- apply(output$samples_posterior, MARGIN = 2, FUN = median)

  names <- colnames(output$sigma)
  names_bycol <- matrix(rep(names, each = ncol(output$sigma)), ncol = ncol(output$sigma))
  names_byrow <- matrix(rep(names, each = ncol(output$sigma)), ncol = ncol(output$sigma), byrow = T)
  names_comb <- matrix(paste0(names_byrow, "-", names_bycol), ncol = ncol(output$sigma))
  index <- names_comb[upper.tri(names_comb)]

  posterior <- cbind(data.frame(posterior_medians, row.names = NULL),
                     data.frame(t(hdi_intervals), row.names = NULL), index)
  colnames(posterior) <- c("posterior_medians", "lower", "upper", "names")
  posterior <- posterior[order(posterior$posterior_medians, decreasing = F),]
  posterior$names <- factor(posterior$names, levels = posterior$names)


  ggplot2::ggplot(data = posterior, aes(x = names, y = posterior_medians, ymin = lower,
                                        ymax = upper)) +
    geom_pointrange(position=position_dodge(width=c(0.3)), size = .3) +
    theme_bw() +
    coord_flip() +
    ylab("Highest Density Interval of Parameter")+
    xlab("") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1.3) +
    theme(axis.text=element_text(size=8), panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.2, "cm"),
          axis.ticks = element_line(size= .8),
          axis.title.x = element_text(size=16,face="bold"), plot.title = element_text(size = 18, face = "bold"))
}


# ---------------------------------------------------------------------------------------------------------------
# Centrality plot

#' Plot centrality measures and 95% highest density interval
#'
#' @param output Output object from the bgm_extract function
#' @param measure Centrality measures that should be plotted. Users can choose "all" or a subsection of the list: "Strength", "Closeness", "Betweenness", or "ExpectedInfluence"
#'
#' @export
#'

plot_centrality <- function(output, measure = "Strength"){
  cent_samples <- output$centrality
  p <- ncol(output$sigma)
  rownames(cent_samples) <- NULL
  # Creating summary statistics

  centrality_means <- cent_samples %>%
    as_tibble() %>%
    group_by(Centrality) %>%
    group_modify(~ as.data.frame(colMeans(.x)))
  centrality_means <- cbind(centrality_means, rep(colnames(output$sigma), 4))
  colnames(centrality_means)[2:3] <- c("value", "node")
  centrality_means <- centrality_means[order(centrality_means$Centrality, centrality_means$node), ]
  centrality_hdi <- cent_samples %>%
    as_tibble() %>%
    group_by(Centrality) %>%
    group_modify(~ as.data.frame(hdi(.x, allowSplit = F)))
  centrality_hdi <- centrality_hdi %>%
    gather(node, value, stress:uncreative) %>%
    add_column(interval = rep(c("lower", "upper"), p*4)) %>%
    spread(interval, value)

  centrality_summary <- merge(centrality_hdi, centrality_means, all = T)

  measure_options <- c("all", "Betweenness", "Closeness", "ExpectedInfluence", "Strength")
  if((measure %in% measure_options) == FALSE) {
    stop("This centrality measure cannot be plotted. Please choose one or several of the following measures: Betweenness, Closeness, ExpectedInfluence, Strength.")
  }
  if(measure == "all"){
    measure <- c("Betweenness", "Closeness", "ExpectedInfluence", "Strength")
  }
  centrality_summary %>%
    filter(Centrality %in% measure) %>%
    ggplot(aes(x = node, y=value, group = Centrality))+
    geom_line()+
    geom_point()+
    geom_errorbar(aes(y= value, ymin =lower, ymax = upper), size = .5, width = 0.4)+
    facet_wrap(.~ Centrality, ncol = 4, scales = "free_x") +
    coord_flip() +
    ylab("Value") +
    xlab("Nodes")
}

