#' @title Plot Posterior Structure Probabilities
#'
#' @description Plots the posterior structure probabilities of all visited structures, sorted from the most to the least probable.
#'
#' @name structure_probs
#'
#' @param output Output object from the easybgm function. Supports also objects from the bgm function of the `bgms` package.
#' @param as_BF If TRUE plots the y-axis as Bayes factors instead of posterior structure probability. Default is FALSE.
#' @param ... Additional arguments passed onto `ggplot2`
#' 
#' @return Returns a plot
#' 
#' @export
#' @importFrom dplyr group_by summarise mutate group_modify filter
#'
#' @examples
#' \donttest{
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#' fit <- easybgm(data, type = "ordinal", save = TRUE, edge_selection = TRUE,
#'                 iter = 1000  # for demonstration only (> 5e4 recommended)
#'                 )
#'
#' plot_structure_probabilities(fit)
#' }

plot_structure_probabilities <- function(output, as_BF = FALSE, ...) {
  if(any(any(class(output) == "easybgm"), any(class(output) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }

  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }
  UseMethod("plot_structure_probabilities", output)

}

# ---------------------------------------------------------------------------------------------------------------

#'
#' @title Plot posterior complexity probabilities
#'
#' @description Plots the posterior complexity probabilities of all visited structures, where complexity comprises the network density.
#'
#' @name complexity_probs
#'
#' @param output Output object from the easybgm function. Supports also objects from the bgm function of the `bgms` package.
#' @param ... Additional arguments passed onto `ggplot2`
#' 
#' @return Returns a plot
#' 
#' @export
#' @import ggplot2
#'
#' @examples
#' \donttest{
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#' fit <- easybgm(data, type = "ordinal", save = TRUE, edge_selection = TRUE,
#'                 iter = 1000  # for demonstration only (> 5e4 recommended)
#'                 )
#'
#' plot_complexity_probabilities(fit)
#' }

plot_complexity_probabilities <- function(output, ...) {
  if(any(any(class(output) == "easybgm"), any(class(output) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }

  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }
  UseMethod("plot_complexity_probabilities", output)

}

# ---------------------------------------------------------------------------------------------------------------

#' @title Edge evidence plot
#'
#' @description The edge evidence plot colors edges according to their hypothesis testing results: blue for included, red for excluded, and gray for inconclusive. This plot can be used to visualize the hypothesis testing results whether edge presence or absence. The edge evidence plot can aid researchers in deciding which edges provide robust inferential conclusions
#'
#' @name edgeevidence
#'
#' @param output Output object from the easybgm function. Supports also objects from the bgm function of the `bgms` package.
#' @param evidence_thresh Bayes Factor which will be considered sufficient evidence for in-/exclusion, default is 10.
#' @param split if TRUE, plot is split in included and excluded edges. Note that by default separate plots are shown and appear after each other in the plot window. To show the plots side-by-side specify par(mfrow = c(1, 2)).
#' @param show specifies which edges should be shown, indicated by "all", "included", "inconclusive", "excluded".
#' @param ... Additional arguments passed onto `qgraph`.
#'
#' @return Returns a plot
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#' fit <- easybgm(data, type = "continuous",
#'                 iter = 1000  # for demonstration only (> 5e4 recommended)
#'                 )
#'
#' plot_edgeevidence(fit)
#'
#' oldpar <- par(mfrow = c(1,1))
#' 
#' par(mfrow = c(1, 2))
#' plot_edgeevidence(fit, split = TRUE)
#'
#' #' par(mfrow = c(1, 3))
#' plot_edgeevidence(fit, show = "included")
#' plot_edgeevidence(fit, show = "inconclusive")
#' plot_edgeevidence(fit, show = "excluded")
#' 
#' par(oldpar)
#' }


plot_edgeevidence <- function(output, evidence_thresh = 10, split = FALSE, show = "all",...) {
  if(any(any(class(output) == "easybgm"), any(class(output) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }

  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }

  UseMethod("plot_edgeevidence", output)

}

# ---------------------------------------------------------------------------------------------------------------

#' @title Network plot
#'
#' @description The network plot visualizes the strength of interactions between two nodes, the partial associations. Solely edges with a posterior inclusion probability larger than the `exc_prob` argument (default = 0.5) are shown. Edge thickness and saturation represent the strength of the association; the thicker the edge, the stronger the association. Red edges indicate negative relations and blue edges indicate positive associations.
#'
#' @name network
#'
#' @param output Output object from the easybgm function. Supports also objects from the bgm function of the `bgms` package.
#' @param exc_prob The threshold for excluding edges. All edges with a lower inclusion probability will not be shown. The default is set to 0.5 in line with the median probability plot.

#' @param dashed A binary parameter indicating whether edges with inconclusive evidence should be dashed. Default is FALSE
#' @param evidence_thresh If dashed = TRUE, users can specify the threshold for sufficient evidence for inclusion. All edges with evidence lower than `evidence_tresh` are dashed.
#' @param ... Additional arguments passed onto `qgraph`.
#' 
#' @return Returns a plot
#'
#' @export
#' @examples
#'
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#' fit <- easybgm(data, type = "continuous",
#'                 iter = 1000  # for demonstration only (> 5e4 recommended)
#'                 )
#'
#' plot_network(fit)
#'
#' # Shows all edges with an inclusion probability larger than 0.1
#' plot_network(fit, exc_prob = 0.1)
#'
#' # Indicate which edges have insufficient evidence for inclusion through a dashed line
#' plot_network(fit, dashed = TRUE, evidence_thresh = 10)
#' 

plot_network <- function(output, exc_prob = .5, evidence_thresh = 10, dashed = FALSE, ...) {

  if(any(any(class(output) == "easybgm"), any(class(output) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }

  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }
  UseMethod("plot_network", output)

}


# -------------------------------------------------

#' @title Structure plot
#'
#' @description The plot shows the resulting graph structure, i.e. all edges with some evidence of inclusion (i.e., inclusion Bayes factor greater than 1).
#'
#' @name structure
#'
#' @param output Output object from the easybgm function. Supports also objects from the bgm function of the `bgms` package.
#' @param ... Additional arguments passed onto `qgraph`
#'
#' @return Returns a plot
#' 
#' @export
#'
#' @import qgraph
#'
#' @examples
#' \donttest{
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#' fit <- easybgm(data, type = "ordinal",
#'                 iter = 1000  # for demonstration only (> 5e4 recommended)
#'                )
#'
#' plot_structure(fit)
#' }

plot_structure <- function(output, ...) {
  if(any(any(class(output) == "easybgm"), any(class(output) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }

  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }
  UseMethod("plot_structure", output)

}



# ---------------------------------------------------------------------------------------------------------------

#' @title Plot of interaction parameters and their 95% highest density intervals
#'
#' @description Plots the 95% highest density interval of the posterior distribution of the parameter estimates. The plot can be used to visualize the uncertainty of the partial association estimates. The x-axis indicates the strength of the partial association. The y-axis indicates the edge between nodes $i$ and $j$. The farther the posterior estimates (i.e., the points in the plot) are from zero, the stronger the partial association of the edge. The wider the highest density intervals (i.e., the error bar around the point), the less certain we are about the strength of the association.
#'
#' @name HDI
#'
#' @param output Output object from the easybgm function. Supports also objects from the bgm function of the `bgms` package.
#' @param ... Additional arguments passed onto `ggplot2`
#'
#' @return Returns a plot
#' 
#' @export
#' @import ggplot2 HDInterval
#' @importFrom stats median
#' @examples
#' \donttest{
#'
#' library(easybgm)
#' library(bgms)
#'
#'
#' data <- na.omit(Wenchuan)
#' fit <- easybgm(data, type = "ordinal",
#'               iter = 1000,  # for demonstration only (> 5e4 recommended)
#'               save = TRUE, edge_selection = TRUE)
#' plot_parameterHDI(fit)
#' }



plot_parameterHDI <- function(output, ...) {
  if(any(any(class(output) == "easybgm"), any(class(output) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }

  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }
  UseMethod("plot_parameterHDI", output)

}

# ---------------------------------------------------------------------------------------------------------------
# Centrality plot

#' @title Plot strength centralities and 95% highest density interval
#'
#' @description Visualize the strength centralities and their uncertainties. The centrality estimate can be obtained for each sample of the posterior distribution of the association parameters to obtain an estimate of the uncertainty of the strength centrality estimate.
#'
#' @name centrality
#'
#' @param output One output object or a list of several output objects from the easybgm function. Supports also objects from the bgm function of the `bgms` package.
#' @param group_names Specifying the group names, when providing a list of output objects. Needs to be a vector in the same length as the provided number of objects.
#' @param ... Additional arguments passed onto `ggplot2`
#'
#' @return Returns a plot
#' 
#' @importFrom dplyr arrange
#' @importFrom ggplot2 .data
#' @export
#'
#' @examples
#' \donttest{
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)
#' fit <- easybgm(data, type = "ordinal",
#'                 iter = 1000,  # for demonstration only (> 5e4 recommended)
#'                 save = TRUE, edge_selection = TRUE,
#'                 centrality = TRUE)
#'
#' plot_centrality(fit)
#' }

plot_centrality <- function(output, group_names = NULL, ...){
  if(any(any(class(output) == "easybgm"), any(class(output) == "list"), any(class(output) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }

  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }
  UseMethod("plot_centrality", output)

}


#' Prior sensitivity plot
#' @title Plot sensitivity to edge inclusion prior setting
#' @description For a given list of easybgm outputs with different prior edge inclusion probabilities, the function
#'  plots the percentage of edges that are included, excluded, and inconclusive.
#' @name prior_sensitivity 
#' @param output A list of easybgm outputs with different prior edge inclusion probabilities
#' @param ... Additional arguments passed onto ggplot2.
#'
#' @return Returns a plot
#' 
#' @export 
#'
#' @examples
#' \donttest{
#'
#' library(easybgm)
#' library(bgms)
#'
#' #data <- na.omit(Wenchuan)
#' #fit1 <- easybgm(data, type = "ordinal",
#' #               iter = 1000  # for demonstration only (> 5e4 recommended),
#' #                inclusion_probability = .1
#' #               )
#' #fit2 <- easybgm(data, type = "ordinal",
#' #                  iter = 1000,
#' #                  inclusion_probability = .5
#' #             )
#' #fit3 <- easybgm(data, type = "ordinal",
#' #                iter = 1000, inclusion_probability = .9)              
#' 
#' #plot_prior_sensitivity(list(fit1, fit2, fit3))
#' }

plot_prior_sensitivity <- function(output, ...) {
  if (!is.list(output))
    stop("Wrong input provided. Please provide a list of outputs of the easybgm function.")
  if(any(any(class(output[[1]]) == "easybgm"), any(class(output[[1]]) == "bgms")) == FALSE){
    stop("Wrong input provided. The function requires as input the output of the easybgm or bgm function.")
  }
  
  if(any(class(output)=="bgms") & (packageVersion("bgms") < "0.1.1")){
    stop("The fit of this version of bgms is not compatible with the plot. Please install the latest package version and refit the data.")
  }
  
  UseMethod("plot_prior_sensitivity", output)
}



