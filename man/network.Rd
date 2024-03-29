% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgm_plot_class.R
\name{network}
\alias{network}
\alias{plot_network}
\title{Network plot}
\usage{
plot_network(output, exc_prob = 0.5, evidence_thresh = 10, dashed = FALSE, ...)
}
\arguments{
\item{output}{Output object from the easybgm function. Supports also objects from the bgm function of the \code{bgms} package.}

\item{exc_prob}{The threshold for excluding edges. All edges with a lower inclusion probability will not be shown. The default is set to 0.5 in line with the median probability plot.}

\item{evidence_thresh}{If dashed = TRUE, users can specify the threshold for sufficient evidence for inclusion. All edges with evidence lower than \code{evidence_tresh} are dashed.}

\item{dashed}{A binary parameter indicating whether edges with inconclusive evidence should be dashed. Default is FALSE}

\item{...}{Additional arguments passed onto \code{qgraph}.}
}
\value{
Returns a plot
}
\description{
The network plot visualizes the strength of interactions between two nodes, the partial associations. Solely edges with a posterior inclusion probability larger than the \code{exc_prob} argument (default = 0.5) are shown. Edge thickness and saturation represent the strength of the association; the thicker the edge, the stronger the association. Red edges indicate negative relations and blue edges indicate positive associations.
}
\examples{


library(easybgm)
library(bgms)

data <- na.omit(Wenchuan)
fit <- easybgm(data, type = "continuous",
                iter = 1000  # for demonstration only (> 5e4 recommended)
                )

plot_network(fit)

# Shows all edges with an inclusion probability larger than 0.1
plot_network(fit, exc_prob = 0.1)

# Indicate which edges have insufficient evidence for inclusion through a dashed line
plot_network(fit, dashed = TRUE, evidence_thresh = 10)

}
