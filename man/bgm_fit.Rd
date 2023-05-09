% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgm_fit.R, R/bgm_fit_class.R
\name{bgm_fit}
\alias{bgm_fit}
\title{Fit a Bayesian analysis of networks}
\usage{
bgm_fit(x, ...)

bgm_fit(x, ...)
}
\arguments{
\item{x}{Object with a particular class that will dispatch to the respective package functions}

\item{...}{Additional arguments to be passed onto the respective fitting functions}

\item{data}{An n x p matrix or dataframe containing the variables for n independent observations on p variables.}

\item{type}{What is the data type? Options: continuos, mixed, ordinal, binary}

\item{package}{The R-package that should be used for fitting the network model}

\item{not.cont}{If data-type is mixed, a vector of length p, specifying the not-continuous
variables (1 = not continuous, 0 = continuous)}

\item{iter}{number of iterations for the sampler}

\item{save}{Logical. Should the posterior samples be obtained (default = FALSE)?}

\item{centrality}{Logical. Should the centrality measures be extracted (default = FALSE)?}

\item{progress}{Logical. Should a progress bar be shown (default = TRUE)?}
}
\description{
Fit a Bayesian analysis of networks

Fit a Bayesian analysis of networks
}