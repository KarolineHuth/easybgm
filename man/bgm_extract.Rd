% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bgm_extract.R
\name{bgm_extract}
\alias{bgm_extract}
\title{Extract results from a Bayesian analysis of networks}
\usage{
bgm_extract(
  fit,
  method,
  edge.prior = 0.5,
  package = "BDgraph",
  posterior_samples = F,
  not.cont = NULL,
  data = NULL,
  centrality = F
)
}
\arguments{
\item{fit}{fit object of the respective package used. Note for objects from the package 'BGGM', the package requires the input from explore(data)}

\item{method}{type of model estimated, e.g., ggm, gcgm, dgm-binary, Ising}

\item{package}{package used to obtain the fit object}

\item{posterior_samples}{if TRUE, the posterior samples will be extracted. Note will significantly increase the computation time for 'BDgraph'.}

\item{not.cont}{only if method = "gcgm" vector indicating the not-continuous variables}

\item{data}{if posterior_samples = T, provide the raw data used to estimate the network}

\item{centrality}{if TRUE, the centrality samples will be extracted. Note will significantly increase the computation time.}
}
\description{
Extract results from a Bayesian analysis of networks
}
