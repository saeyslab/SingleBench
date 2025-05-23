% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_Visual_CompositionMaps.R
\name{PlotCompositionMap}
\alias{PlotCompositionMap}
\title{Create cluster composition map}
\usage{
PlotCompositionMap(benchmark, idx.subpipeline, idcs.n_param = NULL)
}
\arguments{
\item{benchmark}{an object of class \code{Benchmark}, as generated by the constructor \code{Benchmark} and evaluated using \code{Evaluate}}

\item{idx.subpipeline}{integer value: index of sub-pipeline that includes a clustering step}

\item{idcs.n_param}{integer or vector of integers: indices of \emph{n}-parameter iteration(s) of interest. Default value is \code{NULL} (all of them, or \emph{n}-parameter was not used)}
}
\value{
list of plots, one per each n-parameter iteration
}
\description{
Creates a plot showing the composition of each cluster in terms of manually annotated populations.
You need to specify a sub-pipeline and \emph{n}-parameter iteration(s) by index.
If multiple \emph{n}-parameter iterations are chosen, a list of plots is returned.
}
