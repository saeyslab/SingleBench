% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Layout.R
\name{AddLayout}
\alias{AddLayout}
\title{Add 2-dimensional layout to a Benchmark object}
\usage{
AddLayout(
  benchmark,
  idx.subpipeline,
  idx.n_param,
  method,
  idcs_training,
  coordinates,
  ask_overwrite,
  seed,
  verbose
)
}
\arguments{
\item{benchmark}{object of type \code{Benchmark} object created using the \code{Benchmark} constructor and successfully evaluated using \code{Evaluate}}

\item{idx.subpipeline}{optional integer: sub-pipeline index to identify existing latent-space projection to use. Default value is \code{NULL}}

\item{idx.n_param}{optional integer: \emph{n}-parameter iteration index to identify existing latent-space projection to use. Default value is \code{NULL}}

\item{method}{optional object of class \code{WrapperWithParameters}, generated by the function \code{Fix}, that is set to produce a 2-dimensional layout of input data associated with the benchmark object}

\item{idcs_training}{optional vector of integers: if input dataset consists of multiple expression matrices (corresponding to multiple samples), which of them should be used for training the dimension-reduction model? Default value is \code{NULL}, which causes all of them to be used}

\item{coordinates}{optional numeric matrix: row-wise coordinates of a 2-dimensional layout to use. If specified, this parameter overrides \code{idx.subpipeline}, \code{idx.n_param}, \code{wrapper} and \code{params}}

\item{ask_overwrite}{logical: if a 2-dimensional layout is already available, should the user be asked whether to overwrite it? Default value is \code{TRUE}}

\item{seed}{optional numeric value: value random seed to be used prior to the \code{wrapper} call (or \code{NULL} to avoid setting a seed). Default value is \code{1}}

\item{verbose}{logical: should progress messages be printed during evaluation? Default value is \code{TRUE}}
}
\description{
Applies a projection tool to produce a 2-dimensional layout of all input data, for visualisation purposes downstream.
}
\details{
During the evaluation, a 2-dimensional latent-space projection of the input dataset might have been created already.
In that case, you can specify, by sub-pipeline and \emph{n}-parameter iteration index, the projection to use as layout (instead of generating a new one).
Moreover, you can specify a more-than-2-dimensional projection, in which case the first two column of its coordinate matrix are used (this might be sensible when using PCA or diffusion maps).

Finally, if you want to use an existing coordinate matrix instead, you can pass it as a parameter to this function and it will be used as layout.
}
\seealso{
\itemize{
\item \strong{\code{Plot}}: lets you produce a plot of your choosing, visualising results of a previously evaluated benchmark pipeline.
}
}
