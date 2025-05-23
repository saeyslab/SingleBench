% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_Jacc.R
\name{Jacc}
\alias{Jacc}
\title{Match clusters across two clusterings and plot results}
\usage{
Jacc(
  c1,
  c2,
  obj = "f1",
  title = "Jaccard heatmap",
  unassigned = NULL,
  generate_plot = TRUE,
  c1_name = "c1",
  c2_name = "c2",
  scoring_matrix = NULL,
  verbose = FALSE
)
}
\arguments{
\item{c1}{factor, numeric or character vector: assignment of each data point to a cluster or otherwise defined population}

\item{c2}{factor, numeric or character vector: assignment of each data point to a cluster}

\item{obj}{string: evaluation metric used for matching groups in \code{c1} and \code{c2}; one of \code{f1} (default), \code{precision} and \code{recall}}

\item{title}{string: tile of Jaccard similarity heatmap plot (default value is '\code{Jaccard heatmap}')}

\item{unassigned}{optional string vector: names of levels of \code{c1} denoting unlabelled data points}

\item{generate_plot}{logical: whether a Jaccard heatmap-style plot should be generated (default value is \code{TRUE})}

\item{c1_name}{optional string: name of the \code{c1} vector to be used in text of the plot (default value is '\code{c1}')}

\item{c2_name}{optional string: name of the \code{c2} vector to be used in text of the plot (default value is '\code{c2}')}

\item{scoring_matrix}{optional numeric matrix: scoring matrix for hierarchical penalties (see function \code{Benchmark}). Default value is \code{NULL}}

\item{verbose}{logical: indicates whether to display progress messages (default value is \code{FALSE})}
}
\value{
list of results for evaluation approach (i) \code{Results.Bijective}, approach (ii) \code{Results.FixedC1} and approach \code{Results.FixedC2}, as well as a Jaccard similarity heatmap diagram \code{Plot} (if produced)
}
\description{
Compares two labels of cluster assignment per data point (or a vector of ground-truth labels and a clustering vector) \code{c1} and \code{c2}, matching groups in each vector to each other while maximising the value of an evaluation metric \code{obj}.
The evaluation metric \code{obj} is either \code{f1} (default), \code{precision} or \code{recall}.
}
\details{
Three approaches are used to solve the cluster-cluster (or label-cluster) matching problem. All of them seek to maximise the total value of \code{obj}.
Approach (i) gives 1-to-1 matches, whereby each group in \code{c1} is matched to a (different) group in \code{c2}. (In the special case where the number of groups in \code{c1} is equal to the number of groups in \code{c2}, this guarantees no unmatched groups.)
Approach (ii) uses a relaxed fixed-\code{c1} matching, whereby each group in \code{c1} is matched to the group in \code{c2} that maximises \code{obj} value of the match. This can result in 1-to-many matches.
Approach (iii) uses a relaxed fixed-\code{c2} matching, which mirrors approach (ii).

If \code{c1} is in fact a vector of ground-truth labels (or manual annotation of each data point), there may be de-facto unlabelled data points in the original data.
\code{unassigned} is an optional vector of the labels given to data points which don't belong to an annotated population.
If specified, the unassigned groups in \code{c1} are left out of the evaluation: points that are unassigned are ignored in constructing the contingency tables for each match and groups in \code{c2} may not be matched to these unassigned points.

In addition to evaluation results, a heatmap showing agreement between \code{c1} and \code{c2} and agreement between the different matching approaches is produced by default.
}
