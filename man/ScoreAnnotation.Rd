% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_Scoring_Annotation_.R
\name{ScoreAnnotation}
\alias{ScoreAnnotation}
\title{Apply unsupervised evaluation metrics to score manual annotation of data}
\usage{
ScoreAnnotation(exprs, annotation, unassigned_labels = c())
}
\arguments{
\item{exprs}{numeric matrix: expression data, with measured markers in columns and events in rows}

\item{annotation}{factor vector: manual labels per row of \code{exprs}. Can be extracted via \code{GetAnnotation(benchmark, concatenate = TRUE)}}

\item{unassigned_labels}{character vector: names of label(s) given to cells that are left as unassigned by the annotation strategy. Default value is \code{c()} (empty vector)}
}
\description{
Apply unsupervised evaluation metrics to score manual annotation of data
}
