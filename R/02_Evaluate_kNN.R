
SavekNNMatrix <- function(
  benchmark,
  knn,
  verbose
) {
  if (verbose) .msg('Saving k-NN matrix...')
  .h5writekNNMatrix(benchmark, knn)
  if (verbose) .msg_alt_good(' done\n')
}

SaveDistanceMatrix <- function(
  benchmark,
  knn,
  verbose
) {
  if (verbose) .msg('Saving distance matrix...')
  .h5writeDistanceMatrix(benchmark, knn)
  if (verbose) .msg_alt_good(' done\n')
}


Evaluate_ComputekNNMatrix <- function(
  benchmark,
  verbose
) {
  exprs <- GetExpressionMatrix(benchmark, concatenate = TRUE)
  if (!is.null(benchmark$rel_idcs.knn_features))
    exprs <- exprs[, benchmark$rel_idcs.knn_features]
  if (verbose) .msg('Computing k-NN matrix...')
  
  systime <- NA
  res <- ComputekNNMatrix(exprs, benchmark$knn.k, out.systime = systime)
  
  if (verbose) .msg_alt_good(' done in ', round(systime['elapsed'], 2), ' seconds\n')
  res
}

#' Compute a \code{k}-nearest-neighbour matrix for expression data
#'
#' Finds \code{k} closest neighbours to each point in high-dimensional space and the distances to these points.
#' The approximate k-NN search algorith from \code{RcppHNSW} is used.
#' 
#' @param exprs numeric matrix: a coordinate matrix of biological expression data (columns correspond to markers, rows correspond to cells)
#' @param k integer: number of nearest neighbours to find for each point (not including self)
#' @param metric string: distance metric to pass to \code{RcppHNSW::hnsw_knn}. Default value is \code{'euclidean'}
#' @param out.systime optional out-variable: if an object is passed as \code{out.systime}, a side-effect of executing this function is that this object will be assigned elapsed time (in seconds) needed to complete the \code{k}-NN search
#'
#' @return list with two slots: \code{Indices} contains a matrix of nearest neighbours to each point (per row) and \code{Distances} contains a matrix of corresponding Euclidean distances
#'
#' @export
ComputekNNMatrix <- function(
  exprs,
  k,
  metric = 'euclidean',
  out.systime = NULL
) {
  systime <- system.time(
    knn <- RcppHNSW::hnsw_knn(exprs, k = k+1, distance = metric)
  )
  
  res <- list('Indices' = knn$idx[, 2:(k+1)], 'Distances' = knn$dist[, 2:(k+1)])
  if (!is.null(out.systime))
    eval.parent(substitute(out.systime <- systime['elapsed']))
  res
}
  