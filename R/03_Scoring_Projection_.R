
ScoreProjection <- function(
  exprs, result, exprs_distance_matrix, knn, knn.k, knn.algorithm, knn.distance, verbose
) {
  
  if (is.list(exprs))
    exprs <- do.call(rbind, exprs)
  if (is.list(result$Projection))
    result$Projection <- do.call(rbind, result$Projection)
  
  n <- nrow(exprs)
  k <- projection_neighbourhood
  
  if (verbose)
    .msg('\t-> computing projection distance matrix... ')
  systime <- system.time({
    proj_distance_matrix <- coRanking:::euclidean_C(result$Projection)
  })
  if (verbose) {
    .msg_alt_good('done in ', round(systime['elapsed'], 2), ' seconds\n')
  }
  
  exprs_rankmatrix <- coRanking::rankmatrix(exprs_distance_matrix, input = 'dist')
  proj_rankmatrix <- coRanking::rankmatrix(proj_distance_matrix, input = 'dist')
  
  Q <- coRanking:::coranking_C(exprs_rankmatrix, proj_rankmatrix)
  lcmc <- coRanking::LCMC(Q, K = as.integer(k))
  
  Gk <- ifelse(k < n / 2, n*k*(2*n-3*k-1), n*(n-k)*(n-k-1))
  
  LL <- Q[(k+1):nrow(Q), 1:k] # ext (lower-left quadrant)
  UR <- Q[1:k, (k+1):ncol(Q)] # int (upper-right quadrant)
  
  eT <- 1-2/Gk*sum(LL * 1:nrow(LL))
  eC <- 1-2/Gk*sum(t(t(UR) * 1:ncol(UR)))
  
  Un <- 1 / (k * n) * sum(Q[1:k, 1:k][upper.tri(Q[1:k, 1:k])])
  Ux <- 1 / (k * n) * sum(Q[1:k, 1:k][lower.tri(Q[1:k, 1:k])])
  Bnx <- Un - Ux # Bnx > 0 <=> intrusive, Bnx < 0 <=> extrusive
  
  knn_res <- list(Indices = NA, Distances = NA)
  Q <- as.matrix(Q)
  class(Q) <- 'matrix'
  
  list(
    'Layout k-NNG' = knn_res,
    'Collapsed' = FALSE,
    'Co-Ranking Matrix' = Q,
    'Local Continuity Meta-Criterion' = lcmc,
    'Relative Intrusiveness' = Bnx,
    'Trustworthiness' = eT,
    'Continuity' = eC,
    'Projection Neighbourhood' = knn.k
  )
}
