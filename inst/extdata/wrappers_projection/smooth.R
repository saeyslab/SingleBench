
wrapper.projection.smooth <- WrapTool(
  name = 'smooth',
  type = 'projection',
  r_packages      = c(),  # no non-base R packages needed
  use_knn_graph   = TRUE, # this triggers pre-computation of a k-NNG
  fun.build_model = function(
    input, knn, k = 50, n_iter = 1, lam = 0.5
  ) {
    if (n_iter == 0) return(input)
    for (iter in seq_len(n_iter)) {
      
      targets <- t(apply(
        X = knn$Indices[, seq_len(k)],
        MARGIN = 1,
        FUN = function(idcs) colMeans(input[idcs, ])
      ))
      shifts <- (targets-input)*lam
      input <- input+shifts
    }
    input
  }
)
