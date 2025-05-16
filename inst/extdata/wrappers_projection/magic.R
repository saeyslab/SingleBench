
wrapper.projection.magic <- WrapTool(
  name = 'magic',
  type = 'projection',
  r_packages = c(),
  python_modules = c('magic'),
  use_knn_graph = FALSE,
  random_seed_argument = 'random_state',
  fun.build_model = function(
    input,
    k = 5,
    k_max = NULL,
    decay = 1,
    t = 3,
    n_pca = 100,
    solver = 'exact',
    knn_dist = 'euclidean',
    n_jobs = -2,
    random_state = 0
  ) {
    
    magic <- reticulate::import('magic')
    
    k <- as.integer(k)
    decay <- as.integer(decay)
    if (is.null(k_max)) {
      k_max <- k*3
    }
    k_max <- as.integer(k_max)
    t <- as.integer(t)
    n_pca <- as.integer(n_pca)
    n_jobs <- as.integer(n_jobs)
    if (!is.null(random_state)) {
      random_state <- as.integer(random_state)
    }
    
    
    reticulate::py_capture_output(model <- magic$MAGIC(
      knn = k, knn_max = k_max, decay = decay, t = t, n_pca = n_pca, solver = solver,
      knn_dist = knn_dist, n_jobs = n_jobs, random_state = random_state, verbose = 0L
    ))
    reticulate::py_capture_output(input <- model$fit_transform(input))
    
    input
  }
)
