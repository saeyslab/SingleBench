
#' Evaluate a benchmark
#'
#' This function evaluates a benchmark pipeline, specified by an object of type \code{Benchmark}.
#' This means that all the projection, clustering or projection->clustering subpipelines that were set up when creating the benchmark object are executed, and their performance is scored.
#' Both the benchmark object and its auxiliary HDF5 file (created when the \code{Benchmark} constructor was called) are needed for this.
#' 
#' # Optional scoring of projection steps
#' 
#' Optionally, results of projection steps (if included) can be scored using evaluation metrics designed to measure the quality of dimension reduction (preservation of information versus original high-dimensional data).
#' This makes sense for methods that reduce dimensionality of the original data for the purposes of visualisation or to make the data more amenable to clustering.
#' To turn on scoring of projection steps, set parameter \code{score_projection} to \code{TRUE}.
#' Beware, this is only tractable for small datasets, due to the necessity to compute a co-ranking matrix (quadratic complexity with size of dataset).
#' The local continuity meta-criterion (LCMC), trustworthiness and continuity are computed.
#' 
#' # Parallelisation
#' 
#' For stability analysis of clustering tools, repeated runs of the tool can be run in parallel (unless this is forbidden in the tool wrapper).
#' The \code{parallel}, \code{doParallel} and \code{foreach} packages are needed for this.
#' To do this, specify the parameter \code{n_cores}.
#' Use \code{parallel::detectCores()} to check how many CPU cores are available.
#' 
#' @param benchmark object of class \code{Benchmark}, as generated by the constructor \code{Benchmark}
#' @param score_projections logical: whether results of projection steps should be scored. Default value is \code{FALSE}
#' @param projection_neighbourhood integer: size of neighbourhood considered as 'local' for scoring of projections. Default value is \code{100}
#' @param n_cores optional integer: number of CPU threads to use for parallelisation of repeated runs of clustering for stability analysis. Default value is \code{NULL} (no parallelisation)
#' @param which_python optional string: path to Python if Python needs to be used via \code{reticulate}. Default value is \code{NULL} (\code{reticulate} uses its default Python configuration)
#' @param seed.projection optional numeric value: value random seed to be used prior to each deployment of a projection method. (Use \code{NULL} to avoid setting a seed.) Default value is \code{1}
#' @param seed.clustering optional numeric value: value random seed to be used prior to each deployment of a clustering method. (Use \code{NULL} to avoid setting a seed.) Default value is \code{1}
#' @param ask_overwrite logical: if \code{benchmark} was evaluated before, should the user be asked prior to overwriting the previous evaluation results? Default value is \code{TRUE}
#' @param verbose logical: should progress messages be printed during evaluation? Default value is \code{TRUE}
#'
#' @seealso 
#' 
#' * **\code{AddLayout}**: allows you to add a separate 2-dimensional layout of the input dataset or to use an existing projection (produced in the evaluation) as a visualisation layout.
#'
#' @export
Evaluate <- function(benchmark, score_projections, projection_neighbourhood, n_cores, which_python, seed.projection, seed.clustering, ask_overwrite, verbose) UseMethod('Evaluate', benchmark)

#' Evaluate a benchmark
#'
#' @export
Evaluate.Benchmark <- function(
  benchmark,
  score_projections        = FALSE,
  projection_neighbourhood = 100,
  n_cores                  = NULL,
  which_python             = NULL,
  seed.projection          = 1,
  seed.clustering          = 1,
  ask_overwrite            = TRUE,
  verbose                  = TRUE
) {
  
  if (!file.exists(benchmark$h5_path))
    stop(paste0('Auxiliary HDF5 file ', benchmark$h5_path, ' not found'))
  
  if (benchmark$evaluated_previously && ask_overwrite) {
    cat(crayon::bgRed('     ')); .msg(' (?) '); cat(crayon::bgRed('     \n'))
    response <- readline(prompt = paste0('Benchmark was evaluated previously. Overwrite evaluation results? (Y/n) '))
    if (response != 'Y') {
      .msg_alt_bad('Aborting Benchmark evaluation, returning NULL\n')
      return(NULL)
    }
  }
  
  if (benchmark$uses_python) {
    if (!is.null(which_python)) {
      if (verbose) { .msg_python('Configuring reticulate') }
      reticulate::use_python(which_python, required = TRUE)
    }
  }
  
  if (score_projections) {
    benchmark$compute_dist <- TRUE
  }
  
  if (verbose) {
    .msg('Starting evaluation of '); .msg_name(benchmark$name);
    .msg(', time stamp: '); .msg_alt(as.character(Sys.time()), '\n')
  }
  
  if (benchmark$compute_knn && !benchmark$knn_available) {
    knn <- Evaluate_ComputekNNMatrix(benchmark, verbose)
    SavekNNMatrix(benchmark, knn, verbose)
    benchmark$knn_available <- TRUE
  }
  
  if (benchmark$compute_dist && !benchmark$dist_available) {
    if (verbose) .msg('Computing distance matrix... ')
    systime <- system.time({
      d <- coRanking:::euclidean_C(GetExpressionMatrix(benchmark, concatenate = TRUE))
    })
    if (verbose) {
      .msg_alt_good('done in ', round(systime['elapsed'], 2), ' seconds\n')
    }
    SaveDistanceMatrix(benchmark, d, verbose)
    if (verbose)
      .msg_alt_good('done\n')
  }
  
  benchmark$n_cores                  <- n_cores
  benchmark$seed.projection          <- seed.projection
  benchmark$seed.clustering          <- seed.clustering
  benchmark$score_projections        <- score_projections
  
  HDF5_InitialiseEvaluationResults(benchmark)
  
  EvalProjection(benchmark, verbose)
  EvalClustering(benchmark, verbose)
  
  benchmark$evaluated_previously <- TRUE
  if (verbose) { .msg_alt_good('Evaluation complete'); .msg(', time stamp: '); .msg_alt(as.character(Sys.time())) }
  
  gc(verbose = FALSE)
  invisible(benchmark)
}
