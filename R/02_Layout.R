
#' Add 2-dimensional layout to a Benchmark object
#'
#' Applies a projection tool to produce a 2-dimensional layout of all input data, for visualisation purposes downstream. 
#'
#' During the evaluation, a 2-dimensional latent-space projection of the input dataset might have been created already.
#' In that case, you can specify, by sub-pipeline and *n*-parameter iteration index, the projection to use as layout (instead of generating a new one).
#' Moreover, you can specify a more-than-2-dimensional projection, in which case the first two column of its coordinate matrix are used (this might be sensible when using PCA or diffusion maps).
#' 
#' Finally, if you want to use an existing coordinate matrix instead, you can pass it as a parameter to this function and it will be used as layout.
#'
#' @param benchmark object of type \code{Benchmark} object created using the \code{Benchmark} constructor and successfully evaluated using \code{Evaluate}
#' @param idx.subpipeline optional integer: sub-pipeline index to identify existing latent-space projection to use. Default value is \code{NULL}
#' @param idx.n_param optional integer: *n*-parameter iteration index to identify existing latent-space projection to use. Default value is \code{NULL}
#' @param method optional object of class \code{WrapperWithParameters}, generated by the function \code{Fix}, that is set to produce a 2-dimensional layout of input data associated with the benchmark object
#' @param idcs_training optional vector of integers: if input dataset consists of multiple expression matrices (corresponding to multiple samples), which of them should be used for training the dimension-reduction model? Default value is \code{NULL}, which causes all of them to be used
#' @param coordinates optional numeric matrix: row-wise coordinates of a 2-dimensional layout to use. If specified, this parameter overrides \code{idx.subpipeline}, \code{idx.n_param}, \code{wrapper} and \code{params}
#' @param ask_overwrite logical: if a 2-dimensional layout is already available, should the user be asked whether to overwrite it? Default value is \code{TRUE}
#' @param seed optional numeric value: value random seed to be used prior to the \code{wrapper} call (or \code{NULL} to avoid setting a seed). Default value is \code{1}
#' @param verbose logical: should progress messages be printed during evaluation? Default value is \code{TRUE}
#'
#' @seealso
#' 
#' * **\code{Plot}**: lets you produce a plot of your choosing, visualising results of a previously evaluated benchmark pipeline.
#'
#' @export
AddLayout <- function(benchmark, idx.subpipeline, idx.n_param, method, idcs_training, coordinates, ask_overwrite, seed, verbose) UseMethod('AddLayout', benchmark)

#' Add 2-dimensional layout to a Benchmark object
#'
#' @export
AddLayout.Benchmark <- function(
  benchmark,
  idx.subpipeline          = NULL,
  idx.n_param              = NULL,
  method                   = NULL,
  idcs_training            = NULL,
  coordinates              = NULL,
  ask_overwrite            = TRUE,
  seed                     = 1,
  verbose                  = TRUE
) {
  
  .AddLayout.ValidityChecks(environment())
  
  if (!is.null(coordinates)) {
      if (ncol(coordinates) != 2)
        stop('If specified, "coordinates" must have 2 columns')
      if (nrow(coordinates) != benchmark$row_count)
        stop('If specified, "coordinates" must have the same number of rows as the input expression data')
    if (verbose) .msg('Saving pre-computed 2-dimensional layout...')
    .h5write(obj = NA, file = benchmark$h5_path, name = 'Input/Layout/IsReferenceTo')
    .h5write(obj = coordinates, file = benchmark$h5_path, name = 'Input/Layout/Coordinates')
    benchmark$layout_available <- TRUE
    if (verbose) .msg_alt_good(' done\n')
    return(invisible(benchmark))
  }
  
  if (!is.null(idx.subpipeline) && !is.null(idx.n_param)) {
    if (verbose) {
      name <- GetNParameterIterationName_DR(benchmark, idx.subpipeline, idx.n_param)
      .msg('Using projection from '); .msg_alt(name); .msg(' for 2-dimensional layout\n')
    }
    .h5write(obj = c(idx.subpipeline, idx.n_param), file = benchmark$h5_path, name = 'Input/Layout/IsReferenceTo')
    benchmark$layout_available <- TRUE
    return(invisible(benchmark))
  }
  
  exprs <- GetExpressionMatrix(benchmark)
  
  if (verbose) {
    .msg('Generating layout using '); .msg_alt(method$name)
  }
  
  train   <- fTrain.WrapperWithParameters(method)
  extract <- fExtract.WrapperWithParameters(method)
  map     <- fMap.WrapperWithParameters(method)
  
  input <- GetExpressionMatrix(benchmark)
  knn   <- if (method$wrapper$uses_knn_graph) GetkNNMatrix(benchmark) else NULL
  
  set.seed(seed)
  res <- 
    train(
      input = if (is.null(idcs_training)) input else input[idcs_training],
      knn   = knn,
      exprs = input
    )
  
  res <-
    if (is.null(idcs_training))
      extract(res)
    else
      map(res, input)
  if (verbose) {
    .msg_alt_good(' done\n')
    .msg('Saving layout...')
  }
  
  .h5write(obj = NA, file = benchmark$h5_path, name = 'Input/Layout/IsReferenceTo')
  .h5write(obj = res[, 1:2], file = benchmark$h5_path, name = 'Input/Layout/Coordinates')
  if (verbose)
    .msg_alt_good(' done\n')
  
  benchmark$layout_available <- TRUE
  invisible(benchmark)
}
