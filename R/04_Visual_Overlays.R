
PlotOverlay <- function(
    layout,
    labels,
    labels_name,
    show_labels,
    labelsize,
    raster_threshold,
    svg_pointsize,
    raster_pointsize,
    plot_title,
    plot_subtitle,
    fix_aspect_ratio
) {
  
  
  d <- cbind(layout, 'PointLabel' = labels)
  
  points <-
    if (nrow(layout) > raster_threshold) {
      scattermore::geom_scattermore(pointsize = raster_pointsize)
    } else {
      geom_point(size = svg_pointsize, alpha = 0.85)
    }
  
  res <- ggplot(
    data = d,
    mapping = aes(x = Component1, y = Component2, col = PointLabel)
  ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      col = labels_name
    ) +
    theme_dark() + theme(plot.margin = unit(c(1, 1, 1, 1), 'cm')) +
    points
  if (fix_aspect_ratio) res <- res + coord_fixed()
  
  if (show_labels) {
    res <- res + geom_label(
      data = d %>%
        dplyr::group_by(PointLabel) %>%
        dplyr::mutate(Center1 = median(Component1), Center2 = median(Component2)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(PointLabel, .keep_all = TRUE),
      aes(x = Center1, y = Center2, label = PointLabel),
      colour = 'black',
      size = labelsize,
      alpha = 0.6,
      fill = 'white'
    )
  }
  res
}

#' Create a labels-overlay plot
#'
#' Creates a plot showing separation of manually annotated populations using previously generated 2-dimensional layout of input expression data.
#' For large numbers of data points, this function uses raster graphics to visualise the layout faster.
#'
#' @param benchmark object of type \code{Benchmark}
#' @param exclude_unassigned logical: if \code{TRUE}, data points that are considered unassigned per manual annotation are omitted. Default value is \code{FALSE}
#' @param raster_threshold integer: maximum number of data points for which vector graphics should be used. Default value is \code{5000}
#' @param plot_title string: title of the plot. Default value is '*Clustering*'
#' @param population_labels logical: whether to label populations within the plot. Default value is \code{FALSE}
#' @param labelsize numeric: size of the population labels within the plot, if shown. Default value is \code{3}
#' @param palette string vector: palette of distinctive colours for each cluster. Defaults to a palette of 32 colours
#' @param svg_pointsize numeric: size of the data points in the plot when using vector graphics. Default value is \code{0.5}
#' @param raster_pointsize numeric: size of the data points in the plot when using raster graphics. Default value is \code{0.005}
#' @param fix_aspect_ratio logical: whether to maintain a 1:1 aspect ratio of the plot. Default value is \code{FALSE}#'
#' @seealso
#' 
#' * **\code{Benchmark}**: constructs a \code{Benchmark} object (benchmark pipeline set-up)\
#' 
#' * **\code{Evaluate}**: runs all benchmark sub-pipelines and scores the performance of each tool
#' 
#' * **\code{AddLayout}**: allows you to add a separate 2-dimensional layout of the input dataset or to use an existing projection (produced in the evaluation) as a visualisation layout.
#'
#' @export
PlotLabelsOverlay <- function(
  benchmark,
  exclude_unassigned = FALSE,
  raster_threshold = 5000L,
  plot_title = 'Cell annotation',
  population_labels = FALSE,
  labelsize = 3,
  palette = c(
    RColorBrewer::brewer.pal(8, 'Dark2'),
    RColorBrewer::brewer.pal(12,'Paired'),
    RColorBrewer::brewer.pal(12,'Set3')
  ),
  svg_pointsize = 0.5,
  raster_pointsize = 0.005,
  fix_aspect_ratio = FALSE
) {
  if (class(benchmark) != 'Benchmark') stop('"benchmark" is not of class "Benchmark"')
  if (!benchmark$layout_available) stop('No 2-dimensional layout available')
  
  layout <- as.data.frame(GetLayout(benchmark), concatenate = TRUE)
  annotation <- GetAnnotation(benchmark, concatenate = TRUE)
  
  mask <- rep(TRUE, nrow(layout))
  if (exclude_unassigned) {
    levels(annotation)[levels(annotation) %in% benchmark$unassigned_labels] <- NA
    mask <- !is.na(annotation)
    layout <- layout[mask, ]
  }
  
  PlotOverlay(
    layout, labels = annotation, labels_name = 'Population', show_labels = population_labels,
    labelsize = labelsize, raster_threshold = raster_threshold, svg_pointsize = svg_pointsize,
    raster_pointsize = raster_pointsize, plot_title = plot_title, plot_subtitle = NULL,
    fix_aspect_ratio = fix_aspect_ratio
  )
}

#' Create a clusters-overlay plot
#'
#' Creates a plot showing separation of generated clusters using previously generated 2-dimensional layout of input expression data.
#' For large numbers of data points, this function uses raster graphics to visualise the layout faster.
#' 
#' You need to specify a single sub-pipeline and one or more *n*-parameter iteration by index.
#' If multiple *n*-parameter iterations are chosen, a list of plots will be generated.
#' 
#' You can choose a custom title for your plot.
#' The plot subtitle is always the *n*-parameter iteration name.
#'
#' @param benchmark object of type \code{Benchmark}
#' @param idx.subpipeline integer: index of sub-pipeline that includes a clustering step
#' @param idx.n_param integer: index of *n*-parameter of interest. Default value is \code{NULL}
#' @param idx.run integer: index of repeated run of interest. Default value is 1
#' @param exclude_unassigned logical: if \code{TRUE}, data points that are considered unassigned per manual annotation are omitted. Default value is \code{FALSE}
#' @param raster_threshold integer: maximum number of data points for which vector graphics should be used. Default value is \code{5000}
#' @param plot_title string: title of the plot. Default value is '*Clustering*'
#' @param cluster_labels logical: whether to label clusters within the plot. Default value is \code{FALSE}
#' @param labelsize numeric: size of the cluster labels within the plot, if shown. Default value is \code{3}
#' @param palette string vector: palette of distinctive colours for each cluster. Defaults to a palette of 32 colours
#' @param svg_pointsize numeric: size of the data points in the plot when using vector graphics. Default value is \code{0.5}
#' @param raster_pointsize numeric: size of the data points in the plot when using raster graphics. Default value is \code{0.005}
#' @param fix_aspect_ratio logical: whether to maintain a 1:1 aspect ratio of the plot. Default value is \code{FALSE}
#'
#' @export
PlotClustersOverlay <- function(
  benchmark,
  idx.subpipeline,
  idx.n_param = NULL,
  idx.run = 1,
  exclude_unassigned = FALSE,
  raster_threshold = 5000L,
  plot_title = 'Clustering',
  cluster_labels = FALSE,
  labelsize = 3,
  palette = c(
    RColorBrewer::brewer.pal(8, 'Dark2'),
    RColorBrewer::brewer.pal(12,'Paired'),
    RColorBrewer::brewer.pal(12,'Set3')
  ),
  svg_pointsize = 0.5,
  raster_pointsize = 0.005,
  fix_aspect_ratio = FALSE
) {
  if (class(benchmark) != 'Benchmark') stop('"benchmark" is not of class "Benchmark"')
  if (!benchmark$layout_available) stop('No 2-dimensional layout available')
  if (!benchmark$evaluated_previously) stop('Benchmark pipeline has not been evaluated')
  
  layout <- as.data.frame(GetLayout(benchmark))
  
  mask <- rep(TRUE, nrow(layout))
  if (exclude_unassigned) {
    annotation <- GetAnnotation(benchmark, concatenate = TRUE)
    levels(annotation)[levels(annotation) %in% benchmark$unassigned_labels] <- NA
    mask <- !is.na(annotation)
    layout <- layout[mask, ]
  }
  
  subpipeline_name <- GetNParameterIterationName(benchmark, idx.subpipeline, idx.n_param)
  clustering <- GetClustering(benchmark, idx.subpipeline, idx.n_param, idx.run)[mask]
  
  PlotOverlay(
    layout, labels = clustering, labels_name = 'Cluster', show_labels = cluster_labels,
    labelsize = labelsize, raster_threshold = raster_threshold, svg_pointsize = svg_pointsize,
    raster_pointsize = raster_pointsize, plot_title = plot_title, plot_subtitle = subpipeline_name,
    fix_aspect_ratio = fix_aspect_ratio
  )
}