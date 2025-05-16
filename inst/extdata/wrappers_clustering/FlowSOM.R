
wrapper.clustering.FlowSOM <- WrapTool(
  name = 'FlowSOM',
  type = 'clustering',
  r_packages = c('FlowSOM'),
  fun.build_model =
    function(
      input, n_clusters, grid_width = 10, grid_height = 10, max_metaclusters = 90
    ) {
      ## Clustering using self-organising map (SOM)
      som <- FlowSOM::SOM(
        data = input, xdim = grid_width, ydim = grid_height, silent = TRUE
      )
      clus <- som$mapping[, 1]
      
      ## Merge clusters into metaclusters
      meta <-
        if (n_clusters == 0)
          FlowSOM::MetaClustering(
            som$codes,
            method = 'metaClustering_consensus',
            max = max_metaclusters
          )
        else
          FlowSOM::metaClustering_consensus(som$codes, k = n_clusters)
      list(
        som            = som,  # trained model
        clustering     = clus, # cluster assignments
        metaclustering = meta  # metacluster assignments
      )
    },
  fun.extract = function(model)
    model$metaclustering[model$clustering], # metacluster assignments per cell
  fun.apply_model = function(model, input)
    model$meta[ # metacluster assignment per cell for unseen data
      FlowSOM::MapDataToCodes(codes = model$som$codes, newdata = input)
    ]
)
