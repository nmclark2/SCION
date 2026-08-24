#' Cluster genes by dynamic time warping (temporal data)
#'
#' For each gene, compares its normalized profile to every other
#' not-yet-clustered gene using dynamic time warping; genes whose warped
#' alignment matches sufficiently closely are placed in the same cluster.
#'
#' @param clustering_data data frame or matrix, rows = genes, columns = samples/timepoints.
#' @param threshold proportion of matching aligned indices (0-1) required for two genes to
#'   be clustered together. Values above 1 are clamped to 1.
#' @return a data frame: normalized `clustering_data` plus a `clusters` column, row names
#'   preserved from `clustering_data`.
#' @keywords internal
dtw_clustering <- function(clustering_data, threshold) {
  if (threshold > 1) {
    threshold <- 1
  }
  mydata <- clustering_data
  normmatrix <- (mydata[, ] - rowMeans(mydata) * matrix(1, nrow = dim(mydata)[1], ncol = dim(mydata)[2])) /
    apply(mydata[, ], 1, stats::sd)
  normmatrix[is.na(rowSums(normmatrix)), ] <- rep(1, dim(normmatrix)[2])

  clusters <- matrix(0, nrow = dim(normmatrix)[1], ncol = 1)
  ind <- 1
  # for each gene, compare its profile to the rest of the genes
  # any gene with a sufficiently matching profile is clustered together
  for (j in seq_len(dim(normmatrix)[1])) {
    # first check if gene j has already been clustered -- if it has, we skip
    # it; if it hasn't, we initialize the cluster #
    if (clusters[j] != 0) {
      next
    } else {
      clusters[j] <- ind
      ind <- ind + 1
    }
    for (k in (j + 1):dim(normmatrix)[1] - 1) {
      # first check if gene k has already been clustered -- if it has, we skip it
      if (clusters[k] != 0) {
        next
      }
      # otherwise perform dtw on the two genes -- DTW expects the time
      # series as a vector
      ts1 <- as.vector(t(normmatrix[j, ]))
      ts2 <- as.vector(t(normmatrix[k, ]))

      dtwresults <- dtw::dtw(ts1, ts2)

      # DTW gives us the indices that matched; if the time series matches
      # perfectly, then x1=y1, x2=y2, etc. Test via proportion of matches.
      xinds <- dtwresults$index1
      yinds <- dtwresults$index2
      propmatches <- sum(xinds == yinds) / length(xinds)

      # if propmatches >= threshold, cluster gene k with gene j; otherwise,
      # move on
      if (propmatches >= threshold) {
        clusters[k] <- clusters[j]
      }
    }
  }
  message(sprintf("Found %d clusters", length(unique(clusters))))
  # put together gene names, normalized expression, and clusters
  results <- data.frame(normmatrix, clusters, row.names = row.names(mydata))

  # prepare data for plotting
  # plotdata <- as.data.frame(t(results[, seq_len(ncol(normmatrix))]))
  # stacked <- utils::stack(plotdata)
  # stacked[, 3] <- rep(colnames(results)[seq_len(ncol(normmatrix))], ncol(plotdata))
  # stacked[, 4] <- rep(clusters, each = ncol(normmatrix))
  # colnames(stacked) <- c("Norm.Intensity", "gene", "group", "cluster")
  # # add means of each cluster as reference lines
  # reflines <- by(results[, seq_len(ncol(normmatrix))], results$clusters, colMeans)
  # reflinedata <- as.data.frame(do.call(cbind, reflines))
  # reflinestacked <- utils::stack(reflinedata)
  # reflinestacked[, 2] <- rep(colnames(results)[seq_len(ncol(normmatrix))], max(clusters))
  # reflinestacked[, 3] <- rep(seq_len(max(clusters)), each = ncol(normmatrix))
  # colnames(reflinestacked) <- c("Norm.Intensity", "group", "cluster")
  # # plot one cluster per file
  # for (i in seq_len(max(results$clusters))) {
  #   g <- ggplot2::ggplot(data = stacked[stacked$cluster == i, ],
  #                         mapping = ggplot2::aes(x = group, y = Norm.Intensity,
  #                                                 colour = as.factor(cluster), group = gene)) +
  #     ggplot2::geom_line() + ggplot2::theme(legend.position = "none") +
  #     ggplot2::scale_x_discrete(limits = unique(stacked$group), labels = colnames(mydata))
  #   g <- g + ggplot2::geom_line(data = reflinestacked[reflinestacked$cluster == i, ],
  #                                ggplot2::aes(x = group, y = Norm.Intensity, group = cluster),
  #                                colour = "black")
  #   plotly::ggplotly(g)
  #   ggplot2::ggsave(paste0("Cluster Plots/cluster", i, ".png"), device = "png",
  #                    width = 5, height = 3, dpi = 600)
  # }

  results
}
