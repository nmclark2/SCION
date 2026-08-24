#' Cluster genes by independent component analysis (non-temporal data)
#'
#' Original author Mitch Elmore; adapted for SCION by Natalie Clark.
#' See <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181195>.
#'
#' @param clustering_data data frame or matrix, rows = genes, columns = samples.
#' @param k Mojena cutoff multiplier for the hierarchical clustering of independent
#'   components. Smaller `k` produces more, tighter clusters.
#' @return a data frame: centered/scaled `clustering_data` plus a `clusters` column.
#' @keywords internal
ica_clustering <- function(clustering_data, k) {
  # center and scale expression matrix
  normmatrix <- t(scale(t(clustering_data), scale = TRUE, center = TRUE))

  x_ica <- fastICA::fastICA(normmatrix, n.comp = ncol(normmatrix), alg.typ = "parallel",
                             fun = "logcosh", alpha = 1.0, method = "C", row.norm = FALSE,
                             maxit = 5000, tol = 1e-03, verbose = FALSE)
  hc_ica <- stats::hclust(stats::dist(x_ica$S), method = "ward.D", members = NULL)
  mojena <- mean(hc_ica$height) + k * stats::sd(hc_ica$height)
  cluster_num <- length(hc_ica$height[hc_ica$height > mojena]) + 1
  clusters <- stats::cutree(hc_ica, k = cluster_num)
  message(sprintf("Found %d clusters", cluster_num))

  # convert matrix to dataframe
  results <- as.data.frame(cbind(normmatrix, clusters))

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
  #     ggplot2::scale_x_discrete(limits = unique(stacked$group), labels = colnames(normmatrix))
  #   g <- g + ggplot2::geom_line(data = reflinestacked[reflinestacked$cluster == i, ],
  #                                ggplot2::aes(x = group, y = Norm.Intensity, group = cluster),
  #                                colour = "black")
  #   plotly::ggplotly(g)
  #   ggplot2::ggsave(paste0("Cluster Plots/cluster", i, ".png"), device = "png",
  #                    width = 5, height = 3, dpi = 600)
  # }

  results
}
