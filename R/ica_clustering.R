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
  normmatrix <- t(scale(t(clustering_data), scale = TRUE, center = TRUE))

  x_ica <- fastICA::fastICA(normmatrix, n.comp = ncol(normmatrix), alg.typ = "parallel",
                             fun = "logcosh", alpha = 1.0, method = "C", row.norm = FALSE,
                             maxit = 5000, tol = 1e-03, verbose = TRUE)
  hc_ica <- stats::hclust(stats::dist(x_ica$S), method = "ward.D", members = NULL)
  mojena <- mean(hc_ica$height) + k * stats::sd(hc_ica$height)
  cluster_num <- length(hc_ica$height[hc_ica$height > mojena]) + 1
  clusters <- stats::cutree(hc_ica, k = cluster_num)

  results <- as.data.frame(cbind(normmatrix, clusters))
  results
}
