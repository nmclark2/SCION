#' Cluster genes by k-means (non-temporal data)
#'
#' Tests several values of k around `kmid` and picks the one with the largest
#' average silhouette width.
#'
#' @param clustering_data data frame or matrix, rows = genes, columns = samples.
#' @param kmid starting/center value of k to test around. Candidate k values
#'   derived from this are clamped to `[2, nrow(clustering_data) - 1]`, since
#'   kmeans/silhouette need k < the number of genes being clustered.
#' @return a data frame: centered/scaled `clustering_data` plus a `clusters` column.
#' @keywords internal
kmeans_clustering <- function(clustering_data, kmid) {
  normmatrix <- t(scale(t(clustering_data), scale = TRUE, center = TRUE))

  n <- nrow(normmatrix)
  max_k <- n - 1 # kmeans/silhouette require k < number of distinct points
  if (max_k < 2) {
    stop("kmeans_clustering() needs at least 3 genes to cluster; got ", n, ".")
  }

  # test different clustering configurations, and use the silhouette index to pick the best one
  dist_mat <- stats::dist(normmatrix)
  threshold <- NULL
  s_index <- NULL
  row <- 1
  kint <- floor(kmid * 0.1) # how much to iterate k by
  kint <- ifelse(kint < 8, 8, kint) # min kmid is 20 so this works
  # clamp candidate k's to what this data can actually support (kmid is chosen from
  # the regulator matrix's dimensions, which can be larger than the number of genes
  # actually being clustered here)
  k_values <- unique(pmin(pmax(seq(kmid - 2 * kint, kmid + 2 * kint, by = kint), 2), max_k))
  for (k in k_values) {
    message(sprintf("Testing k=%d", k))
    km_all <- suppressWarnings(stats::kmeans(normmatrix, k, nstart = 5, iter.max = 20))
    s <- cluster::silhouette(km_all$cluster, dist_mat)
    s_sum <- summary(s)
    threshold[row] <- k
    s_index[row] <- s_sum$avg.width
    row <- row + 1
  }

  # save the best configuration, which is the largest silhouette index (first on ties)
  my_k <- threshold[s_index == max(s_index)][1]
  message(sprintf("Choose k=%s", my_k))
  km_all <- stats::kmeans(normmatrix, my_k, nstart = 50, iter.max = 20)
  clusters <- data.frame(clusters = km_all$cluster)

  as.data.frame(cbind(normmatrix, clusters))
}
