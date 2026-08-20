#' Cluster genes by k-means (non-temporal data)
#'
#' Tests several values of k around `kmid` and picks the one with the largest
#' average silhouette width.
#'
#' @param clustering_data data frame or matrix, rows = genes, columns = samples.
#' @param kmid starting/center value of k to test around.
#' @return a data frame: centered/scaled `clustering_data` plus a `clusters` column.
#' @keywords internal
kmeans_clustering <- function(clustering_data, kmid) {
  normmatrix <- t(scale(t(clustering_data), scale = TRUE, center = TRUE))

  # test different clustering configurations, and use the silhouette index to pick the best one
  dist_mat <- stats::dist(normmatrix)
  threshold <- NULL
  s_index <- NULL
  row <- 1
  kint <- floor(kmid * 0.1) # how much to iterate k by
  kint <- ifelse(kint < 8, 8, kint) # min kmid is 20 so this works
  for (k in seq(kmid - 2 * kint, kmid + 2 * kint, by = kint)) {
    message(sprintf("Testing k=%d", k))
    km_all <- suppressWarnings(stats::kmeans(normmatrix, k, nstart = 5, iter.max = 20))
    s <- cluster::silhouette(km_all$cluster, dist_mat)
    s_sum <- summary(s)
    threshold[row] <- k
    s_index[row] <- s_sum$avg.width
    row <- row + 1
  }

  # save the best configuration, which is the largest silhouette index
  my_k <- threshold[s_index == max(s_index)]
  message(sprintf("Choose k=%s", my_k))
  km_all <- stats::kmeans(normmatrix, my_k, nstart = 50, iter.max = 20)
  clusters <- data.frame(clusters = km_all$cluster)

  as.data.frame(cbind(normmatrix, clusters))
}
