#' Cluster genes prior to network inference
#'
#' Computes a fixed gene-to-cluster assignment used by [infer_network()]. This
#' assignment is computed once from the real, unpermuted data; when running
#' permutations via [permute_network()], the same assignment must be reused
#' rather than recomputed, so that permuted networks are compared against the
#' real network on identical gene groupings.
#'
#' @param clustering_data data frame or matrix of expression data to cluster on
#'   (rows = genes, columns = samples). Ignored when `method = "none"`.
#' @param method clustering method: `"none"` (no clustering, one network for all
#'   genes), `"dtw"` (temporal, dynamic time warping), `"ica"` (non-temporal,
#'   independent component analysis), `"kmeans"` (non-temporal, k-means), or
#'   `"upload"` (use a pre-computed clusters file).
#' @param threshold clustering threshold used by `"dtw"` and `"ica"`; ignored by
#'   `"kmeans"` and `"upload"`.
#' @param clusters_file path to a CSV of pre-computed clusters (first column =
#'   gene names, a `clusters` column with cluster numbers). Required when
#'   `method = "upload"`.
#' @param reg_data the processed regulator matrix (rows = genes, columns =
#'   samples), used only by `method = "kmeans"` to pick a starting `k`.
#'   Historically SCION scales k off the *regulator* matrix's dimensions, not
#'   the clustering matrix's -- preserved here for behavior-compatibility.
#'   Falls back to `clustering_data`'s own dimensions if not supplied.
#' @return a data frame with a `clusters` column (row names = gene names), or
#'   `NULL` when `method = "none"`.
#' @export
cluster_genes <- function(clustering_data, method = c("none", "dtw", "ica", "kmeans", "upload"),
                           threshold = 0.5, clusters_file = NULL, reg_data = NULL) {
  method <- match.arg(method)

  if (method == "none") {
    return(NULL)
  }

  if (method == "upload") {
    clusters <- utils::read.csv(clusters_file, row.names = 1)
    rownames(clusters) <- make.names(rownames(clusters))
    return(clusters)
  }

  if (method == "dtw") {
    return(dtw_clustering(clustering_data, threshold))
  }

  if (method == "ica") {
    return(ica_clustering(clustering_data, threshold))
  }

  # kmeans: scale starting k based on the smallest dimension of the regulator matrix
  dims_source <- if (!is.null(reg_data)) reg_data else clustering_data
  kmid <- min(dim(dims_source)[1], dim(dims_source)[2])
  kmid <- ifelse(kmid < 20, 20, kmid) # if kmid < 20, change to 20 for a better starting number
  kmeans_clustering(clustering_data, kmid)
}
