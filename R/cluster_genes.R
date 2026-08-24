#' Cluster genes prior to network inference
#'
#' Computes a fixed gene-to-cluster assignment used by [infer_network()]. This
#' assignment is computed once from the real, unpermuted data; when running
#' permutations via [permute_network()], the same assignment must be reused
#' rather than recomputed, so that permuted networks are compared against the
#' real network on identical gene groupings.
#'
#' @param clustering_data data frame or matrix of expression data to cluster on
#'   (rows = genes, columns = samples). Ignored when `method = "none"`. If
#'   `NULL` and `method` needs data (`"dtw"`/`"ica"`/`"kmeans"`), it is derived
#'   from `target_data`/`reg_data` -- see Details.
#' @param method clustering method: `"none"` (no clustering, one network for all
#'   genes), `"dtw"` (temporal, dynamic time warping), `"ica"` (non-temporal,
#'   independent component analysis), `"kmeans"` (non-temporal, k-means), or
#'   `"upload"` (use a pre-computed clusters file).
#' @param threshold clustering threshold used by `"dtw"` and `"ica"`; ignored by
#'   `"kmeans"` and `"upload"`.
#' @param clusters_file path to a delimited text file (`.csv`/`.tsv`/`.txt`/`.ssv`,
#'   dispatched by extension) of pre-computed clusters (first column = gene
#'   names, a `clusters` column with cluster numbers). Required when
#'   `method = "upload"`.
#' @param target_data,reg_data the processed target/regulator matrices (rows =
#'   genes, columns = samples). `reg_data` alone is also used by `method =
#'   "kmeans"` to pick a starting `k` (see Details); `target_data` is only used
#'   to derive `clustering_data` when it isn't supplied.
#' @return a data frame with a `clusters` column (row names = gene names), or
#'   `NULL` when `method = "none"`.
#' @details
#' When `clustering_data` is `NULL` and clustering is requested, it defaults to
#' `reg_data` plus whatever rows of `target_data` aren't already in `reg_data`
#' (regulator data takes precedence for genes that are both a target and a
#' regulator) -- i.e. the same genes that will be used for network inference,
#' combined into one matrix. This matches the historical PANOPLY behavior of
#' auto-deriving a clustering matrix from the target/regulator data when no
#' separate clustering file is provided.
#' @export
cluster_genes <- function(clustering_data = NULL, method = c("none", "dtw", "ica", "kmeans", "upload"),
                           threshold = 0.5, clusters_file = NULL, target_data = NULL, reg_data = NULL) {
  method <- match.arg(method)

  if (method == "none") {
    return(NULL)
  }

  message("SCION_STAGE: clustering started")

  result <- if (method == "upload") {
    clusters <- read_delimited_matrix(clusters_file)
    rownames(clusters) <- make.names(rownames(clusters))
    clusters
  } else {
    if (is.null(clustering_data)) {
      if (is.null(target_data) || is.null(reg_data)) {
        stop("clustering_data was not supplied, and target_data/reg_data are required to derive ",
             "it automatically for method = '", method, "'.")
      }
      message("No clustering_data supplied; combining target_data and reg_data for clustering.")
      target_data <- as.data.frame(target_data)
      reg_data <- as.data.frame(reg_data)
      target_only <- target_data[!rownames(target_data) %in% rownames(reg_data), , drop = FALSE]
      clustering_data <- rbind(reg_data, target_only)
    }

    if (method == "dtw") {
      dtw_clustering(clustering_data, threshold)
    } else if (method == "ica") {
      ica_clustering(clustering_data, threshold)
    } else {
      # kmeans: scale starting k based on the smallest dimension of the regulator matrix
      dims_source <- if (!is.null(reg_data)) reg_data else clustering_data
      kmid <- min(dim(dims_source)[1], dim(dims_source)[2])
      kmid <- ifelse(kmid < 20, 20, kmid) # if kmid < 20, change to 20 for a better starting number
      kmeans_clustering(clustering_data, kmid)
    }
  }

  message("SCION_STAGE: clustering complete")
  result
}
