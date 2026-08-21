#' Plot the distribution of edge weights in a network
#'
#' @param network an edge table with a `Weight` column, e.g. the `network`
#'   element of [run_scion()]'s result.
#' @param bins number of histogram bins.
#' @return a `ggplot2` object.
#' @export
plot_weight_distribution <- function(network, bins = 100) {
  ggplot2::ggplot(network, ggplot2::aes(x = Weight)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::labs(x = "Edge weight", y = "Count") +
    ggplot2::theme_minimal()
}

#' Compute each regulator's out-degree (number of edges) in a network
#'
#' @param network an edge table with a `Regulator` column.
#' @return a data frame with `Regulator` and `outdegree` columns, one row per
#'   regulator that has at least one edge, sorted by decreasing out-degree.
#' @export
compute_outdegree <- function(network) {
  counts <- table(network$Regulator)
  result <- data.frame(Regulator = names(counts), outdegree = as.integer(counts),
                        stringsAsFactors = FALSE)
  result[order(-result$outdegree), , drop = FALSE]
}

#' Plot the distribution of regulator out-degrees in a network
#'
#' @param network an edge table with a `Regulator` column.
#' @param bins number of histogram bins.
#' @return a `ggplot2` object.
#' @export
plot_outdegree_distribution <- function(network, bins = 30) {
  outdegree <- compute_outdegree(network)
  ggplot2::ggplot(outdegree, ggplot2::aes(x = outdegree)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::labs(x = "Regulator out-degree (number of target edges)", y = "Count") +
    ggplot2::theme_minimal()
}
