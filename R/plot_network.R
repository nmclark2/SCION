#' Plot a SCION network
#'
#' A quick in-R sanity-check visualization -- for anything beyond a few dozen
#' nodes, or for a publication-quality figure, importing the edge table into
#' Cytoscape (see [write_scion_network()]) will still give better layout
#' control.
#'
#' @param edge_table an edge table with `Regulator`, `Target`, `Weight` columns
#'   (e.g. the `network` element of [run_scion()]'s result, or a
#'   [compute_fdr_threshold()] result's `thresholded_network`).
#' @param interactive if `FALSE` (default), plots via base `igraph` plotting
#'   (static). If `TRUE`, returns an interactive `visNetwork` HTML widget
#'   (requires the optional `visNetwork` package).
#' @param ... additional arguments passed to `igraph::plot.igraph()` (static)
#'   or `visNetwork::visNetwork()` (interactive).
#' @return for `interactive = FALSE`, invisibly returns the `igraph` object
#'   (after plotting it as a side effect); for `interactive = TRUE`, a
#'   `visNetwork` htmlwidget.
#' @export
plot_network <- function(edge_table, interactive = FALSE, ...) {
  g <- igraph::graph_from_data_frame(edge_table[, c("Regulator", "Target", "Weight")], directed = TRUE)
  weights <- igraph::E(g)$Weight

  if (!interactive) {
    edge_width <- rescale_for_plot(weights, 1, 5)
    igraph::plot.igraph(g, edge.width = edge_width, vertex.label.cex = 0.7, vertex.size = 8,
                         edge.arrow.size = 0.4, ...)
    return(invisible(g))
  }

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("interactive = TRUE requires the 'visNetwork' package. Install it with install.packages('visNetwork').")
  }
  nodes <- data.frame(id = igraph::V(g)$name, label = igraph::V(g)$name)
  edges <- igraph::as_data_frame(g, what = "edges") # already has "from"/"to" columns
  edges$value <- rescale_for_plot(edges$Weight, 1, 10)
  visNetwork::visEdges(visNetwork::visNetwork(nodes, edges, ...), arrows = "to")
}

#' @keywords internal
rescale_for_plot <- function(x, lo, hi) {
  if (length(unique(x)) <= 1) {
    return(rep(mean(c(lo, hi)), length(x)))
  }
  lo + (hi - lo) * (x - min(x)) / (max(x) - min(x))
}
