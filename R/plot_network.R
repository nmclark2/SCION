#' Plot a SCION network
#'
#' A quick in-R sanity-check visualization -- for anything beyond a few dozen
#' nodes, or for a publication-quality figure, importing the edge table into
#' Cytoscape (see [write_scion_network()]) will still give better layout
#' control. Both the static and interactive versions use the same visual
#' encoding: regulators are blue squares, targets are orange circles, and edge
#' width is proportional to weight.
#'
#' @param edge_table an edge table with `Regulator`, `Target`, `Weight` columns
#'   (e.g. the `network` element of [run_scion()]'s result, or a
#'   [compute_fdr_threshold()] result's `thresholded_network`).
#' @param interactive if `FALSE` (default), plots via base `igraph` plotting
#'   (static). If `TRUE`, returns an interactive `visNetwork` HTML widget
#'   (requires the optional `visNetwork` package) with on-canvas zoom/pan
#'   controls.
#' @param ... additional arguments passed to `igraph::plot.igraph()` (static)
#'   or `visNetwork::visNetwork()` (interactive).
#' @return for `interactive = FALSE`, invisibly returns the `igraph` object
#'   (after plotting it, with a legend, as a side effect); for
#'   `interactive = TRUE`, a `visNetwork` htmlwidget.
#' @export
plot_network <- function(edge_table, interactive = FALSE, ...) {
  g <- igraph::graph_from_data_frame(edge_table[, c("Regulator", "Target", "Weight")], directed = TRUE)
  weights <- igraph::E(g)$Weight

  node_names <- igraph::V(g)$name
  regulator_names <- unique(edge_table$Regulator)
  node_group <- ifelse(node_names %in% regulator_names, "Regulator", "Target")
  regulator_color <- "#4C72B0"
  target_color <- "#DD8452"
  # stored as a real vertex attribute (not just a plotting parameter) so the
  # returned graph is self-describing -- e.g. igraph::V(g)$group -- and the
  # shape/color mapping used for the plot is inspectable/reproducible after
  igraph::V(g)$group <- node_group

  if (!interactive) {
    edge_width <- rescale_for_plot(weights, 1, 5)
    igraph::V(g)$shape <- ifelse(node_group == "Regulator", "square", "circle")
    igraph::V(g)$color <- ifelse(node_group == "Regulator", regulator_color, target_color)
    igraph::plot.igraph(g, edge.width = edge_width, vertex.label.cex = 0.7, vertex.label.color = "#1a1a1a",
                         vertex.size = 8, vertex.frame.color = "#444444", edge.arrow.size = 0.4, ...)
    graphics::legend("topright", legend = c("Regulator", "Target"), pch = c(15, 19),
                      col = c(regulator_color, target_color), bty = "n", title = "Node type")
    graphics::mtext("Edge width proportional to weight", side = 1, line = 4, cex = 0.8, col = "#555555")
    return(invisible(g))
  }

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("interactive = TRUE requires the 'visNetwork' package. Install it with install.packages('visNetwork').")
  }
  nodes <- data.frame(id = node_names, label = node_names, group = node_group, stringsAsFactors = FALSE)
  edges <- igraph::as_data_frame(g, what = "edges") # already has "from"/"to" columns
  edges$value <- rescale_for_plot(edges$Weight, 1, 10)

  vis <- visNetwork::visNetwork(
    nodes, edges,
    submain = list(text = "Edge width proportional to weight",
                    style = "font-size:13px;color:#555555;font-weight:normal;"),
    ...
  )
  vis <- visNetwork::visNodes(vis, font = list(size = 16, color = "#1a1a1a"))
  vis <- visNetwork::visGroups(vis, groupname = "Regulator", shape = "square",
                                color = list(background = regulator_color, border = "#2c4a75",
                                             highlight = "#6f97d6"))
  vis <- visNetwork::visGroups(vis, groupname = "Target", shape = "dot",
                                color = list(background = target_color, border = "#a8622f",
                                             highlight = "#f0a878"))
  vis <- visNetwork::visLegend(vis, useGroups = TRUE, main = "Node type", width = 0.3)
  vis <- visNetwork::visInteraction(vis, navigationButtons = TRUE, keyboard = TRUE, zoomView = TRUE,
                                     dragView = TRUE)
  visNetwork::visEdges(vis, arrows = "to")
}

#' @keywords internal
rescale_for_plot <- function(x, lo, hi) {
  if (length(unique(x)) <= 1) {
    return(rep(mean(c(lo, hi)), length(x)))
  }
  lo + (hi - lo) * (x - min(x)) / (max(x) - min(x))
}
