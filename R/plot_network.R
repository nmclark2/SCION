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
#' @param legend if `TRUE` (default) and `interactive = TRUE`, attaches the
#'   HTML legend to the returned widget via `htmlwidgets::prependContent()`.
#'   Only takes effect for standalone use (R Markdown, `htmlwidgets::saveWidget()`)
#'   -- Shiny's `renderVisNetwork()` doesn't transmit prepended content to the
#'   client at all (and warns about it), so the Shiny app calls this with
#'   `legend = FALSE` and renders [network_legend_html()] as its own separate
#'   UI element instead (see `tab_visualize.R`). Ignored for `interactive = FALSE`,
#'   which always draws its own legend directly on the plot.
#' @param ... additional arguments passed to `igraph::plot.igraph()` (static)
#'   or `visNetwork::visNetwork()` (interactive).
#' @return for `interactive = FALSE`, invisibly returns the `igraph` object
#'   (after plotting it, with a legend, as a side effect); for
#'   `interactive = TRUE`, a `visNetwork` htmlwidget.
#' @export
plot_network <- function(edge_table, interactive = FALSE, legend = TRUE, ...) {
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
    graphics::mtext("Thicker edge = higher weight", side = 1, line = 4, cex = 0.8, col = "#555555")
    return(invisible(g))
  }

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("interactive = TRUE requires the 'visNetwork' package. Install it with install.packages('visNetwork').")
  }
  nodes <- data.frame(id = node_names, label = node_names, group = node_group, stringsAsFactors = FALSE)
  edges <- igraph::as_data_frame(g, what = "edges") # already has "from"/"to" columns
  edges$value <- rescale_for_plot(edges$Weight, 1, 10)

  vis <- visNetwork::visNetwork(nodes, edges, ...)
  vis <- visNetwork::visNodes(vis, font = list(size = 22, color = "#1a1a1a"))
  vis <- visNetwork::visGroups(vis, groupname = "Regulator", shape = "square",
                                color = list(background = regulator_color, border = "#2c4a75",
                                             highlight = "#6f97d6"))
  vis <- visNetwork::visGroups(vis, groupname = "Target", shape = "dot",
                                color = list(background = target_color, border = "#a8622f",
                                             highlight = "#f0a878"))
  vis <- visNetwork::visInteraction(vis, navigationButtons = TRUE, keyboard = TRUE, zoomView = TRUE,
                                     dragView = TRUE)
  vis <- visNetwork::visEdges(vis, arrows = "to")
  # visLegend() renders its own separate, independently zoomable vis.js canvas --
  # at typical widget sizes that reads as a comically oversized, pointlessly
  # interactive legend. A plain static HTML/CSS caption has no such quirks.
  # NOTE: htmlwidgets::prependContent() correctly attaches this to the widget
  # object (confirmed via htmlwidgets::saveWidget()), but Shiny's
  # renderVisNetwork()/visNetworkOutput() binding does not transmit prepend/
  # append content to the client at all -- it's silently dropped, and as of
  # recent htmlwidgets versions also emits "Ignoring prepended content;
  # prependContent can't be used in a Shiny render call" to the console. It
  # still renders fine for standalone/non-Shiny use (R Markdown,
  # saveWidget()), so it's still applied by default here -- but the Shiny app
  # passes legend = FALSE and instead renders network_legend_html() as its
  # own independent UI element (see tab_visualize.R), so this call never runs
  # (and never warns) from inside the app.
  if (isTRUE(legend)) {
    vis <- htmlwidgets::prependContent(vis, network_legend_html(regulator_color, target_color))
  }
  vis
}

#' @keywords internal
network_legend_html <- function(regulator_color = "#4C72B0", target_color = "#DD8452") {
  htmltools::div(
    style = paste("display:flex; flex-wrap:wrap; align-items:center; gap:16px;",
                  "font-size:13px; color:#333; padding:4px 2px 8px 2px;"),
    htmltools::div(
      style = "display:flex; align-items:center; gap:6px;",
      htmltools::tags$span(style = sprintf(
        "display:inline-block; width:13px; height:13px; background:%s;", regulator_color
      )),
      "Regulator"
    ),
    htmltools::div(
      style = "display:flex; align-items:center; gap:6px;",
      htmltools::tags$span(style = sprintf(
        "display:inline-block; width:13px; height:13px; border-radius:50%%; background:%s;", target_color
      )),
      "Target"
    ),
    htmltools::div(style = "color:#555;", "Edge width -> weight (thicker = higher)")
  )
}

#' @keywords internal
rescale_for_plot <- function(x, lo, hi) {
  if (length(unique(x)) <= 1) {
    return(rep(mean(c(lo, hi)), length(x)))
  }
  # rank-based (percentile), not a raw linear min-max: a SCION network is
  # normally already thresholded (by weight cutoff or FDR), so the *retained*
  # edges' raw weights are naturally clustered in a narrow high band -- a
  # linear rescale of that narrow band squashes nearly all of them toward one
  # end of [lo, hi], making "thicker = higher weight" invisible for all but a
  # few outliers. Scaling by each edge's percentile rank instead spreads the
  # same edges evenly across the full width range regardless of how skewed
  # the raw weights are.
  pct <- (rank(x, ties.method = "average") - 1) / (length(x) - 1)
  lo + (hi - lo) * pct
}
