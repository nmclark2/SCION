#' Plot the distribution of edge weights in a network
#'
#' @param network an edge table with a `Weight` column, e.g. the `network`
#'   element of [run_scion()]'s result. Always pass the full, unthresholded
#'   network here (not an already-filtered one) so `cutoff` is visible in
#'   context against the whole distribution.
#' @param bins number of histogram bins.
#' @param cutoff optional edge-weight cutoff (FDR-based or manual) to mark
#'   with a dashed vertical line. `NULL` (default) or `NA` draws no line.
#' @param title optional plot title -- `NULL` (default) omits it. The
#'   interactive app view skips this (its box header already names the plot);
#'   [save_diagnostic_plots()] and the app's PDF export set one, since an
#'   exported plot has no other title to identify it by.
#' @param label_line if `FALSE` (default), the cutoff line carries its value
#'   as plotly hover text only (via [plotly::ggplotly()]) -- appropriate for
#'   the interactive app. If `TRUE`, also draws the value as permanent on-plot
#'   text, for contexts with no hover available (the static PDF export).
#' @return a `ggplot2` object.
#' @export
plot_weight_distribution <- function(network, bins = 100, cutoff = NULL, title = NULL, label_line = FALSE) {
  p <- ggplot2::ggplot(network, ggplot2::aes(x = Weight)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::labs(x = "Edge weight", y = "Count", title = title) +
    ggplot2::theme_minimal() +
    # x = 0 is always in view even when every edge already clears some
    # weightthreshold applied upstream (at inference time) -- otherwise the
    # axis auto-scales tightly to the data and a cutoff drawn right at its
    # lower edge looks like the plot itself starts there, not at 0.
    ggplot2::expand_limits(x = 0)
  if (!is.null(cutoff) && !is.na(cutoff)) {
    cutoff_label <- paste0("Threshold = ", signif(cutoff, 4))
    p <- hover_line(p, "vertical", cutoff, cutoff_label)
    if (isTRUE(label_line)) {
      p <- p + ggplot2::annotate("text", x = cutoff, y = Inf, label = cutoff_label,
                                  hjust = -0.05, vjust = 1.2, color = "red", size = 3, fontface = "bold")
    }
  }
  p
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
#' @param title optional plot title -- see [plot_weight_distribution()]'s `title` argument.
#' @return a `ggplot2` object.
#' @export
plot_outdegree_distribution <- function(network, bins = 30, title = NULL) {
  outdegree <- compute_outdegree(network)
  ggplot2::ggplot(outdegree, ggplot2::aes(x = outdegree)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::labs(x = "Regulator out-degree (number of target edges)", y = "Count", title = title) +
    ggplot2::theme_minimal()
}
