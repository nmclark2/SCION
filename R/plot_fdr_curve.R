#' Plot diagnostics from a permutation-based FDR result
#'
#' @param fdr_result the output of [compute_fdr_threshold()].
#' @param type `"curve"` (default; weight vs. FDR scatter with the chosen
#'   threshold marked), `"fdr_hist"` (histogram of FDR values across ranks), or
#'   `"weight_comparison"` (real vs. one permuted network's weight distribution
#'   overlay).
#' @param permutation_index which permuted network's weights to show for
#'   `type = "weight_comparison"` (default 1; matches `permuted_networks[[1]]`
#'   used to build `fdr_result`).
#' @param title optional plot title -- see [plot_weight_distribution()]'s `title` argument.
#' @param threshold the edge-weight cutoff to mark on `type = "curve"` (a
#'   dashed vertical line, labeled with its value) -- defaults to
#'   `fdr_result$threshold`, but the caller can pass a different value (e.g. a
#'   manual cutoff that overrides the FDR-based one) so the marked line always
#'   reflects whichever cutoff is actually in effect. `NULL`/`NA` draws no line.
#' @param show_target_fdr_line if `TRUE` (default), also draws a dashed
#'   horizontal line at `fdr_result$target_fdr`. Set `FALSE` when `threshold`
#'   is a manual override, since the FDR target is no longer the criterion
#'   that produced it and showing it would be misleading.
#' @param label_lines if `FALSE` (default), the threshold/target-FDR lines
#'   carry their value as plotly hover text only (via [plotly::ggplotly()]) --
#'   appropriate for the interactive app, where hovering is how you'd read a
#'   value off any other point on the plot too. If `TRUE`, also draws the
#'   value as permanent on-plot text, for contexts with no hover available
#'   (the static PDF export).
#' @return a `ggplot2` object.
#' @export
plot_fdr_curve <- function(fdr_result, type = c("curve", "fdr_hist", "weight_comparison"),
                            permutation_index = 1, title = NULL,
                            threshold = fdr_result$threshold, show_target_fdr_line = TRUE,
                            label_lines = FALSE) {
  type <- match.arg(type)
  curve <- fdr_result$curve

  if (type == "curve") {
    p <- ggplot2::ggplot(curve, ggplot2::aes(x = weight, y = fdr)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Edge weight", y = "FDR", title = title) +
      ggplot2::theme_minimal()
    if (!is.null(threshold) && !is.na(threshold)) {
      threshold_label <- paste0("Threshold = ", signif(threshold, 4))
      p <- hover_line(p, "vertical", threshold, threshold_label)
      if (isTRUE(label_lines)) {
        p <- p + ggplot2::annotate("text", x = threshold, y = max(curve$fdr, na.rm = TRUE),
                                    label = threshold_label, hjust = -0.05, vjust = -0.5, color = "red",
                                    size = 3, fontface = "bold")
      }
    }
    if (isTRUE(show_target_fdr_line)) {
      fdr_label <- paste0("Target FDR = ", fdr_result$target_fdr)
      p <- hover_line(p, "horizontal", fdr_result$target_fdr, fdr_label)
      if (isTRUE(label_lines)) {
        # right-aligned at the rightmost data point -- the left side of this
        # plot is where the real/permuted weights are densest, so a label
        # there reliably overlapped points; the right side is comparatively
        # empty.
        p <- p + ggplot2::annotate("text", x = max(curve$weight, na.rm = TRUE), y = fdr_result$target_fdr,
                                    label = fdr_label, hjust = 1, vjust = -0.5, color = "red",
                                    size = 3, fontface = "bold")
      }
    }
    return(p)
  }

  if (type == "fdr_hist") {
    return(
      ggplot2::ggplot(curve, ggplot2::aes(x = fdr)) +
        ggplot2::geom_histogram(bins = 100) +
        ggplot2::labs(x = "FDR", y = "Count", title = title) +
        ggplot2::theme_minimal()
    )
  }

  # weight_comparison
  real_weights <- curve$weight
  perm_weights <- stats::na.omit(fdr_result$permuted_weights[, permutation_index])
  df <- data.frame(
    weight = c(real_weights, perm_weights),
    source = c(rep("Real network", length(real_weights)), rep("Permutation", length(perm_weights)))
  )
  ggplot2::ggplot(df, ggplot2::aes(x = weight, fill = source)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), alpha = 0.5, bins = 100,
                             position = "identity") +
    ggplot2::labs(x = "Edge weight", y = "Density", fill = NULL, title = title) +
    ggplot2::theme_minimal()
}
