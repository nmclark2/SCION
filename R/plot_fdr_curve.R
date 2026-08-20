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
#' @return a `ggplot2` object.
#' @export
plot_fdr_curve <- function(fdr_result, type = c("curve", "fdr_hist", "weight_comparison"),
                            permutation_index = 1) {
  type <- match.arg(type)
  curve <- fdr_result$curve

  if (type == "curve") {
    p <- ggplot2::ggplot(curve, ggplot2::aes(x = weight, y = fdr)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Edge weight", y = "FDR") +
      ggplot2::theme_minimal()
    if (!is.na(fdr_result$threshold)) {
      p <- p +
        ggplot2::geom_vline(xintercept = fdr_result$threshold, linetype = "dashed", color = "red") +
        ggplot2::geom_hline(yintercept = fdr_result$target_fdr, linetype = "dashed", color = "red")
    }
    return(p)
  }

  if (type == "fdr_hist") {
    return(
      ggplot2::ggplot(curve, ggplot2::aes(x = fdr)) +
        ggplot2::geom_histogram(bins = 100) +
        ggplot2::labs(x = "FDR", y = "Count") +
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
    ggplot2::labs(x = "Edge weight", y = "Density", fill = NULL) +
    ggplot2::theme_minimal()
}
