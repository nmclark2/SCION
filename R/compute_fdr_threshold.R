#' @keywords internal
prep_weights_for_fdr <- function(weights) {
  w <- sort(weights, decreasing = TRUE)
  w <- w[w > 0]
  signif(w, 4)
}

#' Compute a permutation-based FDR and edge-weight cutoff
#'
#' Ports the validated permutation-testing procedure exactly: for each rank
#' position in the sorted, nonzero real network weights, computes an empirical
#' p-value from how often permuted networks exceed the real weight at that
#' rank, BH-adjusts across ranks, and picks the smallest real weight for which
#' FDR is still below `target_fdr`.
#'
#' @param real_network the real network's edge table (must have a `Weight`
#'   column), e.g. the `network` element of [run_scion()]'s result.
#' @param permuted_networks a list of permuted edge tables, as returned by
#'   [permute_network()]. **Must have been produced with the same `engine`** as
#'   `real_network` -- weights from different engines are not on comparable
#'   scales, which would invalidate this rank-matching.
#' @param target_fdr the FDR threshold to select a weight cutoff at (default 0.05).
#' @return a list:
#'   \item{curve}{a data frame with one row per ranked real-network weight:
#'     `weight`, `p_value`, `fdr` -- the input to [plot_fdr_curve()].}
#'   \item{threshold}{the chosen edge-weight cutoff, or `NA` if no rank achieves
#'     `target_fdr`.}
#'   \item{target_fdr}{echoes the input.}
#'   \item{thresholded_network}{`real_network` filtered to `Weight >= threshold`.}
#'   \item{permuted_weights}{the rank x permutation matrix of (padded) permuted
#'     weights, useful for diagnostic plots.}
#' @export
compute_fdr_threshold <- function(real_network, permuted_networks, target_fdr = 0.05) {
  truth <- prep_weights_for_fdr(real_network$Weight)
  n <- length(truth)

  perm_weights <- lapply(permuted_networks, function(net) {
    w <- prep_weights_for_fdr(net$Weight)
    if (length(w) < n) {
      w <- c(w, rep(NA_real_, n - length(w)))
    }
    w[seq_len(n)]
  })
  perm_matrix <- do.call(cbind, perm_weights)

  p_vals <- vapply(seq_len(n), function(i) {
    perm_i <- stats::na.omit(perm_matrix[i, ])
    (sum(perm_i > truth[i]) + 1) / (length(perm_i) + 1)
  }, numeric(1))

  fdr <- stats::p.adjust(p_vals, method = "fdr")

  below <- which(fdr < target_fdr)
  threshold <- if (length(below) > 0) truth[max(below)] else NA_real_

  thresholded_network <- if (!is.na(threshold)) {
    real_network[real_network$Weight >= threshold, , drop = FALSE]
  } else {
    real_network[0, , drop = FALSE]
  }

  list(
    curve = data.frame(weight = truth, p_value = p_vals, fdr = fdr),
    threshold = threshold,
    target_fdr = target_fdr,
    thresholded_network = thresholded_network,
    permuted_weights = perm_matrix
  )
}
