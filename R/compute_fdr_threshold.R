#' @keywords internal
prep_weights_for_fdr <- function(weights) {
  # descending order so rank 1 is the strongest edge; zero-weight edges carry
  # no signal and would just be noise in the rank-by-rank comparison below.
  w <- sort(weights, decreasing = TRUE)
  w <- w[w > 0]
  # round to 4 significant figures so real and permuted weights that are
  # effectively equal don't fail to compare equal due to floating-point noise
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
#'   [permute_network()], produced against the same `real_network`.
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
  # the real network's sorted, nonzero weights -- one "rank" per edge, rank 1
  # being the strongest. This ranking, not the raw weight, is what each
  # permutation's weights get compared against below.
  truth <- prep_weights_for_fdr(real_network$Weight)
  n <- length(truth)

  # each permuted network's own sorted weights, padded with NA up to the real
  # network's length `n` -- a permutation that produced fewer nonzero edges
  # than the real network simply has no weight at the missing ranks, and NA
  # lets those ranks drop out of the comparison at that position instead of
  # recycling/misaligning shorter vectors against `truth`.
  perm_weights <- lapply(permuted_networks, function(net) {
    w <- prep_weights_for_fdr(net$Weight)
    if (length(w) < n) {
      w <- c(w, rep(NA_real_, n - length(w)))
    }
    w[seq_len(n)]
  })
  # rows = rank, columns = permutation, so perm_matrix[i, ] is every
  # permutation's weight at rank i, directly comparable to truth[i].
  perm_matrix <- do.call(cbind, perm_weights)

  # empirical p-value at each rank: the fraction of permutations whose weight
  # at that same rank meets or beats the real weight there (+1/+1 is the
  # standard continuity correction, so a rank can never score exactly p = 0
  # even if no permutation ever beat it).
  p_vals <- vapply(seq_len(n), function(i) {
    perm_i <- stats::na.omit(perm_matrix[i, ])
    (sum(perm_i > truth[i]) + 1) / (length(perm_i) + 1)
  }, numeric(1))

  # Benjamini-Hochberg across all ranks, converting each rank's raw p-value
  # into a false discovery rate for "keep every edge down to this rank."
  fdr <- stats::p.adjust(p_vals, method = "fdr")

  # walk down through ranks (decreasing weight) as far as FDR still stays
  # below target -- the threshold is the weight at the LAST (lowest-weight)
  # rank that still qualifies, since every edge above it also qualifies.
  below <- which(fdr < target_fdr)
  threshold <- if (length(below) > 0) truth[max(below)] else NA_real_

  thresholded_network <- if (!is.na(threshold)) {
    real_network[real_network$Weight >= threshold, , drop = FALSE]
  } else {
    # no rank ever reached the target FDR -- nothing survives the cutoff
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
