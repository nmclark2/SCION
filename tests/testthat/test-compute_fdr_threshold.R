test_that("compute_fdr_threshold's p-values/FDR match the underlying formula directly", {
  real_network <- data.frame(Weight = c(10, 8, 6, 4, 2))
  # every permutation's weights are below every real weight at every rank
  permuted_networks <- replicate(3, data.frame(Weight = rep(1, 5)), simplify = FALSE)

  result <- compute_fdr_threshold(real_network, permuted_networks, target_fdr = 0.05)

  expected_p <- rep((0 + 1) / (3 + 1), 5) # sum(perm > truth) == 0 at every rank
  expected_fdr <- stats::p.adjust(expected_p, method = "fdr")

  expect_equal(result$curve$weight, c(10, 8, 6, 4, 2))
  expect_equal(result$curve$p_value, expected_p)
  expect_equal(result$curve$fdr, expected_fdr)
  # none of the ranks reach FDR < 0.05 here, so no threshold is chosen
  expect_true(all(expected_fdr >= 0.05))
  expect_true(is.na(result$threshold))
  expect_equal(nrow(result$thresholded_network), 0)
})

test_that("compute_fdr_threshold selects a threshold and filters the real network by it", {
  real_network <- data.frame(Weight = c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1))
  set.seed(42)
  permuted_networks <- replicate(
    20, data.frame(Weight = sort(stats::runif(10, 0, 3), decreasing = TRUE)), simplify = FALSE
  )

  result <- compute_fdr_threshold(real_network, permuted_networks, target_fdr = 0.05)

  expect_false(is.na(result$threshold))
  expect_true(all(result$thresholded_network$Weight >= result$threshold))
  expect_true(nrow(result$thresholded_network) <= nrow(real_network))
})

test_that("shorter permutation weight vectors are NA-padded, not recycled", {
  real_network <- data.frame(Weight = c(5, 4, 3, 2, 1))
  permuted_networks <- list(data.frame(Weight = c(0.5, 0.4))) # much shorter than real

  result <- compute_fdr_threshold(real_network, permuted_networks, target_fdr = 0.05)

  expect_equal(nrow(result$curve), 5)
  expect_true(all(is.na(result$permuted_weights[3:5, 1])))
})
