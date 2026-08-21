test_that("kmeans_clustering clamps candidate k's to what the data supports", {
  # Regression test: kmid can come from the regulator matrix's dimensions (per
  # cluster_genes()'s documented behavior) and be much larger than the number
  # of genes actually being clustered here, which used to error with
  # "more cluster centers than distinct data points".
  mats <- make_test_matrices(n_targets = 20, n_regs = 6, n_samples = 10)
  small_clustering_data <- mats$target[1:10, ] # only 10 genes, but kmid will be 20+

  result <- SCION:::kmeans_clustering(small_clustering_data, kmid = 20)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 10)
  expect_true(max(result$clusters) <= 9) # k must stay < nrow
})

test_that("kmeans_clustering errors informatively with too few genes", {
  mats <- make_test_matrices(n_targets = 2, n_regs = 2, n_samples = 5)
  expect_error(SCION:::kmeans_clustering(mats$target, kmid = 20), "at least 3 genes")
})
