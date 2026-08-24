test_that("permute_network is reproducible regardless of num.cores", {
  skip_on_cran()
  mats <- make_test_matrices()

  perms_1core <- permute_network(mats$target, mats$reg, n_permutations = 3, num.cores = 1,
                                  weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                                  nb.trees = 50, trace = FALSE)
  perms_3core <- permute_network(mats$target, mats$reg, n_permutations = 3, num.cores = 3,
                                  weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                                  nb.trees = 50, trace = FALSE)

  expect_identical(perms_1core, perms_3core)
})

test_that("permute_network works with a clustered assignment that has singleton clusters", {
  # Regression test: shuffle_matrix() converts target/reg to plain matrices, and
  # subsetting a matrix to exactly one row without drop = FALSE silently collapses
  # it to a vector, breaking infer_network_clustered()'s dim() guard with
  # "missing value where TRUE/FALSE needed" the first time a cluster had one gene.
  mats <- make_test_matrices(n_targets = 6, n_regs = 4, n_samples = 8)
  # clusters 2 and 3 are singletons -- exactly the case that broke without drop = FALSE
  cluster_assignment <- data.frame(clusters = c(1, 1, 1, 1, 2, 3), row.names = rownames(mats$target))

  expect_no_error(
    permute_network(mats$target, mats$reg, cluster_assignment = cluster_assignment,
                     n_permutations = 2, num.cores = 1, weightthreshold = 0, normalize = FALSE,
                     connect_hubs = FALSE, nb.trees = 30, trace = FALSE)
  )
})

test_that("shuffle_matrix preserves dimnames and only reorders values", {
  mat <- matrix(1:12, nrow = 3, ncol = 4,
                dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3", "s4")))

  shuffled_col <- shuffle_matrix(mat, "col")
  expect_equal(dimnames(shuffled_col), dimnames(mat))
  # each column's set of values is preserved, just reordered across rows
  for (j in seq_len(ncol(mat))) {
    expect_setequal(unname(shuffled_col[, j]), unname(mat[, j]))
  }

  shuffled_row <- shuffle_matrix(mat, "row")
  expect_equal(dimnames(shuffled_row), dimnames(mat))
  for (i in seq_len(nrow(mat))) {
    expect_setequal(unname(shuffled_row[i, ]), unname(mat[i, ]))
  }
})
