test_that("permute_network is reproducible regardless of num.cores", {
  skip_on_cran()
  mats <- make_test_matrices()

  perms_1core <- permute_network(mats$target, mats$reg, n_permutations = 3, num.cores = 1,
                                  weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                                  engine = "randomForest", nb.trees = 50, trace = FALSE)
  perms_3core <- permute_network(mats$target, mats$reg, n_permutations = 3, num.cores = 3,
                                  weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                                  engine = "randomForest", nb.trees = 50, trace = FALSE)

  expect_identical(perms_1core, perms_3core)
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
