test_that("cluster_genes derives clustering_data from target_data + reg_data when not supplied", {
  mats <- make_test_matrices(n_targets = 20, n_regs = 6, n_samples = 10)

  result <- cluster_genes(method = "kmeans", target_data = mats$target, reg_data = mats$reg)

  expect_true(is.data.frame(result))
  expect_setequal(rownames(result), c(rownames(mats$target), rownames(mats$reg)))
})

test_that("derived clustering_data gives regulator data precedence for overlapping genes", {
  mats <- make_test_matrices(n_targets = 20, n_regs = 6, n_samples = 10)
  # make one gene both a target and a regulator, with DIFFERENT values in each
  overlap_gene <- rownames(mats$reg)[1]
  target_with_overlap <- mats$target
  rownames(target_with_overlap)[1] <- overlap_gene
  target_with_overlap[1, ] <- target_with_overlap[1, ] + 1000 # distinguishably different

  result <- cluster_genes(method = "kmeans", target_data = target_with_overlap, reg_data = mats$reg)

  # exactly one row for the overlapping gene, using the regulator matrix's values
  expect_equal(sum(rownames(result) == overlap_gene), 1)
})

test_that("cluster_genes errors informatively when clustering_data, target_data, and reg_data are all missing", {
  expect_error(cluster_genes(method = "kmeans"), "clustering_data")
})

test_that("cluster_genes(method = 'none') ignores missing clustering_data entirely", {
  expect_null(cluster_genes(method = "none"))
})
