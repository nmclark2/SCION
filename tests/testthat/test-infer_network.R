test_that("infer_network (no clustering) is invariant to num.cores given the same seed", {
  skip_on_cran()
  mats <- make_test_matrices()

  net_1core <- infer_network(mats$target, mats$reg, weightthreshold = 0, normalize = FALSE,
                              num.cores = 1, engine = "randomForest", nb.trees = 50, trace = FALSE,
                              seed = 123)
  net_3core <- infer_network(mats$target, mats$reg, weightthreshold = 0, normalize = FALSE,
                              num.cores = 3, engine = "randomForest", nb.trees = 50, trace = FALSE,
                              seed = 123)

  expect_equal(net_1core, net_3core)
})

test_that("infer_network with MULTIPLE clusters is invariant to num.cores given the same seed", {
  skip_on_cran()
  # Regression test: a single RS.Get.Weight.Matrix() call being num.cores-invariant is not
  # enough on its own -- cluster 2+'s own seed draw must not depend on ambient RNG state left
  # behind by cluster 1, which used to differ between forked (num.cores > 2) and serial
  # execution of the *previous* cluster.
  mats <- make_test_matrices(n_targets = 16, n_regs = 6, n_samples = 10)
  cluster_assignment <- data.frame(clusters = rep(1:2, each = 8), row.names = rownames(mats$target))

  net_1core <- infer_network(mats$target, mats$reg, cluster_assignment = cluster_assignment,
                              weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                              num.cores = 1, engine = "randomForest", nb.trees = 50, trace = FALSE,
                              seed = 123)
  net_3core <- infer_network(mats$target, mats$reg, cluster_assignment = cluster_assignment,
                              weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                              num.cores = 3, engine = "randomForest", nb.trees = 50, trace = FALSE,
                              seed = 123)

  expect_equal(net_1core, net_3core)
})

test_that("infer_network skips clusters with exactly one regulator (NA %IncMSE) without warning", {
  skip_on_cran()
  # randomForest's %IncMSE importance is deterministically NA with exactly one
  # regulator, regardless of the data -- a cluster with only 1 regulator must be
  # skipped outright rather than silently producing an all-NA weight matrix
  # (which used to surface as a confusing "no non-missing arguments to min/max"
  # warning and a suspiciously small final network).
  mats <- make_test_matrices(n_targets = 8, n_regs = 6, n_samples = 10)
  cluster_assignment <- data.frame(
    clusters = c(rep(1, 4), rep(2, 4), 1, rep(2, 5)),
    row.names = c(rownames(mats$target), "reg1", paste0("reg", 2:6))
  )

  expect_no_warning(
    net <- infer_network(mats$target, mats$reg, cluster_assignment = cluster_assignment,
                          weightthreshold = 0, normalize = TRUE, connect_hubs = FALSE,
                          num.cores = 1, engine = "randomForest", nb.trees = 20, trace = FALSE,
                          seed = 123)
  )
  cluster1_targets <- rownames(mats$target)[1:4]
  cluster2_targets <- rownames(mats$target)[5:8]
  expect_false(any(net$Target %in% cluster1_targets))
  expect_true(any(net$Target %in% cluster2_targets))
})

test_that("weight_matrix_to_edges drops edges below threshold and orders target-major", {
  network <- matrix(c(0.1, 0.9, 0.4, 0.2), nrow = 2, byrow = TRUE,
                     dimnames = list(c("t1", "t2"), c("r1", "r2")))
  edges <- weight_matrix_to_edges(network, weightthreshold = 0.3)

  expect_equal(nrow(edges), 2)
  expect_equal(edges$Target, c("t1", "t2"))
  expect_equal(edges$Regulator, c("r2", "r1"))
  expect_equal(edges$Weight, c(0.9, 0.4))
})

test_that("weight_matrix_to_edges returns a valid 0-row edge table when nothing clears the threshold", {
  # e.g. weightthreshold = 0 with raw, unnormalized (possibly negative) importances --
  # the exact scenario permutation testing now always uses. A bare "regulates" scalar
  # doesn't recycle against zero-length columns in data.frame(), so this used to error
  # with "arguments imply differing number of rows: 0, 1".
  network <- matrix(c(-0.1, -0.2, -0.3, -0.4), nrow = 2, byrow = TRUE,
                     dimnames = list(c("t1", "t2"), c("r1", "r2")))
  edges <- weight_matrix_to_edges(network, weightthreshold = 0)

  expect_equal(nrow(edges), 0)
  expect_equal(names(edges), c("Regulator", "Interaction", "Target", "Weight"))
})

test_that("pick_hub_genes returns all tied max-out-degree regulators", {
  edges <- data.frame(Regulator = c("a", "a", "b", "b"), stringsAsFactors = FALSE)
  expect_setequal(pick_hub_genes(edges), c("a", "b"))
})

test_that("pick_hub_genes returns character(0) with no warning for a 0-row edge table", {
  edges <- data.frame(Regulator = character(0), stringsAsFactors = FALSE)
  expect_no_warning(result <- pick_hub_genes(edges))
  expect_equal(result, character(0))
})
