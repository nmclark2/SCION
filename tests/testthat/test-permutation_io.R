test_that("permute_network(indices=) runs a single permutation matching the equivalent full run", {
  skip_on_cran()
  mats <- make_test_matrices()

  full_run <- permute_network(mats$target, mats$reg, n_permutations = 3, num.cores = 1,
                               weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                               engine = "randomForest", nb.trees = 50, trace = FALSE)
  one_task <- permute_network(mats$target, mats$reg, indices = 2, num.cores = 1,
                               weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                               engine = "randomForest", nb.trees = 50, trace = FALSE)

  expect_equal(names(one_task), "2")
  expect_equal(one_task[["2"]], full_run[["2"]])
})

test_that("save_permutation()/load_permutations() round-trip an HPC-array-style workflow", {
  skip_on_cran()
  mats <- make_test_matrices()
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  for (task_id in c(3, 1, 2)) { # out of order, as array tasks might finish
    perm <- permute_network(mats$target, mats$reg, indices = task_id, num.cores = 1,
                             weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                             engine = "randomForest", nb.trees = 50, trace = FALSE)[[1]]
    save_permutation(perm, task_id, dir = tmp_dir)
  }

  loaded <- load_permutations(tmp_dir)
  expect_equal(names(loaded), c("1", "2", "3")) # sorted by index, not save order

  direct <- permute_network(mats$target, mats$reg, indices = 1:3, num.cores = 1,
                             weightthreshold = 0, normalize = FALSE, connect_hubs = FALSE,
                             engine = "randomForest", nb.trees = 50, trace = FALSE)
  expect_equal(loaded, direct)
})
