write_test_csvs <- function(mats, dir) {
  target_file <- file.path(dir, "target.csv")
  reg_file <- file.path(dir, "reg.csv")
  utils::write.csv(mats$target, target_file)
  utils::write.csv(mats$reg, reg_file)
  list(target_file = target_file, reg_file = reg_file)
}

test_that("run_scion(permute = TRUE) runs the real network, permutations, and FDR threshold in one call", {
  skip_on_cran()
  mats <- make_test_matrices()
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  files <- write_test_csvs(mats, tmp_dir)

  result <- run_scion(files$target_file, files$reg_file, weightthreshold = 0, normalize = FALSE,
                       num.cores = 1, engine = "randomForest", seed = 1, permute = TRUE,
                       n_permutations = 3, base_seed = 0, target_fdr = 0.05,
                       nb.trees = 50, trace = FALSE)

  expect_true(is.data.frame(result$network))
  expect_length(result$permuted_networks, 3)
  expect_true(all(c("curve", "threshold", "thresholded_network") %in% names(result$fdr_result)))
  expect_identical(result$network_thresholded, result$fdr_result$thresholded_network)
  if (!is.na(result$fdr_result$threshold)) {
    expect_true(all(result$network_thresholded$Weight >= result$fdr_result$threshold))
  }
})

test_that("run_scion(permute = TRUE)'s permutations match a standalone permute_network() call", {
  skip_on_cran()
  mats <- make_test_matrices()
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  files <- write_test_csvs(mats, tmp_dir)

  result <- run_scion(files$target_file, files$reg_file, weightthreshold = 0, normalize = FALSE,
                       num.cores = 1, engine = "randomForest", seed = 1, permute = TRUE,
                       n_permutations = 3, base_seed = 0, nb.trees = 50, trace = FALSE)

  standalone <- permute_network(result$target, result$reg, n_permutations = 3, base_seed = 0,
                                 num.cores = 1, weightthreshold = 0, normalize = FALSE,
                                 engine = "randomForest", nb.trees = 50, trace = FALSE)

  expect_identical(result$permuted_networks, standalone)
})

test_that("run_scion() without permute = TRUE does not compute FDR/permutation results", {
  skip_on_cran()
  mats <- make_test_matrices()
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  files <- write_test_csvs(mats, tmp_dir)

  result <- run_scion(files$target_file, files$reg_file, weightthreshold = 0, normalize = FALSE,
                       num.cores = 1, engine = "randomForest", nb.trees = 50, trace = FALSE)

  expect_false("fdr_result" %in% names(result))
  expect_false("permuted_networks" %in% names(result))
  expect_false("network_thresholded" %in% names(result))
})

test_that("output_file writes the thresholded network when permute = TRUE", {
  skip_on_cran()
  mats <- make_test_matrices()
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  files <- write_test_csvs(mats, tmp_dir)
  out_file <- file.path(tmp_dir, "out.txt")

  result <- run_scion(files$target_file, files$reg_file, weightthreshold = 0, normalize = FALSE,
                       num.cores = 1, engine = "randomForest", permute = TRUE, n_permutations = 3,
                       nb.trees = 50, trace = FALSE, output_file = out_file)

  written <- utils::read.delim(out_file)
  expect_equal(nrow(written), nrow(result$network_thresholded))
})
