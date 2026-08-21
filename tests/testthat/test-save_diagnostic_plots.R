test_that("save_diagnostic_plots() writes weight/outdegree PDFs, and FDR PDFs only when present", {
  result <- list(network = data.frame(Regulator = c("r1", "r1", "r2"), Target = c("t1", "t2", "t3"),
                                        Weight = c(0.5, 0.9, 0.2), stringsAsFactors = FALSE))
  tmp_dir <- tempfile("scion-plots-")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  paths <- save_diagnostic_plots(result, tmp_dir, prefix = "myrun")

  expect_true(file.exists(file.path(tmp_dir, "myrun_weight_distribution.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "myrun_outdegree_distribution.pdf")))
  expect_false(file.exists(file.path(tmp_dir, "myrun_fdr_curve.pdf")))
  expect_equal(length(paths), 2)
})

test_that("save_diagnostic_plots() also writes FDR curve/weight-comparison PDFs when fdr_result is present", {
  real_network <- data.frame(Regulator = "r1", Target = paste0("t", 1:5),
                              Weight = c(10, 8, 6, 4, 2), stringsAsFactors = FALSE)
  permuted_networks <- replicate(3, data.frame(Weight = rep(1, 5)), simplify = FALSE)
  fdr_result <- compute_fdr_threshold(real_network, permuted_networks)
  result <- list(network = real_network, fdr_result = fdr_result)

  tmp_dir <- tempfile("scion-plots-")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  paths <- save_diagnostic_plots(result, tmp_dir)

  expect_true(file.exists(file.path(tmp_dir, "diagnostics_fdr_curve.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "diagnostics_weight_comparison.pdf")))
  expect_equal(length(paths), 4)
})

test_that("save_diagnostic_plots() creates the target directory if it doesn't exist", {
  result <- list(network = data.frame(Regulator = "r1", Target = "t1", Weight = 0.5,
                                        stringsAsFactors = FALSE))
  tmp_dir <- file.path(tempfile("scion-plots-"), "nested", "dir")
  on.exit(unlink(dirname(dirname(tmp_dir)), recursive = TRUE), add = TRUE)

  expect_false(dir.exists(tmp_dir))
  save_diagnostic_plots(result, tmp_dir)
  expect_true(dir.exists(tmp_dir))
})
