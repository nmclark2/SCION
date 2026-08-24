test_that("detect_scion_format() dispatches on file extension", {
  expect_equal(detect_scion_format("target.csv"), "delimited")
  expect_equal(detect_scion_format("target.CSV"), "delimited")
  expect_equal(detect_scion_format("target.gct"), "gct")
  expect_equal(detect_scion_format("target.GCTX"), "gct")
  expect_equal(detect_scion_format("target.tsv"), "delimited") # not gct -> falls back to delimited-text reading
  expect_equal(detect_scion_format("target"), "delimited") # no extension at all
})

test_that("read_scion_inputs(format = 'auto') routes a .gct file through the GCT reader", {
  skip_if_not_installed("cmapR")

  mats <- make_test_matrices(n_targets = 4, n_regs = 4, n_samples = 5)
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  target_gct <- cmapR::GCT(mat = as.matrix(mats$target))
  reg_gct <- cmapR::GCT(mat = as.matrix(mats$reg))
  target_file <- file.path(tmp_dir, "target")
  reg_file <- file.path(tmp_dir, "reg")
  cmapR::write_gct(target_gct, target_file)
  cmapR::write_gct(reg_gct, reg_file)

  result <- read_scion_inputs(paste0(target_file, "_n5x4.gct"), paste0(reg_file, "_n5x4.gct"),
                               format = "auto")

  expect_equal(sort(rownames(result$target)), sort(make.names(rownames(mats$target))))
  expect_equal(sort(rownames(result$reg)), sort(make.names(rownames(mats$reg))))
})

test_that("read_scion_inputs(format = 'auto') defaults to CSV for a .csv file", {
  mats <- make_test_matrices(n_targets = 4, n_regs = 4, n_samples = 5)
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  target_file <- file.path(tmp_dir, "target.csv")
  reg_file <- file.path(tmp_dir, "reg.csv")
  utils::write.csv(mats$target, target_file)
  utils::write.csv(mats$reg, reg_file)

  result <- read_scion_inputs(target_file, reg_file, format = "auto")

  expect_equal(sort(rownames(result$target)), sort(rownames(mats$target)))
})
