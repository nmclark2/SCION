test_that("bundled Arabidopsis example data has both protein and phospho regulator variants", {
  example_dir <- system.file("extdata", "arabidopsis", package = "SCION")
  skip_if(!nzchar(example_dir), "package not installed")

  expected_files <- c(
    "target_mat_RNA.csv", "target_list_RNA.csv",
    "reg_mat_protein.csv", "reg_list_protein.csv",
    "reg_mat_phospho.csv", "reg_list_phospho.csv"
  )
  expect_true(all(file.exists(file.path(example_dir, expected_files))))
})
