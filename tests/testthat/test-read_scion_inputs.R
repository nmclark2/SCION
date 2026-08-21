test_that("gene_list_header = FALSE reads a headerless gene list without dropping the first gene", {
  mats <- make_test_matrices(n_targets = 4, n_regs = 4, n_samples = 5)
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  target_file <- file.path(tmp_dir, "target.csv")
  reg_file <- file.path(tmp_dir, "reg.csv")
  utils::write.csv(mats$target, target_file)
  utils::write.csv(mats$reg, reg_file)

  # a plain one-gene-per-line file with NO header, e.g. PANOPLY's TF_file convention
  reg_genes_file <- file.path(tmp_dir, "tf_list.txt")
  writeLines(rownames(mats$reg), reg_genes_file) # includes "reg1", the first gene

  result <- read_scion_inputs(target_file, reg_file, reg_genes_file = reg_genes_file,
                               gene_list_header = FALSE)

  expect_true("reg1" %in% rownames(result$reg))
  expect_equal(sort(rownames(result$reg)), sort(rownames(mats$reg)))
})

test_that("gene_list_header = TRUE (default) matches plain read.csv behavior", {
  mats <- make_test_matrices(n_targets = 4, n_regs = 4, n_samples = 5)
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  target_file <- file.path(tmp_dir, "target.csv")
  reg_file <- file.path(tmp_dir, "reg.csv")
  utils::write.csv(mats$target, target_file)
  utils::write.csv(mats$reg, reg_file)

  reg_genes_file <- file.path(tmp_dir, "reg_genes.csv")
  utils::write.csv(data.frame(gene = rownames(mats$reg)), reg_genes_file, row.names = FALSE)

  result <- read_scion_inputs(target_file, reg_file, reg_genes_file = reg_genes_file)

  expect_equal(sort(rownames(result$reg)), sort(rownames(mats$reg)))
})
