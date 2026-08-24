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

test_that("gene-list filtering matches on raw names before make.names() is applied, so genes
           whose names make.names() would alter are not silently dropped", {
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # names make.names() actually changes: a leading digit gets an "X" prefix,
  # a hyphen becomes a dot
  target_genes <- c("0610005C13Rik", "HLA-A", "NormalGene1")
  reg_genes <- c("TF-1", "2900026A02Rik", "NormalTF1")
  samples <- paste0("S", 1:5)

  target_mat <- matrix(stats::rnorm(length(target_genes) * length(samples)),
                        nrow = length(target_genes), dimnames = list(target_genes, samples))
  reg_mat <- matrix(stats::rnorm(length(reg_genes) * length(samples)),
                     nrow = length(reg_genes), dimnames = list(reg_genes, samples))

  target_file <- file.path(tmp_dir, "target.csv")
  reg_file <- file.path(tmp_dir, "reg.csv")
  utils::write.csv(target_mat, target_file)
  utils::write.csv(reg_mat, reg_file)

  # gene list files use the RAW gene names, as a real user would provide --
  # not the make.names()-escaped versions
  target_genes_file <- file.path(tmp_dir, "target_genes.csv")
  reg_genes_file <- file.path(tmp_dir, "reg_genes.csv")
  utils::write.csv(data.frame(gene = target_genes), target_genes_file, row.names = FALSE)
  utils::write.csv(data.frame(gene = reg_genes), reg_genes_file, row.names = FALSE)

  result <- read_scion_inputs(target_file, reg_file, target_genes_file = target_genes_file,
                               reg_genes_file = reg_genes_file)

  expect_equal(nrow(result$target), length(target_genes))
  expect_equal(nrow(result$reg), length(reg_genes))
  expect_setequal(rownames(result$target), make.names(target_genes))
  expect_setequal(rownames(result$reg), make.names(reg_genes))
})

test_that("clustering_data_file filtering is also make.names()-safe for tricky gene names", {
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  genes <- c("0610005C13Rik", "HLA-A", "NormalGene1", "NormalGene2")
  samples <- paste0("S", 1:5)
  mat <- matrix(stats::rnorm(length(genes) * length(samples)),
                nrow = length(genes), dimnames = list(genes, samples))

  target_file <- file.path(tmp_dir, "target.csv")
  reg_file <- file.path(tmp_dir, "reg.csv")
  clustering_file <- file.path(tmp_dir, "clustering.csv")
  utils::write.csv(mat, target_file)
  utils::write.csv(mat, reg_file)
  utils::write.csv(mat, clustering_file)

  genes_file <- file.path(tmp_dir, "genes.csv")
  utils::write.csv(data.frame(gene = genes), genes_file, row.names = FALSE)

  result <- read_scion_inputs(target_file, reg_file, target_genes_file = genes_file,
                               clustering_data_file = clustering_file)

  expect_equal(nrow(result$cluster_data), length(genes))
  expect_setequal(rownames(result$cluster_data), make.names(genes))
})
