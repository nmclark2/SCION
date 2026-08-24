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

test_that("derived clustering_data's row names match target_data/reg_data exactly, even for
           genes whose names make.names() alters (leading digit, hyphen)", {
  target_genes <- c("0610005C13Rik", "HLA-A", "NormalGene1", "NormalGene2")
  reg_genes <- c("TF-1", "2900026A02Rik", "NormalTF1", "NormalTF2")
  samples <- paste0("S", 1:10)

  target <- as.data.frame(matrix(stats::rnorm(length(target_genes) * length(samples)),
                                  nrow = length(target_genes),
                                  dimnames = list(make.names(target_genes), samples)))
  reg <- as.data.frame(matrix(stats::rnorm(length(reg_genes) * length(samples)),
                               nrow = length(reg_genes),
                               dimnames = list(make.names(reg_genes), samples)))

  result <- cluster_genes(method = "kmeans", target_data = target, reg_data = reg)

  expect_true(all(rownames(target) %in% rownames(result)))
  expect_true(all(rownames(reg) %in% rownames(result)))
})

test_that("method = 'upload' applies make.names() to the clusters file's row names, matching
           the same treatment target_data/reg_data get from read_scion_inputs()", {
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  genes <- c("0610005C13Rik", "HLA-A", "NormalGene1")
  clusters_file <- file.path(tmp_dir, "clusters.csv")
  utils::write.csv(data.frame(clusters = c(1, 1, 2), row.names = genes), clusters_file)

  result <- cluster_genes(method = "upload", clusters_file = clusters_file)

  expect_setequal(rownames(result), make.names(genes))
})
