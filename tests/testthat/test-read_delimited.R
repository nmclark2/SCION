test_that("read_delimited() reads .csv, .tsv, .txt, and .ssv identically (modulo delimiter)", {
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  df <- data.frame(gene = c("g1", "g2", "g3"), value = c(1.5, 2.5, 3.5))

  csv_file <- file.path(tmp_dir, "data.csv")
  tsv_file <- file.path(tmp_dir, "data.tsv")
  txt_file <- file.path(tmp_dir, "data.txt")
  ssv_file <- file.path(tmp_dir, "data.ssv")
  utils::write.csv(df, csv_file, row.names = FALSE)
  utils::write.table(df, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(df, txt_file, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(df, ssv_file, sep = ";", row.names = FALSE, quote = FALSE)

  expected <- df
  for (f in c(csv_file, tsv_file, txt_file, ssv_file)) {
    result <- read_delimited(f)
    expect_equal(result$gene, expected$gene, info = f)
    expect_equal(result$value, expected$value, info = f)
  }
})

test_that("read_delimited(header = FALSE) does not drop the first row", {
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  genes_file <- file.path(tmp_dir, "genes.txt")
  writeLines(c("gene1", "gene2", "gene3"), genes_file)

  result <- read_delimited(genes_file, header = FALSE)
  expect_true("gene1" %in% result[[1]])
})

test_that("read_delimited() errors informatively on an unsupported extension", {
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  bad_file <- file.path(tmp_dir, "data.xlsx")
  writeLines("gene,value", bad_file)

  expect_error(read_delimited(bad_file), "Unsupported file extension")
})

test_that("read_delimited_matrix() uses the first column as row names, matching
           utils::read.csv(path, row.names = 1)", {
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  mat <- matrix(1:6, nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  csv_file <- file.path(tmp_dir, "mat.csv")
  tsv_file <- file.path(tmp_dir, "mat.tsv")
  utils::write.csv(mat, csv_file)
  utils::write.table(mat, tsv_file, sep = "\t", col.names = NA, quote = FALSE)

  expected <- utils::read.csv(csv_file, row.names = 1)

  result_csv <- read_delimited_matrix(csv_file)
  result_tsv <- read_delimited_matrix(tsv_file)

  expect_equal(result_csv[sort(rownames(result_csv)), ], expected[sort(rownames(expected)), ])
  expect_equal(result_tsv[sort(rownames(result_tsv)), ], expected[sort(rownames(expected)), ],
               ignore_attr = TRUE)
})

test_that("read_scion_inputs() accepts .tsv target/regulator matrices, not just .csv", {
  mats <- make_test_matrices(n_targets = 4, n_regs = 4, n_samples = 5)
  tmp_dir <- tempfile("scion-test-")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  target_file <- file.path(tmp_dir, "target.tsv")
  reg_file <- file.path(tmp_dir, "reg.tsv")
  utils::write.table(mats$target, target_file, sep = "\t", col.names = NA, quote = FALSE)
  utils::write.table(mats$reg, reg_file, sep = "\t", col.names = NA, quote = FALSE)

  result <- read_scion_inputs(target_file, reg_file)

  expect_setequal(rownames(result$target), rownames(mats$target))
  expect_setequal(rownames(result$reg), rownames(mats$reg))
})
