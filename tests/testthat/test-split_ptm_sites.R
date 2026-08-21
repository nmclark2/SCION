test_that("split_ptm_sites splits gene+site regulators and adds a Site column", {
  network <- data.frame(
    Regulator = c("SOX2.S35", "SOX2.T90", "OCT4"),
    Interaction = "regulates",
    Target = c("t1", "t2", "t3"),
    Weight = c(0.5, 0.3, 0.8),
    stringsAsFactors = FALSE
  )

  result <- split_ptm_sites(network, ptm_sep = ".")

  expect_equal(result$Regulator, c("SOX2", "SOX2", "OCT4"))
  expect_equal(result$Site, c("S35", "T90", NA_character_))
  expect_equal(names(result), c("Regulator", "Site", "Interaction", "Target", "Weight"))
  # SOX2.S35 and SOX2.T90 are now recognized as the same regulator node
  expect_equal(sum(result$Regulator == "SOX2"), 2)
})

test_that("split_ptm_sites leaves non-PTM networks completely unchanged (no Site column added)", {
  network <- data.frame(Regulator = c("AT1G01260", "AT2G36080"), Interaction = "regulates",
                         Target = c("t1", "t2"), Weight = c(0.5, 0.3), stringsAsFactors = FALSE)

  result <- split_ptm_sites(network, ptm_sep = ".")

  expect_identical(result, network)
})

test_that("split_ptm_sites handles underscore separator and NULL/empty input", {
  network <- data.frame(Regulator = "MEF2C_S453s", Interaction = "regulates", Target = "t1",
                         Weight = 0.5, stringsAsFactors = FALSE)
  result <- split_ptm_sites(network, ptm_sep = "_")
  expect_equal(result$Regulator, "MEF2C")
  expect_equal(result$Site, "S453s")

  expect_null(split_ptm_sites(NULL, ptm_sep = "."))
  empty <- network[0, ]
  expect_identical(split_ptm_sites(empty, ptm_sep = "."), empty)
})

test_that("infer_network splits PTM regulators in its returned edge table", {
  set.seed(1)
  n_samples <- 8
  target <- as.data.frame(matrix(stats::rnorm(n_samples), nrow = 1, dimnames = list("SOX2", NULL)))
  reg <- as.data.frame(matrix(stats::rnorm(2 * n_samples), nrow = 2,
                               dimnames = list(c("SOX2.S35", "OCT4"), NULL)))
  colnames(target) <- colnames(reg) <- paste0("sample", seq_len(n_samples))

  net <- infer_network(target, reg, weightthreshold = -Inf, normalize = FALSE, num.cores = 1,
                        engine = "randomForest", ptm_sep = ".", nb.trees = 10, trace = FALSE)

  expect_true("Site" %in% names(net))
  expect_true(all(net$Regulator %in% c("SOX2", "OCT4")))
})
