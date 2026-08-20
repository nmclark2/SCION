test_that("plot_network returns an igraph object with the right node/edge counts", {
  edges <- data.frame(Regulator = c("r1", "r2", "r1"), Target = c("t1", "t1", "t2"),
                       Weight = c(0.5, 0.9, 0.2), stringsAsFactors = FALSE)

  plot_file <- tempfile(fileext = ".png")
  grDevices::png(plot_file)
  g <- plot_network(edges)
  grDevices::dev.off()
  unlink(plot_file)

  expect_s3_class(g, "igraph")
  expect_equal(igraph::ecount(g), 3)
  expect_equal(igraph::vcount(g), 4) # r1, r2, t1, t2
})

test_that("plot_fdr_curve returns a ggplot object for each type", {
  real_network <- data.frame(Weight = c(10, 8, 6, 4, 2))
  permuted_networks <- replicate(3, data.frame(Weight = c(1, 1, 1, 1, 1)), simplify = FALSE)
  fdr_result <- compute_fdr_threshold(real_network, permuted_networks)

  expect_s3_class(plot_fdr_curve(fdr_result, type = "curve"), "ggplot")
  expect_s3_class(plot_fdr_curve(fdr_result, type = "fdr_hist"), "ggplot")
  expect_s3_class(plot_fdr_curve(fdr_result, type = "weight_comparison"), "ggplot")
})
