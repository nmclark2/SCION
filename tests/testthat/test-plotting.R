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

test_that("plot_network(interactive = TRUE) groups nodes by Regulator/Target with a legend", {
  skip_if_not_installed("visNetwork")
  edges <- data.frame(Regulator = c("r1", "r2", "r1"), Target = c("t1", "t1", "t2"),
                       Weight = c(0.5, 0.9, 0.2), stringsAsFactors = FALSE)

  vis <- plot_network(edges, interactive = TRUE)

  nodes <- vis$x$nodes
  expect_equal(nodes$group[nodes$id == "r1"], "Regulator")
  expect_equal(nodes$group[nodes$id == "r2"], "Regulator")
  expect_equal(nodes$group[nodes$id == "t1"], "Target")
  expect_equal(nodes$group[nodes$id == "t2"], "Target")
  expect_true(all(!is.na(nodes$label)))
  # a plain static HTML/CSS legend, not visNetwork::visLegend() -- that
  # renders its own independently zoomable vis.js canvas, which reads as a
  # comically oversized, pointlessly interactive legend at normal widget sizes
  legend_html <- as.character(htmltools::doRenderTags(vis$prepend[[1]]))
  expect_match(legend_html, "Regulator")
  expect_match(legend_html, "Target")
  expect_match(legend_html, "weight", ignore.case = TRUE)
  expect_true(isTRUE(vis$x$options$interaction$navigationButtons))
})

test_that("plot_network() static and interactive both encode Regulator/Target by shape/color", {
  edges <- data.frame(Regulator = c("r1", "r2", "r1"), Target = c("t1", "t1", "t2"),
                       Weight = c(0.5, 0.9, 0.2), stringsAsFactors = FALSE)

  plot_file <- tempfile(fileext = ".png")
  grDevices::png(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  })
  g <- plot_network(edges, interactive = FALSE)

  is_regulator <- igraph::V(g)$name %in% c("r1", "r2")
  expect_true(all(igraph::V(g)$shape[is_regulator] == "square"))
  expect_true(all(igraph::V(g)$shape[!is_regulator] == "circle"))
})

test_that("plot_fdr_curve returns a ggplot object for each type", {
  real_network <- data.frame(Weight = c(10, 8, 6, 4, 2))
  permuted_networks <- replicate(3, data.frame(Weight = c(1, 1, 1, 1, 1)), simplify = FALSE)
  fdr_result <- compute_fdr_threshold(real_network, permuted_networks)

  expect_s3_class(plot_fdr_curve(fdr_result, type = "curve"), "ggplot")
  expect_s3_class(plot_fdr_curve(fdr_result, type = "fdr_hist"), "ggplot")
  expect_s3_class(plot_fdr_curve(fdr_result, type = "weight_comparison"), "ggplot")
})
