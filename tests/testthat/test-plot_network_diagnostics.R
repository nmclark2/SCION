test_that("compute_outdegree() counts edges per regulator, sorted descending", {
  network <- data.frame(Regulator = c("r1", "r1", "r2", "r1"), Target = paste0("t", 1:4),
                         stringsAsFactors = FALSE)

  result <- compute_outdegree(network)

  expect_equal(result$Regulator[1], "r1")
  expect_equal(result$outdegree[result$Regulator == "r1"], 3)
  expect_equal(result$outdegree[result$Regulator == "r2"], 1)
  expect_true(all(diff(result$outdegree) <= 0)) # sorted descending
})

test_that("plot_weight_distribution() and plot_outdegree_distribution() return ggplot objects", {
  network <- data.frame(Regulator = c("r1", "r1", "r2"), Target = c("t1", "t2", "t3"),
                         Weight = c(0.5, 0.9, 0.2), stringsAsFactors = FALSE)

  expect_s3_class(plot_weight_distribution(network), "ggplot")
  expect_s3_class(plot_outdegree_distribution(network), "ggplot")
})

test_that("plot_weight_distribution() adds a cutoff line only when cutoff is given", {
  network <- data.frame(Regulator = "r1", Target = c("t1", "t2"), Weight = c(0.5, 0.9),
                         stringsAsFactors = FALSE)

  # the cutoff is drawn as a dense geom_line() (not geom_vline()) so it stays
  # hoverable along its whole length in the interactive app -- see hover_line()
  is_vline_layer <- function(p) {
    any(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))
  }

  expect_false(is_vline_layer(plot_weight_distribution(network)))
  expect_false(is_vline_layer(plot_weight_distribution(network, cutoff = NA)))
  expect_true(is_vline_layer(plot_weight_distribution(network, cutoff = 0.6)))
})
