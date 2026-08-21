test_that("RS.Get.Weight.Matrix returns NULL (no warning) when every importance value is NA", {
  # %IncMSE can come back NA for every target/regulator pair in a small or
  # low-variance cluster -- mock the per-target fit directly so this is
  # deterministic, rather than depending on randomForest's actual stochastic
  # NA behavior (real trigger: e.g. a 2-target, 3-regulator DTW cluster from
  # the bundled Arabidopsis phospho example data).
  testthat::local_mocked_bindings(
    rsgwm2_randomforest = function(target.gene.name, ...) {
      stats::setNames(c(NA_real_, NA_real_), c("r1", "r2"))
    },
    .package = "SCION"
  )

  target <- data.frame(t1 = c(1, 2, 3), t2 = c(2, 3, 1))
  reg <- data.frame(r1 = c(1, 1, 2), r2 = c(2, 1, 1))

  expect_no_warning(
    net <- RS.Get.Weight.Matrix(target, reg, num.cores = 1, engine = "randomForest",
                                 normalize = TRUE, trace = FALSE)
  )
  expect_null(net)
})
