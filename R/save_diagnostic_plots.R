#' Save SCION's diagnostic plots as publication-quality PDFs
#'
#' Writes the edge weight distribution and regulator out-degree distribution
#' (and, if permutations were run, the FDR curve and real-vs-permuted weight
#' comparison) as separate single-plot PDFs via `ggplot2::ggsave()` -- vector
#' output, so resolution is never a concern. This is what [run_scion()] calls
#' automatically whenever `output_file` is set (i.e. CLI/scripted use, where
#' there's no app to view the interactive versions in); call it directly for
#' any other `run_scion()` result you want plots from.
#'
#' @param result a [run_scion()] result list.
#' @param dir directory to save into (created, recursively, if it doesn't exist).
#' @param prefix filename prefix for each saved plot.
#' @param width,height,units passed to `ggplot2::ggsave()`.
#' @return invisibly, a character vector of the file paths written.
#' @export
save_diagnostic_plots <- function(result, dir, prefix = "diagnostics",
                                   width = 7, height = 5, units = "in") {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  cutoff <- if (!is.null(result$fdr_result)) result$fdr_result$threshold else NULL

  paths <- character(0)
  save_one <- function(plot, name) {
    path <- file.path(dir, sprintf("%s_%s.pdf", prefix, name))
    ggplot2::ggsave(path, plot, device = "pdf", width = width, height = height, units = units)
    paths[[length(paths) + 1]] <<- path
  }

  save_one(plot_weight_distribution(result$network, cutoff = cutoff), "weight_distribution")
  save_one(plot_outdegree_distribution(result$network), "outdegree_distribution")
  if (!is.null(result$fdr_result)) {
    save_one(plot_fdr_curve(result$fdr_result, type = "curve"), "fdr_curve")
    save_one(plot_fdr_curve(result$fdr_result, type = "weight_comparison"), "weight_comparison")
  }

  invisible(paths)
}
