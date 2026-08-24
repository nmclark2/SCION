#' Save one permutation's network to disk
#'
#' Companion to [load_permutations()] for an HPC job-array workflow: each array
#' task runs one permutation (via `permute_network(..., indices = task_id)`)
#' and saves its result with this function; a final job reads them all back
#' with [load_permutations()] and passes them to [compute_fdr_threshold()].
#'
#' @param network a single permuted network (edge table), e.g.
#'   `permute_network(..., indices = task_id)[[1]]`.
#' @param index the permutation index this network corresponds to (e.g. the
#'   HPC job array task ID).
#' @param dir directory to save into (created if it doesn't exist).
#' @param prefix filename prefix; the file is saved as `<prefix>_<index>.rds`.
#' @return the path written to, invisibly.
#' @export
save_permutation <- function(network, index, dir = ".", prefix = "permutation") {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- file.path(dir, sprintf("%s_%d.rds", prefix, as.integer(index)))
  saveRDS(network, path)
  invisible(path)
}

#' Read back permutation networks saved by `save_permutation()`
#'
#' @param dir directory to read from.
#' @param prefix filename prefix used when saving (must match).
#' @return a list of edge tables, named by permutation index (as a string) and
#'   ordered by increasing index -- directly usable as the `permuted_networks`
#'   argument to [compute_fdr_threshold()].
#' @export
load_permutations <- function(dir = ".", prefix = "permutation") {
  pattern <- sprintf("^%s_([0-9]+)\\.rds$", prefix)
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("No files matching '%s_<index>.rds' found in '%s'.", prefix, dir))
  }
  indices <- as.integer(sub(pattern, "\\1", basename(files)))
  ord <- order(indices)
  results <- lapply(files[ord], readRDS)
  stats::setNames(results, as.character(indices[ord]))
}
