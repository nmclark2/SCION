#' Small synthetic target/regulator matrices for fast tests.
#' Genes as rows, samples as columns (matches read_scion_inputs() output).
make_test_matrices <- function(n_targets = 8, n_regs = 6, n_samples = 10, seed = 1) {
  set.seed(seed)
  target <- matrix(stats::rnorm(n_targets * n_samples), nrow = n_targets, ncol = n_samples)
  reg <- matrix(stats::rnorm(n_regs * n_samples), nrow = n_regs, ncol = n_samples)
  rownames(target) <- paste0("target", seq_len(n_targets))
  rownames(reg) <- paste0("reg", seq_len(n_regs))
  colnames(target) <- colnames(reg) <- paste0("sample", seq_len(n_samples))
  list(target = as.data.frame(target), reg = as.data.frame(reg))
}
