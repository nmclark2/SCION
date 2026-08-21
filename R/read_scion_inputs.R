#' Read SCION input matrices from CSV or GCT files
#'
#' Reads target and regulator expression matrices (and, optionally, a separate
#' clustering matrix), restricts each to a gene list if provided, and drops any
#' row containing a missing value -- SCION does not support missing values, so
#' unlike the private-version `na.max` proportional filter, this is always an
#' all-or-nothing drop.
#'
#' @param target_data_file,reg_data_file path to target/regulator expression
#'   matrices. CSV: first column = gene names, remaining columns = samples.
#'   GCT: a GCT(x) file, read via `cmapR::parse_gctx()` (requires the optional
#'   `cmapR` package -- install with `BiocManager::install("cmapR")`).
#' @param target_genes_file,reg_genes_file optional path to a file listing which
#'   genes to keep as targets/regulators (first column = gene names). `NULL`
#'   (default) keeps every gene present in the corresponding data file.
#' @param gene_list_header whether `target_genes_file`/`reg_genes_file` have a
#'   header row. Default `TRUE` (matches this package's own tutorial data).
#'   Set to `FALSE` for a plain one-gene-symbol-per-line file with no header
#'   (e.g. PANOPLY's `TF_file` convention) -- with the default `TRUE`, such a
#'   file would silently have its first gene misread as a column header and
#'   dropped.
#' @param clustering_data_file optional path to a CSV clustering matrix (rows =
#'   genes, columns = samples), used by [cluster_genes()]. Restricted to genes
#'   that appear in `target_genes_file` or `reg_genes_file` when those are given.
#' @param format `"csv"` (default) or `"gct"`. Regulator gene names may use
#'   either a dot (`SOX2.S35`) or underscore (`MEF2C_S453s`) PTM-site
#'   convention; both are left as-is here (see the `ptm_sep` argument of
#'   [infer_network()] for where the convention matters downstream).
#' @return a list with `target` and `reg` data frames (genes as rows, samples as
#'   columns, row names made syntactically valid via [make.names()]), and
#'   `cluster_data` (or `NULL` if `clustering_data_file` was not given).
#' @export
read_scion_inputs <- function(target_data_file, reg_data_file, target_genes_file = NULL,
                               reg_genes_file = NULL, gene_list_header = TRUE,
                               clustering_data_file = NULL, format = c("csv", "gct")) {
  format <- match.arg(format)

  if (format == "gct") {
    if (!requireNamespace("cmapR", quietly = TRUE)) {
      stop("Reading GCT files requires the 'cmapR' package. Install it with BiocManager::install('cmapR').")
    }
    target_data <- as.data.frame(cmapR::parse_gctx(target_data_file)@mat)
    reg_data <- as.data.frame(cmapR::parse_gctx(reg_data_file)@mat)
  } else {
    target_data <- utils::read.csv(target_data_file, row.names = 1)
    reg_data <- utils::read.csv(reg_data_file, row.names = 1)
  }

  target_genes <- if (!is.null(target_genes_file)) {
    utils::read.csv(target_genes_file, header = gene_list_header, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  reg_genes <- if (!is.null(reg_genes_file)) {
    utils::read.csv(reg_genes_file, header = gene_list_header, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  if (!is.null(target_genes)) {
    target_data <- target_data[row.names(target_data) %in% target_genes[, 1], ]
  }
  if (!is.null(reg_genes)) {
    reg_data <- reg_data[row.names(reg_data) %in% reg_genes[, 1], ]
  }

  rownames(target_data) <- make.names(rownames(target_data))
  rownames(reg_data) <- make.names(rownames(reg_data))

  if (sum(is.na(target_data), is.na(reg_data)) > 0) {
    message("Missing values detected in target and/or regulator matrix. SCION cannot use missing ",
            "values. Features with missing values have been removed. To retain these features, ",
            "please impute missing values.")
    target_data <- na.omit(target_data)
    reg_data <- na.omit(reg_data)
  }

  cluster_data <- NULL
  if (!is.null(clustering_data_file)) {
    cluster_data <- utils::read.csv(clustering_data_file, row.names = 1)
    keep <- rep(TRUE, nrow(cluster_data))
    if (!is.null(target_genes) || !is.null(reg_genes)) {
      keep <- row.names(cluster_data) %in% c(
        if (!is.null(target_genes)) target_genes[, 1] else character(0),
        if (!is.null(reg_genes)) reg_genes[, 1] else character(0)
      )
    }
    cluster_data <- cluster_data[keep, ]
    rownames(cluster_data) <- make.names(rownames(cluster_data))
    if (sum(is.na(cluster_data)) > 0) {
      message("Missing values detected in clustering matrix. SCION cannot use missing values. ",
              "Features with missing values have been removed. To retain these features, please ",
              "impute missing values.")
      cluster_data <- na.omit(cluster_data)
    }
  }

  list(target = target_data, reg = reg_data, cluster_data = cluster_data)
}
