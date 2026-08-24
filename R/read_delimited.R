#' Read a delimited text file, dispatching on its extension
#'
#' Mirrors protigy-v2's file-upload convention: any of CSV, TSV, plain-text
#' (tab-delimited), or SSV (semicolon-delimited) is accepted, detected from
#' the file's own extension, rather than assuming every input is a CSV.
#'
#' @param path path to the file.
#' @param header whether the file has a header row.
#' @return a data frame (never a tibble -- callers throughout this package
#'   use base-R `[`/`row.names<-` semantics).
#' @keywords internal
read_delimited <- function(path, header = TRUE) {
  ext <- tolower(tools::file_ext(path))
  reader <- switch(ext,
    csv = readr::read_csv,
    tsv = ,
    txt = readr::read_tsv,
    ssv = function(p, ...) readr::read_delim(p, delim = ";", ...),
    stop("Unsupported file extension '.", ext, "' for '", path,
         "' -- expected .csv, .tsv, .txt, or .ssv.")
  )
  # suppressMessages(): readr reports it when it has to repair a blank/duplicate
  # column name (e.g. the empty header above a row-names column) -- read.csv()
  # handles the same input silently, and this should too.
  as.data.frame(suppressMessages(reader(path, col_names = header, show_col_types = FALSE)),
                stringsAsFactors = FALSE)
}

#' Read a delimited text file as a row-named matrix (first column = row names)
#'
#' @inheritParams read_delimited
#' @return a data frame with the first column's values as row names, that
#'   column itself dropped -- the same shape [utils::read.csv()] with
#'   `row.names = 1` produces.
#' @keywords internal
read_delimited_matrix <- function(path) {
  df <- read_delimited(path, header = TRUE)
  row.names(df) <- df[[1]]
  df[[1]] <- NULL
  df
}
