#' Merge strand-filtered CpG tables from multiple pipelines
#'
#' @param filtered A named list of strand-filtered data.tables (output of
#'   [filter_by_strand()]), one per pipeline. Each must have `chr`, `pos`,
#'   `TRead`, `MRead`, `ML`, `n`.
#'
#' @return A wide data.table with `chr`, `pos`, and each pipeline's own
#'   `<name>_TR`, `<name>_MR`, `<name>_ML`, `<name>_StrN` columns.
#' @export
merge_filtered_pipelines <- function(filtered) {
  renamed <- Map(function(dt, name) {
    dt <- data.table::copy(dt)
    data.table::setnames(dt, c("TRead", "MRead", "ML", "n"),
                         paste0(name, c("_TR", "_MR", "_ML", "_StrN")))
    dt
  }, filtered, names(filtered))

  Reduce(function(x, y) merge(x, y, by = c("chr", "pos"), all = TRUE), renamed)
}
