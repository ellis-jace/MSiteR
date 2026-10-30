#' Filter a pipeline's CpG calls using strand-based thresholds
#'
#' Keeps sites where double-strand (S2) support meets `cutoff2`, or
#' single-strand (S1) support meets `cutoff1`, using the thresholds
#' calculated by [calculate_strand_thresholds()].
#'
#' @param dt A data.table for one pipeline with `n` (strand count: 1 or 2)
#'   and a read-depth column (`TRead` by default).
#' @param pipeline_name The pipeline's name, used to look up its row in
#'   `thresholds`.
#' @param thresholds The threshold table returned by
#'   [calculate_strand_thresholds()].
#' @param reads_col Name of the read-depth column in `dt`. Default `"TRead"`.
#'
#' @return The filtered data.table.
#' @export
filter_by_strand <- function(dt, pipeline_name, thresholds, reads_col = "TRead") {
  row <- thresholds[Pipeline == pipeline_name]
  if (nrow(dt) == 0L) {
    return(dt)
  }
  if (nrow(row) != 1L || anyNA(row[, .(cutoff1, cutoff2)])) {
    stop("Expected exactly one valid strand-threshold row for '", pipeline_name, "'.")
  }
  if (!all(c("n", reads_col) %in% names(dt))) {
    stop("`dt` must contain `n` and `", reads_col, "`.")
  }
  dt[(n == 2L & get(reads_col) >= row$cutoff2) |
       (n == 1L & get(reads_col) >= row$cutoff1)]
}
