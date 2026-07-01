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
  cutoff1 <- thresholds[Pipeline == pipeline_name, cutoff1]
  cutoff2 <- thresholds[Pipeline == pipeline_name, cutoff2]
  dt[n == 2 & get(reads_col) >= cutoff2 | (n == 1 & get(reads_col) >= cutoff1)]
}
