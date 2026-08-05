#' Calculate cross-pipeline consensus statistics
#'
#' Adds `n_methods` (how many pipelines called each site), and median read
#' depth / methylation level across whichever pipelines called it, to a
#' wide table produced by [merge_filtered_pipelines()].
#'
#' @param dt A wide data.table with one or more `<pipeline>_TR` and
#'   `<pipeline>_ML` columns (output of [merge_filtered_pipelines()]).
#'
#' @return `dt` with four new columns added: `n_methods`, `median_Tread`,
#'   `median_MethyL`, `median_Mread`.
#' @export
calculate_consensus_stats <- function(dt) {
  read_cols <- grep("_TR$", names(dt), value = TRUE)
  ml_cols   <- grep("_ML$", names(dt), value = TRUE)

  if (length(read_cols) == 0) {
    stop("No `_TR` columns found in `dt`. Expected output of merge_filtered_pipelines().")
  }

  dt <- data.table::copy(dt)
  dt[, n_methods := rowSums(!is.na(.SD) & .SD > 0), .SDcols = read_cols]
  dt[, median_Tread := ceiling(matrixStats::rowMedians(as.matrix(.SD), na.rm = TRUE)), .SDcols = read_cols]
  dt[, median_MethyL := matrixStats::rowMedians(as.matrix(.SD), na.rm = TRUE), .SDcols = ml_cols]
  dt[, median_Mread := round(median_Tread * median_MethyL)]
  dt[]
}
