#' Calculate cross-pipeline consensus statistics
#'
#' Adds `n_methods` (how many pipelines called each site), and median read
#' depth / methylation level across whichever pipelines called it, to a
#' wide table produced by [merge_filtered_pipelines()].
#' `median_Mread` is the rounded *product* of median depth and median
#' methylation level, not the median of the observed `_MR` columns.
#'
#' @param dt A wide data.table with one or more `<pipeline>_TR` and
#'   `<pipeline>_ML` columns (output of [merge_filtered_pipelines()]).
#'
#' @return `dt` with four new columns added: `n_methods`, `median_Tread`,
#'   `median_MethyL`, `median_Mread`.
#' @export
calculate_consensus_stats <- function(dt) {
  read_cols <- grep("_TR$", names(dt), value = TRUE)

  if (length(read_cols) == 0) {
    stop("No `_TR` columns found in `dt`. Expected output of merge_filtered_pipelines().")
  }
  tools <- sub("_TR$", "", read_cols)
  ml_cols <- paste0(tools, "_ML")
  if (!all(ml_cols %in% names(dt))) {
    stop("Missing `_ML` columns for: ",
         paste(tools[!ml_cols %in% names(dt)], collapse = ", "))
  }

  dt <- data.table::copy(dt)
  depths <- as.matrix(dt[, ..read_cols])
  depths[!is.finite(depths) | depths <= 0] <- NA_real_
  present <- !is.na(depths)
  dt[, n_methods := rowSums(present)]
  dt[, median_Tread := ceiling(matrixStats::rowMedians(depths, na.rm = TRUE))]

  methylation <- as.matrix(dt[, ..ml_cols])
  methylation[!present] <- NA_real_
  dt[, median_MethyL := matrixStats::rowMedians(methylation, na.rm = TRUE)]
  dt[, median_Mread := round(median_Tread * median_MethyL)]
  dt[]
}
