#' Filter CpG sites using tiered consensus thresholds
#'
#' Keeps sites where at least one pipeline's read depth meets the threshold
#' for that site's `n_methods` tier (see [calculate_consensus_thresholds()]).
#'
#' @param dt A data.table with `n_methods` and `<pipeline>_TR` columns.
#' @param thresholds The threshold table from
#'   [calculate_consensus_thresholds()].
#'
#' @return The filtered data.table.
#' @export
filter_by_consensus <- function(dt, thresholds) {
  read_cols <- grep("_TR$", names(dt), value = TRUE)
  tools <- sub("_TR$", "", read_cols)
  if (!length(read_cols) || !"n_methods" %in% names(dt)) {
    stop("`dt` needs `n_methods` and pipeline `_TR` columns.")
  }
  if (length(tools) > 4L) {
    stop("Consensus thresholding currently supports at most four pipelines.")
  }
  needed <- paste0("Threshold_", 1:4, "M")
  if (!all(c("Tool", needed) %in% names(thresholds)) ||
      anyDuplicated(thresholds$Tool) || !all(tools %in% thresholds$Tool)) {
    stop("`thresholds` must contain one complete row for every pipeline.")
  }
  if (anyNA(thresholds[, ..needed])) {
    stop("Consensus thresholds cannot contain missing cutoffs.")
  }
  if (anyNA(dt$n_methods) || any(dt$n_methods < 0 | dt$n_methods > length(tools))) {
    stop("Invalid `n_methods` values.")
  }

  # Apply all thresholds to one logical mask; avoid four successive table copies.
  keep <- rep(FALSE, nrow(dt))
  cutoffs <- thresholds[match(tools, thresholds$Tool)]
  for (k in seq_along(tools)) {
    rows <- which(dt$n_methods == k)
    if (!length(rows)) next
    cutoff <- cutoffs[[paste0("Threshold_", k, "M")]]
    for (j in seq_along(read_cols)) {
      values <- dt[[read_cols[j]]][rows]
      keep[rows] <- keep[rows] | (!is.na(values) & values >= cutoff[j])
    }
  }
  dt[keep]
}
