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

  dt <- data.table::copy(dt)
  dt[, n_methods := as.integer(n_methods)]

  kept <- list()
  remaining <- dt

  for (k in 1:4) {
    col <- paste0("Threshold_", k, "M")
    group <- remaining[n_methods == k]

    if (nrow(group) > 0) {
      keep <- Reduce(`|`, lapply(tools, function(t) {
        tr_col <- paste0(t, "_TR")
        cutoff <- thresholds[Tool == t, get(col)]
        !is.na(group[[tr_col]]) & group[[tr_col]] >= cutoff
      }))
      kept[[k]] <- group[keep]
    }

    remaining <- remaining[n_methods != k]
  }

  data.table::rbindlist(kept)
}
