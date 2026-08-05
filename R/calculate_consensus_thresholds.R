#' Calculate tiered consensus read-depth thresholds
#'
#' For sites called by only one pipeline (`n_methods == 1`), computes each
#' pipeline's 95th-percentile read depth (floored at 3) as its threshold.
#' Sites called by 2, 3, or 4 pipelines get progressively looser thresholds
#' (half, a third, a quarter of the 1-method threshold, each clamped to a
#' fixed floor/ceiling), reflecting that more pipeline agreement requires
#' less depth evidence to trust.
#'
#' @param dt A data.table with `n_methods` and `<pipeline>_TR` columns
#'   (output of [calculate_consensus_stats()]).
#'
#' @return A data.table with one row per pipeline: `Tool`, `Threshold_1M`,
#'   `Threshold_2M`, `Threshold_3M`, `Threshold_4M`.
#' @export
calculate_consensus_thresholds <- function(dt) {
  read_cols <- grep("_TR$", names(dt), value = TRUE)
  tools <- sub("_TR$", "", read_cols)

  dt_n1 <- dt[n_methods == 1]
  cutoffs_n1 <- sapply(read_cols, function(col) {
    vals <- dt_n1[[col]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(3)
    max(quantile(vals, 0.95, na.rm = TRUE), 3)
  })
  names(cutoffs_n1) <- tools

  cutoffs_n2 <- pmax(5, pmin(10, ceiling(cutoffs_n1 / 2)))
  cutoffs_n3 <- pmax(4, pmin(8,  ceiling(cutoffs_n1 / 3)))
  cutoffs_n4 <- pmax(3, pmin(6,  ceiling(cutoffs_n1 / 4)))

  data.table::data.table(
    Tool = tools,
    Threshold_1M = as.numeric(cutoffs_n1),
    Threshold_2M = as.numeric(cutoffs_n2),
    Threshold_3M = as.numeric(cutoffs_n3),
    Threshold_4M = as.numeric(cutoffs_n4)
  )
}
