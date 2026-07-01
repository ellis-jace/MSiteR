#' Calculate per-pipeline strand-based read-depth thresholds
#'
#' For each pipeline, computes a single-strand (S1) read-depth cutoff based
#' on the 80th percentile of S1 reads (clamped between 5 and 12), and a
#' double-strand (S2) cutoff equal to half of the S1 cutoff, rounded up.
#'
#' @param dt A long-format data.table with columns `Pipeline`, `Strand`
#'   (values `"S1"`/`"S2"`), and `reads`.
#'
#' @return A data.table with one row per pipeline: `Pipeline`, `median_val`,
#'   `quantile80`, `cutoff1` (S1 threshold), `cutoff2` (S2 threshold).
#' @export
calculate_strand_thresholds <- function(dt) {
  s1_data <- dt[Strand == "S1"]
  s1_data[, .(
    median_val = as.numeric(median(reads)),
    quantile80 = as.numeric(quantile(reads, 0.8)),
    cutoff1 = min(12, max(5, quantile(reads, 0.8))),
    cutoff2 = ceiling(min(12, max(5, quantile(reads, 0.8))) / 2)
  ), by = Pipeline][order(Pipeline)]
}
