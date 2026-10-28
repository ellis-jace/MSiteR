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
  required <- c("Pipeline", "Strand", "reads")
  if (!all(required %in% names(dt))) {
    stop("`dt` must contain: ", paste(required, collapse = ", "))
  }
  if (nrow(dt) == 0L) stop("Cannot calculate thresholds from empty data.")

  # Include pipelines with no S1 calls. In that case use the documented
  # minimum S1 cutoff of 5 (and therefore S2 cutoff of 3).
  dt[, {
    values <- reads[Strand == "S1" & is.finite(reads)]
    median_val <- if (length(values)) median(values) else NA_real_
    quantile80 <- if (length(values)) {
      as.numeric(stats::quantile(values, 0.8, names = FALSE))
    } else NA_real_
    cutoff1 <- if (length(values)) min(12, max(5, quantile80)) else 5
    .(median_val = median_val, quantile80 = quantile80,
      cutoff1 = cutoff1, cutoff2 = ceiling(cutoff1 / 2))
  }, by = Pipeline][order(Pipeline)]
}
