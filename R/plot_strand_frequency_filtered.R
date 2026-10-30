#' Frequency distribution: S1 and S2 (filtered)
#'
#' Plot read-depth frequencies after strand filtering. The function reconstructs
#' a plotting table from the wide `merged` output of [run_strand_filtering()].
#'
#' @param result Output from [run_strand_filtering()].
#' @param outdir Optional directory in which to save a PDF.
#' @param sample_name Sample ID used in the PDF filename.
#' @param reads_max Maximum read depth displayed.
#' @param show_thresholds Add lines for the actual filtering cutoffs.
#' @return A ggplot object (invisibly).
#' @export
plot_strand_frequency_filtered <- function(result,
                                           outdir = NULL,
                                           sample_name = "sample",
                                           reads_max = 75,
                                           show_thresholds = FALSE) {
  plot_strand_frequency_common(result, filtered = TRUE, outdir = outdir,
                               sample_name = sample_name, reads_max = reads_max,
                               show_thresholds = show_thresholds)
}
