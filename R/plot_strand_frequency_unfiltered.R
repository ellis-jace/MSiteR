#' Frequency distribution: S1 and S2 (unfiltered)
#'
#' Plot read-depth frequencies before strand filtering. Optional dashed lines
#' show the actual strand-filtering cutoffs in [run_strand_filtering()].
#'
#' @param result Output from [run_strand_filtering()].
#' @param outdir Optional directory in which to save a PDF.
#' @param sample_name Sample ID used in the PDF filename.
#' @param reads_max Maximum read depth displayed.
#' @param show_thresholds Add lines for the actual filtering cutoffs.
#' @return A ggplot object (invisibly).
#' @export
plot_strand_frequency_unfiltered <- function(result,
                                             outdir = NULL,
                                             sample_name = "sample",
                                             reads_max = 75,
                                             show_thresholds = FALSE) {
  plot_strand_frequency_common(result, filtered = FALSE, outdir = outdir,
                               sample_name = sample_name, reads_max = reads_max,
                               show_thresholds = show_thresholds)
}
