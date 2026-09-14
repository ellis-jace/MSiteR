#' Boxplot for Strand 1 vs Strand 2 reads
#'
#' Shows read depth distribution for single-strand (S1) vs double-strand (S2)
#' methylation sites across pipelines. Includes Wilcoxon test statistics.
#'
#' @param result Output from [prepare_filtered_cpg_table()]
#' @param outdir Optional directory to save PDF. If NULL, plot is not saved.
#' @param sample_name Sample ID for filename (default "sample")
#'
#' @return ggplot object (invisibly if saved to file)
#' @export
#'
#' @examples
#' \dontrun{
#' result <- prepare_filtered_cpg_table(pipelines)
#' plot_strand_boxplot(result, outdir = "figs/", sample_name = "Lab1")
#' }
plot_strand_boxplot <- function(result,
                                  outdir = NULL,
                                  sample_name = "sample") {

  # Validate input
  unfiltered <- result$unfiltered
  if (!all(c("Pipeline", "Strand", "reads") %in% names(unfiltered))) {
    stop("result$unfiltered must contain columns: Pipeline, Strand, reads")
  }

  # Prepare data
  plot_dt <- data.table::copy(unfiltered)
  plot_dt[, Pipeline := factor(Pipeline,
                               levels = unique(result$thresholds$Pipeline))]
  plot_dt[, Strand := factor(Strand, levels = c("S1", "S2"))]

  # Downsample for speed
  set.seed(123)
  plot_sub <- plot_dt[, .SD[sample(.N, min(.N, 1000000))],
                      by = .(Pipeline, Strand)]

  # Build plot
  p_box <- ggplot2::ggplot(plot_sub, ggplot2::aes(x = Pipeline, y = reads, fill = Strand)) +
    ggplot2::geom_boxplot(outlier.shape = NA,
                          notch = TRUE,
                          width = 0.6,
                          size = 0.7,
                          color = "black") +
    ggplot2::scale_fill_manual(
      values = c("S1" = "#BDBDBD", "S2" = "#E64B35"),
      labels = c("S1" = "Single-strand", "S2" = "Double-strand")) +
    ggpubr::stat_compare_means(ggplot2::aes(group = Strand),
                               label = "p.format",
                               method = "wilcox.test",
                               label.y = 75,
                               size = 4,
                               fontface = "italic") +
    ggpubr::stat_compare_means(ggplot2::aes(group = Strand),
                               label = "p.signif",
                               method = "wilcox.test",
                               label.y = 71,
                               size = 5) +
    ggplot2::coord_cartesian(ylim = c(0, 80)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 80, 20),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggpubr::theme_pubr(base_size = 14, legend = "top") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linetype = "dashed")) +
    ggplot2::labs(x = "Pipelines",
                  y = "Total Reads Count",
                  title = "Total Reads count distribution of S1 and S2",
                  fill = "Strand Status")

  # Save if outdir provided
  if (!is.null(outdir)) {
    filename <- file.path(outdir, paste0(sample_name, "_01Fig_Boxplot_S1_S2.pdf"))
    ggplot2::ggsave(filename, p_box, width = 8, height = 6)
    message("Saved: ", filename)
  }

  return(invisible(p_box))
}
