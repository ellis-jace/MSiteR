#' Frequency distribution: S1 and S2 (filtered)
#'
#' Shows the distribution of read depths AFTER filtering thresholds are applied.
#' Demonstrates the effect of the filtering step. Optionally display threshold
#' cutoff lines to visualize where filtering occurred.
#'
#' @param result Output from [prepare_filtered_cpg_table()]
#' @param outdir Optional directory to save PDF. If NULL, plot is not saved.
#' @param sample_name Sample ID for filename (default "sample")
#' @param reads_max Maximum read depth to display (default 75)
#' @param show_thresholds Logical: add threshold cutoff lines? (default FALSE)
#'
#' @return ggplot object
#' @export
#'
#' @details
#' When `show_thresholds = TRUE`, dashed vertical lines mark the minimum reads
#' per strand in the filtered dataset, showing exactly where filtering boundaries are.
plot_strand_frequency_filtered <- function(result,
                                           outdir = NULL,
                                           sample_name = "sample",
                                           reads_max = 75,
                                           show_thresholds = FALSE) {

  # Validate input
  merged <- result$merged
  if (!all(c("Pipeline", "Strand", "reads") %in% names(merged))) {
    stop("result$merged must contain columns: Pipeline, Strand, reads")
  }

  # Calculate frequency counts on filtered data
  freq_dt <- merged[reads >= 1 & reads <= reads_max,
                    .N,
                    by = .(Pipeline, Strand, reads)]
  data.table::setnames(freq_dt, "N", "Frequency")

  # Set factor levels
  freq_dt[, Pipeline := factor(Pipeline,
                               levels = unique(result$thresholds$Pipeline))]
  freq_dt[, Strand := factor(Strand, levels = c("S1", "S2"))]

  # If showing thresholds, extract minimum reads per strand from filtered data
  if (show_thresholds) {
    s1_marks <- merged[Strand == "S1", .(Min_Reads = min(reads)), by = Pipeline]
    s2_marks <- merged[Strand == "S2", .(Min_Reads = min(reads)), by = Pipeline]
  }

  # Build base plot
  pub_colors <- c("S1" = "#6d6c6cff", "S2" = "#E64B35")

  p_freq <- ggplot2::ggplot(freq_dt,
                            ggplot2::aes(x = reads, y = Frequency,
                                         color = Strand, fill = Strand)) +
    ggplot2::geom_area(position = "identity", alpha = 0.2, size = 0.8) +
    ggplot2::geom_line(size = 1)

  # Conditionally add threshold lines
  if (show_thresholds) {
    p_freq <- p_freq +
      ggplot2::geom_vline(data = s1_marks,
                          ggplot2::aes(xintercept = Min_Reads),
                          color = "#757575", linetype = "dashed", size = 0.7) +
      ggplot2::geom_vline(data = s2_marks,
                          ggplot2::aes(xintercept = Min_Reads),
                          color = "#E64B35", linetype = "dashed", size = 0.7) +
      ggplot2::geom_text(data = s1_marks,
                         ggplot2::aes(x = Min_Reads, y = -0.01, label = Min_Reads),
                         inherit.aes = FALSE, color = "#757575",
                         size = 4, fontface = "bold", vjust = 1.5) +
      ggplot2::geom_text(data = s2_marks,
                         ggplot2::aes(x = Min_Reads, y = -0.02, label = Min_Reads),
                         inherit.aes = FALSE, color = "#E64B35",
                         size = 4, fontface = "bold", vjust = 1.5)
  }

  # Add remaining layers
  p_freq <- p_freq +
    ggplot2::facet_wrap(~Pipeline, scales = "free_y", nrow = 1) +
    ggplot2::scale_x_continuous(breaks = c(1, 25, 50, 75, 100)) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::scale_color_manual(values = pub_colors) +
    ggplot2::scale_fill_manual(values = pub_colors) +
    ggplot2::labs(
      x = "Total Reads Count",
      y = "Frequency (Number of Sites)",
      title = "Single-strand and Double-strand methylation sites across pipelines (AFTER FILTERING)",
      subtitle = if (show_thresholds)
        "S1: Single-strand | S2: Double-strand | Dashed lines: filtering thresholds"
      else
        "S1: Single-strand | S2: Double-strand") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "gray95"),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "top",
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = if (show_thresholds)
        ggplot2::element_text(margin = ggplot2::margin(t = 10))
      else
        ggplot2::element_text())

  # Save if outdir provided
  if (!is.null(outdir)) {
    suffix <- if (show_thresholds) "_04Fig_Frequency_Filtered_with_Thresholds.pdf"
    else "_04Fig_Frequency_Filtered.pdf"
    filename <- file.path(outdir, paste0(sample_name, suffix))
    ggplot2::ggsave(filename, p_freq, width = 12, height = 5)
    message("Saved: ", filename)
  }

  return(invisible(p_freq))
}
