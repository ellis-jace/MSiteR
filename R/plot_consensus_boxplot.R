#' Boxplot: Read depth distribution by consensus tier
#'
#' Compares median read depth distributions across consensus tiers (1M, 2M, 3M, 4M)
#' using boxplots. Shows that higher consensus tiers typically have higher read depths.
#'
#' @param result Output from [run_consensus_filtering()]
#' @param outdir Optional directory to save PDF. If NULL, plot is not saved.
#' @param sample_name Sample ID for filename (default "sample")
#' @param show_filtered Logical: show filtered data? If FALSE (default), shows
#'   unfiltered data. If TRUE, shows filtered data.
#'
#' @return ggplot object (invisibly if saved)
#' @export
plot_consensus_boxplot <- function(result,
                                   outdir = NULL,
                                   sample_name = "sample",
                                   show_filtered = FALSE) {

  # Validate input
  if (!all(c("merged", "unfiltered") %in% names(result))) {
    stop("result must be output from run_consensus_filtering()")
  }

  # Validate that consensus metrics exist
  dt <- if (show_filtered) result$merged else result$unfiltered
  if (!all(c("n_methods", "median_Tread") %in% names(dt))) {
    stop("result must have n_methods and median_Tread columns")
  }

  # Prepare data: extract median_Tread and n_methods
  plot_dt <- dt[, .(median_Tread, n_methods)]
  plot_dt[, dataset := factor(n_methods, levels = c(1, 2, 3, 4), labels = c("1M", "2M", "3M", "4M"))]

  # Publication color scheme
  pub_colors <- c("1M" = "#D32F2F", "2M" = "#1976D2", "3M" = "#388E3C", "4M" = "#7B1FA2")

  title_suffix <- if (show_filtered) " (After Consensus Filtering)" else " (Before Consensus Filtering)"

  # Build boxplot
  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = dataset, y = median_Tread, fill = dataset)) +
    ggplot2::geom_boxplot(outlier.shape = NA,
                          width = 0.5,
                          color = "black",
                          size = 0.6,
                          notch = TRUE) +
    ggplot2::scale_fill_manual(values = pub_colors) +
    ggplot2::coord_cartesian(ylim = c(0, 80)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 10)) +
    ggplot2::labs(
      x = "Consensus Tier (Pipelines Overlap)",
      y = "Median Read Depth",
      title = paste0("Methylation Distribution by Consensus Tier", title_suffix)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"))

  # Save if outdir provided
  if (!is.null(outdir)) {
    suffix <- if (show_filtered) "_Fig_Boxplot_Consensus_filtered.pdf"
    else "_Fig_Boxplot_Consensus_unfiltered.pdf"
    filename <- file.path(outdir, paste0(sample_name, suffix))
    ggplot2::ggsave(filename, p, width = 6, height = 6, device = "pdf")
    message("Saved: ", filename)
    return(invisible(p))
  }

  return(invisible(p))
}
