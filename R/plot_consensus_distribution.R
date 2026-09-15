#' Distribution of read depths by consensus tier
#'
#' Shows the distribution of median read depths across sites grouped by
#' consensus tier (1M, 2M, 3M, 4M = how many pipelines detected each site).
#' Displays as area + line plot, one series per tier.
#'
#' @param result Output from [run_consensus_filtering()]
#' @param outdir Optional directory to save PDF. If NULL, plot is not saved.
#' @param sample_name Sample ID for filename (default "sample")
#' @param show_filtered Logical: show filtered data? If FALSE (default), shows
#'   unfiltered data. If TRUE, shows filtered data.
#' @param reads_max Maximum read depth to display (default 100)
#'
#' @return ggplot object (invisibly if saved)
#' @export
plot_consensus_distribution <- function(result,
                                        outdir = NULL,
                                        sample_name = "sample",
                                        show_filtered = FALSE,
                                        reads_max = 100) {

  # Validate input
  if (!all(c("merged", "unfiltered") %in% names(result))) {
    stop("result must be output from run_consensus_filtering()")
  }

  # Validate that consensus metrics exist
  dt <- if (show_filtered) result$merged else result$unfiltered
  if (!all(c("n_methods", "median_Tread") %in% names(dt))) {
    stop("result must have n_methods and median_Tread columns")
  }

  # Split by consensus tier and extract read depth
  dt1 <- dt[n_methods == 1, .(median_Tread)]
  dt2 <- dt[n_methods == 2, .(median_Tread)]
  dt3 <- dt[n_methods == 3, .(median_Tread)]
  dt4 <- dt[n_methods == 4, .(median_Tread)]

  # Add tier labels
  dt1[, dataset := "1M"]
  dt2[, dataset := "2M"]
  dt3[, dataset := "3M"]
  dt4[, dataset := "4M"]

  # Combine all tiers
  combined_data <- rbind(dt1, dt2, dt3, dt4)

  # Aggregate: Count how many sites at each read depth within each tier
  final_summary <- combined_data[median_Tread <= reads_max,
                                 .N,
                                 by = .(dataset, median_Tread)][order(dataset, median_Tread)]
  data.table::setnames(final_summary, "N", "frequency")

  # Publication color scheme
  pub_colors <- c("1M" = "#D32F2F", "2M" = "#1976D2", "3M" = "#388E3C", "4M" = "#7B1FA2")

  title_suffix <- if (show_filtered) " (After Consensus Filtering)" else " (Before Consensus Filtering)"

  # Plot distribution as area + line
  p <- ggplot2::ggplot(final_summary,
                       ggplot2::aes(x = median_Tread, y = frequency,
                                    color = dataset, fill = dataset)) +
    ggplot2::geom_area(position = "identity", alpha = 0.15, size = 0.8) +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_x_continuous(breaks = seq(0, reads_max, by = 10)) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_color_manual(values = pub_colors) +
    ggplot2::scale_fill_manual(values = pub_colors) +
    ggplot2::labs(
      x = "Median Read Depth",
      y = "Number of Sites",
      title = paste0("Methylation Signal Distribution by Consensus Tier", title_suffix),
      color = "Consensus",
      fill = "Consensus") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = c(0.85, 0.85),
      legend.background = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

  # Save if outdir provided
  if (!is.null(outdir)) {
    suffix <- if (show_filtered) "_Fig_Distribution_Consensus_filtered.pdf"
    else "_Fig_Distribution_Consensus_unfiltered.pdf"
    filename <- file.path(outdir, paste0(sample_name, suffix))
    ggplot2::ggsave(filename, p, width = 8, height = 6, device = "pdf")
    message("Saved: ", filename)
    return(invisible(p))
  }

  return(invisible(p))
}
