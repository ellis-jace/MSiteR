# Internal helper shared by the filtered and unfiltered frequency plots.
# The public API remains plot_strand_frequency_filtered/unfiltered().
plot_strand_frequency_common <- function(result, filtered, outdir,
                                         sample_name, reads_max,
                                         show_thresholds) {
  if (!all(c("merged", "unfiltered", "thresholds") %in% names(result))) {
    stop("`result` must be output from run_strand_filtering().")
  }

  if (filtered) {
    wide <- result$merged
    tr_cols <- grep("_TR$", names(wide), value = TRUE)
    if (!length(tr_cols)) stop("No pipeline `_TR` columns in result$merged.")

    # Reconstruct the plotting-only long table from the documented wide output.
    pieces <- lapply(tr_cols, function(tr_col) {
      name <- sub("_TR$", "", tr_col)
      str_col <- paste0(name, "_StrN")
      if (!str_col %in% names(wide)) {
        stop("Missing strand-count column: ", str_col)
      }
      reads <- wide[[tr_col]]
      strands <- wide[[str_col]]
      valid <- !is.na(reads) & !is.na(strands) & strands %in% c(1L, 2L)
      data.table::data.table(
        reads = reads[valid], Pipeline = name,
        Strand = ifelse(strands[valid] == 1L, "S1", "S2")
      )
    })
    plot_dt <- data.table::rbindlist(pieces, use.names = TRUE)
  } else {
    plot_dt <- data.table::copy(result$unfiltered)
    if (!all(c("Pipeline", "Strand", "reads") %in% names(plot_dt))) {
      stop("result$unfiltered must contain Pipeline, Strand and reads.")
    }
  }

  frequency <- plot_dt[is.finite(reads) & reads >= 1 & reads <= reads_max,
                       .N, by = .(Pipeline, Strand, reads)]
  data.table::setnames(frequency, "N", "Frequency")
  frequency[, Pipeline := factor(Pipeline, levels = result$thresholds$Pipeline)]
  frequency[, Strand := factor(Strand, levels = c("S1", "S2"))]

  colors <- c(S1 = "#6d6c6cff", S2 = "#E64B35")
  p <- ggplot2::ggplot(frequency, ggplot2::aes(
    x = reads, y = Frequency, color = Strand, fill = Strand
  )) +
    ggplot2::geom_area(position = "identity", alpha = 0.2, linewidth = 0.8) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~Pipeline, scales = "free_y", nrow = 1) +
    ggplot2::scale_x_continuous(breaks = pretty(c(0, reads_max), n = 4)) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      x = "Total Reads Count", y = "Frequency (Number of Sites)",
      title = paste0("Single-strand and double-strand sites (",
                     if (filtered) "AFTER" else "BEFORE", " FILTERING)"),
      subtitle = if (show_thresholds) "Dashed lines: actual S1 and S2 cutoffs" else NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())

  if (show_thresholds) {
    thresholds <- result$thresholds
    marks <- data.table::rbindlist(list(
      data.table::data.table(Pipeline = thresholds$Pipeline,
                             Strand = "S1", cutoff = thresholds$cutoff1),
      data.table::data.table(Pipeline = thresholds$Pipeline,
                             Strand = "S2", cutoff = thresholds$cutoff2)
    ))
    marks[, Pipeline := factor(Pipeline, levels = thresholds$Pipeline)]
    marks[, Strand := factor(Strand, levels = c("S1", "S2"))]
    p <- p + ggplot2::geom_vline(
      data = marks, ggplot2::aes(xintercept = cutoff, color = Strand),
      linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE
    )
  }

  if (!is.null(outdir)) {
    suffix <- if (filtered) "_04Fig_Frequency_Filtered" else "_02Fig_Frequency_Unfiltered"
    if (show_thresholds) suffix <- paste0(suffix, "_with_Thresholds")
    filename <- file.path(outdir, paste0(sample_name, suffix, ".pdf"))
    ggplot2::ggsave(filename, p, width = 12, height = 5)
    message("Saved: ", filename)
  }
  invisible(p)
}
