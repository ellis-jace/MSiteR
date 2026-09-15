#' UpSet plot: Pipeline consensus overlaps
#'
#' Shows which combinations of pipelines detected each site using an UpSet plot.
#' UpSet is more interpretable than Venn diagrams for 4-way comparisons.
#'
#' @param result Output from [run_consensus_filtering()]
#' @param outdir Optional directory to save PDF. If NULL, plot is not saved.
#' @param sample_name Sample ID for filename (default "sample")
#' @param show_filtered Logical: show filtered data? If FALSE (default), shows
#'   unfiltered data (before consensus filtering). If TRUE, shows filtered data.
#'
#' @return ggplot object (invisibly if saved)
#' @export
plot_consensus_upset <- function(result,
                                 outdir = NULL,
                                 sample_name = "sample",
                                 show_filtered = FALSE) {

  # Validate input
  if (!all(c("merged", "unfiltered") %in% names(result))) {
    stop("result must be output from run_consensus_filtering()")
  }

  # Select data based on show_filtered
  dt <- if (show_filtered) result$merged else result$unfiltered

  # Create logical matrix: TRUE if pipeline detected this site
  upset_dt <- dt[, .(
    Bismark = !is.na(Bismark_TR),
    Bwameth = !is.na(Bwameth_TR),
    Biscuit = !is.na(Biscuit_TR),
    ENCODE = !is.na(Encode_TR)
  )]

  plot_df <- as.data.frame(upset_dt)
  setDT(plot_df)
  cols <- c("Bismark", "Bwameth", "Biscuit", "ENCODE")

  # Convert to integer (1/0) for UpSet
  for (col in cols) {
    data.table::set(plot_df, j = col, value = as.integer(plot_df[[col]]))
  }

  # Color scheme for publication
  main_color <- "#285291"  # Solid deep blue
  sets_color <- "#555555"  # Neutral gray

  # Generate UpSet plot
  pdf_file <- tempfile(fileext = ".pdf")
  pdf(pdf_file, width = 10, height = 7, onefile = FALSE)

  UpSetR::upset(plot_df,
                sets = cols,
                order.by = "freq",
                decreasing = TRUE,
                nintersects = 15,
                text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.3),
                main.bar.color = main_color,
                sets.bar.color = sets_color,
                matrix.color = main_color,
                shade.color = "#E0E0E0",
                mb.ratio = c(0.65, 0.35),
                point.size = 3.5,
                line.size = 1,
                mainbar.y.label = "Intersection Size (CpG Sites)",
                sets.x.label = "Total Sites per Method")
  dev.off()

  # Save if outdir provided
  if (!is.null(outdir)) {
    suffix <- if (show_filtered) "_Fig_UpSet_Consensus_filtered.pdf"
    else "_Fig_UpSet_Consensus_unfiltered.pdf"
    filename <- file.path(outdir, paste0(sample_name, suffix))
    file.copy(pdf_file, filename)
    file.remove(pdf_file)
    message("Saved: ", filename)
    return(invisible(NULL))
  }

  file.remove(pdf_file)
  return(invisible(NULL))
}
