#' Venn diagram: Pipeline consensus overlaps
#'
#' Shows 4-way overlap of pipeline detection using a Venn diagram.
#' Displays percentages of sites in each region.
#'
#' @param result Output from [run_consensus_filtering()]
#' @param outdir Optional directory to save PDF. If NULL, plot is not saved.
#' @param sample_name Sample ID for filename (default "sample")
#' @param show_filtered Logical: show filtered data? If FALSE (default), shows
#'   unfiltered data. If TRUE, shows filtered data.
#'
#' @return ggplot object (invisibly if saved)
#' @export
plot_consensus_venn <- function(result,
                                outdir = NULL,
                                sample_name = "sample",
                                show_filtered = FALSE) {

  # Validate input
  if (!all(c("merged", "unfiltered") %in% names(result))) {
    stop("result must be output from run_consensus_filtering()")
  }

  # Select data based on show_filtered
  dt <- if (show_filtered) result$merged else result$unfiltered

  # Create logical matrix
  upset_dt <- dt[, .(
    Bismark = !is.na(Bismark_TR),
    Bwameth = !is.na(Bwameth_TR),
    Biscuit = !is.na(Biscuit_TR),
    ENCODE = !is.na(Encode_TR)
  )]

  # Convert to list of indices for Venn diagram
  venn_list <- list(
    Bismark = which(upset_dt$Bismark == TRUE),
    Bwameth = which(upset_dt$Bwameth == TRUE),
    Biscuit = which(upset_dt$Biscuit == TRUE),
    ENCODE = which(upset_dt$ENCODE == TRUE)
  )

  # Define colors
  venn_colors <- c("#D32F2F", "#1976D2", "#388E3C", "#7B1FA2")

  # Build Venn diagram
  title_suffix <- if (show_filtered) " (After Consensus Filtering)" else " (Before Consensus Filtering)"

  pvenn <- ggvenn::ggvenn(
    venn_list,
    columns = c("Bismark", "Bwameth", "Biscuit", "ENCODE"),
    fill_color = venn_colors,
    fill_alpha = 0.8,
    stroke_size = 0,
    stroke_color = "transparent",
    set_name_size = 8,
    set_name_color = "black",
    text_size = 6,
    text_color = "white",
    show_percentage = TRUE,
    digits = 2) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold", color = "black"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)) +
    ggplot2::labs(title = paste0("Consensus of Methylation Sites Across 4 Pipelines", title_suffix))

  # Save if outdir provided
  if (!is.null(outdir)) {
    suffix <- if (show_filtered) "_Fig_Venn_Consensus_filtered.pdf"
    else "_Fig_Venn_Consensus_unfiltered.pdf"
    filename <- file.path(outdir, paste0(sample_name, suffix))
    ggplot2::ggsave(filename, pvenn, width = 10, height = 10, device = "pdf")
    message("Saved: ", filename)
    return(invisible(pvenn))
  }

  return(invisible(pvenn))
}
