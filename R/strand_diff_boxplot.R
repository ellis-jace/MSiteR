#' Boxplot for Strand 1 vs Strand 2.
#' @param result Output from [prepare_filtered_cpg_table()]
#' @param outdir Where to save PDFs
#' @export
plot_strand_filtering <- function(result, outdir = "."){
  unfiltered <- result$unfiltered

  # Use unfiltered data.table for before plot
  plot_dt <- unfiltered[, Pipeline := factor(Pipeline, levels = result$thresholds['Pipeline', ])]
  plot_dt[, Strand := factor(Strand, levels = c("S1", "S2"))]

  # Optimization: Downsample large datasets to 1M rows per group
  # for faster statistical tests and plotting.
  plot_sub <- plot_dt[, .SD[sample(.N, min(.N, 1000000))], by = .(Pipeline, Strand)]

  # Generate boxplot with statistical annotations
  p_box <- ggplot2::ggplot(plot_sub, aes(x = Pipeline, y = reads, fill = Strand)) +
    ggplot2::geom_boxplot(outlier.shape = NA,    # Hide individual outlier points)
                          notch = T,
                          width = 0.6,
                          size = 0.7,
                          color = 'black') +
    # Color scheme: S1 light gray, S2 warm red
    scale_fill_manual(values = c("S1" = "#BDBDBD", "S2" = "#E64B35"),
                      labels = c("S1" = "Single-strand", "S2" = "Double-strand")) +
    # Add Wilcoxon test p-values
    ggpubr::stat_compare_means(aes(group = Strand),
                               label = "p.format",
                               method = "wilcox.test",
                               label.y = 75,
                               size = 4,
                               fontface = 'italic') +
    # Add significance stars (* p < 0.05, ** p < 0.01, etc.)
    ggpubr::stat_compare_means(aes(group = Strand),
                               label = 'p.signif',
                               method = 'wilcox.test',
                               label.y = 71,
                               size = 5) +
    # Focus on core range (0-80 read covers 99% of data)
    ggplot2::coord_cartesian(ylim = c(0,80)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 80, 20),
                                expand = ggplot2::expansion(mult = c(0,0.0,5))) +
    # Publication-quality theme
    ggpubr::theme_pubr(base_size = 14, legend = 'top') +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = 'bold'),
      axis.text = ggplot2::element_text(color = 'black'),
      plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold', size = 16),
      panel.grid.major.y = ggplot2::element_line(color = 'gray90', linetype = 'dashed')) +
    ggplot2::labs(x = 'Pipelines',
                  y = 'Total Reads Count',
                  title = 'Total Reads count distribution of S1 and S2',
                  fill = 'Strand Status')

    # Save figure to provided out-directory
    ggplot2::ggsave(paste0(outdir,'/',sample,'_01Fig_Boxplot_S1_S2.pdf'), p_box,
                    width = 8, height = 6)
}
