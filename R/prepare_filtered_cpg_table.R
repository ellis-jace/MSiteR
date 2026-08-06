#' Collapse strand, calculate thresholds, filter, and merge pipelines
#'
#' Convenience wrapper running the full first stage of the pipeline in one
#' call: collapses +/- strand calls per pipeline, computes per-pipeline
#' strand-depth thresholds, filters each pipeline using those thresholds,
#' then joins the filtered tables by chr/pos into one wide table.
#'
#' Each underlying step ([collapse_cpg_strand()], [prepare_unfiltered_dt()],
#' [calculate_strand_thresholds()], [filter_by_strand()],
#' [merge_filtered_pipelines()], which are called from [run_pipeline_on_chunk()])
#' remains independently callable if you need to inspect or re-plot an
#' intermediate stage. Use [chunk_by_chromosome()] to partition data by
#' chromosome before custom workflows.
#'
#' @param pipelines A named list of raw file paths or data.tables, one per
#'   pipeline (e.g. Bismark, Bwameth, Biscuit, Encode).
#' @param strand_reference File path or data.table providing strand info for
#'   pipelines whose raw files lack a `strand` column, matched against
#'   `strand_reference_for`.
#' @param strand_reference_for Character vector of pipeline names needing
#'   `strand_reference` (default `"Biscuit"`).
#' @param chunk_by_chromosome Logical; if TRUE, process each chromosome
#'   separately to optimize memory usage (default FALSE).
#'
#' @return A list with `merged` (the joined, strand-filtered wide
#'   data.table), `thresholds` (the per-pipeline threshold table), and
#'   `unfiltered` (the long-format pre-filter table, for plotting/QC).
#' @export
prepare_filtered_cpg_table <- function(pipelines,
                                       strand_reference = NULL,
                                       strand_reference_for = "Biscuit",
                                       chunk_by_chromosome = FALSE) {

  if (!chunk_by_chromosome) {
    # Original behavior: process full dataset in one pass
    result <- run_pipeline_on_chunk(pipelines,
                                    strand_reference,
                                    strand_reference_for)
    thresholds <- calculate_strand_thresholds(result$unfiltered)

    list(merged = result$merged, thresholds = thresholds, unfiltered = result$unfiltered)

  } else {
    # Chunked behavior: process one chromosome at a time
    chunks <- chunk_by_chromosome(pipelines)

    all_merged <- list()
    all_unfiltered <- list()

    for (chr in names(chunks)) {
      result <- run_pipeline_on_chunk(chunks[[chr]],
                                      strand_reference,
                                      strand_reference_for)
      all_merged[[chr]] <- result$merged
      all_unfiltered[[chr]] <- result$unfiltered

      rm(result)
      gc(verbose = FALSE)
    }

    # Combine and recalculate global thresholds
    unfiltered <- data.table::rbindlist(all_unfiltered)
    merged <- data.table::rbindlist(all_merged)
    thresholds <- calculate_strand_thresholds(unfiltered)

    list(merged = merged, thresholds = thresholds, unfiltered = unfiltered)
  }
}
