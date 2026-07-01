#' Collapse strand, calculate thresholds, filter, and merge pipelines
#'
#' Convenience wrapper running the full first stage of the pipeline in one
#' call: collapses +/- strand calls per pipeline, computes per-pipeline
#' strand-depth thresholds, filters each pipeline using those thresholds,
#' then joins the filtered tables by chr/pos into one wide table.
#'
#' Each underlying step ([collapse_cpg_strand()], [prepare_unfiltered_dt()],
#' [calculate_strand_thresholds()], [filter_by_strand()],
#' [merge_filtered_pipelines()]) remains independently callable if you need
#' to inspect or re-plot an intermediate stage (e.g. thresholds, or
#' unfiltered vs. filtered distributions).
#'
#' @param pipelines A named list of raw file paths or data.tables, one per
#'   pipeline (e.g. Bismark, Bwameth, Biscuit, Encode).
#' @param strand_reference File path or data.table providing strand info for
#'   pipelines whose raw files lack a `strand` column, matched against
#'   `strand_reference_for`.
#' @param strand_reference_for Character vector of pipeline names needing
#'   `strand_reference` (default `"Biscuit"`).
#'
#' @return A list with `merged` (the joined, strand-filtered wide
#'   data.table), `thresholds` (the per-pipeline threshold table), and
#'   `unfiltered` (the long-format pre-filter table, for plotting/QC).
#' @export
prepare_filtered_cpg_table <- function(pipelines,
                                       strand_reference = NULL,
                                       strand_reference_for = "Biscuit") {

  # 1. collapse strand, per pipeline
  collapsed <- Map(function(x, name) {
    ref <- if (name %in% strand_reference_for) strand_reference else NULL
    collapse_cpg_strand(x, strand_reference = ref)
  }, pipelines, names(pipelines))

  # 2. long-format table + thresholds
  unfiltered <- prepare_unfiltered_dt(collapsed)
  thresholds <- calculate_strand_thresholds(unfiltered)

  # 3. filter each pipeline, then merge
  filtered <- Map(function(dt, name) filter_by_strand(dt, name, thresholds),
                  collapsed, names(collapsed))

  merged <- merge_filtered_pipelines(filtered)

  list(merged = merged, thresholds = thresholds, unfiltered = unfiltered)
}
