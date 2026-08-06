#' Helper to run the six-function pipeline on a single chromosome's pipelines
#'
#' Run each data pipeline on pre-chunked chromosomes from the output of
#' [chunk_by_chromosome()], or from a single-chromosome data.table. The underlying
#' steps ([collapse_cpg_strand()], [prepare_unfiltered_dt()],
#' [filter_by_strand()]) remain independently callable for personal inspection
#' plotting.
#'
#' @param chr_pipelines A named list of pipeline data.tables subsetted to a single
#' chromosome (e.g. from [chunk_by_chromosome()], or manually created).
#' @param strand_reference File path or data.table providing strand info for
#'   pipelines whose raw files lack a `strand` column, matched against
#'   `strand_reference_for`.
#' @param strand_reference_for Character vector of pipeline names needing
#'   `strand_reference` (default `"Biscuit"`).
#'
#' @return A list with `merged` (the merged, strand-filtered
#'   data.table), and `unfiltered` (raw, chromosome-specific data.table)
#' @export
run_pipeline_on_chunk <- function(chr_pipelines,
                                  strand_reference = NULL,
                                  strand_reference_for = 'Biscuit') {
  chr_collapsed <- Map(function(x, name) {
    ref <- if (name %in% strand_reference_for) strand_reference else NULL
    collapse_cpg_strand(x, strand_reference = ref)
  }, chr_pipelines, names(chr_pipelines))

  chr_unfiltered <- prepare_unfiltered_dt(chr_collapsed)
  chr_thresholds <- calculate_strand_thresholds(chr_unfiltered)
  chr_filtered <- Map(function(dt, name) filter_by_strand(dt, name, chr_thresholds),
                      chr_collapsed, names(chr_collapsed))
  chr_merged <- merge_filtered_pipelines(chr_filtered)

  list(merged = chr_merged, unfiltered = chr_unfiltered)
}
