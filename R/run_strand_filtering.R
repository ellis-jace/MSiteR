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
run_strand_filtering <- function(pipelines,
                                       strand_reference = NULL,
                                       strand_reference_for = "Biscuit",
                                       chunk_by_chromosome = FALSE) {

  # Validate and normalize consistently, including already-loaded data.tables.
  pipelines <- read_pipeline_files(pipelines)

  # Parse the reference once. A table passes through unchanged.
  if (!is.null(strand_reference)) {
    strand_reference <- resolve_strand_reference(strand_reference)
  }

  if (!chunk_by_chromosome) {
    # Original behavior: process full dataset in one pass
    result <- run_pipeline_on_chunk(pipelines,
                                    strand_reference,
                                    strand_reference_for)
    thresholds <- calculate_strand_thresholds(result$unfiltered)

    list(merged = result$merged, thresholds = thresholds, unfiltered = result$unfiltered)

  } else {
    # Do not construct an in-memory list containing all chromosome subsets.
    chromosomes <- sort(unique(unlist(lapply(pipelines, function(x) x$chr),
                                      use.names = FALSE)))
    if (length(chromosomes) == 0L) stop("No chromosomes in the input pipelines.")

    get_chunk <- function(chr_name) {
      lapply(pipelines, function(x) x[x$chr == chr_name, ])
    }
    get_reference <- function(chr_name) {
      if (is.null(strand_reference)) return(NULL)
      strand_reference[strand_reference$chr == chr_name, ]
    }

    # Pass 1: collect pre-filter data to calculate ONE set of global cutoffs.
    all_unfiltered <- lapply(chromosomes, function(chr_name) {
      run_pipeline_on_chunk(get_chunk(chr_name),
                            strand_reference = get_reference(chr_name),
                            strand_reference_for = strand_reference_for,
                            unfiltered_only = TRUE)
    })
    unfiltered <- data.table::rbindlist(all_unfiltered, use.names = TRUE)
    rm(all_unfiltered)
    thresholds <- calculate_strand_thresholds(unfiltered)

    # Pass 2: filter each chromosome using the SAME global cutoffs.
    all_merged <- lapply(chromosomes, function(chr_name) {
      run_pipeline_on_chunk(get_chunk(chr_name),
                            strand_reference = get_reference(chr_name),
                            strand_reference_for = strand_reference_for,
                            thresholds = thresholds)$merged
    })
    merged <- data.table::rbindlist(all_merged, use.names = TRUE, fill = TRUE)

    list(merged = merged, thresholds = thresholds, unfiltered = unfiltered)
  }
}
