#' Run consensus filtering across pipelines
#'
#' Convenience wrapper for the second stage of the pipeline. Takes output from
#' [run_strand_filtering()], adds consensus metrics (how many pipelines detected
#' each site), calculates per-tier filtering thresholds, and filters sites based
#' on cross-pipeline agreement.
#'
#' Each underlying step ([calculate_consensus_stats()],
#' [calculate_consensus_thresholds()], [filter_by_consensus()]) remains
#' independently callable if you need to inspect or re-run an intermediate stage.
#'
#' @param result Output from [run_strand_filtering()]
#'
#' @return A list with components:
#'   - `merged`: Filtered consensus table (wide format, chr+pos+pipelines)
#'   - `unfiltered`: Same table before consensus filtering (for QC/plotting)
#'   - `thresholds`: Data.table with per-tool, per-tier (1M/2M/3M/4M) thresholds
#'   - `filtering_stats`: Data.table with overall before/after counts and retention %
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Stage 1: Strand filtering
#' result_strand <- run_strand_filtering(
#'   pipelines = list(
#'     Bismark = bismark_dt,
#'     Bwameth = bwameth_dt,
#'     Biscuit = biscuit_dt,
#'     Encode = encode_dt
#'   )
#' )
#'
#' # Stage 2: Consensus filtering
#' result_consensus <- run_consensus_filtering(result_strand)
#'
#' # Access results
#' final_cpgs <- result_consensus$merged
#' stats <- result_consensus$filtering_stats
#' thresholds <- result_consensus$thresholds
#' }
run_consensus_filtering <- function(result) {

  # Validate Input
  if (!is.list(result)) {
    stop("result must be a list (output from run_strand_filtering())")
  }

  required_components <- c("merged", "unfiltered", "thresholds")
  missing <- setdiff(required_components, names(result))

  if (length(missing) > 0) {
    stop("result missing required components: ", paste(missing, collapse = ", "),
         "\nExpected output from run_strand_filtering()")
  }

  # Validate that merged and unfiltered are data.tables
  if (!data.table::is.data.table(result$merged)) {
    stop("result$merged must be a data.table")
  }

  if (!data.table::is.data.table(result$unfiltered)) {
    stop("result$unfiltered must be a data.table")
  }

  # STEP 1: Extract merged table and add consensus metrics
  message("Adding consensus metrics...")
  dt <- data.table::copy(result$merged)
  dt <- calculate_consensus_stats(dt)

  # Validate that n_methods was added
  if (!("n_methods" %in% names(dt))) {
    stop("Failed to add n_methods column. Check calculate_consensus_stats().")
  }

  message("  - Added: n_methods, median_Tread, median_MethyL, median_Mread")

  # STEP 2: Calculate consensus-tier thresholds
  message("Calculating consensus thresholds...")
  consensus_thresholds <- calculate_consensus_thresholds(dt)

  if (!data.table::is.data.table(consensus_thresholds)) {
    stop("calculate_consensus_thresholds() failed to return a data.table")
  }

  # STEP 3: Filter by consensus thresholds
  message("Applying consensus filtering...")
  dt_filtered <- filter_by_consensus(dt, consensus_thresholds)

  if (!data.table::is.data.table(dt_filtered)) {
    stop("filter_by_consensus() failed to return a data.table")
  }

  # STEP 4: Collect filtering statistics (original: total + 1M tier)
  # Overall
  n_before <- nrow(dt)
  n_after <- nrow(dt_filtered)

  # 1M sites specifically (strictest filtering)
  dt_n1 <- dt[n_methods == 1]
  n_1m_before <- nrow(dt_n1)
  n_1m_after <- nrow(dt_filtered[n_methods == 1])
  n_1m_retain_pct <- if (n_1m_before > 0) round(n_1m_after / n_1m_before * 100, 2) else 0

  filtering_stats <- data.table::data.table(
    metric = c("sites_before", "sites_after", "1m_sites_before", "1m_sites_after", "1m_retain_pct"),
    value = c(n_before, n_after, n_1m_before, n_1m_after, n_1m_retain_pct)
  )

  # LOG SUMMARY
  message("Consensus filtering complete:")
  message("  Methylation sites before: ", n_before)
  message("  Methylation sites after: ", n_after)
  message("  Methylation sites (1M) before: ", n_1m_before)
  message("  Methylation sites (1M) after: ", n_1m_after)
  message("  Retain rate (1M): ", scales::percent(n_1m_after / n_1m_before))

  # RETURN RESULTS
  list(
    merged = dt_filtered,
    unfiltered = dt,
    thresholds = consensus_thresholds,
    filtering_stats = filtering_stats
  )
}

