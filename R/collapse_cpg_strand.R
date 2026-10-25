#' Collapse +/- strand CpG calls into symmetric per-site totals
#'
#' Collapses complementary +/- strand rows for the same CpG into one site,
#' summing read counts and recomputing methylation level. Input must already
#' have standard column names (see [read_pipeline_files()]). If `input` has no
#' `strand` column (e.g. Biscuit output), strand is joined in from
#' `strand_reference`.
#'
#' @param input A data.table with columns `chr, pos, TRead, MRead, ML`, and
#'   optionally `strand`.
#' @param strand_reference Optional data.table with columns `chr`, `CpG_pos`,
#'   and `strand` (e.g. from [make_cpg_table()]). Required when `input` has no
#'   `strand` column.
#'
#' @return A data.table with columns `chr, pos, TRead, MRead, ML, n`, one
#'   row per symmetric CpG site. `pos` is the CpG's C position.
#' @export
collapse_cpg_strand <- function(input, strand_reference = NULL) {
  dt <- data.table::copy(data.table::as.data.table(input))

  if (!"strand" %in% names(dt)) {
    if (is.null(strand_reference)) {
      stop("`input` has no `strand` column and no `strand_reference` was provided.")
    }

    ref <- data.table::as.data.table(strand_reference)
    missing_cols <- setdiff(c("chr", "CpG_pos", "strand"), names(ref))
    if (length(missing_cols) > 0) {
      stop("`strand_reference` is missing column(s): ",
           paste(missing_cols, collapse = ", "))
    }
    ref <- ref[, .(chr, CpG_pos, strand)]

    # Join on the call's position. Calls with no reference match are dropped.
    dt[, CpG_pos := pos]
    dt <- ref[dt, on = .(chr, CpG_pos), nomatch = NULL]
  }

  # Map reverse-strand calls (at the G) back to the CpG's C position
  dt[, CpG_pos := ifelse(strand == "-", pos - 1L, pos)]

  merged <- dt[, .(
    pos    = CpG_pos[1],
    TRead  = sum(TRead),
    MRead  = sum(MRead),
    n      = .N
  ), by = .(chr, CpG_pos)]

  merged[, ML := MRead / TRead]
  data.table::setorder(merged, chr, pos)
  merged[, .(chr, pos, TRead, MRead, ML, n)]
}
