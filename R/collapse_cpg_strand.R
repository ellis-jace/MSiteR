#' Collapse CpG Calls Across Complementary Strands
#'
#' Combines methylation calls from the forward (+) and reverse (-)
#' strands into a single CpG site. Reverse-strand positions are
#' shifted by one base pair to align with the forward-strand cytosine.
#'
#' Total and methylated read counts are summed across observations,
#' and methylation levels are recalculated from the combined counts.
#'
#' @param input A data.frame, data.table, or headerless file path
#'   containing CpG methylation calls. Required columns are
#'   `chr`, `pos`, `TRead`, and `MRead`. A `strand` column
#'   is required unless `strand_reference` is provided.
#'
#' @param strand_reference Optional reference containing
#'   `chr`, `CpG_pos`, and `strand`. Can be a data.frame,
#'   data.table, FASTA file path, or delimited table path.
#'   Used when the input lacks strand information.
#'
#' @return A data.table containing:
#' \describe{
#'   \item{chr}{Chromosome identifier.}
#'   \item{pos}{Forward-strand CpG position.}
#'   \item{TRead}{Combined total read count.}
#'   \item{MRead}{Combined methylated read count.}
#'   \item{ML}{Recalculated methylation level.}
#'   \item{n}{Number of distinct strands supporting the site (1 or 2).}
#' }
#'
#' @details
#' Input observations without a matching reference position are
#' discarded with a warning. Duplicate observations from the same
#' strand contribute to read counts but count as one strand.
#'
#' @export
collapse_cpg_strand <- function(input, strand_reference = NULL) {
  if (is.character(input) && length(input) == 1L) {
    dt <- read_pipeline_files(list(input = input))[[1L]]
  } else {
    dt <- data.table::copy(data.table::as.data.table(input))
  }

  needed <- c("chr", "pos", "TRead", "MRead")
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols)) {
    stop("`input` is missing column(s): ", paste(missing_cols, collapse = ", "))
  }

  empty_result <- function() data.table::data.table(
    chr = character(), pos = integer(), TRead = integer(),
    MRead = integer(), ML = numeric(), n = integer()
  )
  if (nrow(dt) == 0L) return(empty_result())

  if (!"strand" %in% names(dt)) {
    if (is.null(strand_reference)) {
      stop("`input` has no `strand` column and no `strand_reference` was provided.")
    }

    ref <- resolve_strand_reference(strand_reference)
    if (anyDuplicated(ref, by = c("chr", "CpG_pos"))) {
      stop("`strand_reference` has duplicate chr/CpG_pos keys.")
    }

    # Join on the call's position. Calls with no reference match are dropped.
    before <- nrow(dt)
    dt[, CpG_pos := pos]
    dt <- ref[dt, on = .(chr, CpG_pos), nomatch = NULL]
    if (nrow(dt) < before) {
      warning(before - nrow(dt), " CpG calls had no strand-reference match and were dropped.")
    }
  }

  # Empty input: return a typed empty table so downstream merges still work
  if (nrow(dt) == 0) {
    return(empty_result())
  }

  if (anyNA(dt$strand) || any(!dt$strand %in% c("+", "-"))) {
    stop("Strand must be '+' or '-' for every matched call.")
  }

  # Map reverse-strand calls (at the G) back to the CpG's C position.
  # fifelse() keeps the integer type; ifelse() can change it.
  dt[, CpG_pos := data.table::fifelse(strand == "-", as.integer(pos) - 1L, as.integer(pos))]

  merged <- dt[, .(
    pos   = CpG_pos[1],
    TRead = sum(TRead),
    MRead = sum(MRead),
    # Count distinct strands, not input records (which may contain duplicates).
    n     = data.table::uniqueN(strand)
  ), by = .(chr, CpG_pos)]

  merged[, ML := data.table::fifelse(TRead > 0, MRead / TRead, NA_real_)]
  data.table::setorder(merged, chr, pos)
  merged[, .(chr, pos, TRead, MRead, ML, n)]
}
