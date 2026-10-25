#' Collapse +/- strand CpG calls into symmetric per-site totals
#'
#' Reads a single pipeline's raw per-strand methylation calls and collapses
#' complementary +/- strand rows for the same CpG dinucleotide into one
#' symmetric site, summing read counts and recomputing methylation level.
#' If the input lacks a `strand` column (e.g. Biscuit output), strand is
#' joined in from an external reference table.
#'
#' @param input Either a file path (character) to a tab-delimited file with
#'   columns `chr, pos, TRead, MRead, ML[, strand]` (no header), or an
#'   already-loaded data.table with those columns.
#' @param strand_reference Optional file path or data.table providing
#'   `chr, CpG_pos, strand` when `input` doesn't already contain a `strand`
#'   column (required for Biscuit-style input).
#'
#' @return A data.table with columns `chr, pos, TRead, MRead, ML, n`, one
#'   row per symmetric CpG site.
#' @export
#'
collapse_cpg_strand <- function(input, strand_reference = NULL) {
  dt <- if (is.character(input)) {
    data.table::fread(input, header = FALSE)
  } else {
    data.table::copy(input)
  }

  # Strand is present if there are 6 columns (positional) or a column named "strand"
  has_strand <- ncol(dt) == 6L || "strand" %in% names(dt)

  if (!has_strand) {
    if (is.null(strand_reference)) {
      stop("`input` has no `strand` column and no `strand_reference` was provided.")
    }
    data.table::setnames(dt, c("chr", "pos", "TRead", "MRead", "ML"))

    ref <- if (is.character(strand_reference)) {
      data.table::fread(strand_reference)
    } else {
      data.table::copy(strand_reference)
    }

    # Case-insensitive column matching
    data.table::setnames(ref, tolower(names(ref)))
    needed <- c("chr", "cpg_pos", "strand")
    miss <- setdiff(needed, names(ref))
    if (length(miss)) {
      stop("`strand_reference` is missing column(s): ", paste(miss, collapse = ", "),
           ". It has: ", paste(names(ref), collapse = ", "))
    }
    ref <- ref[, .(chr = chr, CpG_pos = cpg_pos, strand = strand)]

    dt <- ref[dt, on = .(chr, CpG_pos = pos)]
    dt <- dt[!is.na(strand)]
    data.table::setnames(dt, c("chr", "pos", "strand", "TRead", "MRead", "ML"))
  } else {
    data.table::setnames(dt, c("chr", "pos", "TRead", "MRead", "ML", "strand"))
  }

  dt[, CpG_pos := ifelse(strand == "-", pos - 1, pos)]
  merged <- dt[, {
    if (.N == 2 && all(c("+", "-") %in% strand)) {
      list(pos = pos[strand == "+"][1], TRead = sum(TRead), MRead = sum(MRead), n = .N)
    } else {
      list(pos = pos[1], TRead = sum(TRead), MRead = sum(MRead), n = .N)
    }
  }, by = .(chr, CpG_pos)]

  merged[, ML := MRead / TRead]
  data.table::setorder(merged, chr, pos)
  merged[, .(chr, pos, TRead, MRead, ML, n)]
}
