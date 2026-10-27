collapse_cpg_strand <- function(input, strand_reference = NULL) {
  dt <- data.table::copy(data.table::as.data.table(input))

  if (!"strand" %in% names(dt)) {
    if (is.null(strand_reference)) {
      stop("`input` has no `strand` column and no `strand_reference` was provided.")
    }

    ref <- resolve_strand_reference(strand_reference)

    # Join on the call's position. Calls with no reference match are dropped.
    dt[, CpG_pos := pos]
    dt <- ref[dt, on = .(chr, CpG_pos), nomatch = NULL]
  }

  # Empty input: return a typed empty table so downstream merges still work
  if (nrow(dt) == 0) {
    return(data.table::data.table(
      chr   = character(),
      pos   = integer(),
      TRead = integer(),
      MRead = integer(),
      ML    = numeric(),
      n     = integer()
    ))
  }

  # Map reverse-strand calls (at the G) back to the CpG's C position.
  # fifelse() keeps the integer type; ifelse() can change it.
  dt[, CpG_pos := data.table::fifelse(strand == "-", as.integer(pos) - 1L, as.integer(pos))]

  merged <- dt[, .(
    pos   = CpG_pos[1],
    TRead = sum(TRead),
    MRead = sum(MRead),
    n     = .N
  ), by = .(chr, CpG_pos)]

  merged[, ML := MRead / TRead]
  data.table::setorder(merged, chr, pos)
  merged[, .(chr, pos, TRead, MRead, ML, n)]
}
