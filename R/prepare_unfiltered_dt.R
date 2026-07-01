#' Reshape collapsed CpG tables into long format for thresholding/plotting
#'
#' Converts a named list of strand-collapsed (but not yet filtered)
#' per-pipeline data.tables into one long-format table with `reads`,
#' `Strand` (S1/S2, derived from strand-support count `n`), and `Pipeline`.
#'
#' @param collapsed A named list of data.tables, each the output of
#'   [collapse_cpg_strand()] (must contain `TRead` and `n`).
#'
#' @return A long-format data.table with columns `reads`, `Pipeline`,
#'   `Strand`.
#' @export
prepare_unfiltered_dt <- function(collapsed) {
  data.table::rbindlist(Map(function(dt, name) {
    data.table::data.table(
      reads = dt$TRead,
      Pipeline = name,
      Strand = ifelse(dt$n == 1, "S1", "S2")
    )
  }, collapsed, names(collapsed)))
}
