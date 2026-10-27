#' Resolve a strand reference into a data.table
#'
#' Accepts a genome FASTA path (`.fna`, `.fa`, `.fasta`, or any file whose
#' first line starts with `>`), a delimited table path, or a data.frame/data.table.
#' Returns a data.table with columns `chr`, `CpG_pos`, `strand`.
#'
#' @param strand_reference A FASTA file path, a table file path, or a
#'   data.frame/data.table.
#'
#' @return A data.table with columns `chr`, `CpG_pos`, `strand`.
#'
#' @export
resolve_strand_reference <- function(strand_reference) {
  required <- c("chr", "CpG_pos", "strand")

  if (is.character(strand_reference)) {
    if (length(strand_reference) != 1) {
      stop("`strand_reference` path must be a single string.")
    }
    path <- path.expand(strand_reference)
    if (!file.exists(path)) {
      stop("`strand_reference` file not found: ", path)
    }

    first_line <- readLines(path, n = 1, warn = FALSE)
    is_fasta <- grepl("\\.(fna|fa|fasta)(\\.gz)?$", path, ignore.case = TRUE) ||
      (length(first_line) == 1 && startsWith(first_line, ">"))

    if (is_fasta) {
      ref <- make_cpg_table(path)
      if (is.null(ref)) {
        stop("No CpG sites found in FASTA file: ", path)
      }
    } else {
      ref <- data.table::fread(path)
    }
  } else if (is.data.frame(strand_reference)) {
    ref <- data.table::copy(data.table::as.data.table(strand_reference))
  } else {
    stop("`strand_reference` must be a file path (FASTA or table) or a data.frame/data.table.")
  }

  ref <- data.table::as.data.table(ref)
  missing_cols <- setdiff(required, names(ref))
  if (length(missing_cols) > 0) {
    stop("`strand_reference` is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". It has: ", paste(names(ref), collapse = ", "))
  }
  ref[, .(chr, CpG_pos, strand)]
}
