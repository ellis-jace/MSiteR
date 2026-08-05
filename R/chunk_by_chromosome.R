#' Partition pipelines by chromosome
#'
#' Subset each pipeline's data.table to individual chromosomes. Useful for
#' memory-efficient processing or custom per-chromosome workflows.
#'
#' @param pipelines A named list of raw file paths or data.tables, one per
#'   pipeline (e.g. Bismark, Bwameth, Biscuit, Encode).
#'
#' @return A list of lists: `list(chr1 = list(pipeline1 = ..., pipeline2 = ...),
#'   chr2 = ..., ...)` where each chromosome's element contains the subsetted
#'   pipelines for that chromosome.
#' @export
chunk_by_chromosome <- function(pipelines) {
  # Read file paths into data.tables
  read_if_path <- function(x) {
    if (is.character(x)) data.table::fread(x) else x
  }
  all_data <- lapply(pipelines, read_if_path)

  # Extract unique chromosomes and sort
  chromosomes <- sort(unique(do.call(c, lapply(all_data, function(x) x$chr))))

  # Subset each pipeline to each chromosome
  chunks <- lapply(setNames(chromosomes, chromosomes), function(chr) {
    lapply(all_data, function(x) x[x$chr == chr, ])
  })

  chunks
}
