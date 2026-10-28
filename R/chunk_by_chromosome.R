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

  pipelines <- read_pipeline_files(pipelines)

  # Extract unique chromosomes and sort
  chromosomes <- sort(unique(do.call(c, lapply(pipelines, function(x) x$chr))))

  # Subset each pipeline to each chromosome
  chunks <- lapply(setNames(chromosomes, chromosomes), function(chr_name) {
    lapply(pipelines, function(x) x[x$chr == chr_name, ])
  })

  chunks
}
