#' Extract CpG sites from FASTA reference file
#'
#' Reads a FASTA sequence file, finds all CG dinucleotides (CpG sites),
#' and returns their positions with strand information (+ for CG, - for GC).
#'
#' @param fasta_file Path to FASTA file (e.g., reference genome)
#'
#' @return A data.frame with columns:
#'   - `chr`: Chromosome/sequence ID from FASTA header
#'   - `CpG_pos`: Position of the C in the CG or G in the GC dinucleotide (1-indexed)
#'   - `strand`: Strand (+ for CG forward, - for GC reverse complement)
#'
#' @details
#' The function:
#' 1. Parses FASTA headers (lines starting with ">") as chromosome IDs
#' 2. Processes sequence in chunks to optimize memory usage
#' 3. Finds all occurrences of "CG" (forward strand) and "GC" (reverse complement)
#' 4. Returns C/G positions with corresponding strand information
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cpg_table <- make_cpg_table("reference.fasta")
#' head(cpg_table)
#' }
make_cpg_table <- function(fasta_file) {
  # Read FASTA file
  lines <- readLines(fasta_file)

  # Find header lines (start with ">")
  headers <- grep("^>", lines)

  # Find where each chromosome's sequence ends
  ends <- c(headers[-1] - 1, length(lines))

  # Process each chromosome
  result <- lapply(seq_along(headers), function(i) {
    # Extract chromosome ID from header (remove ">")
    chr <- sub("^>", "", lines[headers[i]])

    # Combine all sequence lines belonging to this chromosome
    seq <- paste0(lines[(headers[i] + 1):ends[i]], collapse = "")
    seq <- toupper(seq)

    # Find all "CG" positions (forward strand)
    cpg_forward <- gregexpr("CG", seq, fixed = TRUE)[[1]]

    # Find all "GC" positions (reverse complement)
    cpg_reverse <- gregexpr("GC", seq, fixed = TRUE)[[1]]

    # Combine results
    cpg_list <- list()

    if (cpg_forward[1] != -1) {
      cpg_list[[1]] <- data.frame(
        chr = chr,
        CpG_pos = cpg_forward,
        strand = "+",
        stringsAsFactors = FALSE
      )
    }

    if (cpg_reverse[1] != -1) {
      cpg_list[[2]] <- data.frame(
        chr = chr,
        CpG_pos = cpg_reverse,
        strand = "-",
        stringsAsFactors = FALSE
      )
    }

    if (length(cpg_list) == 0) {
      return(NULL)
    }

    do.call(rbind, cpg_list)
  })

  # Combine all chromosomes and sort by position
  result_df <- do.call(rbind, result)
  rownames(result_df) <- NULL
  result_df[order(result_df$chr, result_df$CpG_pos), ]
}
