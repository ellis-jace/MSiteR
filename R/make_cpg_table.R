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
  genome <- Biostrings::readDNAStringSet(fasta_file)
  # keep only the first token of the header as the chromosome ID
  names(genome) <- sub("\\s.*$", "", names(genome))

  hits   <- Biostrings::vmatchPattern("CG", genome)
  starts <- IRanges::start(hits)              # one integer vector per sequence

  chr <- rep(names(genome), lengths(starts))
  pos <- unlist(starts, use.names = FALSE)

  out <- data.table::data.table(
    chr     = c(chr, chr),
    CpG_pos = c(pos, pos + 1L),
    strand  = rep(c("+", "-"), each = length(pos)),
    stringsAsFactors = FALSE
  )
  out[order(out$chr, out$CpG_pos), ]
}
