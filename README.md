# MSiteR

Tools for filtering, comparing, and merging CpG methylation calls across multiple alignment/calling pipelines (Bismark, BWA-meth, Biscuit, ENCODE), including strand-level read-depth thresholding and strand-symmetric CpG merging.

This package grew out of a manual, script-based workflow for comparing CpG methylation calls across four pipelines on sheep/cattle genomic references. It replaces four near-duplicate per-pipeline scripts with a small set of tested, reusable functions.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("ellis-jace/MSiteR")
```

Or, for local development:

``` r
devtools::load_all("path/to/MSiteR")
```

## Pipeline overview

The package currently covers **Stage 1** of the methylation-calling workflow: turning raw, per-strand pipeline output into one strand-filtered, cross-pipeline table.

```         
raw pipeline files (+ optional strand reference)
        │
        ▼
chunk_by_chromosome()       # if necessary, separate chromosomes for individual pipeline calls
        │
        ▼
collapse_cpg_strand()        # collapse +/- strand pairs into symmetric sites, per pipeline
        │
        ▼
prepare_unfiltered_dt()      # reshape collapsed tables into long format for thresholding/plotting
        │
        ▼
calculate_strand_thresholds()# per-pipeline S1 (single-strand) / S2 (double-strand) read-depth cutoffs
        │
        ▼
filter_by_strand()           # apply thresholds, per pipeline
        │
        ▼
merge_filtered_pipelines()   # join filtered pipelines by chr + pos into one wide table
```

`prepare_filtered_cpg_table()` is a convenience wrapper that runs this entire chain in one call, while still returning the intermediate `thresholds` and `unfiltered` tables so nothing is hidden.

Every step above is also available and testable on its own — useful for inspecting an intermediate result, re-plotting a distribution, or re-running just one pipeline after a change.

## Quick start

``` r
library(MSiteR)

# Run the full first stage in one call
result <- prepare_filtered_cpg_table(
  pipelines = list(
    Bismark = "Bismark.chr.txt",
    Bwameth = "bwameth.chr.txt",
    Encode  = "ENCODE.chr.txt",
    Biscuit = "biscuit.chr.txt"
  ),
  strand_reference = "sheep_cpg_all.txt",
  strand_reference_for = "Biscuit",
  chunk_by_chromosome = TRUE
)

result$thresholds   # per-pipeline S1/S2 read-depth cutoffs used
result$unfiltered    # long-format table, pre-filtering (for QC/plotting)
result$merged        # final wide table: chr, pos, and each pipeline's
                      # <name>_TR / <name>_MR / <name>_ML / <name>_StrN
```

## Example

Here's a complete example using synthetic methylation data:

``` r
library(MSiteR)

# Create sample data for two pipelines, two chromosomes
bismark <- data.table::data.table(
  chr = c("chr1", "chr1", "chr1", "chr2", "chr2"),
  pos = c(100, 101, 102, 1000, 1001),
  TRead = c(20, 22, 18, 25, 30),
  MRead = c(10, 11, 9, 12, 15),
  ML = c(0.50, 0.50, 0.50, 0.48, 0.50),
  strand = c("+", "-", "+", "+", "-")
)

bwameth <- data.table::data.table(
  chr = c("chr1", "chr1", "chr2", "chr2"),
  pos = c(100, 101, 1000, 1001),
  TRead = c(18, 20, 28, 32),
  MRead = c(9, 10, 14, 16),
  ML = c(0.50, 0.50, 0.50, 0.50),
  strand = c("+", "-", "+", "-")
)

# Run the full pipeline
result <- prepare_filtered_cpg_table(
  pipelines = list(Bismark = bismark, Bwameth = bwameth),
  chunk_by_chromosome = FALSE
)

# View strand-specific thresholds per pipeline
result$thresholds

# View the merged, filtered output
head(result$merged)
#       chr  pos Bismark_TR Bismark_MR Bismark_ML Bwameth_TR Bwameth_MR Bwameth_ML
# 1:  chr1  100         20         10       0.50         18          9       0.50
# 2:  chr1  101         22         11       0.50         20         10       0.50
# 3:  chr2 1000         25         12       0.48         28         14       0.50
# 4:  chr2 1001         30         15       0.50         32         16       0.50
```

### Or call each step manually, for full control

#### `chunk_by_chromosome` == TRUE:

``` r
# 1. Divide Chromosomes into chunks for separate pipeline calls
chunks <- chunk_by_chromosome(pipelines = list(
    Bismark = "Bismark.chr.txt",
    Bwameth = "bwameth.chr.txt",
    Encode  = "ENCODE.chr.txt",
    Biscuit = "biscuit.chr.txt"
  )
)

# 2. Iterate through chromosomes and call [run_pipeline_on_chunk()]
all_merged <- list()      
all_unfiltered <- list()

for (chr in names(chunks)) {
  result <- run_pipeline_on_chunk(chunks[[chr]])
  all_merged[[chr]] <- result$merged
  all_unfiltered[[chr]] <- result$unfiltered

  rm(result)
  gc(verbose = FALSE)
}

# 3. Combine and recalculate global thresholds
unfiltered <- data.table::rbindlist(all_unfiltered)
merged <- data.table::rbindlist(all_merged)
thresholds <- calculate_strand_thresholds(unfiltered)

# 4. Filter each pipeline using its own thresholds
filtered <- Map(function(dt, name) filter_by_strand(dt, name, thresholds),
                 collapsed, names(collapsed))

# 5. Merge into one wide table
merged <- merge_filtered_pipelines(filtered)
```

#### `chunk_by_chromosome` == FALSE:

``` r
# 1. Collapse +/- strand pairs, per pipeline
collapsed <- list(
  Bismark = collapse_cpg_strand("Bismark.chr.txt"),
  Bwameth = collapse_cpg_strand("bwameth.chr.txt"),
  Encode  = collapse_cpg_strand("ENCODE.chr.txt"),
  Biscuit = collapse_cpg_strand("biscuit.chr.txt", strand_reference = "sheep_cpg_all.txt")
)

# 2. Reshape to long format and compute thresholds
unfiltered <- prepare_unfiltered_dt(collapsed)
thresholds <- calculate_strand_thresholds(unfiltered)

# 3. Filter each pipeline using its own thresholds
filtered <- Map(function(dt, name) filter_by_strand(dt, name, thresholds),
                 collapsed, names(collapsed))

# 4. Merge into one wide table
merged <- merge_filtered_pipelines(filtered)
```

## Function reference

| Function | Purpose |
|------------------------------------|------------------------------------|
| `chunk_by_chromosome()` | Subset each pipeline's data.table to individual chromosomes. |
| `collapse_cpg_strand()` | Collapses complementary +/- strand CpG calls for one pipeline into symmetric per-site totals. Joins in strand from a reference file if the input lacks a `strand` column (e.g. Biscuit). |
| `prepare_unfiltered_dt()` | Reshapes a named list of collapsed pipeline tables into one long-format table (`reads`, `Pipeline`, `Strand`) for thresholding and plotting. |
| `calculate_strand_thresholds()` | Computes per-pipeline single-strand (S1) and double-strand (S2) read-depth cutoffs from the long-format table. |
| `filter_by_strand()` | Filters one pipeline's collapsed table using its strand-based thresholds. |
| `run_pipeline_on_chunk()` | Wrapper function that combines [collapse_cpg_strand()], [prepare_unfiltered_dt()], and [filter_by_strand()] to call data pipelines for a selected subset of data. |
| `merge_filtered_pipelines()` | Joins strand-filtered tables from multiple pipelines by `chr`/`pos` into one wide table, with pipeline-prefixed column names. |
| `prepare_filtered_cpg_table()` | Convenience wrapper chaining all of the above; returns `merged`, `thresholds`, and `unfiltered`. |

Full argument/return documentation is available via `?function_name` once the package is loaded, e.g. `?collapse_cpg_strand`.

## Input format

Raw pipeline files are expected as tab-delimited, no header, with columns:

```         
chr  pos  TRead  MRead  ML  [strand]
```

`strand` is optional in the file itself — if absent, pass a `strand_reference` (a `chr`, `CpG_pos`, `strand` table) so `collapse_cpg_strand()` can join it in.

### Status

Stage 1 (strand collapsing → strand filtering → cross-pipeline merge) is implemented and covered by unit tests (100+ passing as of this writing).

Stage 2 (cross-pipeline consensus statistics, tiered consensus filtering, and final strand-symmetric merge) and plotting functions are planned next.

## Author

Jace Ellis ([ellisjacem\@gmail.com](mailto:ellisjacem@gmail.com)), data science student, University of Missouri–Columbia. Developed for CpG methylation pipeline comparison work under the guidance of Shangqian Xie (<https://scholar.google.com/citations?user=HZ8VFAsAAAAJ&hl=zh-CN>).
