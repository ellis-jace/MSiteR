# MSiteR
 
Tools for filtering, comparing, and merging CpG methylation calls across multiple alignment/calling pipelines (Bismark, BWA-meth, Biscuit, ENCODE), including strand-level read-depth thresholding, cross-pipeline consensus analysis, and comprehensive visualization.
 
This package grew out of a manual, script-based workflow for comparing CpG methylation calls across four pipelines on sheep/cattle genomic references. It replaces four near-duplicate per-pipeline scripts with a small set of tested, reusable functions organized into two filtering stages.
 
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
 
The package covers **2 stages** of methylation-calling workflow: turning raw, per-strand pipeline output into a final consensus-filtered, cross-pipeline table.
 
### Stage 1: Strand filtering
 
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
 
`run_strand_filtering()` is a convenience wrapper that runs this entire chain in one call, while still returning the intermediate `thresholds` and `unfiltered` tables so nothing is hidden.
 
### Stage 2: Consensus filtering
 
```
strand-filtered wide table
        │
        ▼
calculate_consensus_stats()    # add n_methods (pipeline agreement) and median metrics
        │
        ▼
calculate_consensus_thresholds()# per-tool, per-tier (1M/2M/3M/4M) read-depth cutoffs
        │
        ▼
filter_by_consensus()           # apply tiered thresholds based on pipeline agreement
        │
        ▼
final consensus-filtered table
```
 
`run_consensus_filtering()` orchestrates Stage 2, taking output from Stage 1 and returning both filtered and unfiltered versions for QC and visualization.
 
Every step in both stages is independently callable — useful for inspecting an intermediate result, re-plotting a distribution, or re-running just one pipeline after a change.
 
## Quick start
 
``` r
library(MSiteR)
 
# Stage 1: Strand filtering
result_strand <- run_strand_filtering(
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
 
result_strand$thresholds   # per-pipeline S1/S2 read-depth cutoffs used
result_strand$unfiltered   # long-format table, pre-filtering (for QC/plotting)
result_strand$merged       # final wide table: chr, pos, and each pipeline's
                            # <name>_TR / <name>_MR / <name>_ML
 
# Stage 2: Consensus filtering
result_consensus <- run_consensus_filtering(result_strand)
 
result_consensus$merged       # final consensus-filtered table (all pipelines merged)
result_consensus$unfiltered   # same table before consensus filtering (for before/after plots)
result_consensus$thresholds   # per-tool, per-tier (1M/2M/3M/4M) thresholds
result_consensus$filtering_stats  # filtering summary statistics
```
 
## Visualization
 
MSiteR includes plotting functions for both stages to visualize filtering effects.
 
### Strand-level plots
 
All strand plots take `result_strand` from `run_strand_filtering()` and support a `show_thresholds` parameter to toggle threshold visualization:
 
``` r
# Boxplot: S1 vs S2 read depths
plot_strand_boxplot(result_strand, outdir = "figs/", sample_name = "sample1")
 
# Frequency distribution before filtering
plot_strand_frequency_unfiltered(result_strand, outdir = "figs/", 
                                 sample_name = "sample1", show_thresholds = TRUE)
 
# Frequency distribution after filtering
plot_strand_frequency_filtered(result_strand, outdir = "figs/", 
                               sample_name = "sample1", show_thresholds = TRUE)
```
 
### Consensus-level plots
 
All consensus plots take `result_consensus` from `run_consensus_filtering()` and support a `show_filtered` parameter to toggle before/after consensus filtering visualization:
 
``` r
# UpSet plot: pipeline detection overlaps
plot_consensus_upset(result_consensus, outdir = "figs/", sample_name = "sample1")
plot_consensus_upset(result_consensus, outdir = "figs/", sample_name = "sample1", 
                     show_filtered = TRUE)
 
# Venn diagram: 4-way pipeline overlaps
plot_consensus_venn(result_consensus, outdir = "figs/", sample_name = "sample1")
 
# Distribution: read depths by consensus tier (1M, 2M, 3M, 4M)
plot_consensus_distribution(result_consensus, outdir = "figs/", sample_name = "sample1")
 
# Boxplot: read depths by consensus tier
plot_consensus_boxplot(result_consensus, outdir = "figs/", sample_name = "sample1")
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
 
# Stage 1: Strand filtering
result_strand <- run_strand_filtering(
  pipelines = list(Bismark = bismark, Bwameth = bwameth),
  chunk_by_chromosome = FALSE
)
 
# View strand-specific thresholds per pipeline
result_strand$thresholds
 
# View the strand-filtered merged output
head(result_strand$merged)
#       chr  pos Bismark_TR Bismark_MR Bismark_ML Bwameth_TR Bwameth_MR Bwameth_ML
# 1:  chr1  100         20         10       0.50         18          9       0.50
# 2:  chr1  101         22         11       0.50         20         10       0.50
# 3:  chr2 1000         25         12       0.48         28         14       0.50
# 4:  chr2 1001         30         15       0.50         32         16       0.50
 
# Stage 2: Consensus filtering
result_consensus <- run_consensus_filtering(result_strand)
 
# View consensus-filtered output with pipeline agreement counts
head(result_consensus$merged)
#       chr  pos Bismark_TR ... n_methods median_Tread median_MethyL
# 1:  chr1  100         20 ...         2           19             0.50
# 2:  chr1  101         22 ...         2           21             0.50
# 3:  chr2 1000         25 ...         2           27             0.49
# 4:  chr2 1001         30 ...         2           31             0.50
 
# View consensus-tier thresholds
result_consensus$thresholds
```
 
## Function reference
 
### Stage 1: Strand filtering
 
| Function | Purpose |
|----------|---------|
| `chunk_by_chromosome()` | Subset each pipeline's data.table to individual chromosomes. |
| `collapse_cpg_strand()` | Collapse complementary +/- strand CpG calls into symmetric per-site totals. Joins strand from a reference file if absent (e.g., Biscuit). |
| `prepare_unfiltered_dt()` | Reshape collapsed pipeline tables into long-format (`reads`, `Pipeline`, `Strand`) for thresholding and plotting. |
| `calculate_strand_thresholds()` | Compute per-pipeline single-strand (S1) and double-strand (S2) read-depth cutoffs (80th percentile, clamped 5–12). |
| `filter_by_strand()` | Filter one pipeline's collapsed table using strand-based thresholds. |
| `run_pipeline_on_chunk()` | Wrapper combining collapse, reshape, and filtering for a data subset. |
| `merge_filtered_pipelines()` | Join strand-filtered tables by `chr`/`pos` into one wide table. |
| `run_strand_filtering()` | **Stage 1 wrapper:** orchestrates all strand-filtering steps; returns `merged`, `thresholds`, and `unfiltered`. |
 
### Stage 2: Consensus filtering
 
| Function | Purpose |
|----------|---------|
| `calculate_consensus_stats()` | Add consensus metrics: `n_methods` (pipeline agreement count), median read depths/methylation. |
| `calculate_consensus_thresholds()` | Compute per-tool, per-tier read-depth cutoffs (95th percentile for 1M, scaled for 2M–4M). |
| `filter_by_consensus()` | Filter sites based on tiered thresholds matching their consensus tier. |
| `run_consensus_filtering()` | **Stage 2 wrapper:** adds consensus metrics, calculates thresholds, filters; returns `merged`, `unfiltered`, `thresholds`, and `filtering_stats`. |
 
### Strand-level plotting
 
| Function | Purpose |
|----------|---------|
| `plot_strand_boxplot()` | S1 vs. S2 read-depth boxplot with Wilcoxon test statistics. |
| `plot_strand_frequency_unfiltered()` | Frequency distribution of read depths before strand filtering. |
| `plot_strand_frequency_filtered()` | Frequency distribution of read depths after strand filtering. |
 
All support `show_thresholds = TRUE` to overlay threshold lines and labels.
 
### Consensus-level plotting
 
| Function | Purpose |
|----------|---------|
| `plot_consensus_upset()` | UpSet plot showing pipeline detection overlaps. |
| `plot_consensus_venn()` | Venn diagram of 4-way pipeline overlaps. |
| `plot_consensus_distribution()` | Area plot: read-depth distribution by consensus tier (1M, 2M, 3M, 4M). |
| `plot_consensus_boxplot()` | Boxplot: read depths grouped by consensus tier. |
 
All support `show_filtered = TRUE` to visualize before/after consensus filtering.
 
Full argument/return documentation is available via `?function_name` once the package is loaded, e.g. `?run_strand_filtering` or `?plot_consensus_upset`.
 
## Input format
 
Raw pipeline files are expected as tab-delimited, no header, with columns:
 
```
chr  pos  TRead  MRead  ML  [strand]
```
 
`strand` is optional — if absent, pass a `strand_reference` (a `chr`, `CpG_pos`, `strand` table) to `run_strand_filtering()` so it can join strand info in.
 
## Status
 
**Stages 1 and 2 are production-ready:**
- Strand filtering (Stage 1): strand collapsing → per-pipeline thresholding → cross-pipeline merge. Covered by 100+ unit tests.
- Consensus filtering (Stage 2): pipeline-agreement metrics → tiered thresholding → final filtering. Covered by 40+ unit tests.
- Plotting functions (both stages): 7 strand/consensus visualization functions with 70+ tests.
**Total test coverage:** 150+ passing unit tests across all functions and edge cases.
 
## Author
 
Jace Ellis ([ellisjacem\@gmail.com](mailto:ellisjacem@gmail.com)), data science student, University of Missouri–Columbia. Developed for CpG methylation pipeline comparison work under the guidance of Shangqian Xie (<https://scholar.google.com/citations?user=HZ8VFAsAAAAJ&hl=zh-CN>).
 
