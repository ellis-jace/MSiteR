test_that("chunking uses the same GLOBAL cutoffs as a full-genome run", {
  # Each chromosome has a different S1 depth distribution. The original
  # implementation retained chr1 in chunked mode but not in full mode.
  x <- data.table::data.table(
    chr = rep(c("chr1", "chr2"), each = 5),
    pos = rep(1:5, 2),
    TRead = c(rep(5L, 5), rep(12L, 5)),
    MRead = c(rep(2L, 5), rep(6L, 5)),
    ML = c(rep(0.4, 5), rep(0.5, 5)),
    strand = "+"
  )
  full <- run_strand_filtering(list(Bismark = x), chunk_by_chromosome = FALSE)
  chunks <- run_strand_filtering(list(Bismark = x), chunk_by_chromosome = TRUE)

  sorted <- function(x) data.table::setorder(data.table::copy(x), chr, pos)
  expect_equal(sorted(full$merged), sorted(chunks$merged))
  expect_equal(full$thresholds, chunks$thresholds)
  expect_true(all(full$merged$chr == "chr2"))
})

test_that("duplicate same-strand records still count as S1", {
  x <- data.table::data.table(
    chr = "chr1", pos = c(100L, 100L), TRead = c(5L, 7L),
    MRead = c(2L, 3L), ML = c(0.4, 3/7), strand = c("+", "+")
  )
  result <- collapse_cpg_strand(x)
  expect_equal(result$n, 1L)
  expect_equal(result$TRead, 12)
})

test_that("a pipeline with only S2 calls still gets thresholds", {
  x <- data.table::data.table(
    Pipeline = "Bismark", Strand = "S2", reads = c(6L, 7L)
  )
  thresholds <- calculate_strand_thresholds(x)
  expect_equal(thresholds$cutoff1, 5)
  expect_equal(thresholds$cutoff2, 3)
})

test_that("both strand-frequency functions accept real wrapper output", {
  x <- data.table::data.table(
    chr = "chr1", pos = 1:5, TRead = 5:9,
    MRead = 2:6, ML = (2:6) / (5:9), strand = "+"
  )
  result <- run_strand_filtering(list(Bismark = x))
  expect_s3_class(plot_strand_frequency_filtered(result), "ggplot")
  expect_s3_class(plot_strand_frequency_filtered(result, show_thresholds = TRUE), "ggplot")
  expect_s3_class(plot_strand_frequency_unfiltered(result, show_thresholds = TRUE), "ggplot")
})

test_that("duplicate pipeline merge coordinates are rejected", {
  x <- data.table::data.table(
    chr = "chr1", pos = c(100L, 100L), TRead = c(5L, 7L),
    MRead = c(2L, 3L), ML = c(0.4, 3/7), n = c(1L, 1L)
  )
  expect_error(merge_filtered_pipelines(list(Bismark = x)), "duplicate chr/pos")
})

test_that("input normalization does not mutate an already-loaded table", {
  x <- data.table::data.table(
    chr = "chr1", pos = 100L, TRead = 9L,
    MRead = 5L, ml = 5/9, strand = "+"
  )
  original <- data.table::copy(x)
  result <- read_pipeline_files(list(Bismark = x))
  expect_true("ML" %in% names(result$Bismark))
  expect_equal(x, original)
})

test_that("a missing reference match is reported", {
  x <- data.table::data.table(
    chr = "chr1", pos = c(100L, 999L), TRead = c(9L, 9L),
    MRead = c(5L, 5L), ML = c(5/9, 5/9)
  )
  ref <- data.table::data.table(chr = "chr1", CpG_pos = 100L, strand = "+")
  expect_warning(collapse_cpg_strand(x, ref), "no strand-reference match")
})

test_that("zero-depth noncalls do not lower 1M percentile thresholds", {
  x <- data.table::data.table(
    n_methods = c(1L, 1L, 1L, 1L),
    A_TR = c(20, 21, 22, 0),
    B_TR = c(0, 0, 0, 8)
  )
  result <- calculate_consensus_thresholds(x)
  expect_true(result[Tool == "A", Threshold_1M] > 20)
  expect_equal(result[Tool == "B", Threshold_1M], 8)
})
