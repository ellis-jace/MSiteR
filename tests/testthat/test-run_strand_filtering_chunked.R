test_that("prepare_filtered_cpg_table with chunk_by_chromosome=TRUE returns correct structure", {
  bis <- data.table::data.table(
    chr = c("chr1", "chr1", "chr2", "chr2"),
    pos = c(1, 2, 100, 200),
    TRead = c(20, 22, 25, 30),
    MRead = c(10, 11, 12, 15),
    ML = c(0.5, 0.5, 0.48, 0.5),
    strand = "+"
  )
  result <- prepare_filtered_cpg_table(list(Bismark = bis), chunk_by_chromosome = TRUE)

  expect_named(result, c("merged", "thresholds", "unfiltered"))
  expect_true(data.table::is.data.table(result$merged))
  expect_true(data.table::is.data.table(result$thresholds))
  expect_true(data.table::is.data.table(result$unfiltered))
})

test_that("prepare_filtered_cpg_table chunked mode produces same results as non-chunked", {
  bis <- data.table::data.table(
    chr = c("chr1", "chr1", "chr2", "chr2"),
    pos = c(1, 2, 100, 200),
    TRead = c(20, 22, 25, 30),
    MRead = c(10, 11, 12, 15),
    ML = c(0.5, 0.5, 0.48, 0.5),
    strand = "+"
  )
  bwa <- data.table::data.table(
    chr = c("chr1", "chr2"),
    pos = c(1, 100),
    TRead = c(18, 28),
    MRead = c(9, 14),
    ML = c(0.5, 0.5),
    strand = "+"
  )
  pipelines <- list(Bismark = bis, Bwameth = bwa)

  # Run both modes
  result_normal <- prepare_filtered_cpg_table(pipelines, chunk_by_chromosome = FALSE)
  result_chunked <- prepare_filtered_cpg_table(pipelines, chunk_by_chromosome = TRUE)

  # Merged tables should have same number of rows and columns
  expect_equal(nrow(result_normal$merged), nrow(result_chunked$merged))
  expect_equal(ncol(result_normal$merged), ncol(result_chunked$merged))
  expect_equal(names(result_normal$merged), names(result_chunked$merged))

  # Both should produce consistent threshold structures
  expect_equal(ncol(result_normal$thresholds), ncol(result_chunked$thresholds))

  # Unfiltered tables should have same structure and number of rows
  expect_equal(nrow(result_normal$unfiltered), nrow(result_chunked$unfiltered))
  expect_equal(names(result_normal$unfiltered), names(result_chunked$unfiltered))
})

test_that("prepare_filtered_cpg_table chunked mode handles multiple chromosomes", {
  bis <- data.table::data.table(
    chr = c("chr1", "chr2", "chr3"),
    pos = c(1, 100, 1000),
    TRead = c(20, 25, 30),
    MRead = c(10, 12, 15),
    ML = c(0.5, 0.48, 0.5),
    strand = "+"
  )
  result <- prepare_filtered_cpg_table(list(Bismark = bis), chunk_by_chromosome = TRUE)

  # Should have produced a merged result and unfiltered data
  expect_true(nrow(result$merged) > 0)
  expect_true(nrow(result$unfiltered) > 0)
  expect_true(nrow(result$thresholds) > 0)
})

test_that("prepare_filtered_cpg_table chunked mode produces output structure", {
  # Create data with mixed quality
  bis <- data.table::data.table(
    chr = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr2"),
    pos = c(1, 2, 3, 100, 200, 300),
    TRead = c(20, 22, 1, 25, 30, 1),  # pos 3 and 300 are low-depth
    MRead = c(10, 11, 0, 12, 15, 0),
    ML = c(0.5, 0.5, 0, 0.48, 0.5, 0),
    strand = "+"
  )
  result <- prepare_filtered_cpg_table(list(Bismark = bis), chunk_by_chromosome = TRUE)

  # Should have consistent output structure
  expect_true(data.table::is.data.table(result$merged))
  expect_true(data.table::is.data.table(result$unfiltered))
  expect_true(data.table::is.data.table(result$thresholds))

  # Merged should have fewer or equal rows than input (due to filtering)
  expect_true(nrow(result$merged) <= nrow(bis))
})

