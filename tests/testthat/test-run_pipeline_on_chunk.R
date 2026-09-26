# Note: run_pipeline_on_chunk is a helper function nested inside run_strand_filtering
# It is not directly exported and cannot be tested independently.
# Instead, we test its behavior through run_strand_filtering's behavior.

test_that("run_strand_filtering produces merged and unfiltered output", {
  bis <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3),
    TRead = c(15, 20, 18), MRead = c(7, 10, 9), ML = c(0.47, 0.5, 0.5),
    strand = "+"
  )
  result <- run_strand_filtering(list(Bismark = bis), chunk_by_chromosome = FALSE)

  # Wrapper should produce merged and unfiltered (plus thresholds)
  expect_named(result, c("merged", "thresholds", "unfiltered"))
  expect_true(data.table::is.data.table(result$merged))
  expect_true(data.table::is.data.table(result$unfiltered))
  expect_true(data.table::is.data.table(result$thresholds))
})

test_that("run_strand_filtering respects strand_reference", {
  biscuit <- data.table::data.table(
    chr = "chr1", pos = c(50, 100),
    TRead = c(15, 18), MRead = c(7, 9), ML = c(0.47, 0.5)
  )
  ref <- data.table::data.table(
    chr = "chr1", CpG_pos = c(50, 100), strand = c("+", "-")
  )
  result <- run_strand_filtering(
    list(Biscuit = biscuit),
    strand_reference = ref,
    strand_reference_for = "Biscuit",
    chunk_by_chromosome = FALSE
  )

  expect_true("Biscuit_TR" %in% names(result$merged))
})

test_that("run_strand_filtering chunked and non-chunked share same logic", {
  bis <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3, 4),
    TRead = c(20, 22, 25, 30),
    MRead = c(10, 11, 12, 15),
    ML = c(0.5, 0.5, 0.48, 0.5),
    strand = "+"
  )

  # Both modes should produce output with same structure
  normal_result <- run_strand_filtering(list(Bismark = bis), chunk_by_chromosome = FALSE)
  chunked_result <- run_strand_filtering(list(Bismark = bis), chunk_by_chromosome = TRUE)

  # Should have same column names
  expect_equal(names(normal_result$merged), names(chunked_result$merged))
  expect_equal(names(normal_result$unfiltered), names(chunked_result$unfiltered))
})

