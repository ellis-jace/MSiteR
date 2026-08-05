test_that("calculate_consensus_stats computes n_methods correctly", {
  dt <- data.table::data.table(
    chr = "chr1", pos = c(1, 2),
    Bismark_TR = c(10, NA), Bismark_ML = c(0.5, NA),
    Bwameth_TR = c(8, 12),  Bwameth_ML = c(0.5, 0.4)
  )
  result <- calculate_consensus_stats(dt)

  expect_equal(result[pos == 1, n_methods], 2)
  expect_equal(result[pos == 2, n_methods], 1)
})

test_that("calculate_consensus_stats computes median read depth and methylation", {
  dt <- data.table::data.table(
    chr = "chr1", pos = 1,
    Bismark_TR = 10, Bismark_ML = 0.5,
    Bwameth_TR = 20, Bwameth_ML = 0.5
  )
  result <- calculate_consensus_stats(dt)

  expect_equal(result$median_Tread, 15)
  expect_equal(result$median_MethyL, 0.5)
  expect_equal(result$median_Mread, round(15 * 0.5))
})

test_that("calculate_consensus_stats errors with no _TR columns", {
  dt <- data.table::data.table(chr = "chr1", pos = 1, foo = 5)
  expect_error(calculate_consensus_stats(dt), "_TR")
})
