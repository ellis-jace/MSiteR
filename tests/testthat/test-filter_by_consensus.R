test_that("filter_by_consensus keeps n=1 sites meeting the 1-method threshold", {
  dt <- data.table::data.table(
    chr = "chr1", pos = c(1, 2), n_methods = c(1, 1),
    Bismark_TR = c(2, 20)
  )
  thresholds <- data.table::data.table(
    Tool = "Bismark", Threshold_1M = 10, Threshold_2M = 5, Threshold_3M = 4, Threshold_4M = 3
  )
  result <- filter_by_consensus(dt, thresholds)
  expect_equal(result$pos, 2)
})

test_that("filter_by_consensus keeps a site if ANY pipeline meets its tier threshold", {
  dt <- data.table::data.table(
    chr = "chr1", pos = 1, n_methods = 2,
    Bismark_TR = 2, Bwameth_TR = 20   # Bwameth alone should be enough
  )
  thresholds <- data.table::data.table(
    Tool = c("Bismark", "Bwameth"),
    Threshold_1M = c(10, 10), Threshold_2M = c(15, 15),
    Threshold_3M = c(4, 4), Threshold_4M = c(3, 3)
  )
  result <- filter_by_consensus(dt, thresholds)
  expect_equal(nrow(result), 1)
})

test_that("filter_by_consensus drops a site if no pipeline meets its tier threshold", {
  dt <- data.table::data.table(
    chr = "chr1", pos = 1, n_methods = 2,
    Bismark_TR = 2, Bwameth_TR = 3
  )
  thresholds <- data.table::data.table(
    Tool = c("Bismark", "Bwameth"),
    Threshold_1M = c(10, 10), Threshold_2M = c(15, 15),
    Threshold_3M = c(4, 4), Threshold_4M = c(3, 3)
  )
  result <- filter_by_consensus(dt, thresholds)
  expect_equal(nrow(result), 0)
})

test_that("filter_by_consensus preserves rows across all four tiers", {
  dt <- data.table::data.table(
    chr = "chr1", pos = 1:4, n_methods = 1:4,
    Bismark_TR = c(20, 20, 20, 20)
  )
  thresholds <- data.table::data.table(
    Tool = "Bismark", Threshold_1M = 5, Threshold_2M = 5, Threshold_3M = 5, Threshold_4M = 5
  )
  result <- filter_by_consensus(dt, thresholds)
  expect_equal(nrow(result), 4)
})
