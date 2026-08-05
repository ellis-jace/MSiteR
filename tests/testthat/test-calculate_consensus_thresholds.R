test_that("calculate_consensus_thresholds computes tier 1 as 95th pctile floored at 3", {
  dt <- data.table::data.table(
    n_methods = rep(1, 10),
    Bismark_TR = 1:10
  )
  result <- calculate_consensus_thresholds(dt)
  expect_equal(result$Threshold_1M, as.numeric(max(quantile(1:10, 0.95), 3)))
})

test_that("tier 2/3/4 thresholds are derived from tier 1 with correct clamps", {
  dt <- data.table::data.table(n_methods = rep(1, 10), Bismark_TR = rep(40, 10))
  result <- calculate_consensus_thresholds(dt)

  expect_equal(result$Threshold_1M, 40)
  expect_equal(result$Threshold_2M, 10)  # ceiling(40/2)=20, clamped to max 10
  expect_equal(result$Threshold_3M, 8)   # ceiling(40/3)=14, clamped to max 8
  expect_equal(result$Threshold_4M, 6)   # ceiling(40/4)=10, clamped to max 6
})

test_that("falls back to floor of 3 when no n_methods==1 sites exist", {
  dt <- data.table::data.table(n_methods = rep(2, 5), Bismark_TR = 1:5)
  result <- calculate_consensus_thresholds(dt)
  expect_equal(result$Threshold_1M, 3)
})

test_that("computes thresholds independently per pipeline", {
  dt <- data.table::data.table(
    n_methods = rep(1, 5),
    Bismark_TR = c(1, 2, 3, 4, 5),
    Bwameth_TR = c(50, 60, 70, 80, 90)
  )
  result <- calculate_consensus_thresholds(dt)
  expect_true(result[Tool == "Bwameth", Threshold_1M] > result[Tool == "Bismark", Threshold_1M])
})
