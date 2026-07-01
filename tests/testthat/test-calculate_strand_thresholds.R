test_that("calculate_strand_thresholds computes cutoff1 as clamped 80th percentile of S1 reads", {
  toy <- data.table::data.table(
    Pipeline = rep("Bismark", 6),
    Strand   = c("S1", "S1", "S1", "S2", "S2", "S2"),
    reads    = c(3, 8, 20, 5, 10, 15)
  )
  result <- calculate_strand_thresholds(toy)

  expected_q80 <- quantile(c(3, 8, 20), 0.8)
  expected_cutoff1 <- min(12, max(5, expected_q80))

  expect_equal(result$cutoff1, as.numeric(expected_cutoff1))
})

test_that("cutoff2 is always half of cutoff1, rounded up", {
  toy <- data.table::data.table(
    Pipeline = rep("Bismark", 4),
    Strand   = c("S1", "S1", "S2", "S2"),
    reads    = c(1, 2, 5, 6)
  )
  result <- calculate_strand_thresholds(toy)
  expect_equal(result$cutoff2, ceiling(result$cutoff1 / 2))
})

test_that("cutoff1 is floored at 5 even when the 80th percentile is lower", {
  toy <- data.table::data.table(
    Pipeline = rep("Bismark", 3),
    Strand   = c("S1", "S1", "S1"),
    reads    = c(1, 2, 3)   # 80th percentile well below 5
  )
  result <- calculate_strand_thresholds(toy)
  expect_equal(result$cutoff1, 5)
})

test_that("cutoff1 is capped at 12 even when the 80th percentile is higher", {
  toy <- data.table::data.table(
    Pipeline = rep("Bismark", 3),
    Strand   = c("S1", "S1", "S1"),
    reads    = c(50, 60, 70)   # 80th percentile well above 12
  )
  result <- calculate_strand_thresholds(toy)
  expect_equal(result$cutoff1, 12)
})

test_that("thresholds are computed independently per pipeline", {
  toy <- data.table::data.table(
    Pipeline = c("Bismark", "Bismark", "Bwameth", "Bwameth"),
    Strand   = c("S1", "S1", "S1", "S1"),
    reads    = c(1, 2, 50, 60)
  )
  result <- calculate_strand_thresholds(toy)

  expect_equal(nrow(result), 2)
  expect_equal(result[Pipeline == "Bismark"]$cutoff1, 5)
  expect_equal(result[Pipeline == "Bwameth"]$cutoff1, 12)
})

test_that("results are ordered alphabetically by Pipeline", {
  toy <- data.table::data.table(
    Pipeline = c("Encode", "Bismark"),
    Strand   = c("S1", "S1"),
    reads    = c(10, 10)
  )
  result <- calculate_strand_thresholds(toy)
  expect_equal(result$Pipeline, c("Bismark", "Encode"))
})
