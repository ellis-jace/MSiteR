test_that("filter_by_strand keeps S1 sites meeting cutoff1", {
  dt <- data.table::data.table(chr = "chr1", pos = c(1, 2), TRead = c(4, 10), n = c(1, 1))
  thresholds <- data.table::data.table(Pipeline = "Bismark", cutoff1 = 5, cutoff2 = 3)

  result <- filter_by_strand(dt, "Bismark", thresholds)

  expect_equal(nrow(result), 1)
  expect_equal(result$pos, 2)   # TRead 10 >= cutoff1 5; TRead 4 < 5, dropped
})

test_that("filter_by_strand keeps S2 sites meeting cutoff2", {
  dt <- data.table::data.table(chr = "chr1", pos = c(1, 2), TRead = c(2, 8), n = c(2, 2))
  thresholds <- data.table::data.table(Pipeline = "Bismark", cutoff1 = 5, cutoff2 = 3)

  result <- filter_by_strand(dt, "Bismark", thresholds)

  expect_equal(nrow(result), 1)
  expect_equal(result$pos, 2)   # TRead 8 >= cutoff2 3; TRead 2 < 3, dropped
})

test_that("filter_by_strand applies the correct cutoff per strand type in mixed data", {
  dt <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3, 4),
    TRead = c(4, 6, 2, 4),
    n = c(1, 1, 2, 2)
  )
  thresholds <- data.table::data.table(Pipeline = "Bismark", cutoff1 = 5, cutoff2 = 3)

  result <- filter_by_strand(dt, "Bismark", thresholds)

  # pos 1: S1, TRead 4 < 5 -> dropped
  # pos 2: S1, TRead 6 >= 5 -> kept
  # pos 3: S2, TRead 2 < 3 -> dropped
  # pos 4: S2, TRead 4 >= 3 -> kept
  expect_setequal(result$pos, c(2, 4))
})

test_that("filter_by_strand looks up thresholds for the correct pipeline only", {
  dt <- data.table::data.table(chr = "chr1", pos = 1, TRead = 6, n = 1)
  thresholds <- data.table::data.table(
    Pipeline = c("Bismark", "Bwameth"),
    cutoff1  = c(10, 3),
    cutoff2  = c(5, 2)
  )
  # Using Bwameth's (lower) threshold should keep this site
  result <- filter_by_strand(dt, "Bwameth", thresholds)
  expect_equal(nrow(result), 1)

  # Using Bismark's (higher) threshold should drop it
  result2 <- filter_by_strand(dt, "Bismark", thresholds)
  expect_equal(nrow(result2), 0)
})

test_that("filter_by_strand respects a custom reads_col", {
  dt <- data.table::data.table(chr = "chr1", pos = 1, CustomReads = 10, n = 1)
  thresholds <- data.table::data.table(Pipeline = "Bismark", cutoff1 = 5, cutoff2 = 3)

  result <- filter_by_strand(dt, "Bismark", thresholds, reads_col = "CustomReads")
  expect_equal(nrow(result), 1)
})
