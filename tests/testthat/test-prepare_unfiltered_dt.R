test_that("prepare_unfiltered_dt stacks pipelines into long format with correct labels", {
  bis <- data.table::data.table(chr = "chr1", pos = c(1, 2), TRead = c(3, 20), n = c(1, 2))
  bwa <- data.table::data.table(chr = "chr1", pos = c(1, 2), TRead = c(4, 22), n = c(1, 2))

  result <- prepare_unfiltered_dt(list(Bismark = bis, Bwameth = bwa))

  expect_equal(nrow(result), 4)
  expect_setequal(names(result), c("reads", "Pipeline", "Strand"))
  expect_setequal(unique(result$Pipeline), c("Bismark", "Bwameth"))
})

test_that("prepare_unfiltered_dt correctly derives Strand from n", {
  bis <- data.table::data.table(chr = "chr1", pos = c(1, 2), TRead = c(5, 15), n = c(1, 2))
  result <- prepare_unfiltered_dt(list(Bismark = bis))

  expect_equal(result[reads == 5]$Strand, "S1")
  expect_equal(result[reads == 15]$Strand, "S2")
})

test_that("prepare_unfiltered_dt handles a single pipeline", {
  bis <- data.table::data.table(chr = "chr1", pos = 1, TRead = 10, n = 1)
  result <- prepare_unfiltered_dt(list(Bismark = bis))
  expect_equal(nrow(result), 1)
  expect_equal(result$Pipeline, "Bismark")
})
