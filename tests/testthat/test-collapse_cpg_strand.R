test_that("collapse_cpg_strand merges a +/- pair into one symmetric site", {
  # + strand at pos 100, - strand at pos 101 -> should collapse to one row
  dt <- data.table::data.table(
    V1 = "chr1", V2 = c(100, 101),
    V3 = c(10, 12), V4 = c(5, 6), V5 = c(0.5, 0.5),
    V6 = c("+", "-")
  )
  data.table::setnames(dt, c("chr", "pos", "TRead", "MRead", "ML", "strand"))

  result <- collapse_cpg_strand(dt)

  expect_equal(nrow(result), 1)
  expect_equal(result$pos, 100)             # keeps the + strand position
  expect_equal(result$TRead, 22)            # 10 + 12
  expect_equal(result$MRead, 11)            # 5 + 6
  expect_equal(result$ML, 11 / 22)          # recalculated, not averaged
  expect_equal(result$n, 2)
})

test_that("collapse_cpg_strand keeps an unpaired site as-is", {
  dt <- data.table::data.table(
    chr = "chr1", pos = 200, TRead = 8, MRead = 4, ML = 0.5, strand = "+"
  )
  result <- collapse_cpg_strand(dt)

  expect_equal(nrow(result), 1)
  expect_equal(result$pos, 200)
  expect_equal(result$TRead, 8)
  expect_equal(result$n, 1)
})

test_that("collapse_cpg_strand joins strand info from a reference when missing", {
  # Biscuit-style input: no strand column
  biscuit_raw <- data.table::data.table(
    chr = "chr1", pos = c(100, 101), TRead = c(10, 12), MRead = c(5, 6), ML = c(0.5, 0.5)
  )
  ref <- data.table::data.table(chr = "chr1", CpG_pos = c(100, 101), strand = c("+", "-"))

  result <- collapse_cpg_strand(biscuit_raw, strand_reference = ref)

  expect_equal(nrow(result), 1)
  expect_equal(result$TRead, 22)
  expect_equal(result$n, 2)
})

test_that("collapse_cpg_strand errors when strand is missing and no reference given", {
  dt <- data.table::data.table(chr = "chr1", pos = 100, TRead = 10, MRead = 5, ML = 0.5)
  expect_error(collapse_cpg_strand(dt), "strand_reference")
})

test_that("collapse_cpg_strand accepts a file path", {
  tmp <- tempfile(fileext = ".txt")
  data.table::fwrite(
    list(rep("chr1", 2), c(100, 101), c(10, 12), c(5, 6), c(0.5, 0.5), c("+", "-")),
    tmp, sep = "\t", col.names = FALSE
  )

  result <- collapse_cpg_strand(tmp)
  expect_equal(nrow(result), 1)
  expect_equal(result$TRead, 22)

  unlink(tmp)
})
