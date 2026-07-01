test_that("prepare_filtered_cpg_table runs end-to-end and returns all three pieces", {
  bis <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3, 4, 5),
    TRead = c(3, 15, 8, 20, 6), MRead = c(1, 7, 4, 10, 3), ML = c(0.33, 0.47, 0.5, 0.5, 0.5),
    strand = c("+", "+", "+", "+", "+")   # all unpaired -> all S1 (n = 1)
  )
  bwa <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3),
    TRead = c(4, 16, 9), MRead = c(2, 8, 4), ML = c(0.5, 0.5, 0.44),
    strand = c("+", "+", "+")
  )

  result <- prepare_filtered_cpg_table(list(Bismark = bis, Bwameth = bwa))

  expect_named(result, c("merged", "thresholds", "unfiltered"))
  expect_true(data.table::is.data.table(result$merged))
  expect_true(data.table::is.data.table(result$thresholds))
  expect_true(data.table::is.data.table(result$unfiltered))
})

test_that("prepare_filtered_cpg_table's merged output has expected columns", {
  bis <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3, 4, 5),
    TRead = c(10, 20, 8, 15, 6), MRead = c(5, 10, 4, 7, 3), ML = c(0.5, 0.5, 0.5, 0.47, 0.5),
    strand = "+"
  )
  bwa <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3, 4),
    TRead = c(9, 22, 11, 14), MRead = c(4, 11, 5, 7), ML = c(0.44, 0.5, 0.45, 0.5),
    strand = "+"
  )

  result <- prepare_filtered_cpg_table(list(Bismark = bis, Bwameth = bwa))

  expect_true(all(c("Bismark_TR", "Bwameth_TR") %in% names(result$merged)))
})

test_that("prepare_filtered_cpg_table actually filters out low-depth sites", {
  # 5 sites: 4 with decent depth to establish a real cutoff1, 1 clearly too low to survive
  bis <- data.table::data.table(
    chr = "chr1", pos = c(1, 2, 3, 4, 5),
    TRead = c(20, 22, 25, 30, 1),      # site 5 (pos=5, TRead=1) should fail
    MRead = c(10, 11, 12, 15, 0),
    ML    = c(0.5, 0.5, 0.48, 0.5, 0),
    strand = "+"
  )

  result <- prepare_filtered_cpg_table(list(Bismark = bis))

  expect_false(5 %in% result$merged$pos)
})

test_that("prepare_filtered_cpg_table passes strand_reference through to the correct pipeline", {
  biscuit_raw <- data.table::data.table(
    chr = "chr1", pos = c(50, 100, 101),
    TRead = c(15, 15, 18), MRead = c(7, 7, 9), ML = c(0.47, 0.47, 0.5)
    # pos 50: unpaired -> S1 after collapse (gives a valid cutoff1)
    # pos 100/101: +/- pair -> collapses to one S2 site
  )
  ref <- data.table::data.table(chr = "chr1", CpG_pos = c(50, 100, 101), strand = c("+", "+", "-"))

  result <- prepare_filtered_cpg_table(
    list(Biscuit = biscuit_raw),
    strand_reference = ref,
    strand_reference_for = "Biscuit"
  )

  expect_true("Biscuit_TR" %in% names(result$merged))
  expect_equal(nrow(result$merged), 2)  # one S1 site + one collapsed S2 site
})
