test_that("merge_filtered_pipelines joins by chr/pos and renames columns per pipeline", {
  bis <- data.table::data.table(chr = "chr1", pos = c(1, 2), TRead = c(10, 20),
                                MRead = c(5, 10), ML = c(0.5, 0.5), n = c(1, 2))
  bwa <- data.table::data.table(chr = "chr1", pos = c(1, 3), TRead = c(8, 15),
                                MRead = c(4, 7), ML = c(0.5, 0.47), n = c(1, 1))

  result <- merge_filtered_pipelines(list(Bismark = bis, Bwameth = bwa))

  expect_setequal(names(result),
                  c("chr", "pos", "Bismark_TR", "Bismark_MR", "Bismark_ML", "Bismark_StrN",
                    "Bwameth_TR", "Bwameth_MR", "Bwameth_ML", "Bwameth_StrN"))
})

test_that("merge_filtered_pipelines performs a full outer join (keeps unmatched sites)", {
  bis <- data.table::data.table(chr = "chr1", pos = c(1, 2), TRead = c(10, 20),
                                MRead = c(5, 10), ML = c(0.5, 0.5), n = c(1, 2))
  bwa <- data.table::data.table(chr = "chr1", pos = c(1, 3), TRead = c(8, 15),
                                MRead = c(4, 7), ML = c(0.5, 0.47), n = c(1, 1))

  result <- merge_filtered_pipelines(list(Bismark = bis, Bwameth = bwa))

  # pos 1: in both; pos 2: only Bismark; pos 3: only Bwameth -> 3 rows total
  expect_equal(nrow(result), 3)
  expect_true(is.na(result[pos == 2]$Bwameth_TR))
  expect_true(is.na(result[pos == 3]$Bismark_TR))
})

test_that("merge_filtered_pipelines works with more than two pipelines", {
  a <- data.table::data.table(chr = "chr1", pos = 1, TRead = 10, MRead = 5, ML = 0.5, n = 1)
  b <- data.table::data.table(chr = "chr1", pos = 1, TRead = 12, MRead = 6, ML = 0.5, n = 1)
  c <- data.table::data.table(chr = "chr1", pos = 1, TRead = 9,  MRead = 4, ML = 0.44, n = 2)

  result <- merge_filtered_pipelines(list(Bismark = a, Bwameth = b, Encode = c))

  expect_equal(nrow(result), 1)
  expect_true(all(c("Bismark_TR", "Bwameth_TR", "Encode_TR") %in% names(result)))
})

test_that("merge_filtered_pipelines does not mutate the caller's original tables", {
  bis <- data.table::data.table(chr = "chr1", pos = 1, TRead = 10, MRead = 5, ML = 0.5, n = 1)
  original_names <- names(bis)

  merge_filtered_pipelines(list(Bismark = bis))

  expect_equal(names(bis), original_names)  # TRead should still be TRead, not Bismark_TR
})
