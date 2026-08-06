test_that("chunk_by_chromosome partitions data correctly by chromosome", {
  bis <- data.table::data.table(
    chr = c("chr1", "chr1", "chr2", "chr2"),
    pos = c(1, 2, 100, 200),
    TRead = c(10, 15, 20, 25),
    MRead = c(5, 7, 10, 12),
    ML = c(0.5, 0.47, 0.5, 0.48),
    strand = "+"
  )
  bwa <- data.table::data.table(
    chr = c("chr1", "chr2"),
    pos = c(1, 100),
    TRead = c(12, 22),
    MRead = c(6, 11),
    ML = c(0.5, 0.5),
    strand = "+"
  )
  chunks <- chunk_by_chromosome(list(Bismark = bis, Bwameth = bwa))

  # Should have two chromosomes
  expect_equal(length(chunks), 2)
  expect_named(chunks, c("chr1", "chr2"))

  # chr1 chunk should have partitioned data
  expect_true(all(chunks$chr1$Bismark$chr == "chr1"))
  expect_true(all(chunks$chr1$Bwameth$chr == "chr1"))

  # chr2 chunk should have partitioned data
  expect_true(all(chunks$chr2$Bismark$chr == "chr2"))
  expect_true(all(chunks$chr2$Bwameth$chr == "chr2"))

  # Each chromosome should only contain its own rows
  expect_equal(nrow(chunks$chr1$Bismark) + nrow(chunks$chr2$Bismark), nrow(bis))
  expect_equal(nrow(chunks$chr1$Bwameth) + nrow(chunks$chr2$Bwameth), nrow(bwa))
})

test_that("chunk_by_chromosome preserves pipeline structure", {
  bis <- data.table::data.table(
    chr = c("chr1", "chr2"), pos = c(1, 100),
    TRead = c(10, 20), MRead = c(5, 10), ML = c(0.5, 0.5), strand = "+"
  )
  bwa <- data.table::data.table(
    chr = c("chr1", "chr2"), pos = c(1, 100),
    TRead = c(12, 22), MRead = c(6, 11), ML = c(0.5, 0.5), strand = "+"
  )
  chunks <- chunk_by_chromosome(list(Bismark = bis, Bwameth = bwa))

  # Each chromosome chunk should have the same pipeline names
  expect_named(chunks$chr1, c("Bismark", "Bwameth"))
  expect_named(chunks$chr2, c("Bismark", "Bwameth"))
})

test_that("chunk_by_chromosome sorts chromosomes alphabetically", {
  bis <- data.table::data.table(
    chr = c("chr10", "chr2", "chr1"),
    pos = c(1, 2, 3),
    TRead = c(10, 20, 30),
    MRead = c(5, 10, 15),
    ML = c(0.5, 0.5, 0.5),
    strand = "+"
  )
  chunks <- chunk_by_chromosome(list(Bismark = bis))

  # Chromosomes should be in sorted order (chr1, chr10, chr2)
  expect_equal(names(chunks), c("chr1", "chr10", "chr2"))
})

