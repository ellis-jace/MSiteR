# Helper: write a tab-delimited, headerless methylation file and return its path
write_pipeline_file <- function(df) {
  path <- tempfile(fileext = ".txt")
  data.table::fwrite(df, path, sep = "\t", col.names = FALSE)
  path
}

six_col_df <- data.frame(
  chr = c("chr1", "chr1", "chr2"),
  pos = c(100L, 101L, 1000L),
  TRead = c(20L, 22L, 25L),
  MRead = c(10L, 11L, 12L),
  ML = c(0.5, 0.5, 0.48),
  strand = c("+", "-", "+")
)

five_col_df <- six_col_df[, 1:5]

# ============================================================================
# Input validation
# ============================================================================

test_that("errors when input is not a list", {
  expect_error(
    read_pipeline_files("file.txt"),
    "must be a named list"
  )
})

test_that("errors when list is unnamed", {
  f <- write_pipeline_file(six_col_df)
  expect_error(
    read_pipeline_files(list(f)),
    "must be a named list"
  )
})

test_that("errors when any name is empty", {
  f <- write_pipeline_file(six_col_df)
  expect_error(
    read_pipeline_files(setNames(list(f, f), c("Bismark", ""))),
    "must be a named list"
  )
})

test_that("errors when a file path does not exist", {
  expect_error(
    read_pipeline_files(list(Bismark = "does_not_exist.txt")),
    "File not found for 'Bismark'"
  )
})

test_that("errors when an element is neither a path nor a table", {
  expect_error(
    read_pipeline_files(list(Bismark = 42)),
    "must be a single file path or a data.frame/data.table"
  )
})

test_that("errors when a character element has length > 1", {
  f <- write_pipeline_file(six_col_df)
  expect_error(
    read_pipeline_files(list(Bismark = c(f, f))),
    "must be a single file path or a data.frame/data.table"
  )
})

test_that("errors and names the pipeline when a file has the wrong column count", {
  bad <- data.frame(chr = "chr1", pos = 1L, TRead = 5L)
  f <- write_pipeline_file(bad)
  expect_error(
    read_pipeline_files(list(Biscuit = f)),
    "Input for 'Biscuit' has 3 columns"
  )
})

test_that("errors when an in-memory table has the wrong column count", {
  bad <- data.frame(chr = "chr1", pos = 1L, TRead = 5L)
  expect_error(
    read_pipeline_files(list(Biscuit = bad)),
    "Input for 'Biscuit' has 3 columns"
  )
})

# ============================================================================
# Accepting data.tables already in memory
# ============================================================================

test_that("accepts a data.table already in memory", {
  dt <- data.table::as.data.table(six_col_df)
  result <- read_pipeline_files(list(Bismark = dt))

  expect_s3_class(result$Bismark, "data.table")
  expect_equal(nrow(result$Bismark), 3)
})

test_that("accepts a data.frame and converts it to data.table", {
  result <- read_pipeline_files(list(Bismark = six_col_df))

  expect_s3_class(result$Bismark, "data.table")
  expect_equal(nrow(result$Bismark), 3)
})

test_that("in-memory tables get standard column names", {
  dt <- data.table::as.data.table(five_col_df)
  data.table::setnames(dt, c("V1", "V2", "V3", "V4", "V5"))
  result <- read_pipeline_files(list(Biscuit = dt))

  expect_equal(names(result$Biscuit), c("chr", "pos", "TRead", "MRead", "ML"))
})

test_that("in-memory 6-column table gets lowercase strand name", {
  dt <- data.table::as.data.table(six_col_df)
  data.table::setnames(dt, "strand", "Strand")
  result <- read_pipeline_files(list(Bismark = dt))

  expect_equal(names(result$Bismark)[6], "strand")
})

test_that("does not modify the caller's in-memory table", {
  dt <- data.table::as.data.table(six_col_df)
  original_names <- names(dt)
  data.table::setnames(dt, "strand", "Strand")
  before <- names(dt)

  read_pipeline_files(list(Bismark = dt))

  expect_equal(names(dt), before)
  expect_false(identical(names(dt), original_names))
})

# ============================================================================
# Mixed inputs (paths and tables in one list)
# ============================================================================

test_that("accepts a mix of file paths and in-memory tables", {
  f <- write_pipeline_file(six_col_df)
  result <- read_pipeline_files(list(
    Biscuit = f,
    Bismark = data.table::as.data.table(six_col_df)
  ))

  expect_equal(names(result), c("Biscuit", "Bismark"))
  expect_true(all(vapply(result, data.table::is.data.table, logical(1))))
  expect_equal(nrow(result$Biscuit), nrow(result$Bismark))
})

test_that("accepts a mix of 5-column file and 6-column table", {
  f <- write_pipeline_file(five_col_df)
  result <- read_pipeline_files(list(
    Biscuit = f,
    Bismark = data.table::as.data.table(six_col_df)
  ))

  expect_equal(ncol(result$Biscuit), 5)
  expect_equal(ncol(result$Bismark), 6)
})

# ============================================================================
# File-path behaviour (unchanged)
# ============================================================================

test_that("returns a named list with names matching input", {
  f1 <- write_pipeline_file(six_col_df)
  f2 <- write_pipeline_file(six_col_df)
  result <- read_pipeline_files(list(Bismark = f1, Encode = f2))

  expect_type(result, "list")
  expect_equal(names(result), c("Bismark", "Encode"))
})

test_that("preserves names with hyphens", {
  f <- write_pipeline_file(six_col_df)
  result <- read_pipeline_files(list(`Bwa-meth` = f))

  expect_equal(names(result), "Bwa-meth")
})

test_that("6-column file gets chr, pos, TRead, MRead, ML, strand", {
  f <- write_pipeline_file(six_col_df)
  result <- read_pipeline_files(list(Bismark = f))

  expect_equal(names(result$Bismark), c("chr", "pos", "TRead", "MRead", "ML", "strand"))
})

test_that("5-column file gets chr, pos, TRead, MRead, ML", {
  f <- write_pipeline_file(five_col_df)
  result <- read_pipeline_files(list(Biscuit = f))

  expect_equal(names(result$Biscuit), c("chr", "pos", "TRead", "MRead", "ML"))
})

test_that("file rows are not treated as a header", {
  f <- write_pipeline_file(six_col_df)
  result <- read_pipeline_files(list(Bismark = f))

  expect_equal(nrow(result$Bismark), 3)
})

test_that("row values match the file contents", {
  f <- write_pipeline_file(six_col_df)
  result <- read_pipeline_files(list(Bismark = f))

  expect_equal(result$Bismark$chr, c("chr1", "chr1", "chr2"))
  expect_equal(result$Bismark$pos, c(100, 101, 1000))
  expect_equal(result$Bismark$strand, c("+", "-", "+"))
})
