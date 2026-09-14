context("Strand plotting functions (refactored)")

# Setup: Create mock result object matching prepare_filtered_cpg_table() output
create_mock_result <- function(n_unfiltered = 1000, n_filtered = 800) {
  # Create unfiltered (long-format) data
  unfiltered <- data.table::data.table(
    Pipeline = rep(c("Bismark", "Bwameth", "Biscuit", "Encode"),
                   each = n_unfiltered / 4),
    Strand = rep(c("S1", "S2"), n_unfiltered / 2),
    reads = c(
      rnorm(n_unfiltered / 2, mean = 10, sd = 5),  # S1 tends lower
      rnorm(n_unfiltered / 2, mean = 20, sd = 8)   # S2 tends higher
    )
  )
  unfiltered[reads < 0, reads := 0L]
  unfiltered[, reads := as.integer(reads)]

  # Create merged (filtered) data - subset with higher read counts
  merged <- unfiltered[reads >= 5]

  # Create thresholds table
  thresholds <- data.table::data.table(
    Pipeline = c("Bismark", "Bwameth", "Biscuit", "Encode"),
    Threshold_1M = c(8L, 8L, 7L, 9L),
    Threshold_2M = c(4L, 4L, 3L, 4L),
    Threshold_3M = c(2L, 2L, 2L, 3L),
    Threshold_4M = c(1L, 1L, 1L, 2L)
  )

  list(unfiltered = unfiltered, merged = merged, thresholds = thresholds)
}

# ============================================================================
# Tests for plot_strand_boxplot()
# ============================================================================

test_that("plot_strand_boxplot returns a ggplot object", {
  result <- create_mock_result()
  p <- plot_strand_boxplot(result, outdir = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("plot_strand_boxplot validates required columns", {
  result <- create_mock_result()
  result$unfiltered[, reads := NULL]

  expect_error(
    plot_strand_boxplot(result),
    "result\\$unfiltered must contain columns"
  )
})

test_that("plot_strand_boxplot saves to file when outdir provided", {
  skip_on_cran()

  result <- create_mock_result()
  tmpdir <- tempdir()

  plot_strand_boxplot(result, outdir = tmpdir, sample_name = "test_sample")

  expected_file <- file.path(tmpdir, "test_sample_01Fig_Boxplot_S1_S2.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_strand_boxplot doesn't save when outdir is NULL", {
  result <- create_mock_result()
  tmpdir <- tempdir()

  old_files <- list.files(tmpdir, pattern = ".*_01Fig_Boxplot.*")
  if (length(old_files) > 0) file.remove(file.path(tmpdir, old_files))

  p <- plot_strand_boxplot(result, outdir = NULL, sample_name = "test_no_save")

  new_files <- list.files(tmpdir, pattern = "test_no_save.*_01Fig_Boxplot")
  expect_length(new_files, 0)
})

# ============================================================================
# Tests for plot_strand_frequency_unfiltered()
# ============================================================================

test_that("plot_strand_frequency_unfiltered returns a ggplot object", {
  result <- create_mock_result()
  p <- plot_strand_frequency_unfiltered(result, outdir = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("plot_strand_frequency_unfiltered validates required columns", {
  result <- create_mock_result()
  result$unfiltered[, reads := NULL]

  expect_error(
    plot_strand_frequency_unfiltered(result),
    "result\\$unfiltered must contain columns"
  )
})

test_that("plot_strand_frequency_unfiltered without thresholds has no vlines", {
  result <- create_mock_result()
  p <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = FALSE)

  geom_types <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_false(any(grepl("GeomVline", geom_types)))
})

test_that("plot_strand_frequency_unfiltered with thresholds has vlines", {
  result <- create_mock_result()
  p <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = TRUE)

  geom_types <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_true(any(grepl("GeomVline", geom_types)))
})

test_that("plot_strand_frequency_unfiltered with thresholds has text labels", {
  result <- create_mock_result()
  p <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = TRUE)

  geom_types <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_true(any(grepl("GeomText", geom_types)))
})

test_that("plot_strand_frequency_unfiltered subtitle changes with show_thresholds", {
  result <- create_mock_result()

  p_no_thresh <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = FALSE)
  p_with_thresh <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = TRUE)

  subtitle_no_thresh <- p_no_thresh$labels$subtitle
  subtitle_with_thresh <- p_with_thresh$labels$subtitle

  expect_not_equal(subtitle_no_thresh, subtitle_with_thresh)
  expect_true(grepl("thresholds", subtitle_with_thresh, ignore.case = TRUE))
})

test_that("plot_strand_frequency_unfiltered filters by reads_max", {
  result <- create_mock_result()
  p <- plot_strand_frequency_unfiltered(result, outdir = NULL, reads_max = 50)

  plot_data <- ggplot2::layer_data(p)
  expect_true(all(plot_data$x <= 50, na.rm = TRUE))
})

test_that("plot_strand_frequency_unfiltered saves with correct filename", {
  skip_on_cran()

  result <- create_mock_result()
  tmpdir <- tempdir()

  plot_strand_frequency_unfiltered(result,
                                   outdir = tmpdir,
                                   sample_name = "test_unfiltered",
                                   show_thresholds = FALSE)

  expected_file <- file.path(tmpdir, "test_unfiltered_02Fig_Frequency_Unfiltered.pdf")
  expect_true(file.exists(expected_file))
  file.remove(expected_file)
})

test_that("plot_strand_frequency_unfiltered saves with threshold filename", {
  skip_on_cran()

  result <- create_mock_result()
  tmpdir <- tempdir()

  plot_strand_frequency_unfiltered(result,
                                   outdir = tmpdir,
                                   sample_name = "test_unfiltered_thresh",
                                   show_thresholds = TRUE)

  expected_file <- file.path(tmpdir, "test_unfiltered_thresh_02Fig_Frequency_Unfiltered_with_Thresholds.pdf")
  expect_true(file.exists(expected_file))
  file.remove(expected_file)
})

# ============================================================================
# Tests for plot_strand_frequency_filtered()
# ============================================================================

test_that("plot_strand_frequency_filtered returns a ggplot object", {
  result <- create_mock_result()
  p <- plot_strand_frequency_filtered(result, outdir = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("plot_strand_frequency_filtered validates required columns in merged", {
  result <- create_mock_result()
  result$merged[, reads := NULL]

  expect_error(
    plot_strand_frequency_filtered(result),
    "result\\$merged must contain columns"
  )
})

test_that("plot_strand_frequency_filtered without thresholds has no vlines", {
  result <- create_mock_result()
  p <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = FALSE)

  geom_types <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_false(any(grepl("GeomVline", geom_types)))
})

test_that("plot_strand_frequency_filtered with thresholds has vlines", {
  result <- create_mock_result()
  p <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = TRUE)

  geom_types <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_true(any(grepl("GeomVline", geom_types)))
})

test_that("plot_strand_frequency_filtered with thresholds has text labels", {
  result <- create_mock_result()
  p <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = TRUE)

  geom_types <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_true(any(grepl("GeomText", geom_types)))
})

test_that("plot_strand_frequency_filtered title indicates filtered data", {
  result <- create_mock_result()
  p <- plot_strand_frequency_filtered(result, outdir = NULL)

  title <- p$labels$title
  expect_true(grepl("AFTER FILTERING", title, ignore.case = TRUE))
})

test_that("plot_strand_frequency_filtered has different subtitle with thresholds", {
  result <- create_mock_result()

  p_no_thresh <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = FALSE)
  p_with_thresh <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = TRUE)

  subtitle_no_thresh <- p_no_thresh$labels$subtitle
  subtitle_with_thresh <- p_with_thresh$labels$subtitle

  expect_not_equal(subtitle_no_thresh, subtitle_with_thresh)
  expect_true(grepl("thresholds", subtitle_with_thresh, ignore.case = TRUE))
})

test_that("plot_strand_frequency_filtered shows filtered (subset) data", {
  result <- create_mock_result()

  n_unfiltered <- nrow(result$unfiltered)
  n_filtered <- nrow(result$merged)

  expect_lt(n_filtered, n_unfiltered)
})

test_that("plot_strand_frequency_filtered saves with correct filename", {
  skip_on_cran()

  result <- create_mock_result()
  tmpdir <- tempdir()

  plot_strand_frequency_filtered(result,
                                 outdir = tmpdir,
                                 sample_name = "test_filtered",
                                 show_thresholds = FALSE)

  expected_file <- file.path(tmpdir, "test_filtered_04Fig_Frequency_Filtered.pdf")
  expect_true(file.exists(expected_file))
  file.remove(expected_file)
})

test_that("plot_strand_frequency_filtered saves with threshold filename", {
  skip_on_cran()

  result <- create_mock_result()
  tmpdir <- tempdir()

  plot_strand_frequency_filtered(result,
                                 outdir = tmpdir,
                                 sample_name = "test_filtered_thresh",
                                 show_thresholds = TRUE)

  expected_file <- file.path(tmpdir, "test_filtered_thresh_04Fig_Frequency_Filtered_with_Thresholds.pdf")
  expect_true(file.exists(expected_file))
  file.remove(expected_file)
})

# ============================================================================
# Cross-function integration tests
# ============================================================================

test_that("all plot functions work with same result object", {
  result <- create_mock_result()

  p1 <- plot_strand_boxplot(result, outdir = NULL)
  p2 <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = FALSE)
  p3 <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = FALSE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("all plot combinations work together", {
  result <- create_mock_result()

  # All combinations of show_thresholds
  p1 <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = FALSE)
  p2 <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = TRUE)
  p3 <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = FALSE)
  p4 <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = TRUE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("all plots can save to same directory without conflicts", {
  skip_on_cran()

  result <- create_mock_result()
  tmpdir <- tempdir()
  sample_id <- "test_all_plots"

  plot_strand_boxplot(result, outdir = tmpdir, sample_name = sample_id)
  plot_strand_frequency_unfiltered(result, outdir = tmpdir, sample_name = sample_id,
                                   show_thresholds = FALSE)
  plot_strand_frequency_unfiltered(result, outdir = tmpdir, sample_name = sample_id,
                                   show_thresholds = TRUE)
  plot_strand_frequency_filtered(result, outdir = tmpdir, sample_name = sample_id,
                                 show_thresholds = FALSE)
  plot_strand_frequency_filtered(result, outdir = tmpdir, sample_name = sample_id,
                                 show_thresholds = TRUE)

  # Check all 5 files exist
  files <- list.files(tmpdir, pattern = paste0(sample_id, "_.*Fig"))
  expect_length(files, 5)

  # Cleanup
  file.remove(file.path(tmpdir, files))
})

test_that("plots handle edge case: very small dataset", {
  result <- create_mock_result(n_unfiltered = 10, n_filtered = 5)

  p1 <- plot_strand_boxplot(result, outdir = NULL)
  p2 <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = FALSE)
  p3 <- plot_strand_frequency_unfiltered(result, outdir = NULL, show_thresholds = TRUE)
  p4 <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = FALSE)
  p5 <- plot_strand_frequency_filtered(result, outdir = NULL, show_thresholds = TRUE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
  expect_s3_class(p5, "ggplot")
})

test_that("all functions return invisible plots", {
  result <- create_mock_result()

  out1 <- capture.output({ p1 <- plot_strand_boxplot(result, outdir = NULL) })
  out2 <- capture.output({ p2 <- plot_strand_frequency_unfiltered(result, outdir = NULL) })
  out3 <- capture.output({ p3 <- plot_strand_frequency_filtered(result, outdir = NULL) })

  expect_length(out1, 0)
  expect_length(out2, 0)
  expect_length(out3, 0)
})
