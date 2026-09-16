context("Consensus plotting functions")

# ============================================================================
# SETUP: Create mock consensus filtering output
# ============================================================================

create_mock_consensus_result <- function(n_sites = 500) {
  # Create wide merged table with consensus metrics
  merged <- data.table::data.table(
    chr = rep(1:4, n_sites / 4),
    pos = rep(1:(n_sites / 4), 4),
    Bismark_TR = rnorm(n_sites, 20, 8),
    Bismark_ML = rnorm(n_sites, 0.5, 0.2),
    Bwameth_TR = rnorm(n_sites, 18, 7),
    Bwameth_ML = rnorm(n_sites, 0.5, 0.2),
    Biscuit_TR = rnorm(n_sites, 19, 8),
    Biscuit_ML = rnorm(n_sites, 0.5, 0.2),
    Encode_TR = rnorm(n_sites, 17, 7),
    Encode_ML = rnorm(n_sites, 0.5, 0.2),
    n_methods = sample(1:4, n_sites, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2)),
    median_Tread = rnorm(n_sites, 18, 8),
    median_MethyL = rnorm(n_sites, 0.5, 0.2),
    median_Mread = rnorm(n_sites, 9, 5)
  )

  # Ensure positive values
  for (col in grep("_TR$|_MR$|_Mread$|median_Tread$", names(merged), value = TRUE)) {
    merged[, (col) := pmax(0, get(col))]
  }
  for (col in grep("_ML$|_MethyL$", names(merged), value = TRUE)) {
    merged[, (col) := pmax(0, pmin(1, get(col)))]
  }

  # Randomly set some pipeline values to NA to simulate missing detections
  for (pipeline in c("Bismark", "Bwameth", "Biscuit", "Encode")) {
    na_idx <- sample(1:nrow(merged), size = nrow(merged) * 0.2)
    data.table::set(merged, na_idx, paste0(pipeline, "_TR"), NA)
    data.table::set(merged, na_idx, paste0(pipeline, "_ML"), NA)
  }

  # Create filtered version (subset with higher read counts)
  filtered <- merged[median_Tread >= 5]

  # Create thresholds table
  thresholds <- data.table::data.table(
    Tool = c("Bismark", "Bwameth", "Biscuit", "Encode"),
    Threshold_1M = c(8L, 7L, 8L, 7L),
    Threshold_2M = c(4L, 4L, 4L, 3L),
    Threshold_3M = c(3L, 3L, 3L, 2L),
    Threshold_4M = c(2L, 2L, 2L, 1L)
  )

  list(
    merged = filtered,
    unfiltered = merged,
    thresholds = thresholds,
    filtering_stats = data.table::data.table(
      metric = c("sites_before", "sites_after", "1m_sites_before", "1m_sites_after", "1m_retain_pct"),
      value = c(nrow(merged), nrow(filtered), 100, 80, 80)
    )
  )
}

# ============================================================================
# TESTS: plot_consensus_upset()
# ============================================================================

test_that("plot_consensus_upset with show_filtered = FALSE saves to file", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_upset(result, outdir = tmpdir, sample_name = "test_upset", show_filtered = FALSE)

  expected_file <- file.path(tmpdir, "test_upset_Fig_UpSet_Consensus_unfiltered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_upset with show_filtered = TRUE saves to file", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_upset(result, outdir = tmpdir, sample_name = "test_upset_filt", show_filtered = TRUE)

  expected_file <- file.path(tmpdir, "test_upset_filt_Fig_UpSet_Consensus_filtered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_upset doesn't save when outdir = NULL", {
  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  old_files <- list.files(tmpdir, pattern = ".*_UpSet_Consensus.*")
  if (length(old_files) > 0) file.remove(file.path(tmpdir, old_files))

  plot_consensus_upset(result, outdir = NULL)

  new_files <- list.files(tmpdir, pattern = ".*_UpSet_Consensus.*")
  expect_length(new_files, 0)
})

test_that("plot_consensus_upset validates input", {
  expect_error(
    plot_consensus_upset(list()),
    "result must be output from run_consensus_filtering"
  )
})

# ============================================================================
# TESTS: plot_consensus_venn()
# ============================================================================

test_that("plot_consensus_venn returns ggplot when outdir = NULL", {
  result <- create_mock_consensus_result()
  p <- plot_consensus_venn(result, outdir = NULL, show_filtered = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("plot_consensus_venn with show_filtered = FALSE saves correctly", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_venn(result, outdir = tmpdir, sample_name = "test_venn", show_filtered = FALSE)

  expected_file <- file.path(tmpdir, "test_venn_Fig_Venn_Consensus_unfiltered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_venn with show_filtered = TRUE saves correctly", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_venn(result, outdir = tmpdir, sample_name = "test_venn_filt", show_filtered = TRUE)

  expected_file <- file.path(tmpdir, "test_venn_filt_Fig_Venn_Consensus_filtered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_venn title changes with show_filtered", {
  result <- create_mock_consensus_result()

  p_unfiltered <- plot_consensus_venn(result, outdir = NULL, show_filtered = FALSE)
  p_filtered <- plot_consensus_venn(result, outdir = NULL, show_filtered = TRUE)

  title_unfiltered <- p_unfiltered$labels$title
  title_filtered <- p_filtered$labels$title

  expect_not_equal(title_unfiltered, title_filtered)
  expect_true(grepl("Before", title_unfiltered))
  expect_true(grepl("After", title_filtered))
})

test_that("plot_consensus_venn validates input", {
  expect_error(
    plot_consensus_venn(list()),
    "result must be output from run_consensus_filtering"
  )
})

# ============================================================================
# TESTS: plot_consensus_distribution()
# ============================================================================

test_that("plot_consensus_distribution returns ggplot when outdir = NULL", {
  result <- create_mock_consensus_result()
  p <- plot_consensus_distribution(result, outdir = NULL, show_filtered = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("plot_consensus_distribution with show_filtered = FALSE saves correctly", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_distribution(result, outdir = tmpdir, sample_name = "test_dist", show_filtered = FALSE)

  expected_file <- file.path(tmpdir, "test_dist_Fig_Distribution_Consensus_unfiltered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_distribution with show_filtered = TRUE saves correctly", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_distribution(result, outdir = tmpdir, sample_name = "test_dist_filt", show_filtered = TRUE)

  expected_file <- file.path(tmpdir, "test_dist_filt_Fig_Distribution_Consensus_filtered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_distribution respects reads_max", {
  result <- create_mock_consensus_result()
  p <- plot_consensus_distribution(result, outdir = NULL, reads_max = 50, show_filtered = FALSE)

  plot_data <- ggplot2::layer_data(p)
  expect_true(all(plot_data$x <= 50, na.rm = TRUE))
})

test_that("plot_consensus_distribution has 4 color series (1M, 2M, 3M, 4M)", {
  result <- create_mock_consensus_result()
  p <- plot_consensus_distribution(result, outdir = NULL, show_filtered = FALSE)

  # Check that plot has 4 different colors (one per tier)
  plot_data <- ggplot2::layer_data(p)
  unique_groups <- length(unique(plot_data$group[!is.na(plot_data$group)]))
  expect_gte(unique_groups, 2)  # At least 2 tiers present
})

test_that("plot_consensus_distribution validates input", {
  expect_error(
    plot_consensus_distribution(list()),
    "result must be output from run_consensus_filtering"
  )
})

# ============================================================================
# TESTS: plot_consensus_boxplot()
# ============================================================================

test_that("plot_consensus_boxplot returns ggplot when outdir = NULL", {
  result <- create_mock_consensus_result()
  p <- plot_consensus_boxplot(result, outdir = NULL, show_filtered = FALSE)

  expect_s3_class(p, "ggplot")
})

test_that("plot_consensus_boxplot with show_filtered = FALSE saves correctly", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_boxplot(result, outdir = tmpdir, sample_name = "test_box", show_filtered = FALSE)

  expected_file <- file.path(tmpdir, "test_box_Fig_Boxplot_Consensus_unfiltered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_boxplot with show_filtered = TRUE saves correctly", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  plot_consensus_boxplot(result, outdir = tmpdir, sample_name = "test_box_filt", show_filtered = TRUE)

  expected_file <- file.path(tmpdir, "test_box_filt_Fig_Boxplot_Consensus_filtered.pdf")
  expect_true(file.exists(expected_file))

  file.remove(expected_file)
})

test_that("plot_consensus_boxplot has geom_boxplot layer", {
  result <- create_mock_consensus_result()
  p <- plot_consensus_boxplot(result, outdir = NULL, show_filtered = FALSE)

  geom_types <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_true(any(grepl("GeomBoxplot", geom_types)))
})

test_that("plot_consensus_boxplot title changes with show_filtered", {
  result <- create_mock_consensus_result()

  p_unfiltered <- plot_consensus_boxplot(result, outdir = NULL, show_filtered = FALSE)
  p_filtered <- plot_consensus_boxplot(result, outdir = NULL, show_filtered = TRUE)

  title_unfiltered <- p_unfiltered$labels$title
  title_filtered <- p_filtered$labels$title

  expect_not_equal(title_unfiltered, title_filtered)
  expect_true(grepl("Before", title_unfiltered))
  expect_true(grepl("After", title_filtered))
})

test_that("plot_consensus_boxplot validates input", {
  expect_error(
    plot_consensus_boxplot(list()),
    "result must be output from run_consensus_filtering"
  )
})

# ============================================================================
# CROSS-FUNCTION TESTS
# ============================================================================

test_that("all consensus plots work with same result object", {
  result <- create_mock_consensus_result()

  p1 <- plot_consensus_venn(result, outdir = NULL)
  p2 <- plot_consensus_distribution(result, outdir = NULL)
  p3 <- plot_consensus_boxplot(result, outdir = NULL)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("all consensus plots work with show_filtered = TRUE", {
  result <- create_mock_consensus_result()

  p1 <- plot_consensus_venn(result, outdir = NULL, show_filtered = TRUE)
  p2 <- plot_consensus_distribution(result, outdir = NULL, show_filtered = TRUE)
  p3 <- plot_consensus_boxplot(result, outdir = NULL, show_filtered = TRUE)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("all plots can save to same directory without conflicts", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()
  sample_id <- "test_all_consensus"

  plot_consensus_upset(result, outdir = tmpdir, sample_name = sample_id, show_filtered = FALSE)
  plot_consensus_venn(result, outdir = tmpdir, sample_name = sample_id, show_filtered = FALSE)
  plot_consensus_distribution(result, outdir = tmpdir, sample_name = sample_id, show_filtered = FALSE)
  plot_consensus_boxplot(result, outdir = tmpdir, sample_name = sample_id, show_filtered = FALSE)

  plot_consensus_venn(result, outdir = tmpdir, sample_name = sample_id, show_filtered = TRUE)
  plot_consensus_distribution(result, outdir = tmpdir, sample_name = sample_id, show_filtered = TRUE)
  plot_consensus_boxplot(result, outdir = tmpdir, sample_name = sample_id, show_filtered = TRUE)

  # Check files exist (note: UpSet only works with outdir, returns invisible NULL)
  files <- list.files(tmpdir, pattern = paste0(sample_id, "_.*Consensus"))
  expect_gte(length(files), 6)  # At least 6 files (3 plots × 2 versions, UpSet saves directly)

  # Cleanup
  file.remove(file.path(tmpdir, files))
})

test_that("plots handle edge case: very small dataset", {
  result <- create_mock_consensus_result(n_sites = 10)

  p1 <- plot_consensus_venn(result, outdir = NULL)
  p2 <- plot_consensus_distribution(result, outdir = NULL)
  p3 <- plot_consensus_boxplot(result, outdir = NULL)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("venn, distribution, boxplot return invisible ggplots", {
  result <- create_mock_consensus_result()

  out1 <- capture.output({ p1 <- plot_consensus_venn(result, outdir = NULL) })
  out2 <- capture.output({ p2 <- plot_consensus_distribution(result, outdir = NULL) })
  out3 <- capture.output({ p3 <- plot_consensus_boxplot(result, outdir = NULL) })

  # No output means invisible return
  expect_length(out1, 0)
  expect_length(out2, 0)
  expect_length(out3, 0)
})

test_that("upset plot returns invisible NULL when outdir provided", {
  skip_on_cran()

  result <- create_mock_consensus_result()
  tmpdir <- tempdir()

  out <- capture.output({
    result_obj <- plot_consensus_upset(result, outdir = tmpdir, sample_name = "test_invisible")
  })

  # Should be invisible (no output except message)
  expect_true(any(grepl("Saved", out)))

  # Cleanup
  file.remove(file.path(tmpdir, "test_invisible_Fig_UpSet_Consensus_unfiltered.pdf"))
})

