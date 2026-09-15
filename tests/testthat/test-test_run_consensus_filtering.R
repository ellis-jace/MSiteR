context("run_consensus_filtering()")

# ============================================================================
# SETUP: Create mock strand filtering output
# ============================================================================

create_mock_strand_result <- function(n_sites = 1000) {
  # Create merged wide table (chr, pos, Bismark_TR, Bismark_ML, ..., Bwameth_TR, ...)
  # Simulate sites detected by different combinations of pipelines

  merged <- data.table::data.table(
    chr = rep(1:4, n_sites / 4),
    pos = rep(1:(n_sites / 4), 4),
    # Bismark: detected in most sites
    Bismark_TR = c(rnorm(n_sites * 0.9, 20, 8), rep(NA, n_sites * 0.1)),
    Bismark_MR = c(rnorm(n_sites * 0.9, 10, 5), rep(NA, n_sites * 0.1)),
    Bismark_ML = c(rnorm(n_sites * 0.9, 0.5, 0.2), rep(NA, n_sites * 0.1)),
    # Bwameth: detected in ~80% of sites
    Bwameth_TR = c(rnorm(n_sites * 0.8, 18, 7), rep(NA, n_sites * 0.2)),
    Bwameth_MR = c(rnorm(n_sites * 0.8, 9, 4), rep(NA, n_sites * 0.2)),
    Bwameth_ML = c(rnorm(n_sites * 0.8, 0.5, 0.2), rep(NA, n_sites * 0.2)),
    # Biscuit: detected in ~75% of sites
    Biscuit_TR = c(rnorm(n_sites * 0.75, 19, 8), rep(NA, n_sites * 0.25)),
    Biscuit_MR = c(rnorm(n_sites * 0.75, 9.5, 4), rep(NA, n_sites * 0.25)),
    Biscuit_ML = c(rnorm(n_sites * 0.75, 0.5, 0.2), rep(NA, n_sites * 0.25)),
    # Encode: detected in ~70% of sites
    Encode_TR = c(rnorm(n_sites * 0.7, 17, 7), rep(NA, n_sites * 0.3)),
    Encode_MR = c(rnorm(n_sites * 0.7, 8.5, 4), rep(NA, n_sites * 0.3)),
    Encode_ML = c(rnorm(n_sites * 0.7, 0.5, 0.2), rep(NA, n_sites * 0.3))
  )

  # Ensure positive values
  for (col in grep("_TR$|_MR$", names(merged), value = TRUE)) {
    merged[, (col) := pmax(0, get(col))]
  }
  for (col in grep("_ML$", names(merged), value = TRUE)) {
    merged[, (col) := pmax(0, pmin(1, get(col)))]
  }

  # Convert to integers where appropriate
  for (col in grep("_TR$|_MR$", names(merged), value = TRUE)) {
    merged[, (col) := as.integer(get(col))]
  }

  # Create strand thresholds (mock)
  strand_thresholds <- data.table::data.table(
    Pipeline = c("Bismark", "Bwameth", "Biscuit", "Encode"),
    cutoff1 = c(8L, 7L, 8L, 7L),
    cutoff2 = c(4L, 4L, 4L, 3L)
  )

  # Create unfiltered (same as merged for this mock)
  unfiltered <- data.table::copy(merged)

  list(
    merged = merged,
    unfiltered = unfiltered,
    thresholds = strand_thresholds
  )
}

# ============================================================================
# TESTS: Input validation
# ============================================================================

test_that("run_consensus_filtering rejects non-list input", {
  expect_error(
    run_consensus_filtering(data.frame()),
    "result must be a list"
  )
})

test_that("run_consensus_filtering requires merged component", {
  result <- create_mock_strand_result()
  result$merged <- NULL

  expect_error(
    run_consensus_filtering(result),
    "result missing required components"
  )
})

test_that("run_consensus_filtering requires unfiltered component", {
  result <- create_mock_strand_result()
  result$unfiltered <- NULL

  expect_error(
    run_consensus_filtering(result),
    "result missing required components"
  )
})

test_that("run_consensus_filtering requires thresholds component", {
  result <- create_mock_strand_result()
  result$thresholds <- NULL

  expect_error(
    run_consensus_filtering(result),
    "result missing required components"
  )
})

test_that("run_consensus_filtering rejects non-data.table merged", {
  result <- create_mock_strand_result()
  result$merged <- as.data.frame(result$merged)

  expect_error(
    run_consensus_filtering(result),
    "result\\$merged must be a data.table"
  )
})

# ============================================================================
# TESTS: Output structure
# ============================================================================

test_that("run_consensus_filtering returns correct list structure", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  expect_true(is.list(result_consensus))
  expect_equal(names(result_consensus),
               c("merged", "unfiltered", "thresholds", "filtering_stats"))
})

test_that("run_consensus_filtering returns data.tables", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  expect_s3_class(result_consensus$merged, "data.table")
  expect_s3_class(result_consensus$unfiltered, "data.table")
  expect_s3_class(result_consensus$thresholds, "data.table")
  expect_s3_class(result_consensus$filtering_stats, "data.table")
})

# ============================================================================
# TESTS: Consensus metrics added
# ============================================================================

test_that("run_consensus_filtering adds n_methods column", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  expect_true("n_methods" %in% names(result_consensus$unfiltered))
})

test_that("run_consensus_filtering adds median_Tread column", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  expect_true("median_Tread" %in% names(result_consensus$unfiltered))
})

test_that("run_consensus_filtering adds median_MethyL column", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  expect_true("median_MethyL" %in% names(result_consensus$unfiltered))
})

test_that("run_consensus_filtering adds median_Mread column", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  expect_true("median_Mread" %in% names(result_consensus$unfiltered))
})

test_that("n_methods is 1-4 or NA", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  n_methods_vals <- result_consensus$unfiltered$n_methods
  expect_true(all(n_methods_vals %in% c(1, 2, 3, 4), na.rm = TRUE))
})

# ============================================================================
# TESTS: Filtering effect
# ============================================================================

test_that("run_consensus_filtering filters sites (fewer after)", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  n_before <- nrow(result_consensus$unfiltered)
  n_after <- nrow(result_consensus$merged)

  # Should remove some sites
  expect_lt(n_after, n_before)
})

test_that("run_consensus_filtering preserves all columns in merged", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  # All columns from unfiltered should exist in merged (it's a subset)
  expect_true(all(names(result_consensus$unfiltered) %in% names(result_consensus$merged)))
})

# ============================================================================
# TESTS: Thresholds table
# ============================================================================

test_that("consensus thresholds has correct structure", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  thresh <- result_consensus$thresholds

  expect_equal(nrow(thresh), 4)  # 4 tools
  expect_true(all(c("Tool", "Threshold_1M", "Threshold_2M", "Threshold_3M", "Threshold_4M")
                  %in% names(thresh)))
})

test_that("consensus thresholds scale down by tier", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  thresh <- result_consensus$thresholds

  # For each tool, thresholds should generally decrease from 1M to 4M
  for (i in 1:nrow(thresh)) {
    t1m <- thresh[i, Threshold_1M]
    t2m <- thresh[i, Threshold_2M]
    t3m <- thresh[i, Threshold_3M]
    t4m <- thresh[i, Threshold_4M]

    # 1M should be highest (95th percentile, strict)
    expect_gte(t1m, t2m)
    expect_gte(t2m, t3m)
    expect_gte(t3m, t4m)
  }
})

# ============================================================================
# TESTS: Filtering statistics
# ============================================================================

test_that("filtering_stats has expected metrics", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  stats <- result_consensus$filtering_stats

  expect_true("sites_before" %in% stats$metric)
  expect_true("sites_after" %in% stats$metric)
  expect_true("1m_sites_before" %in% stats$metric)
  expect_true("1m_sites_after" %in% stats$metric)
  expect_true("1m_retain_pct" %in% stats$metric)
})

test_that("filtering_stats after counts <= before", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  stats <- result_consensus$filtering_stats
  n_before <- stats[metric == "sites_before", value]
  n_after <- stats[metric == "sites_after", value]
  n_1m_before <- stats[metric == "1m_sites_before", value]
  n_1m_after <- stats[metric == "1m_sites_after", value]

  expect_lte(n_after, n_before)
  expect_lte(n_1m_after, n_1m_before)
})

# ============================================================================
# TESTS: Edge cases
# ============================================================================

test_that("run_consensus_filtering handles small dataset", {
  result_strand <- create_mock_strand_result(n_sites = 10)
  result_consensus <- run_consensus_filtering(result_strand)

  expect_true(nrow(result_consensus$merged) <= nrow(result_consensus$unfiltered))
})

test_that("run_consensus_filtering handles all 1M sites", {
  result_strand <- create_mock_strand_result()

  # Keep only single-detection sites
  result_strand$merged <- result_strand$merged[1:50]
  # Blank out all but one pipeline per site
  for (i in 1:nrow(result_strand$merged)) {
    pipelines <- c("Bismark", "Bwameth", "Biscuit", "Encode")
    keep_pipeline <- sample(pipelines, 1)
    for (p in setdiff(pipelines, keep_pipeline)) {
      data.table::set(result_strand$merged, i, paste0(p, "_TR"), NA)
      data.table::set(result_strand$merged, i, paste0(p, "_ML"), NA)
    }
  }
  result_strand$unfiltered <- data.table::copy(result_strand$merged)

  result_consensus <- run_consensus_filtering(result_strand)

  expect_s3_class(result_consensus$merged, "data.table")
  expect_true(nrow(result_consensus$merged) >= 0)
})

test_that("run_consensus_filtering handles all 4M sites", {
  result_strand <- create_mock_strand_result()

  # Keep only complete-detection sites
  result_strand$merged <- result_strand$merged[
    !is.na(Bismark_TR) & !is.na(Bwameth_TR) & !is.na(Biscuit_TR) & !is.na(Encode_TR)
  ]
  result_strand$unfiltered <- data.table::copy(result_strand$merged)

  result_consensus <- run_consensus_filtering(result_strand)

  expect_s3_class(result_consensus$merged, "data.table")
})

# ============================================================================
# TESTS: Integration with plotting
# ============================================================================

test_that("run_consensus_filtering output suitable for plotting", {
  result_strand <- create_mock_strand_result()
  result_consensus <- run_consensus_filtering(result_strand)

  # Plotting functions need:
  # - unfiltered with n_methods for before plot
  # - merged with n_methods for after plot
  # - thresholds table

  expect_true("n_methods" %in% names(result_consensus$unfiltered))
  expect_true("n_methods" %in% names(result_consensus$merged))
  expect_true("Threshold_1M" %in% names(result_consensus$thresholds))
})

test_that("run_consensus_filtering can be called multiple times on same input", {
  result_strand <- create_mock_strand_result()

  result1 <- run_consensus_filtering(result_strand)
  result2 <- run_consensus_filtering(result_strand)

  # Results should be identical
  expect_equal(nrow(result1$merged), nrow(result2$merged))
  expect_equal(nrow(result1$unfiltered), nrow(result2$unfiltered))
})

