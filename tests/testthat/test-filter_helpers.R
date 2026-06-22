# ---- flag_features ----------------------------------------------------------

test_that("flag_features only ever suggests endogenous features", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 6, seed = 1)
  ff <- flag_features(dds)
  expect_equal(nrow(ff), nrow(dds))
  expect_setequal(colnames(ff),
                  c("feature_id", "feature_class", "total_count", "n_detected",
                    "mean_logcounts", "drop_all_zero", "fails_min_count",
                    "fails_min_samples", "fails_filterByExpr", "suggested_drop",
                    "reason"))
  # spike-in / exogenous are exempt: never suggested for removal.
  non_endo <- ff$feature_class != "endogenous"
  expect_true(any(non_endo))
  expect_false(any(ff$suggested_drop[non_endo]))
})

test_that("flag_features always flags all-zero endogenous features", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 2, seed = 2)
  # Zero out two endogenous genes.
  cnt <- SummarizedExperiment::assay(dds, "counts")
  zeroed <- rownames(dds)[SummarizedExperiment::rowData(dds)$feature_class == "endogenous"][1:2]
  cnt[zeroed, ] <- 0L
  SummarizedExperiment::assay(dds, "counts") <- cnt
  ff <- flag_features(dds, use_filter_by_expr = FALSE, min_count = 0)
  expect_true(all(ff$drop_all_zero[match(zeroed, ff$feature_id)]))
  expect_true(all(ff$suggested_drop[match(zeroed, ff$feature_id)]))
  expect_match(ff$reason[match(zeroed, ff$feature_id)], "all-zero")
})

test_that("flag_features manual rule honors min_count and min_samples", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 1, seed = 3)
  total <- rowSums(as.matrix(SummarizedExperiment::assay(dds, "counts")))
  thr <- stats::median(total)
  ff <- flag_features(dds, use_filter_by_expr = FALSE, min_count = thr)
  endo <- ff$feature_class == "endogenous"
  # Endogenous genes below the threshold are flagged; spike/exo never are.
  expect_true(all(ff$fails_min_count[endo] == (ff$total_count[endo] < thr)))
  expect_true(all(ff$suggested_drop[endo] == (ff$total_count[endo] < thr |
                                              ff$total_count[endo] == 0)))
})

test_that("flag_features uses filterByExpr when available", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("edgeR")
  dds <- make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 2, seed = 4)
  ff <- flag_features(dds, use_filter_by_expr = TRUE)
  endo <- ff$feature_class == "endogenous"
  # The filterByExpr column drives the endogenous suggestions (plus all-zero).
  expect_true(all(ff$suggested_drop[endo] ==
                    (ff$fails_filterByExpr[endo] | ff$drop_all_zero[endo])))
  expect_false(any(ff$fails_filterByExpr[!endo]))
})

# ---- flag_samples -----------------------------------------------------------

test_that("flag_samples returns per-reason columns with the metric schema", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 4, n_spike = 2, seed = 5)
  fs <- flag_samples(dds)
  expect_equal(nrow(fs), ncol(dds))
  expect_true(all(c("library_size", "detected", "pct_mito", "within_group_corr",
                    "low_lib_size", "low_detected", "high_mito",
                    "within_group_outlier", "flagged", "reason") %in% colnames(fs)))
  expect_type(fs$flagged, "logical")
  expect_equal(fs$flagged,
               fs$low_lib_size | fs$low_detected | fs$high_mito | fs$within_group_outlier)
})

test_that("flag_samples honors an explicit library-size threshold", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 1, seed = 6)
  lib <- colSums(as.matrix(SummarizedExperiment::assay(dds, "counts")))
  thr <- stats::median(lib)             # ~half the samples fall below the median
  fs <- flag_samples(dds, lib_size_min = thr)
  expect_equal(fs$low_lib_size, fs$library_size < thr)
  expect_true(all(grepl("low library size", fs$reason[fs$low_lib_size])))
})

test_that("flag_samples treats NULL thresholds as disabled (no flag)", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 1, seed = 13)
  # All thresholds NULL -> every rule off -> nothing flagged.
  fs <- flag_samples(dds, lib_size_min = NULL, detected_min = NULL,
                     pct_mito_max = NULL, within_group_z = NULL)
  expect_false(any(fs$low_lib_size))
  expect_false(any(fs$high_mito))
  expect_false(any(fs$within_group_outlier))
  expect_false(any(fs$flagged))
})

test_that("suggest_sample_thresholds returns finite-or-NA fence values", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 4, n_spike = 2, seed = 14)
  th <- suggest_sample_thresholds(dds)
  expect_setequal(names(th), c("lib_size_min", "detected_min", "pct_mito_max"))
  for (v in th) expect_true(is.na(v) || (is.finite(v) && v >= 0))
  # A returned lower fence, when finite, actually drives flag_samples.
  if (is.finite(th$lib_size_min)) {
    fs <- flag_samples(dds, lib_size_min = th$lib_size_min, within_group_z = NULL)
    expect_equal(fs$low_lib_size, fs$library_size < th$lib_size_min)
  }
})

test_that("flag_samples within-group outlier flags an injected bad replicate", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 200, n_per_group = 5, n_spike = 0, seed = 7)
  dds <- ensure_logcounts(dds)
  # Scramble one sample's expression so it correlates poorly with its group.
  lc <- SummarizedExperiment::assay(dds, "logcounts")
  set.seed(99); lc[, 1] <- sample(lc[, 1])
  SummarizedExperiment::assay(dds, "logcounts") <- lc
  fs <- flag_samples(dds, group = "condition", within_group_z = 1)
  expect_true(fs$within_group_outlier[fs$sample == colnames(dds)[1]])
})

# ---- removal_status ---------------------------------------------------------

test_that("removal_status maps the three colour levels", {
  flagged <- c(FALSE, TRUE, TRUE, NA)
  this    <- c(FALSE, TRUE, FALSE, TRUE)
  rs <- removal_status(flagged, this)
  expect_equal(as.character(rs),
               c("pass", "suggested_this", "suggested_other", "pass"))
  expect_equal(levels(rs), c("pass", "suggested_other", "suggested_this"))
  # Without a specific reason, flagged collapses to suggested_other.
  expect_equal(as.character(removal_status(flagged)),
               c("pass", "suggested_other", "suggested_other", "pass"))
})

# ---- drop_features / drop_samples -------------------------------------------

test_that("drop_features subsets, preserves feature_class, and refreshes assays", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 3, seed = 8))
  drop <- rownames(dds)[1:5]
  before_lc <- SummarizedExperiment::assay(dds, "logcounts")
  out <- drop_features(dds, drop)
  expect_equal(nrow(out), nrow(dds) - 5L)
  expect_false(any(drop %in% rownames(out)))
  expect_true("feature_class" %in% colnames(SummarizedExperiment::rowData(out)))
  # logcounts recomputed (library sizes change when features are removed).
  kept <- setdiff(rownames(dds), drop)
  expect_false(isTRUE(all.equal(SummarizedExperiment::assay(out, "logcounts"),
                                before_lc[kept, , drop = FALSE])))
})

test_that("drop_samples subsets and re-estimates size factors only when set", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 2, seed = 9))
  # No size factors set yet -> they stay NULL after dropping.
  out <- drop_samples(dds, colnames(dds)[1])
  expect_equal(ncol(out), ncol(dds) - 1L)
  expect_null(DESeq2::sizeFactors(out))
  # With size factors set, they are recomputed (and stay non-NULL).
  dds2 <- estimate_size_factors_endogenous(dds)
  out2 <- drop_samples(dds2, colnames(dds2)[1])
  expect_false(is.null(DESeq2::sizeFactors(out2)))
  expect_length(DESeq2::sizeFactors(out2), ncol(dds2) - 1L)
})

test_that("drop_samples drops now-empty factor levels", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 1, seed = 12))
  treated <- colnames(dds)[SummarizedExperiment::colData(dds)$condition == "treated"]
  out <- drop_samples(dds, treated)
  # "treated" no longer appears among any remaining sample -> level is dropped.
  expect_false("treated" %in% levels(SummarizedExperiment::colData(out)$condition))
  expect_setequal(levels(SummarizedExperiment::colData(out)$condition), "control")
})

test_that("drop_* are no-ops on empty / unknown ids", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 1, seed = 10)
  expect_equal(dim(drop_features(dds, character(0))), dim(dds))
  expect_equal(dim(drop_samples(dds, "not-a-sample")), dim(dds))
})

# ---- qc_filter_density ------------------------------------------------------

test_that("qc_filter_density returns before/after long data over endogenous", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 4, seed = 11))
  endo_ids <- rownames(dds)[SummarizedExperiment::rowData(dds)$feature_class == "endogenous"]
  keep <- endo_ids[1:20]
  dens <- qc_filter_density(dds, keep)
  expect_setequal(colnames(dens), c("sample", "value", "status"))
  expect_setequal(levels(dens$status), c("before", "after"))
  n_before <- sum(dens$status == "before")
  n_after  <- sum(dens$status == "after")
  expect_equal(n_before, length(endo_ids) * ncol(dds))
  expect_equal(n_after, 20L * ncol(dds))
})
