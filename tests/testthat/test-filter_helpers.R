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

# ---- flag_samples: spike-in criteria (P3e) ----------------------------------

# A synthetic per-sample dose-response summary aligned to a dds, so the spike
# rules can be exercised deterministically without depending on the random fit.
.mock_spike_summary <- function(dds, n_spike_detected = 5L, n_points = 5L,
                                slope = 1, r_squared = 0.95, lod = 1) {
  cn <- colnames(dds); n <- length(cn)
  rep_n <- function(x) if (length(x) == 1L) rep(x, n) else x
  data.frame(sample = cn, pct_spike = NA_real_,
             n_spike_detected = rep_n(n_spike_detected), n_points = rep_n(n_points),
             slope = rep_n(slope), r_squared = rep_n(r_squared), lod = rep_n(lod),
             stringsAsFactors = FALSE)
}

test_that("flag_samples adds the spike columns and is unchanged when spike = NULL", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 4, seed = 21)
  fs <- flag_samples(dds, within_group_z = NULL)        # spike = NULL (default)
  spike_cols <- c("n_spike_detected", "dose_r2", "dose_slope", "lod", "low_spike",
                  "high_spike", "low_spike_detected", "low_dose_r2", "bad_slope",
                  "few_spike_points")
  expect_true(all(spike_cols %in% colnames(fs)))
  expect_true(all(is.na(fs$n_spike_detected)))           # no spike summary -> NA metrics
  expect_false(any(fs$low_spike | fs$high_spike | fs$low_spike_detected |
                   fs$low_dose_r2 | fs$bad_slope | fs$few_spike_points))
})

test_that("flag_samples fires detected / R2 / slope spike rules when thresholds set", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 4, seed = 22)
  sp <- .mock_spike_summary(dds, n_spike_detected = 5L, slope = 1, r_squared = 0.95)
  fs <- flag_samples(dds, within_group_z = NULL, spike = sp,
                     min_spike_detected = 6, dose_r2_min = 0.99,
                     dose_slope_min = 2, dose_slope_max = 3)
  expect_true(all(fs$low_spike_detected))                # 5 < 6
  expect_true(all(fs$low_dose_r2))                       # 0.95 < 0.99
  expect_true(all(fs$bad_slope))                         # slope 1 outside [2, 3]
  expect_true(all(fs$flagged))
  expect_true(all(grepl("few detected spikes", fs$reason)))
})

test_that("flag_samples applies a two-sided % spike-in fence from canonical pct_spike", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 4, seed = 23)
  ps <- flag_samples(dds, within_group_z = NULL)$pct_spike
  skip_if(any(!is.finite(ps)))
  fs_hi <- flag_samples(dds, within_group_z = NULL, pct_spike_max = min(ps) - 1e-9)
  fs_lo <- flag_samples(dds, within_group_z = NULL, pct_spike_min = max(ps) + 1e-9)
  expect_true(all(fs_hi$high_spike))                     # everything over the max
  expect_true(all(fs_lo$low_spike))                      # everything under the min
  expect_true(all(grepl("high % spike-in", fs_hi$reason)))
})

test_that("flag_samples flags < 3-point samples and never trips R2/slope on NA fits", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 4, seed = 24)
  n <- ncol(dds)
  few <- rep(c(2L, 5L), length.out = n)                  # alternate <3 vs ok
  sp <- .mock_spike_summary(dds, n_points = few,
                            slope = ifelse(few < 3L, NA_real_, 1),
                            r_squared = ifelse(few < 3L, NA_real_, 0.95))
  fs <- flag_samples(dds, within_group_z = NULL, spike = sp, dose_r2_min = 0.99)
  expect_equal(fs$few_spike_points, few < 3L)            # dose rule on + <3 points
  expect_false(any(fs$low_dose_r2[is.na(sp$r_squared)])) # NA R2 cannot fire
  expect_true(all(grepl("insufficient spike points", fs$reason[fs$few_spike_points])))
})

test_that("flag_samples ignores a spike summary that shares no sample", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 50, n_per_group = 2, n_spike = 2, seed = 25)
  sp_bad <- data.frame(sample = paste0("nope", 1:3), pct_spike = NA_real_,
                       n_spike_detected = 1L, n_points = 5L, slope = 1,
                       r_squared = 0.99, lod = 1)
  fs <- flag_samples(dds, within_group_z = NULL, spike = sp_bad, min_spike_detected = 99)
  expect_false(any(fs$low_spike_detected))               # non-overlapping -> ignored
  expect_true(all(is.na(fs$n_spike_detected)))
})

test_that("suggest_sample_thresholds adds spike suggestions only when given a spike df", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 4, n_spike = 6, seed = 26)
  th0 <- suggest_sample_thresholds(dds)
  expect_false(any(c("pct_spike_min", "dose_slope_min") %in% names(th0)))
  sp <- spike_dose_response(dds, source = "column")$per_sample
  th <- suggest_sample_thresholds(dds, spike = sp)
  expect_true(all(c("pct_spike_min", "pct_spike_max", "dose_r2_min",
                    "dose_slope_min", "dose_slope_max") %in% names(th)))
  expect_equal(th$dose_slope_min, 0.5)                   # fixed advisory range
  expect_equal(th$dose_slope_max, 1.5)
  expect_false("lod" %in% names(th))                     # no portable lod default
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

# ---- restore_samples / restore_features (reset removal) ---------------------

test_that("restore_features brings back removed features, keeping kept-feature edits", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 2, n_spike = 3, seed = 3))
  keep_edit <- rownames(orig)[10]
  w <- set_feature_class(orig, keep_edit, "exogenous")
  drop <- rownames(orig)[1:5]
  w <- drop_features(w, drop)
  expect_equal(nrow(w), nrow(orig) - 5L)
  out <- restore_features(w, orig)
  expect_setequal(rownames(out), rownames(orig))         # all features back
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(out, "counts"))),
               unname(as.matrix(SummarizedExperiment::assay(orig, "counts"))[, colnames(out)]))
  rd <- SummarizedExperiment::rowData(out)
  expect_equal(as.character(rd[keep_edit, "feature_class"]), "exogenous")  # kept edit
  expect_equal(as.character(rd[drop[1], "feature_class"]),                  # re-added = original
               as.character(SummarizedExperiment::rowData(orig)[drop[1], "feature_class"]))
})

test_that("restore_samples brings back removed samples, keeping edits + feature filtering", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 1, seed = 4))
  w <- add_meta_column(orig, "colData", "rin", "numeric", 8)
  w <- drop_features(w, rownames(w)[1:4])                 # also filter features
  drop_s <- colnames(w)[1:2]
  w <- drop_samples(w, drop_s)
  out <- restore_samples(w, orig)
  expect_equal(ncol(out), ncol(orig))                     # all samples back
  expect_equal(nrow(out), nrow(w))                        # feature filtering preserved
  expect_true("rin" %in% colnames(SummarizedExperiment::colData(out)))     # kept edit
  expect_true(is.na(SummarizedExperiment::colData(out)[drop_s[1], "rin"]))  # re-added -> NA
})

test_that("restore_* are no-ops when nothing was removed; drop |> restore round-trips dims", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 2, n_spike = 1, seed = 5))
  expect_equal(dim(restore_features(orig, orig)), dim(orig))
  expect_equal(dim(restore_samples(orig, orig)), dim(orig))
  rt <- restore_features(drop_features(orig, rownames(orig)[1:3]), orig)
  expect_setequal(rownames(rt), rownames(orig))
})

test_that("restore_features preserves the normalized assay set (not counts-only)", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 1, seed = 6))
  orig <- add_normalized_assays(orig, "CPM")
  w <- drop_features(orig, rownames(orig)[1:5])
  out <- restore_features(w, orig)
  expect_true(all(c("counts", "logcounts", "CPM") %in% SummarizedExperiment::assayNames(out)))
  expect_setequal(SummarizedExperiment::assayNames(out), SummarizedExperiment::assayNames(w))
})

test_that("restore_samples does not error when the design references a post-load column", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 1, seed = 7))
  w <- add_meta_column(orig, "colData", "lane", "factor", "L1")    # post-load factor column
  w <- edit_meta_cell(w, "colData", colnames(w)[2], "lane", "L2")  # give it 2 levels
  DESeq2::design(w) <- ~ lane                                      # re-added rows -> NA for `lane`
  w <- drop_samples(w, colnames(w)[1])
  expect_no_error(out <- restore_samples(w, orig))                 # falls back to ~1, no throw
  expect_equal(ncol(out), ncol(orig))                             # all samples back
})

test_that(".merge_meta widens narrowed factor levels instead of NA-ing re-added rows", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 1, seed = 8))
  # 'condition' has both levels in original; drop_samples droplevels() may narrow them.
  treated <- colnames(orig)[SummarizedExperiment::colData(orig)$condition == "treated"]
  w <- drop_samples(orig, treated)                                 # working keeps only 'control'
  out <- restore_samples(w, orig)
  cd <- SummarizedExperiment::colData(out)
  expect_false(anyNA(cd$condition))                                # re-added 'treated' not NA-ed
  expect_setequal(as.character(cd$condition), c("control", "treated"))
})
