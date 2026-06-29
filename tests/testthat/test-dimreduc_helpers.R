# PCA helpers (R/dimreduc_helpers.R) for Phase 4.

test_that("pca_assay_advice tiers each assay correctly", {
  expect_equal(pca_assay_advice("vst")$tier, "recommended")
  expect_equal(pca_assay_advice("logcounts")$tier, "recommended")
  expect_equal(pca_assay_advice("norm_logcounts")$tier, "recommended")
  expect_false(pca_assay_advice("vst")$recommend_log)
  for (a in c("CPM", "TPM", "FPKM")) {
    ad <- pca_assay_advice(a)
    expect_equal(ad$tier, "log_first")
    expect_true(ad$recommend_log)
    expect_match(ad$msg, "log2")
  }
  expect_equal(pca_assay_advice("counts")$tier, "unsuitable")
  expect_false(pca_assay_advice("counts")$recommend_log)
  expect_match(pca_assay_advice("counts")$msg, "depth")
  expect_equal(pca_assay_advice("some_other")$tier, "log_first")   # unknown -> cautious
})

test_that("pca_input: VST is endogenous-only, returns a matrix (not a DESeqTransform)", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 8, seed = 1))
  res <- pca_input(dds, "vst")
  expect_true(is.matrix(res$mat))
  expect_equal(res$label, "VST (blind)")
  n_endo <- sum(SummarizedExperiment::rowData(dds)$feature_class == "endogenous")
  expect_equal(nrow(res$mat), n_endo)              # endogenous-only, not re-subset
  expect_equal(ncol(res$mat), ncol(dds))
})

test_that("pca_input: VST failure falls back to logcounts with an honest label", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 2, seed = 2))
  # Force qc_vst() to error so the fallback path is exercised.
  local_mocked_bindings(qc_vst = function(...) stop("boom"))
  res <- pca_input(dds, "vst")
  expect_match(res$label, "VST unavailable")
  expect_equal(nrow(res$mat),
               sum(SummarizedExperiment::rowData(dds)$feature_class == "endogenous"))
})

test_that("pca_input: norm_logcounts estimates size factors and log2s normalized counts", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 2, n_spike = 0, seed = 3)
  expect_null(DESeq2::sizeFactors(dds))
  res <- pca_input(dds, "norm_logcounts")
  expect_equal(res$label, "normalized log-counts (log2)")
  expect_true(is.matrix(res$mat) && all(is.finite(res$mat)))
})

test_that("pca_input: stored assay + log_transform applies log2 and labels it", {
  skip_if_not_installed("DESeq2")
  dds <- add_normalized_assays(
    ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 2, seed = 4)), "CPM")
  raw <- pca_input(dds, "CPM", log_transform = FALSE)
  lg  <- pca_input(dds, "CPM", log_transform = TRUE)
  expect_equal(raw$label, "CPM")
  expect_equal(lg$label, "CPM (log2)")
  expect_equal(lg$mat, log2(raw$mat + 1))
})

test_that("top_variable_features drops zero-variance rows and clamps n_top", {
  m <- rbind(g1 = c(1, 5, 2, 9), g2 = c(3, 3, 3, 3), g3 = c(0, 2, 8, 1), g4 = c(7, 1, 4, 2))
  v <- top_variable_features(m, n_top = 10)
  expect_false("g2" %in% v)                 # constant row dropped
  expect_equal(length(v), 3L)               # clamped to available (4 - 1 constant)
  expect_equal(v[1], names(which.max(apply(m[c("g1","g3","g4"), ], 1, var))))
})

test_that("compute_pca matches a hand prcomp, sums variance to 100, deterministic sign", {
  set.seed(42)
  m <- matrix(rnorm(50 * 6), nrow = 50, dimnames = list(paste0("g", 1:50), paste0("s", 1:6)))
  pc <- compute_pca(m, n_top = 30)
  expect_equal(rownames(pc$scores), colnames(m))
  expect_equal(sum(pc$var_pct), 100, tolerance = 1e-6)
  expect_equal(pc$n_genes, 30L)
  # Deterministic sign -> identical across calls (no mirror flip).
  expect_equal(compute_pca(m, n_top = 30)$scores, pc$scores)
  # Variance magnitudes match a plain prcomp on the same top genes.
  top <- top_variable_features(m, 30)
  ref <- stats::prcomp(t(m[top, ]), center = TRUE, scale. = FALSE)
  expect_equal(unname(pc$var_pct), unname(ref$sdev^2 / sum(ref$sdev^2) * 100), tolerance = 1e-6)
})

test_that("compute_pca errors when too few variable features", {
  m <- matrix(1, nrow = 3, ncol = 4, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  expect_error(compute_pca(m), "variable features")     # all constant -> nothing to PCA
})

test_that("expr_transform applies log2/log10 with a pseudocount, none is identity", {
  x <- c(0, 1, 9, 99)
  expect_equal(expr_transform(x, "none"), x)
  expect_equal(expr_transform(x, "log2", 1), log2(x + 1))
  expect_equal(expr_transform(x, "log10", 1), log10(x + 1))
  expect_equal(expr_transform(c(0, 0.5), "log2", 0.5), log2(c(0.5, 1)))
  expect_error(expr_transform(x, "ln"))                 # match.arg guards the set
})
