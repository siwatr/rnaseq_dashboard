test_that("detect_mito_features flags the mitochondrial rows in the fixture", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 100, n_per_group = 2, n_spike = 2, seed = 1)
  mito <- detect_mito_features(dds)
  expect_type(mito, "logical")
  expect_length(mito, nrow(dds))
  # The fixture creates 5 mito features (chromosome "MT", gene_name "mt-Gene*").
  expect_equal(sum(mito), 5L)
  expect_true(all(grepl("^mt-", SummarizedExperiment::rowData(dds)$gene_name[mito])))
})

test_that("qc_per_sample_metrics returns one row per sample with the expected schema", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 100, n_per_group = 3, n_spike = 5, seed = 1)
  m <- qc_per_sample_metrics(dds)
  expect_equal(nrow(m), ncol(dds))
  expect_equal(rownames(m), colnames(dds))
  expect_setequal(colnames(m),
                  c("sample", "library_size", "detected", "pct_mito", "pct_spike"))
  expect_true(all(is.finite(m$library_size)))
  expect_true(all(m$detected >= 0 & m$detected <= nrow(dds)))
  expect_true(all(m$pct_mito  >= 0 & m$pct_mito  <= 100))
  expect_true(all(m$pct_spike >= 0 & m$pct_spike <= 100))
})

test_that("library size and detected match a direct colSums computation", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 2, n_spike = 4, seed = 2)
  m <- qc_per_sample_metrics(dds)
  counts <- as.matrix(SummarizedExperiment::assay(dds, "counts"))
  expect_equal(m$library_size, unname(colSums(counts)))
  expect_equal(m$detected, unname(colSums(counts > 0)))
})

test_that("the scater and base-R metric paths agree", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("scater")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 2, n_spike = 4, seed = 3)
  is_mito  <- detect_mito_features(dds)
  is_spike <- ddsdashboard:::.detect_spike_features(dds)
  s <- ddsdashboard:::.qc_metrics_scater(dds, is_mito, is_spike)
  b <- ddsdashboard:::.qc_metrics_base(dds, is_mito, is_spike)
  expect_equal(unname(s$library_size), unname(b$library_size))
  expect_equal(unname(s$detected), unname(b$detected))
  expect_equal(unname(s$pct_mito), unname(b$pct_mito), tolerance = 1e-8)
  expect_equal(unname(s$pct_spike), unname(b$pct_spike), tolerance = 1e-8)
})

test_that("percentages are zero when a feature subset is empty", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 50, n_per_group = 2, n_spike = 0, seed = 4)
  # Drop the mito signal so both subsets are empty.
  SummarizedExperiment::rowData(dds)$chromosome <- "1"
  SummarizedExperiment::rowData(dds)$gene_name <-
    sub("^mt-", "", SummarizedExperiment::rowData(dds)$gene_name)
  m <- qc_per_sample_metrics(dds)
  expect_true(all(m$pct_mito == 0))
  expect_true(all(m$pct_spike == 0))
})

# ---- Dataset-level diagnostics (P3b) ---------------------------------------

test_that("qc_vst returns a DESeqTransform on the small mock (vst fallback)", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 100, n_per_group = 3, n_spike = 3, seed = 1)
  v <- qc_vst(dds)
  expect_s4_class(v, "DESeqTransform")
  expect_equal(ncol(SummarizedExperiment::assay(v)), ncol(dds))
})

test_that("qc_sample_correlation is a symmetric n x n matrix with unit diagonal", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 2, n_spike = 2, seed = 2))
  cm <- qc_sample_correlation(dds, method = "spearman")
  expect_equal(dim(cm), c(ncol(dds), ncol(dds)))
  expect_equal(unname(diag(cm)), rep(1, ncol(dds)), tolerance = 1e-8)
  expect_equal(cm, t(cm), tolerance = 1e-8)
  # method argument is honored (Pearson differs from Spearman in general).
  cp <- qc_sample_correlation(dds, method = "pearson")
  expect_false(isTRUE(all.equal(cm, cp)))
})

test_that("qc_sample_correlation falls back to log2(counts+1) without logcounts", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 1, seed = 3)
  expect_false("logcounts" %in% SummarizedExperiment::assayNames(dds))
  cm <- qc_sample_correlation(dds)
  expect_equal(dim(cm), c(ncol(dds), ncol(dds)))
})

test_that("qc_rle_matrix is median-centered per gene", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 1, seed = 4))
  rle <- qc_rle_matrix(dds)
  expect_equal(dim(rle), c(nrow(dds), ncol(dds)))
  expect_true(all(abs(apply(rle, 1, stats::median)) < 1e-8))
})

test_that("qc_expression_long has one row per feature x sample", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 2, n_spike = 2, seed = 5))
  long <- qc_expression_long(dds)
  expect_equal(nrow(long), nrow(dds) * ncol(dds))
  expect_setequal(colnames(long), c("sample", "value"))
  expect_setequal(levels(long$sample), colnames(dds))
})

test_that("qc_annotation_colors maps discrete and continuous columns, stably", {
  df <- data.frame(
    condition = factor(c("control", "treated", "control", "treated")),
    score     = c(1.5, 2.0, 3.0, 4.0)
  )
  cols <- qc_annotation_colors(df)
  expect_setequal(names(cols), c("condition", "score"))
  # Discrete -> named colour vector, one entry per level.
  expect_setequal(names(cols$condition), c("control", "treated"))
  expect_type(cols$condition, "character")
  # Continuous -> a colour-mapping function (circlize::colorRamp2).
  skip_if_not_installed("circlize")
  expect_type(cols$score, "closure")
  # Deterministic across calls (the whole point - ComplexHeatmap randomizes).
  expect_identical(cols$condition, qc_annotation_colors(df)$condition)
})

test_that("qc_annotation_colors returns NULL for NULL/empty input", {
  expect_null(qc_annotation_colors(NULL))
  expect_null(qc_annotation_colors(data.frame()))
})
