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
