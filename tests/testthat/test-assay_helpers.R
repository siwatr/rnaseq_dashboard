test_that("has_feature_length requires a complete positive length", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)
  expect_true(has_feature_length(dds))
  SummarizedExperiment::rowData(dds)$feature_length[1] <- NA
  expect_false(has_feature_length(dds))
})

test_that("add_normalized_assays adds CPM always and TPM only with length", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)
  dds2 <- add_normalized_assays(dds, c("CPM", "TPM"))
  expect_true(all(c("CPM", "TPM") %in% SummarizedExperiment::assayNames(dds2)))

  SummarizedExperiment::rowData(dds)$feature_length <- NA_real_
  expect_warning(d3 <- add_normalized_assays(dds, c("CPM", "TPM")), "TPM skipped")
  expect_true("CPM" %in% SummarizedExperiment::assayNames(d3))
  expect_false("TPM" %in% SummarizedExperiment::assayNames(d3))
})

test_that("estimate_size_factors_endogenous sets size factors", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 100, n_per_group = 2, n_spike = 2, seed = 1)
  dds <- estimate_size_factors_endogenous(dds)
  expect_false(is.null(DESeq2::sizeFactors(dds)))
  expect_length(DESeq2::sizeFactors(dds), ncol(dds))
})

test_that("refresh_assays recomputes present assays from counts", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)
  dds <- add_normalized_assays(dds, "CPM")
  SummarizedExperiment::assay(dds, "CPM")[] <- 0
  dds <- refresh_assays(dds)
  expect_gt(max(SummarizedExperiment::assay(dds, "CPM")), 0)
})
