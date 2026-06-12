small_counts <- function() {
  data.frame(id = c("Gene1", "ERCC-00001"), S1 = c(5L, 2L), S2 = c(3L, 7L))
}
small_samples <- function() {
  data.frame(sample = c("S1", "S2"), condition = c("a", "b"))
}

test_that("tabular_to_dds builds a dds and matches sample ids", {
  skip_if_not_installed("DESeq2")
  dds <- tabular_to_dds(small_counts(), small_samples())
  expect_s4_class(dds, "DESeqDataSet")
  expect_equal(dim(dds), c(2L, 2L))
  expect_equal(colnames(dds), c("S1", "S2"))
})

test_that("tabular_to_dds errors (listing mismatches) on id mismatch", {
  skip_if_not_installed("DESeq2")
  bad <- data.frame(sample = c("S1", "SX"), condition = c("a", "b"))
  expect_error(tabular_to_dds(small_counts(), bad), "do not match")
})

test_that("ensure_feature_class adds the column, defaults endogenous, flags ERCC", {
  skip_if_not_installed("DESeq2")
  dds <- tabular_to_dds(small_counts(), small_samples())
  expect_false("feature_class" %in% colnames(SummarizedExperiment::rowData(dds)))
  dds <- ensure_feature_class(dds)
  fc <- as.character(SummarizedExperiment::rowData(dds)$feature_class)
  expect_equal(fc, c("endogenous", "spike_in"))
})

test_that("ensure_logcounts adds a logcounts assay once", {
  skip_if_not_installed("DESeq2")
  dds <- tabular_to_dds(small_counts(), small_samples())
  expect_false("logcounts" %in% SummarizedExperiment::assayNames(dds))
  dds <- ensure_logcounts(dds)
  expect_true("logcounts" %in% SummarizedExperiment::assayNames(dds))
  # idempotent
  dds2 <- ensure_logcounts(dds)
  expect_equal(SummarizedExperiment::assayNames(dds2), SummarizedExperiment::assayNames(dds))
})

test_that("as_input_dds passes a DESeqDataSet through unchanged", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 5)
  expect_identical(as_input_dds(dds), dds)
})

test_that("as_input_dds coerces a small sce and rejects a large one", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SingleCellExperiment")
  m <- matrix(rpois(20 * 5, 5), nrow = 20,
              dimnames = list(paste0("g", 1:20), paste0("c", 1:5)))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  dds <- as_input_dds(sce, max_cells = 10)
  expect_s4_class(dds, "DESeqDataSet")
  expect_equal(ncol(dds), 5L)
  expect_error(as_input_dds(sce, max_cells = 3), "pseudobulk")
  expect_true(input_meta(sce, max_cells = 10)$sce_per_cell)
  expect_equal(input_meta(sce, max_cells = 10)$data_type, "single-cell")
})

test_that("detect_feature_type reads <type>_name then Ensembl ids", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)
  ft <- detect_feature_type(dds)        # rowData has gene_name + ENSMUSG ids
  expect_equal(ft$feature_type, "gene")
  expect_true(ft$confident)

  plain <- tabular_to_dds(small_counts(), small_samples())
  expect_equal(detect_feature_type(plain)$feature_type, "feature")
  expect_false(detect_feature_type(plain)$confident)
})
