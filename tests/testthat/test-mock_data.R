test_that("make_mock_dds builds a valid, convention-following DESeqDataSet", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  dds <- make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 8, seed = 42)

  expect_s4_class(dds, "DESeqDataSet")
  # 120 endogenous + 8 spike-in + 1 exogenous = 129 features; 2 x 3 = 6 samples.
  expect_equal(ncol(dds), 6L)
  expect_equal(nrow(dds), 120L + 8L + 1L)

  # Raw counts are integer and non-negative.
  cts <- SummarizedExperiment::assay(dds, "counts")
  expect_true(is.integer(cts))
  expect_true(all(cts >= 0))

  rd <- SummarizedExperiment::rowData(dds)
  # feature_class always present with all three roles represented.
  expect_true(all(c("endogenous", "spike_in", "exogenous") %in% as.character(rd$feature_class)))
  expect_equal(sum(rd$feature_class == "spike_in"), 8L)
  expect_equal(sum(rd$feature_class == "exogenous"), 1L)

  # Required rowData columns.
  expect_true(all(c("gene_name", "feature_length", "chromosome") %in% colnames(rd)))
  expect_true(is.numeric(rd$feature_length))
  expect_true(all(rd$feature_length[rd$feature_class == "endogenous"] > 0))

  # Mitochondrial genes exist (for % mito QC) and ERCC ids are spike-ins.
  expect_true(any(rd$chromosome == "MT"))
  expect_true(all(grepl("^ERCC-", rownames(dds)[rd$feature_class == "spike_in"])))
  # Spike-in concentrations are known, endogenous ones are NA.
  expect_true(all(!is.na(rd$spike_concentration[rd$feature_class == "spike_in"])))

  # colData carries the design factor and an auto-suggestable pseudobulk column.
  cd <- SummarizedExperiment::colData(dds)
  expect_true(all(c("condition", "group", "bio_rep") %in% colnames(cd)))
  expect_setequal(levels(cd$condition), c("control", "treated"))
})

test_that("make_mock_dds is deterministic for a given seed", {
  skip_if_not_installed("DESeq2")
  a <- SummarizedExperiment::assay(make_mock_dds(seed = 7), "counts")
  b <- SummarizedExperiment::assay(make_mock_dds(seed = 7), "counts")
  expect_identical(a, b)
})
