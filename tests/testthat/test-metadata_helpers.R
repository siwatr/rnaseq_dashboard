mk_dds <- function() {
  counts <- data.frame(id = c("Gene1", "ERCC-00001", "Gene2"),
                       S1 = c(5L, 2L, 9L), S2 = c(3L, 7L, 1L), S3 = c(4L, 4L, 4L))
  samples <- data.frame(sample = c("S1", "S2", "S3"),
                        cond = c("a", "b", "a"), score = c(1.5, 2.5, 3.5))
  tabular_to_dds(counts, samples)
}

test_that("edit_coldata_cell coerces to the column type and rejects bad values", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  dds <- edit_coldata_cell(dds, "S1", "score", "5.5")
  expect_equal(SummarizedExperiment::colData(dds)$score[1], 5.5)
  expect_error(edit_coldata_cell(dds, "S1", "score", "abc"), "number")
  expect_error(edit_coldata_cell(dds, "S1", "nope", "x"), "Unknown colData column")
  expect_error(edit_coldata_cell(dds, "SX", "score", "1"), "Unknown sample row")
})

test_that("edit_coldata_cell adds a new factor level when needed", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)
  dds <- edit_coldata_cell(dds, 1, "condition", "rescue")  # new level
  cond <- SummarizedExperiment::colData(dds)$condition
  expect_true("rescue" %in% levels(cond))
  expect_equal(as.character(cond[1]), "rescue")
})

test_that("merge_sample_metadata joins by id, NA-fills, and reports", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  sheet <- data.frame(sample = c("S1", "S2"), batch = c("A", "B"))
  res <- merge_sample_metadata(dds, sheet)
  cd <- SummarizedExperiment::colData(res$dds)
  expect_true("batch" %in% colnames(cd))
  expect_equal(as.character(cd$batch), c("A", "B", NA))   # S3 absent -> NA
  expect_equal(res$report$matched, 2L)
  expect_equal(res$report$unmatched_in_data, "S3")
})

test_that("merge_sample_metadata supports an explicit id column and errors on no match", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  sheet <- data.frame(libname = c("S2", "S3"), batch = c("X", "Y"))
  res <- merge_sample_metadata(dds, sheet, id_col = "libname")
  expect_equal(as.character(SummarizedExperiment::colData(res$dds)$batch), c(NA, "X", "Y"))
  expect_error(merge_sample_metadata(dds, data.frame(sample = c("Z1", "Z2"), x = 1:2)),
               "match")
})

test_that("set_feature_class tags resolved ids and errors when none match", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  dds <- set_feature_class(dds, "Gene2", "exogenous")
  fc <- SummarizedExperiment::rowData(dds)$feature_class
  names(fc) <- rownames(dds)
  expect_equal(as.character(fc["Gene2"]), "exogenous")
  expect_equal(as.character(fc["ERCC-00001"]), "spike_in")   # auto from ensure_feature_class
  expect_equal(as.character(fc["Gene1"]), "endogenous")
  expect_error(set_feature_class(dds, "nope", "exogenous"), "None of the given")
})
