mk_dds <- function() {
  counts <- data.frame(id = c("Gene1", "ERCC-00001", "Gene2"),
                       S1 = c(5L, 2L, 9L), S2 = c(3L, 7L, 1L), S3 = c(4L, 4L, 4L))
  samples <- data.frame(sample = c("S1", "S2", "S3"),
                        cond = c("a", "b", "a"), score = c(1.5, 2.5, 3.5))
  ensure_feature_class(tabular_to_dds(counts, samples))
}
mk_mock <- function() make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)

test_that("edit_meta_cell (colData) coerces type and rejects bad values", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  dds <- edit_meta_cell(dds, "colData", "S1", "score", "5.5")
  expect_equal(SummarizedExperiment::colData(dds)$score[1], 5.5)
  expect_error(edit_meta_cell(dds, "colData", "S1", "score", "abc"), "number")
  expect_error(edit_meta_cell(dds, "colData", "S1", "nope", "x"), "Unknown colData column")
})

test_that("edit_meta_cell (rowData) validates the feature_class value set", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  dds <- edit_meta_cell(dds, "rowData", "Gene1", "feature_class", "exogenous")
  fc <- SummarizedExperiment::rowData(dds)$feature_class
  names(fc) <- rownames(dds)
  expect_equal(as.character(fc["Gene1"]), "exogenous")
  expect_error(edit_meta_cell(dds, "rowData", "Gene1", "feature_class", "bogus"), "must be one of")
})

test_that("protected_columns differ by slot", {
  skip_if_not_installed("DESeq2")
  m <- mk_mock()
  expect_equal(protected_columns(m, "colData"), "condition")   # design var
  expect_equal(protected_columns(m, "rowData"), "feature_class")
})

test_that("add_meta_column works for both slots", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  dds <- add_meta_column(dds, "colData", "batch", "character", "A")
  dds <- add_meta_column(dds, "rowData", "note", "character", "x")
  expect_true("batch" %in% colnames(SummarizedExperiment::colData(dds)))
  expect_true("note" %in% colnames(SummarizedExperiment::rowData(dds)))
  expect_error(add_meta_column(dds, "colData", "batch"), "already exists")
})

test_that("remove_meta_columns is multi and skips protected, reports unknown", {
  skip_if_not_installed("DESeq2")
  m <- mk_mock()  # colData: condition (design), bio_rep, group
  res <- remove_meta_columns(m, "colData", c("bio_rep", "group", "condition", "nope"))
  expect_setequal(res$removed, c("bio_rep", "group"))
  expect_equal(res$skipped, "condition")
  expect_equal(res$unknown, "nope")
  cd <- SummarizedExperiment::colData(res$dds)
  expect_false(any(c("bio_rep", "group") %in% colnames(cd)))
  expect_true("condition" %in% colnames(cd))
})

test_that("rename_meta_column rewrites the design (colData) and protects feature_class (rowData)", {
  skip_if_not_installed("DESeq2")
  m <- mk_mock()
  m2 <- rename_meta_column(m, "colData", "condition", "treatment")
  expect_true("treatment" %in% colnames(SummarizedExperiment::colData(m2)))
  expect_equal(all.vars(DESeq2::design(m2)), "treatment")
  expect_error(rename_meta_column(m, "rowData", "feature_class", "fc"), "cannot be renamed")
  m3 <- rename_meta_column(add_meta_column(m, "rowData", "note", "character", "x"),
                           "rowData", "note", "annotation")
  expect_true("annotation" %in% colnames(SummarizedExperiment::rowData(m3)))
})

test_that("rename_samples enforces uniqueness and existence", {
  skip_if_not_installed("DESeq2")
  dds <- mk_mock()
  s <- colnames(dds)
  expect_equal(colnames(rename_samples(dds, s[1], "sampleA"))[1], "sampleA")
  expect_error(rename_samples(dds, s[1], s[2]), "unique")
  expect_error(rename_samples(dds, "nope", "x"), "Unknown sample")
})

test_that("merge_sample_metadata joins, NA-fills, reports overwritten", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  sheet <- data.frame(sample = c("S1", "S2"), cond = c("X", "Y"), batch = c("A", "B"))
  res <- merge_sample_metadata(dds, sheet)
  cd <- SummarizedExperiment::colData(res$dds)
  expect_equal(as.character(cd$batch), c("A", "B", NA))    # S3 absent -> NA
  expect_true("cond" %in% res$report$overwritten)
  expect_equal(res$report$unmatched_in_data, "S3")
})

test_that("set_feature_class tags resolved ids", {
  skip_if_not_installed("DESeq2")
  dds <- mk_dds()
  dds <- set_feature_class(dds, "Gene2", "exogenous")
  fc <- SummarizedExperiment::rowData(dds)$feature_class
  names(fc) <- rownames(dds)
  expect_equal(as.character(fc["Gene2"]), "exogenous")
  expect_equal(as.character(fc["ERCC-00001"]), "spike_in")
  expect_error(set_feature_class(dds, "nope", "exogenous"), "None of the given")
})

# ---- reset_metadata_slot (slot-scoped "reset to original") ------------------

test_that("reset_metadata_slot reverts colData to original, keeping current samples", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 1, seed = 1))
  w <- add_meta_column(orig, "colData", "rin", "numeric", 7)
  w <- edit_meta_cell(w, "colData", colnames(w)[3], "rin", "9")
  w <- drop_samples(w, colnames(w)[1:2])                 # filtering must NOT be undone
  out <- reset_metadata_slot(w, orig, "colData")
  expect_equal(ncol(out), ncol(w))                       # samples unchanged
  expect_false("rin" %in% colnames(SummarizedExperiment::colData(out)))  # added col dropped
  expect_setequal(colnames(SummarizedExperiment::colData(out)),
                  colnames(SummarizedExperiment::colData(orig)))
})

test_that("reset_metadata_slot reverts rowData edits and leaves the other slot intact", {
  skip_if_not_installed("DESeq2")
  orig <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 2, seed = 2))
  w <- set_feature_class(orig, rownames(orig)[1], "exogenous")        # rowData edit
  w <- add_meta_column(w, "colData", "batch", "character", "b1")      # unrelated colData edit
  out <- reset_metadata_slot(w, orig, "rowData")
  fc <- as.character(SummarizedExperiment::rowData(out)$feature_class)
  expect_equal(fc[1], as.character(SummarizedExperiment::rowData(orig)$feature_class)[1])
  expect_true("batch" %in% colnames(SummarizedExperiment::colData(out)))  # other slot kept
})
