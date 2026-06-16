test_that("detect_id_type classifies common id shapes", {
  expect_equal(detect_id_type(c("ENSMUSG00000000001", "ENSMUSG00000000002")), "ensembl")
  expect_equal(detect_id_type(c("100", "200", "300")), "entrez")
  expect_equal(detect_id_type(c("Actb", "Gapdh")), "symbol")
})

mk_ens_dds <- function(ids) {
  counts <- data.frame(id = ids,
                       S1 = seq_along(ids), S2 = rev(seq_along(ids)))
  samples <- data.frame(sample = c("S1", "S2"), cond = c("a", "b"))
  ensure_feature_class(tabular_to_dds(counts, samples))
}

test_that("annotate_with_orgdb maps Ensembl ids to symbol/chromosome/description", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Mm.eg.db")
  dds <- mk_ens_dds(c("ENSMUSG00000000001", "ENSMUSG00000000003"))
  dds <- annotate_with_orgdb(dds, "mouse", id_type = "ensembl", feature_type = "gene")
  rd <- SummarizedExperiment::rowData(dds)
  expect_true(all(c("gene_name", "description") %in% colnames(rd)))
  expect_equal(as.character(rd$gene_name[1]), "Gnai3")   # ENSMUSG00000000001
})

test_that("annotate_with_orgdb tolerates unmapped ids", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Mm.eg.db")
  dds <- mk_ens_dds(c("NOTAGENE1", "NOTAGENE2"))
  dds <- annotate_with_orgdb(dds, "mouse", id_type = "ensembl")
  expect_true(all(is.na(SummarizedExperiment::rowData(dds)$gene_name)))
})

test_that("annotate_with_orgdb can flag mapped features", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Mm.eg.db")
  dds <- mk_ens_dds(c("ENSMUSG00000000001", "NOTAGENE2"))
  dds <- annotate_with_orgdb(dds, "mouse", id_type = "ensembl", matched_col = "in_orgdb")
  fl <- SummarizedExperiment::rowData(dds)$in_orgdb
  expect_type(fl, "logical")
  expect_equal(fl, c(TRUE, FALSE))   # first maps (Gnai3), second does not
})

test_that("orgdb_annotation_preview returns a small non-committing mapping table", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Mm.eg.db")
  dds <- mk_ens_dds(c("ENSMUSG00000000001", "NOTAGENE2"))
  pv <- orgdb_annotation_preview(dds, "mouse", id_type = "ensembl", feature_type = "gene", n = 20)
  expect_true(all(c("id", "gene_name", "description", "in_orgdb") %in% colnames(pv)))
  expect_equal(nrow(pv), 2)
  expect_equal(pv$in_orgdb, c(TRUE, FALSE))   # first maps, second does not
  # non-committing: the dds itself is untouched
  expect_false("gene_name" %in% colnames(SummarizedExperiment::rowData(dds)))
})

test_that("annotation_overwrites flags only existing, populated targets", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 6, n_per_group = 2, n_spike = 1, seed = 1)  # gene_name populated
  expect_setequal(annotation_overwrites(dds, c("gene_name", "description", "feature_class")),
                  c("gene_name", "feature_class"))   # description absent -> not flagged
  SummarizedExperiment::rowData(dds)$empty <- NA_character_
  expect_false("empty" %in% annotation_overwrites(dds, "empty"))   # all-NA -> not flagged
})

test_that("annotation_coverage counts endogenous matches", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)  # has gene_name
  cov <- annotation_coverage(dds, "gene_name")
  expect_equal(cov$total, sum(SummarizedExperiment::rowData(dds)$feature_class == "endogenous"))
  expect_true(cov$matched >= 1)
})
