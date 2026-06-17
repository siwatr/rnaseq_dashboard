gtf_path <- function() system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")

mk_gtf_dds <- function(extra_absent = FALSE) {
  ids <- c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")
  if (extra_absent) ids <- c(ids, "ENSG00000000999")
  n <- length(ids)
  counts <- matrix(seq_len(n * 2L), nrow = n, dimnames = list(ids, c("S1", "S2")))
  samples <- data.frame(cond = c("a", "b"), row.names = c("S1", "S2"))
  ensure_feature_class(DESeq2::DESeqDataSetFromMatrix(counts, samples, ~ 1))
}

test_that("gtf_feature_lengths takes the union per group, by chosen type", {
  skip_if_not_installed("rtracklayer"); skip_if_not_installed("GenomicRanges")
  gtf <- import_gtf(gtf_path())
  le <- gtf_feature_lengths(gtf, "exon")
  expect_equal(le[["ENSG00000000001"]], 202)   # 100-250 (151) + 400-450 (51), overlap collapsed
  expect_equal(le[["ENSG00000000002"]], 201)   # 1000-1200 (touching exons merged)
  expect_equal(le[["ENSG00000000003"]], 501)
  lg <- gtf_feature_lengths(gtf, "gene")
  expect_equal(lg[["ENSG00000000001"]], 351)   # whole span, differs from exon union
  expect_error(gtf_feature_lengths(gtf, "CDS"), "No 'CDS'")
})

test_that("gtf_preview returns at most n rows as a data.frame", {
  skip_if_not_installed("rtracklayer")
  gtf <- import_gtf(gtf_path())
  p <- gtf_preview(gtf, n = 5)
  expect_s3_class(p, "data.frame")
  expect_lte(nrow(p), 5)
  expect_true(all(c("seqnames", "type") %in% colnames(p)))
  expect_equal(nrow(gtf_preview(NULL)), 0)
})

test_that("gtf helpers expose feature types and importable columns", {
  skip_if_not_installed("rtracklayer")
  gtf <- import_gtf(gtf_path())
  expect_setequal(gtf_feature_types(gtf), c("exon", "gene", "transcript"))
  expect_true(all(c("seqnames", "gene_id", "gene_name", "gene_biotype") %in%
                    available_gtf_columns(gtf)))
})

test_that("gtf_attribute_table is one row per group with seqnames + attributes", {
  skip_if_not_installed("rtracklayer"); skip_if_not_installed("GenomicRanges")
  gtf <- import_gtf(gtf_path())
  tab <- gtf_attribute_table(gtf, "gene_id")
  expect_setequal(rownames(tab),
                  c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003"))
  expect_true(all(c("seqnames", "gene_name", "gene_biotype", "type") %in% colnames(tab)))
  expect_equal(as.character(tab["ENSG00000000001", "seqnames"]), "chr1")
  expect_equal(as.character(tab["ENSG00000000001", "gene_name"]), "GeneA")
})

test_that("annotate_with_gtf fills names/chromosome/length and overrides OrgDb", {
  skip_if_not_installed("rtracklayer"); skip_if_not_installed("GenomicRanges")
  gtf <- import_gtf(gtf_path())
  dds <- mk_gtf_dds()
  SummarizedExperiment::rowData(dds)$gene_name <- c("STALE", "STALE", "STALE")
  res <- annotate_with_gtf(dds, gtf, import_cols = c("gene_name", "seqnames", "gene_biotype"),
                           compute_length = TRUE, length_type = "exon")
  rd <- SummarizedExperiment::rowData(res$dds)
  expect_equal(as.character(rd$gene_name), c("GeneA", "GeneB", "GeneC"))   # GTF overrode STALE
  expect_equal(as.character(rd$chromosome), c("chr1", "chr2", "chr3"))
  expect_equal(as.character(rd$gene_biotype), c("protein_coding", "protein_coding", "lincRNA"))
  expect_equal(as.numeric(rd$feature_length), c(202, 201, 501))
  expect_true(res$report$length_complete)
  expect_true(has_feature_length(res$dds))
})

test_that("annotate_with_gtf handles partial coverage without wiping existing values", {
  skip_if_not_installed("rtracklayer"); skip_if_not_installed("GenomicRanges")
  gtf <- import_gtf(gtf_path())
  dds <- mk_gtf_dds(extra_absent = TRUE)
  SummarizedExperiment::rowData(dds)$gene_name <- c("a", "b", "c", "keepme")
  res <- annotate_with_gtf(dds, gtf, import_cols = "gene_name", compute_length = TRUE)
  rd <- SummarizedExperiment::rowData(res$dds)
  expect_equal(res$report$matched, 3L)
  expect_equal(res$report$total, 4L)
  expect_true(is.na(rd$feature_length[4]))                 # unmatched feature -> NA length
  expect_equal(as.character(rd$gene_name[4]), "keepme")    # existing value kept, not wiped
  expect_false(has_feature_length(res$dds))                # incomplete -> TPM/FPKM stay off
  expect_false(res$report$length_complete)
})

test_that("annotate_with_gtf can flag matched features", {
  skip_if_not_installed("rtracklayer"); skip_if_not_installed("GenomicRanges")
  gtf <- import_gtf(gtf_path())
  res <- annotate_with_gtf(mk_gtf_dds(extra_absent = TRUE), gtf,
                           import_cols = NULL, matched_col = "in_gtf")
  fl <- SummarizedExperiment::rowData(res$dds)$in_gtf
  expect_type(fl, "logical")
  expect_equal(fl, c(TRUE, TRUE, TRUE, FALSE))   # 4th gene absent from the GTF
})

test_that("gtf_match_count tallies dds features present in the GTF", {
  skip_if_not_installed("rtracklayer"); skip_if_not_installed("GenomicRanges")
  gtf <- import_gtf(gtf_path())
  cnt <- gtf_match_count(mk_gtf_dds(extra_absent = TRUE), gtf)
  expect_equal(cnt$total, 4)
  expect_equal(cnt$matched, 3)   # 4th gene absent from the GTF
})

test_that("set_feature_length_from_column adopts a numeric column, rejects non-numeric", {
  skip_if_not_installed("DESeq2")
  dds <- mk_gtf_dds()
  SummarizedExperiment::rowData(dds)$mylen <- c(100, 200, 300)
  dds2 <- set_feature_length_from_column(dds, "mylen")
  expect_equal(as.numeric(SummarizedExperiment::rowData(dds2)$feature_length), c(100, 200, 300))
  expect_true(has_feature_length(dds2))
  SummarizedExperiment::rowData(dds)$notnum <- c("a", "b", "c")
  expect_error(set_feature_length_from_column(dds, "notnum"), "not numeric")
})

test_that("GTF feature_length unlocks TPM via add_normalized_assays", {
  skip_if_not_installed("rtracklayer"); skip_if_not_installed("GenomicRanges")
  gtf <- import_gtf(gtf_path())
  res <- annotate_with_gtf(mk_gtf_dds(), gtf, compute_length = TRUE)
  out <- add_normalized_assays(res$dds, c("CPM", "TPM"))
  expect_true(all(c("CPM", "TPM") %in% SummarizedExperiment::assayNames(out)))
})
