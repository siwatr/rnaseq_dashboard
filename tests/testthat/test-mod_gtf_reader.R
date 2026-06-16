gtf_demo <- function() system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")

test_that("mod_gtf_reader reads, trims columns/types on Confirm, frees on Confirm/Reset", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    expect_false(is.null(raw()))
    expect_setequal(gtf_feature_types(raw()), c("exon", "gene", "transcript"))

    # Keep only gene rows and the gene_name column.
    session$setInputs(keep_types = "gene", keep_cols = "gene_name", confirm = 1)
    g <- confirmed()
    expect_false(is.null(g))
    expect_setequal(as.character(g$type), "gene")          # rows trimmed to chosen type
    expect_null(raw())                                     # full parse released
    md <- colnames(S4Vectors::mcols(g))
    expect_true(all(c("type", "gene_id", "gene_name", "transcript_id") %in% md))  # protected kept
    expect_false("gene_biotype" %in% md)                   # not selected, not protected -> dropped

    session$setInputs(reset = 1)
    expect_null(confirmed())
  })
})

test_that("mod_gtf_reader keeps all columns/types when selections are empty", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    session$setInputs(keep_types = character(0), keep_cols = character(0), confirm = 1)
    g <- confirmed()
    expect_setequal(as.character(unique(g$type)), c("exon", "gene", "transcript"))
    expect_true("gene_biotype" %in% colnames(S4Vectors::mcols(g)))
  })
})
