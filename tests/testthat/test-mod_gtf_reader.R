gtf_demo <- function() system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")

test_that("mod_gtf_reader reads, trims columns/types on Apply filter, frees the parse", {
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

    session$setInputs(remove = 1)
    expect_null(confirmed())
  })
})

test_that("mod_gtf_reader filtering persists: a second Apply narrows the trimmed object", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    # First apply keeps everything; the parse is freed but filtering stays live.
    session$setInputs(keep_types = character(0), keep_cols = character(0), confirm = 1)
    expect_setequal(as.character(unique(confirmed()$type)), c("exon", "gene", "transcript"))
    expect_null(raw())
    # A second apply works on the already-trimmed object and narrows it further.
    session$setInputs(keep_types = "gene", confirm = 2)
    expect_setequal(as.character(unique(confirmed()$type)), "gene")
  })
})

test_that("reader UI can split controls from the preview output", {
  with_preview    <- as.character(mod_gtf_reader_ui("gtf", preview = TRUE))
  without_preview  <- as.character(mod_gtf_reader_ui("gtf", preview = FALSE))
  preview_only     <- as.character(mod_gtf_reader_preview_ui("gtf"))
  expect_true(any(grepl("gtf-preview_ui", with_preview)))
  expect_false(any(grepl("gtf-preview_ui", without_preview)))  # controls only
  expect_true(any(grepl("gtf-preview_ui", preview_only)))      # preview alone
  expect_true(any(grepl("gtf-read", without_preview)))         # controls still present
})

test_that("mod_gtf_reader Remove frees the confirmed object", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    session$setInputs(keep_types = character(0), keep_cols = character(0), confirm = 1)
    expect_false(is.null(confirmed()))
    session$setInputs(remove = 1)
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
