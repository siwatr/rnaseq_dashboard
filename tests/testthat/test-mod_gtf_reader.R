gtf_demo <- function() system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")

test_that("mod_gtf_reader reads into a usable object and trims on Apply filter", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    # Loaded GTF is immediately available (annotation works before any filtering).
    expect_false(is.null(gtf()))
    expect_setequal(gtf_feature_types(gtf()), c("exon", "gene", "transcript"))

    # Keep only gene rows and the gene_name column.
    session$setInputs(keep_types = "gene", keep_cols = "gene_name", confirm = 1)
    g <- gtf()
    expect_setequal(as.character(g$type), "gene")          # rows trimmed to chosen type
    md <- colnames(S4Vectors::mcols(g))
    expect_true(all(c("type", "gene_id", "gene_name", "transcript_id") %in% md))  # protected kept
    expect_false("gene_biotype" %in% md)                   # not selected, not protected -> dropped

    session$setInputs(remove = 1)
    expect_null(gtf())
  })
})

test_that("mod_gtf_reader filtering persists: a second Apply narrows the trimmed object", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    # First apply keeps everything; filtering stays live afterwards.
    session$setInputs(keep_types = character(0), keep_cols = character(0), confirm = 1)
    expect_setequal(as.character(unique(gtf()$type)), c("exon", "gene", "transcript"))
    # A second apply works on the already-trimmed object and narrows it further.
    session$setInputs(keep_types = "gene", confirm = 2)
    expect_setequal(as.character(unique(gtf()$type)), "gene")
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

test_that("mod_gtf_reader Remove frees the loaded object", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    expect_false(is.null(gtf()))
    session$setInputs(remove = 1)
    expect_null(gtf())
  })
})

test_that("mod_gtf_reader keeps all columns/types when filter selections are empty", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  shiny::testServer(mod_gtf_reader_server, {
    session$setInputs(path = gtf_demo(), read = 1)
    session$setInputs(keep_types = character(0), keep_cols = character(0), confirm = 1)
    g <- gtf()
    expect_setequal(as.character(unique(g$type)), c("exon", "gene", "transcript"))
    expect_true("gene_biotype" %in% colnames(S4Vectors::mcols(g)))
  })
})
