test_that("Feature-info OrgDb annotation populates a name column", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("org.Mm.eg.db")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    counts <- data.frame(id = c("ENSMUSG00000000001", "ENSMUSG00000000003"),
                         S1 = c(5L, 2L), S2 = c(3L, 7L))
    samples <- data.frame(sample = c("S1", "S2"), cond = c("a", "b"))
    state_load(state, ensure_feature_class(tabular_to_dds(counts, samples)),
               source = "tabular", meta = list(feature_type = "gene"))
    session$flushReact()
    v0 <- state$data_version
    session$setInputs(organism = "mouse", id_type = "ensembl")
    session$setInputs(annotate = 1)
    expect_equal(state$data_version, v0 + 1L)
    rd <- SummarizedExperiment::rowData(state$working)
    expect_true("gene_name" %in% colnames(rd))
    expect_equal(as.character(rd$gene_name[1]), "Gnai3")
  })
})

test_that("Feature-info GTF annotation sets feature_length and bumps version", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicRanges")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    ids <- c("ENSG00000000001", "ENSG00000000002", "ENSG00000000003")
    counts <- matrix(c(5L, 2L, 9L, 3L, 7L, 1L), nrow = 3, dimnames = list(ids, c("S1", "S2")))
    samples <- data.frame(cond = c("a", "b"), row.names = c("S1", "S2"))
    dds <- ensure_feature_class(DESeq2::DESeqDataSetFromMatrix(counts, samples, ~ 1))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    # Inject the parsed GTF directly (skip the fileInput upload path).
    gtf_obj(import_gtf(system.file("extdata", "demo_annotation.gtf", package = "ddsdashboard")))
    v0 <- state$data_version
    session$setInputs(gtf_match = "auto", gtf_import = c("gene_name", "seqnames"),
                      gtf_len = TRUE, gtf_len_type = "exon", apply_gtf = 1)
    expect_equal(state$data_version, v0 + 1L)
    rd <- SummarizedExperiment::rowData(state$working)
    expect_equal(as.numeric(rd$feature_length), c(202, 201, 501))
    expect_true(has_feature_length(state$working))
  })
})
