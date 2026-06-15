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
