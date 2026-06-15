test_that("Assay module adds CPM, sets size factors, and bumps data_version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_assay_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 100, n_per_group = 2, n_spike = 2, seed = 1),
               source = "demo")
    session$flushReact()
    v0 <- state$data_version
    session$setInputs(assays = "CPM")
    session$setInputs(apply = 1)
    expect_equal(state$data_version, v0 + 1L)
    expect_true("CPM" %in% SummarizedExperiment::assayNames(state$working))
    expect_false(is.null(DESeq2::sizeFactors(state$working)))
  })
})
