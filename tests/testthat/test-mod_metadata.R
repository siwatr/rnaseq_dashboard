test_that("Sample 'Additional Metadata' binds an uploaded sheet onto the draft until Save", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_metadata_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 10, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$flushReact()
    ids <- colnames(state$working)
    tmp <- tempfile(fileext = ".csv")
    utils::write.csv(
      data.frame(sample = ids, batch = rep(c("A", "B"), length.out = length(ids))),
      tmp, row.names = FALSE)

    session$setInputs(sheet = list(datapath = tmp, name = "sheet.csv"), id_col = "sample")
    v0 <- state$data_version
    session$setInputs(merge = 1)

    # Binds onto the draft, not state$working, until Save.
    expect_equal(state$data_version, v0)
    cd <- SummarizedExperiment::colData(editor$draft())
    expect_true("batch" %in% colnames(cd))
    expect_setequal(as.character(cd$batch), c("A", "B"))
  })
})
