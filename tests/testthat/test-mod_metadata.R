test_that("metadata module edits a colData cell and bumps data_version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_metadata_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    v0 <- state$data_version
    # colData column 1 is "condition"; row 1 is the first sample.
    session$setInputs(coldata_cell_edit = list(row = 1L, col = 1L, value = "treated"))
    expect_equal(state$data_version, v0 + 1L)
    expect_equal(as.character(SummarizedExperiment::colData(state$working)$condition[1]),
                 "treated")
  })
})

test_that("metadata module merges an uploaded sample sheet", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_metadata_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    samples <- colnames(state$working)
    tmp <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(sample = samples, batch = rep(c("A", "B"), length.out = length(samples))),
                     tmp, row.names = FALSE)
    session$setInputs(sheet = list(datapath = tmp, name = "sheet.csv"), id_col = "")
    session$setInputs(merge = 1)
    expect_true("batch" %in% colnames(SummarizedExperiment::colData(state$working)))
  })
})

test_that("metadata module tags feature_class", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_metadata_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    target <- rownames(state$working)[1]
    session$setInputs(feat = target, fclass = "exogenous")
    session$setInputs(tag = 1)
    fc <- SummarizedExperiment::rowData(state$working)$feature_class
    names(fc) <- rownames(state$working)
    expect_equal(as.character(fc[target]), "exogenous")
  })
})
