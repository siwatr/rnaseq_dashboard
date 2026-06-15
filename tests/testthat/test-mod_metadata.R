test_that("Sample-info cell edits stay in the draft until Save", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_metadata_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$flushReact()
    v0 <- state$data_version

    # Edit colData column 1 ("condition") of sample 1 in the draft.
    session$setInputs(coldata_cell_edit = list(row = 1L, col = 1L, value = "treated"))
    expect_equal(state$data_version, v0)   # not committed yet
    expect_equal(as.character(SummarizedExperiment::colData(state$working)$condition[1]),
                 "control")                # working unchanged

    session$setInputs(save = 1)
    expect_equal(state$data_version, v0 + 1L)   # one commit
    expect_equal(as.character(SummarizedExperiment::colData(state$working)$condition[1]),
                 "treated")
  })
})

test_that("adding a column is draft-only until Save", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_metadata_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$flushReact()
    v0 <- state$data_version
    session$setInputs(add_name = "batch", add_type = "character", add_default = "A")
    session$setInputs(add_col = 1)
    expect_false("batch" %in% colnames(SummarizedExperiment::colData(state$working)))
    session$setInputs(save = 1)
    expect_equal(state$data_version, v0 + 1L)
    expect_true("batch" %in% colnames(SummarizedExperiment::colData(state$working)))
  })
})

test_that("Reset to original discards draft edits (Save then no-ops)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_metadata_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$flushReact()
    session$setInputs(coldata_cell_edit = list(row = 1L, col = 1L, value = "treated"))
    session$setInputs(reset_orig = 1)
    v0 <- state$data_version
    session$setInputs(save = 1)
    expect_equal(state$data_version, v0)   # nothing to save after reset
    expect_equal(as.character(SummarizedExperiment::colData(state$working)$condition[1]),
                 "control")
  })
})
