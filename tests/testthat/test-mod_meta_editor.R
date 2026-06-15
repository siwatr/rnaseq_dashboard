sample_opts <- list(slot = "colData", title = "Sample", row_noun = "sample",
                    allow_merge = TRUE, allow_row_rename = TRUE, bulk_class = FALSE)
feature_opts <- list(slot = "rowData", title = "Feature", row_noun = "feature",
                     allow_merge = FALSE, allow_row_rename = FALSE, bulk_class = TRUE)

test_that("meta_editor (samples): cell edits stay in draft until Save", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(meta_editor_server, args = list(state = state, opts = sample_opts), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$flushReact()
    v0 <- state$data_version
    session$setInputs(table_cell_edit = list(row = 1L, col = 1L, value = "treated"))  # condition
    expect_equal(state$data_version, v0)                       # draft only
    session$setInputs(save = 1)
    expect_equal(state$data_version, v0 + 1L)                  # one commit
    expect_equal(as.character(SummarizedExperiment::colData(state$working)$condition[1]), "treated")
  })
})

test_that("meta_editor (samples): multi-column delete skips protected", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(meta_editor_server, args = list(state = state, opts = sample_opts), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$flushReact()
    session$setInputs(rm_cols = c("bio_rep", "group", "condition"))
    session$setInputs(remove_cols = 1)
    session$setInputs(save = 1)
    cd <- SummarizedExperiment::colData(state$working)
    expect_false(any(c("bio_rep", "group") %in% colnames(cd)))
    expect_true("condition" %in% colnames(cd))                 # protected, skipped
  })
})

test_that("meta_editor (features): feature_class edits validate; bulk set on filtered rows", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(meta_editor_server, args = list(state = state, opts = feature_opts), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 2, seed = 1),
               source = "demo")
    session$flushReact()
    # rowData columns: gene_name(1), feature_class(2), ... -> col 2 is feature_class
    session$setInputs(table_cell_edit = list(row = 1L, col = 2L, value = "exogenous"))
    session$setInputs(save = 1)
    expect_equal(as.character(SummarizedExperiment::rowData(state$working)$feature_class[1]),
                 "exogenous")

    # bulk: set all (filtered = all) rows to spike_in
    n <- nrow(state$working)
    session$setInputs(bulk_value = "spike_in", table_rows_all = seq_len(n))
    session$setInputs(bulk_apply = 1)
    session$setInputs(save = 1)
    expect_true(all(as.character(SummarizedExperiment::rowData(state$working)$feature_class) == "spike_in"))
  })
})
