sample_opts <- list(slot = "colData", title = "Sample", row_noun = "sample",
                    allow_row_rename = TRUE, bulk_class = FALSE)
feature_opts <- list(slot = "rowData", title = "Feature", row_noun = "feature",
                     allow_row_rename = FALSE, bulk_class = TRUE)

test_that(".meta_display_df gives the id a unique name, even when it collides", {
  # No collision: id column takes the plain row_noun name and holds the rownames.
  df <- data.frame(condition = c("a", "b"), row.names = c("S1", "S2"))
  d1 <- .meta_display_df(df, "sample")
  expect_identical(names(d1), c("sample", "condition"))
  expect_identical(d1[[1]], c("S1", "S2"))

  # Collision: an input column literally named "sample" keeps its name; the id
  # column is suffixed to "sample.1" and still holds the rownames.
  df2 <- data.frame(sample = c("x", "y"), condition = c("a", "b"),
                    row.names = c("S1", "S2"))
  d2 <- .meta_display_df(df2, "sample")
  expect_identical(names(d2), c("sample.1", "sample", "condition"))
  expect_identical(d2[["sample.1"]], c("S1", "S2"))   # id column = rownames
  expect_identical(d2[["sample"]], c("x", "y"))       # data column preserved
})

test_that("meta_editor: a colData column named like the id still edits correctly", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(meta_editor_server, args = list(state = state, opts = sample_opts), {
    dds <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1)
    dds <- add_meta_column(dds, "colData", "sample", "character", "tag")  # collides with id
    state_load(state, dds, source = "demo")
    session$flushReact()
    # colData order: condition(1), bio_rep(2), group(3), sample(4). The id sits at
    # display col 0, so data-column indices are unchanged: editing col 4 lands on
    # the user "sample" column, not on the read-only id.
    scol <- match("sample", colnames(SummarizedExperiment::colData(state$working)))
    session$setInputs(table_cell_edit = list(row = 1L, col = scol, value = "edited"))
    session$setInputs(save = 1)
    expect_equal(as.character(SummarizedExperiment::colData(state$working)$sample[1]), "edited")
  })
})

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
