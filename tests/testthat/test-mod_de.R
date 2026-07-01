# The DE page server (R/mod_de.R) Design & Contrasts tab for P5b / P5b-2.

test_that("mod_de: add a contrast, Run fits + extracts; a data edit goes stale", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_de_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 4, n_spike = 4, seed = 2)),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()

    session$setInputs(shrink = "none", c_var = "condition", c_test = "treated", c_control = "control")
    session$setInputs(c_add = 1); session$flushReact()
    expect_length((state$de)$contrasts, 1L)

    session$setInputs(run = 1); session$flushReact()
    expect_length((state$de)$results, 1L)                 # fit + auto-extract
    expect_equal(de_status(state), "current")
    df <- (state$de)$results[[1]]
    expect_true(all(c("log2FoldChange", "padj", "log2FoldChange_shrunk") %in% names(df)))

    state_mutate(state, function(d) d[, -1], action = list(action = "drop_sample"))
    expect_equal(de_status(state), "stale")               # fit no longer matches data
  })
})

test_that("mod_de: extraction is reactive (auto on) and gated (auto off)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_de_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 4, n_spike = 0, seed = 3)),
               source = "demo")
    session$flushReact()
    session$setInputs(shrink = "none", c_var = "condition", c_test = "treated", c_control = "control")
    session$setInputs(c_add = 1, run = 1); session$flushReact()
    expect_length((state$de)$results, 1L)

    # auto ON: remove then re-add re-extracts from the current fit without a new Run
    session$setInputs(remove_multi = "condition: treated vs control", remove_sel = 1); session$flushReact()
    expect_length((state$de)$results, 0L)
    session$setInputs(c_add = 2); session$flushReact()
    expect_length((state$de)$results, 1L)                 # seamless re-extraction

    # auto OFF: re-adding does not extract until Update results
    session$setInputs(auto_update = FALSE); session$flushReact()
    session$setInputs(remove_multi = "condition: treated vs control", remove_sel = 2); session$flushReact()
    expect_length((state$de)$results, 0L)
    session$setInputs(c_add = 3); session$flushReact()
    expect_length((state$de)$results, 0L)                 # gated
    session$setInputs(update_results = 1); session$flushReact()
    expect_length((state$de)$results, 1L)                 # extracted on demand
  })
})

test_that("mod_de: invalid contrast is skipped at extraction; remove-all works", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_de_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 0, seed = 8)),
               source = "demo")
    session$flushReact()
    session$setInputs(shrink = "none", c_var = "condition", c_test = "treated", c_control = "control")
    session$setInputs(c_add = 1); session$flushReact()
    # inject a stale spec referencing a level that doesn't exist
    de <- state$de
    de$contrasts <- c(de$contrasts, list(list(var = "condition", test = "ghost",
                                              control = "control", label = "condition: ghost vs control")))
    state$de <- de
    session$setInputs(run = 1); session$flushReact()
    res <- (state$de)$results
    expect_true("condition: treated vs control" %in% names(res))   # extractable one ran
    expect_false("condition: ghost vs control" %in% names(res))    # invalid one skipped

    # remove invalid drops the ghost spec; remove all clears everything
    session$setInputs(remove_invalid = 1); session$flushReact()
    expect_length((state$de)$contrasts, 1L)
    session$setInputs(remove_all = 1); session$flushReact()
    expect_length((state$de)$contrasts, 0L)
    expect_length((state$de)$results, 0L)
  })
})

test_that("mod_de: duplicate contrasts are rejected", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_de_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 5)),
               source = "demo")
    session$flushReact()
    session$setInputs(c_var = "condition", c_test = "treated", c_control = "control")
    session$setInputs(c_add = 1); session$flushReact()
    session$setInputs(c_add = 2); session$flushReact()          # same spec again
    expect_length((state$de)$contrasts, 1L)
  })
})
