# The DE page server (R/mod_de.R) Design & Contrasts tab for P5b.

test_that("mod_de: add a contrast, run DESeq2, store results + current stamp; edits go stale", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_de_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 4, n_spike = 4, seed = 2)),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()

    session$setInputs(c_var = "condition", c_test = "treated", c_control = "control")
    session$setInputs(c_add = 1)
    session$flushReact()
    expect_length((state$de)$contrasts, 1L)

    session$setInputs(shrink = "none", run = 1)
    session$flushReact()
    expect_length((state$de)$results, 1L)
    expect_equal(de_status(state), "current")

    df <- (state$de)$results[[1]]
    expect_true(all(c("log2FoldChange", "padj", "log2FoldChange_shrunk") %in% names(df)))
    expect_equal(nrow(df), nrow(state$working))

    # a data edit marks the DE results stale (design-independent caches would survive)
    state_mutate(state, function(d) d[, -1], action = list(action = "drop_sample"))
    expect_equal(de_status(state), "stale")
  })
})

test_that("mod_de: duplicate contrasts are rejected; remove works", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_de_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 5)),
               source = "demo")
    session$flushReact()

    session$setInputs(c_var = "condition", c_test = "treated", c_control = "control")
    session$setInputs(c_add = 1); session$flushReact()
    session$setInputs(c_add = 2); session$flushReact()          # same spec again
    expect_length((state$de)$contrasts, 1L)                      # de-duplicated

    session$setInputs(active = "condition: treated vs control", remove_active = 1)
    session$flushReact()
    expect_length((state$de)$contrasts, 0L)
  })
})
