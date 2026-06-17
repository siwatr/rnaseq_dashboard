test_that("QC module caches the metric table keyed on data_version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 2, seed = 1),
               source = "demo")
    session$flushReact()
    session$setInputs(metric = "library_size", plot_type = "bar",
                      group = "condition", auto = FALSE)
    session$setInputs(render = 1)

    # The derived cache holds the metric table at the current data_version.
    expect_true(exists("qc_metrics", envir = state$derived, inherits = FALSE))
    cached <- get("qc_metrics", envir = state$derived)
    expect_equal(cached$version, state$data_version)
    expect_equal(nrow(cached$value), ncol(state$working))
    expect_setequal(colnames(cached$value),
                    c("sample", "library_size", "detected", "pct_mito", "pct_spike"))
  })
})

test_that("QC module re-renders when auto-render is enabled (no button press)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 1, seed = 2),
               source = "demo")
    session$flushReact()
    session$setInputs(metric = "detected", plot_type = "box",
                      group = "condition", auto = TRUE)
    session$flushReact()
    expect_false(is.null(output$plot))
  })
})
