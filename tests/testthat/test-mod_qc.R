test_that("QC module caches the metric table keyed on data_version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 2, seed = 1),
               source = "demo")
    session$flushReact()
    session$setInputs(x_axis = "sample", metric = "library_size",
                      group = "condition", sort = "none", auto = FALSE)
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
    session$setInputs(x_axis = "pct_spike", metric = "detected",
                      group = "condition", sort = "none", auto = TRUE)
    session$flushReact()
    expect_false(is.null(output$plot))
  })
})

test_that(".qc_metric_plot builds for sample and metric x-axes, both themes", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 2, seed = 1)
  tbl <- qc_per_sample_metrics(dds)
  tbl$group <- factor(rep(c("a", "b"), length.out = nrow(tbl)))

  # Discrete x (per-sample bar), and numeric x (metric-vs-metric scatter).
  p_bar <- ddsdashboard:::.qc_metric_plot(tbl, "sample", "library_size", "condition")
  p_sc  <- ddsdashboard:::.qc_metric_plot(tbl, "pct_spike", "library_size", "condition",
                                          dark_theme = TRUE)
  expect_s3_class(p_bar, "ggplot")
  expect_s3_class(p_sc, "ggplot")
})

test_that(".qc_metric_plot sorts the discrete x-axis by the metric value", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 1, seed = 5)
  tbl <- qc_per_sample_metrics(dds)
  tbl$group <- factor(rep("all", nrow(tbl)))
  p <- ddsdashboard:::.qc_metric_plot(tbl, "sample", "detected", sort = "decreasing")
  lvls <- levels(p$data$x)
  # Levels should follow detected, descending.
  ord <- tbl$sample[order(tbl$detected, decreasing = TRUE)]
  expect_equal(lvls, ord)
})
