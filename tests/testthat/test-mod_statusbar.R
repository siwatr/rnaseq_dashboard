test_that("status bar interactive switch writes the global plot-engine flag", {
  expect_match(paste(as.character(mod_statusbar_ui("statusbar")), collapse = " "),
               "Interactive plots")
  state <- new_app_state()
  shiny::testServer(mod_statusbar_server, args = list(state = state), {
    expect_false(isTRUE(state$plot_interactive))
    session$setInputs(interactive = TRUE)
    expect_true(state$plot_interactive)
    session$setInputs(interactive = FALSE)
    expect_false(state$plot_interactive)
  })
})

test_that("status bar Undo steps back the last edit; Reset restores original", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_statusbar_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 2,
                                                     n_spike = 1, seed = 1)), source = "demo")
    n_full <- ncol(state$working)
    state_mutate(state, function(d) drop_samples(d, colnames(d)[1]),
                 action = list(action = "filter_samples"))
    expect_equal(ncol(state$working), n_full - 1L)
    expect_false(is.null(output$actions))     # force the actions UI to render
    v <- state$data_version
    session$setInputs(undo = 1)
    expect_equal(ncol(state$working), n_full)          # last edit reversed
    expect_gt(state$data_version, v)
    state_mutate(state, function(d) drop_samples(d, colnames(d)[1]),
                 action = list(action = "filter_samples"))
    session$setInputs(reset_confirm = 1)               # confirm-modal OK
    expect_equal(ncol(state$working), ncol(state$original))
  })
})

test_that("status bar Undo is a no-op on an empty stack", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_statusbar_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 2),
               source = "demo")
    v <- state$data_version
    session$setInputs(undo = 1)
    expect_equal(state$data_version, v)
  })
})

test_that("status bar UI mounts the global-actions output", {
  expect_match(as.character(mod_statusbar_ui("sb")), "sb-actions")
})

test_that("status bar flags spike-in / exogenous features when present", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_statusbar_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 3, seed = 1),
               source = "demo")
    html <- paste(as.character(output$status), collapse = " ")
    expect_match(html, "spike-in")
    expect_match(html, "exogenous")
  })
})
