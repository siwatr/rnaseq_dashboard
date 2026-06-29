# Shared plot engine (R/mod_plot_engine.R) extracted from mod_qc (P4-pre).

test_that("pure helpers behave", {
  expect_s3_class(ddsdashboard:::.plot_msg("hi"), "ggplot")
  expect_s3_class(ddsdashboard:::.plot_dual("x"), "shiny.tag")
  expect_equal(ddsdashboard:::.plotly_max_elements(), 5000L)
  withr::local_options(ddsdashboard.plotly_max_elements = 42L)
  expect_equal(ddsdashboard:::.plotly_max_elements(), 42L)
})

test_that("plot_engine_server returns the toggle/deferred closures; engine off by default", {
  # Minimal host module that exposes the engine pieces for assertions.
  host <- function(id, state) {
    shiny::moduleServer(id, function(input, output, session) {
      plot_engine_server(input, output, session, state)
    })
  }
  state <- new_app_state()
  shiny::testServer(host, args = list(state = state), {
    expect_setequal(names(session$returned),
                    c("use_plotly_base", "dual_plot", "deferred", "stale_note"))
    expect_type(session$returned$dual_plot, "closure")
    expect_type(session$returned$deferred, "closure")
    # No dataset + toggle off -> not using the plotly base engine.
    expect_false(isTRUE(session$returned$use_plotly_base()))
  })
})

test_that("deferred renders on the button and goes stale when the signature changes", {
  host <- function(id, state) {
    shiny::moduleServer(id, function(input, output, session) {
      eng <- plot_engine_server(input, output, session, state)
      eng$deferred("auto", "render", spec = function() input$x %||% 0,
                   sig = function() input$x %||% 0)
    })
  }
  state <- new_app_state()
  shiny::testServer(host, args = list(state = state), {
    d <- session$returned
    session$setInputs(auto = FALSE, x = 1, render = 1); session$flushReact()
    expect_equal(d$value(), 1)
    expect_false(d$stale())
    session$setInputs(x = 2); session$flushReact()   # changed since last render
    expect_true(d$stale())
    session$setInputs(render = 2); session$flushReact()
    expect_equal(d$value(), 2); expect_false(d$stale())
  })
})
