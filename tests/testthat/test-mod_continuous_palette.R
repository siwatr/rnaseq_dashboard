# Reusable continuous-palette control (R/mod_continuous_palette.R), extracted from
# the Palette page for the P7c heatmap ramp (and the Phase-8 shared controller).

test_that("continuous_palette_ui renders the widget cluster with the suffix prefix", {
  ui <- as.character(continuous_palette_ui(shiny::NS("h"), "ramp",
                                           continuous_palette_default("RColorBrewer: RdBu")))
  expect_match(ui, "h-ramp_cname")
  expect_match(ui, "h-ramp_cmin")
  expect_match(ui, "h-ramp_cmax")
  expect_match(ui, "h-ramp_crev")
  expect_match(ui, "h-ramp_cedit")
  expect_match(ui, "h-ramp_creset")
})

test_that("continuous_palette_default fills sensible fields", {
  d <- continuous_palette_default("viridis: viridis")
  expect_equal(d$name, "viridis: viridis")
  expect_false(d$reverse)
  expect_length(d$custom, 2L)                      # never NULL (breaks the picker)
})

test_that("server tracks name / reverse / anchors and Edit + Reset", {
  srv <- function(input, output, session) {
    cfg <- continuous_palette_server(input, output, session, "r",
             continuous_palette_default("viridis: viridis"))
  }
  shiny::testServer(srv, {
    # Prime inputs to the UI defaults first (the observers ignoreInit, so the very
    # first value is the "initial" one the real UI supplies); then change them.
    session$setInputs(r_cname = "viridis: viridis", r_crev = FALSE,
                      r_cmin = "", r_cmax = "")
    session$flushReact()

    session$setInputs(r_cname = "RColorBrewer: RdBu", r_crev = TRUE,
                      r_cmin = "-2", r_cmax = "2")
    session$flushReact()
    expect_equal(cfg()$name, "RColorBrewer: RdBu")
    expect_true(cfg()$reverse)
    expect_equal(cfg()$min, "-2")
    expect_equal(cfg()$max, "2")

    # Edit palette: copy the named ramp into an editable 5-stop Custom ramp.
    session$setInputs(r_cname = "viridis: viridis", r_cedit = 1)
    session$flushReact()
    expect_equal(cfg()$name, "Custom ramp")
    expect_length(cfg()$custom, 5L)
    expect_false(cfg()$reverse)                     # reverse baked in then zeroed

    # Reset: clear anchors + reverse.
    session$setInputs(r_creset = 1)
    session$flushReact()
    expect_equal(cfg()$min, "")
    expect_equal(cfg()$max, "")
    expect_false(cfg()$reverse)
  })
})

test_that("server collects the visible custom-ramp pickers into `custom`", {
  srv <- function(input, output, session) {
    cfg <- continuous_palette_server(input, output, session, "r",
             continuous_palette_default("Custom ramp", custom = c("#FFFFFF", "#000000")))
  }
  shiny::testServer(srv, {
    session$setInputs(r_cname = "Custom ramp", r_cnstops = 2,   # prime (ignoreInit)
                      r_ccol1 = "#FFFFFF", r_ccol2 = "#000000")
    session$flushReact()
    session$setInputs(r_cnstops = 3,
                      r_ccol1 = "#FF0000", r_ccol2 = "#00FF00", r_ccol3 = "#0000FF")
    session$flushReact()
    expect_equal(cfg()$custom, c("#FF0000", "#00FF00", "#0000FF"))
  })
})
