test_that("palette UI mounts Setting/Preview tabs + the four domain pills", {
  ui <- as.character(mod_palette_ui("palette"))
  expect_match(ui, "palette-panels_colData")
  expect_match(ui, "palette-addui_colData")
  expect_match(ui, "palette-panels_assays")
  expect_match(ui, "palette-panels_other")
  expect_match(ui, "Setting")
  expect_match(ui, "Preview")
})

test_that("colData discrete: add seeds Okabe-Ito; editing a level flips to Custom; reset; remove", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 30, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    expect_null(state$palette$colData)

    session$setInputs(addsel_colData = "condition", addbtn_colData = 1)
    cfg <- state$palette$colData$condition
    expect_equal(cfg$name, "Okabe-Ito")
    expect_equal(unname(cfg$colors[["treated"]]), "#56B4E9")

    session$setInputs(pin_colData__condition_2 = "#000000")
    cfg <- state$palette$colData$condition
    expect_equal(unname(cfg$colors[["treated"]]), "#000000")
    expect_equal(cfg$name, "Custom palette")

    session$setInputs(reset_colData__condition = 1)
    expect_equal(state$palette$colData$condition$name, "Okabe-Ito")
    expect_equal(unname(state$palette$colData$condition$colors[["treated"]]), "#56B4E9")

    session$setInputs(remove_colData__condition = 1)
    expect_null(state$palette$colData$condition)
  })
})

test_that("colData discrete: removed mapping can be re-added (stale Remove must not re-fire)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(addsel_colData = "condition", addbtn_colData = 1); session$flushReact()
    session$setInputs(remove_colData__condition = 1); session$flushReact()
    expect_null(state$palette$colData$condition)
    session$setInputs(addsel_colData = "condition", addbtn_colData = 2); session$flushReact()
    expect_false(is.null(state$palette$colData$condition))
  })
})

test_that("assays continuous: add seeds viridis; anchors store; reset", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(addsel_assays = "counts", addbtn_assays = 1)
    cfg <- state$palette$assays$counts
    expect_equal(cfg$name, "viridis: viridis")
    expect_equal(cfg$min, "")
    # Set percentile anchors.
    session$setInputs(cmin_assays__counts = "p5", cmax_assays__counts = "p95")
    cfg <- state$palette$assays$counts
    expect_equal(cfg$min, "p5"); expect_equal(cfg$max, "p95")
    # Reset returns to viridis + blank anchors.
    session$setInputs(creset_assays__counts = 1)
    cfg <- state$palette$assays$counts
    expect_equal(cfg$name, "viridis: viridis"); expect_equal(cfg$min, "")
  })
})

test_that("Other pill: removal_status (discrete) and correlation (continuous) configurable without a dataset", {
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    # No dataset loaded -- the Other maps are app-internal, still addable.
    session$setInputs(addsel_other = "removal_status", addbtn_other = 1)
    expect_equal(state$palette$other$removal_status$name, "Okabe-Ito")
    expect_true(all(c("pass", "suggested_other", "suggested_this") %in%
                      names(state$palette$other$removal_status$colors)))
    session$setInputs(addsel_other = "correlation", addbtn_other = 2)
    expect_equal(state$palette$other$correlation$name, "viridis: viridis")
  })
})
