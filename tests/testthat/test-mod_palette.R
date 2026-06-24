test_that("palette UI mounts the colData panels + add control + Reference tab", {
  ui <- as.character(mod_palette_ui("palette"))
  expect_match(ui, "palette-panels")
  expect_match(ui, "palette-add_ui")
  expect_match(ui, "Reference")
})

test_that("adding a mapping seeds a discrete config; editing a level flips to Custom", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 30, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    expect_null(state$palette$colData)                    # empty until opted in

    # Add a discrete colData column -> a default Qualitative/Okabe-Ito config.
    session$setInputs(add_col = "condition", add_btn = 1)
    cfg <- state$palette$colData$condition
    expect_equal(cfg$type, "Qualitative")
    expect_equal(cfg$name, "Okabe-Ito")
    expect_setequal(names(cfg$colors), c("control", "treated"))
    expect_equal(unname(cfg$colors[["treated"]]), "#56B4E9")   # Okabe-Ito #2

    # Hand-edit `treated` (level 2) -> that colour changes and type flips to Custom.
    session$setInputs(pin_condition_2 = "#000000")
    cfg <- state$palette$colData$condition
    expect_equal(unname(cfg$colors[["treated"]]), "#000000")
    expect_equal(cfg$type, "Custom")
    expect_equal(cfg$name, "Custom palette")

    # Reset -> back to the default palette (treated returns to the Okabe-Ito colour).
    session$setInputs(reset_condition = 1)
    cfg <- state$palette$colData$condition
    expect_equal(cfg$type, "Qualitative")
    expect_equal(unname(cfg$colors[["treated"]]), "#56B4E9")

    # Remove drops the whole mapping.
    session$setInputs(remove_condition = 1)
    expect_null(state$palette$colData$condition)
  })
})

test_that("choosing a palette name regenerates the colours", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("viridisLite")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(add_col = "condition", add_btn = 1)
    session$flushReact()
    # Switch type to Sequential, then pick a viridis palette.
    session$setInputs(type_condition = "Sequential")
    session$setInputs(name_condition = "viridis: viridis")
    session$flushReact()
    cfg <- state$palette$colData$condition
    expect_equal(cfg$type, "Sequential")
    expect_equal(cfg$name, "viridis: viridis")
    expect_equal(unname(cfg$colors), palette_colors("Sequential", "viridis: viridis", 2))
  })
})

test_that("a removed mapping can be re-added (stale Remove button must not re-fire)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(add_col = "condition", add_btn = 1)
    session$flushReact()
    session$setInputs(remove_condition = 1)
    session$flushReact()
    expect_null(state$palette$colData$condition)
    session$setInputs(add_col = "condition", add_btn = 2)
    session$flushReact()
    expect_false(is.null(state$palette$colData$condition))
  })
})

test_that("a surviving config resolves against a reloaded dataset's new levels", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(add_col = "condition", add_btn = 1)
    expect_false(is.null(state$palette$colData$condition))

    # Load a dataset whose `condition` has different levels; the config is a UI
    # preference, so it persists. The per-level observers read levels live, so
    # they can never write under a stale level name (the reviewer's concern).
    dds2 <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 2)
    levels(SummarizedExperiment::colData(dds2)$condition) <- c("wt", "ko")
    state_load(state, dds2, source = "demo")
    session$flushReact()
    expect_false(is.null(state$palette$colData$condition))

    # The resolver fills the new levels from the config's palette with no stale
    # corruption -- this is exactly what the QC plots / heatmap annotations consume.
    cfg <- state$palette$colData$condition
    res <- palette_discrete(c("wt", "ko"), cfg$colors, cfg$type %||% "Qualitative",
                            cfg$name %||% "Okabe-Ito", cfg$custom)
    expect_setequal(names(res), c("wt", "ko"))
    expect_true(all(grepl("^#", res)))
  })
})
