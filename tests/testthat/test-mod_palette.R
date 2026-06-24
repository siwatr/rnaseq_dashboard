test_that("palette UI mounts the colData panels + add control", {
  ui <- as.character(mod_palette_ui("palette"))
  expect_match(ui, "palette-panels")
  expect_match(ui, "palette-add_ui")
})

test_that("adding a mapping populates state$palette, a pin updates it, remove clears it", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 30, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    expect_null(state$palette$colData)                    # empty until opted in

    # Add a discrete colData column -> a default discrete config appears.
    session$setInputs(add_col = "condition", add_btn = 1)
    cfg <- state$palette$colData$condition
    expect_equal(cfg$type, "discrete")
    expect_equal(cfg$palette, "Okabe-Ito")
    expect_length(cfg$pins, 0L)

    # Pin a level: condition levels are control(1)/treated(2); the auto colour for
    # `treated` under Okabe-Ito is #56B4E9, so #000000 differs -> recorded as a pin.
    session$setInputs(pin_condition_2 = "#000000")
    expect_equal(state$palette$colData$condition$pins[["treated"]], "#000000")
    expect_false("control" %in% names(state$palette$colData$condition$pins))

    # A picker value equal to the palette's auto colour is NOT a pin.
    session$setInputs(pin_condition_1 = "#E69F00")        # control's auto colour
    expect_false("control" %in% names(state$palette$colData$condition$pins))

    # Reset clears pins.
    session$setInputs(reset_condition = 1)
    expect_length(state$palette$colData$condition$pins, 0L)

    # Remove drops the whole mapping.
    session$setInputs(remove_condition = 1)
    expect_null(state$palette$colData$condition)
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
    session$setInputs(remove_condition = 1)                  # remove
    session$flushReact()
    expect_null(state$palette$colData$condition)
    # Re-add: the freshly registered Remove observer must ignore the stale
    # remove_condition counter (still 1) instead of immediately deleting again.
    session$setInputs(add_col = "condition", add_btn = 2)
    session$flushReact()
    expect_false(is.null(state$palette$colData$condition))
  })
})

test_that("observers re-register against new levels after a dataset change", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(add_col = "condition", add_btn = 1)
    expect_false(is.null(state$palette$colData$condition))

    # Load a second dataset whose `condition` has *different* levels; the config
    # survives (it is a UI pref). Editing a picker must now write a pin under a
    # CURRENT level (wt/ko), not a stale one (control/treated).
    dds2 <- make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 2)
    levels(SummarizedExperiment::colData(dds2)$condition) <- c("wt", "ko")
    state_load(state, dds2, source = "demo")
    expect_false(is.null(state$palette$colData$condition))   # config persisted
    # Flush so the load re-registers observers against the new levels before we
    # interact (in the app a load and a picker click are separate event cycles;
    # testServer batches them, which would otherwise destroy the live observer in
    # the same flush as the picker event).
    session$flushReact()

    session$setInputs(pin_condition_2 = "#000000")           # level 2 == "ko" now
    expect_equal(state$palette$colData$condition$pins[["ko"]], "#000000")
    expect_false("treated" %in% names(state$palette$colData$condition$pins))
  })
})
