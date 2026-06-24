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
