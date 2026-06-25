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
    rs <- state$palette$other$removal_status
    expect_equal(rs$name, "Custom palette")                 # preset, not Okabe-Ito
    expect_true(all(c("pass", "suggested_other", "suggested_this") %in% names(rs$colors)))
    expect_equal(unname(rs$colors[["pass"]]), "#2CA02C")    # QC green preset
    session$setInputs(addsel_other = "correlation", addbtn_other = 2)
    co <- state$palette$other$correlation
    expect_equal(co$name, "RColorBrewer: RdBu")             # reversed RdBu preset
    expect_true(isTRUE(co$reverse))
    expect_equal(co$min, "-1"); expect_equal(co$max, "1")
  })
})

test_that("discrete attributes above the hard cap are hidden from the add control", {
  skip_if_not_installed("DESeq2")
  withr::local_options(ddsdashboard.palette_max_levels = 1L)   # condition has 2 levels
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    html <- as.character(output$addui_colData$html)
    expect_match(html, "hidden")
    expect_match(html, "Affected fields")
  })
})

test_that("adding a discrete attribute above the warn threshold asks for confirmation", {
  skip_if_not_installed("DESeq2")
  withr::local_options(ddsdashboard.palette_warn_levels = 1L)  # condition(2) > 1 -> warn
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(addsel_colData = "condition", addbtn_colData = 1)
    expect_null(state$palette$colData$condition)        # staged in a modal, not added
    session$setInputs(confirm_add = 1)                  # Proceed
    expect_false(is.null(state$palette$colData$condition))
  })
})

test_that("continuous: reverse + custom-ramp pickers update the config", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(addsel_assays = "counts", addbtn_assays = 1)
    session$setInputs(crev_assays__counts = TRUE)
    expect_true(isTRUE(state$palette$assays$counts$reverse))
    session$setInputs(cname_assays__counts = "Custom ramp",
                      ccol1_assays__counts = "#000000", ccol2_assays__counts = "#ffffff")
    cu <- state$palette$assays$counts$custom
    expect_true(all(c("#000000", "#FFFFFF") %in% cu))
  })
})

test_that("custom-ramp UI has an N-stops selector + valid-seeded pickers (no null crash)", {
  cfg <- list(name = "Custom ramp", min = "", max = "", reverse = FALSE, custom = NULL)
  html <- as.character(ddsdashboard:::.palette_item_panel(
    NS("p"), "assays", "counts", cfg, "continuous", "numeric", character(0), has_picker = TRUE))
  expect_match(html, "cnstops_assays__counts")             # number-of-colours selector
  for (j in 1:5) expect_match(html, paste0("ccol", j, "_assays__counts"))
  expect_match(html, "#FFFFFF", fixed = TRUE)              # white->black default seed, not NULL
})

test_that("custom ramp: changing the stop count resamples via colorRampPalette", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(addsel_assays = "counts", addbtn_assays = 1)
    expect_equal(length(state$palette$assays$counts$custom), 2L)   # white -> black default
    session$setInputs(cname_assays__counts = "Custom ramp", cnstops_assays__counts = "4")
    cu <- state$palette$assays$counts$custom
    expect_equal(length(cu), 4L)
    expect_equal(toupper(cu[1]), "#FFFFFF")        # endpoints preserved
    expect_equal(toupper(cu[4]), "#000000")
  })
})

test_that("continuous reset keeps the palette and clears anchors / reverse / custom", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(addsel_assays = "counts", addbtn_assays = 1)
    session$setInputs(cname_assays__counts = "viridis: magma", cmin_assays__counts = "p5",
                      cmax_assays__counts = "p95", crev_assays__counts = TRUE)
    session$setInputs(creset_assays__counts = 1)
    cfg <- state$palette$assays$counts
    expect_equal(cfg$name, "viridis: magma")     # palette kept (not reverted to viridis)
    expect_equal(cfg$min, ""); expect_equal(cfg$max, "")
    expect_false(isTRUE(cfg$reverse))
  })
})

test_that("continuous Custom-ramp reset returns colours to white -> black", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$setInputs(addsel_assays = "counts", addbtn_assays = 1)
    session$setInputs(cname_assays__counts = "Custom ramp", cnstops_assays__counts = "3")
    expect_equal(length(state$palette$assays$counts$custom), 3L)
    session$setInputs(creset_assays__counts = 1)
    cfg <- state$palette$assays$counts
    expect_equal(cfg$name, "Custom ramp")        # palette kept
    expect_equal(toupper(cfg$custom), c("#FFFFFF", "#000000"))   # colours reset
  })
})

test_that("resetting a preset item restores its preset", {
  # correlation -> reversed RdBu / -1..1
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    session$setInputs(addsel_other = "correlation", addbtn_other = 1)
    session$setInputs(cname_other__correlation = "viridis: viridis",
                      crev_other__correlation = FALSE, cmin_other__correlation = "0")
    session$setInputs(creset_other__correlation = 1)
    co <- state$palette$other$correlation
    expect_equal(co$name, "RColorBrewer: RdBu")
    expect_true(isTRUE(co$reverse)); expect_equal(co$min, "-1"); expect_equal(co$max, "1")
  })
})

test_that("resetting removal_status restores the green/amber/red preset", {
  state <- new_app_state()
  shiny::testServer(mod_palette_server, args = list(state = state), {
    session$setInputs(addsel_other = "removal_status", addbtn_other = 1)
    session$setInputs(pin_other__removal_status_1 = "#123456")   # edit pass -> Custom
    session$setInputs(reset_other__removal_status = 1)
    rs <- state$palette$other$removal_status
    expect_equal(rs$name, "Custom palette")
    expect_equal(unname(rs$colors[["pass"]]), "#2CA02C")         # restored
  })
})
