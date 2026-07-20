# Size-factor normalization: decoupled from assay assignment, config carried on
# the dds, control-gene set (endogenous/spike-in/custom) + estimator type.

test_that("sizefactor_config round-trips on the dds and fills defaults", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 4, seed = 1)
  expect_equal(sizefactor_config(dds), default_sizefactor_config())      # default when absent
  cfg <- list(control = "custom", custom_ids = c("g1", "g2"), type = "poscounts")
  dds2 <- set_sizefactor_config(dds, cfg)
  got <- sizefactor_config(dds2)
  expect_equal(got$control, "custom")
  expect_equal(got$custom_ids, c("g1", "g2"))
  expect_equal(got$type, "poscounts")
  expect_equal(got$provenance, "auto")                                   # filled from default
})

test_that("estimate_size_factors honors the control-gene set + records the config", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 8, seed = 2)
  endo_ids <- rownames(dds)[SummarizedExperiment::rowData(dds)$feature_class == "endogenous"]

  e_endo <- estimate_size_factors(dds, list(control = "endogenous"))
  expect_false(is.null(DESeq2::sizeFactors(e_endo)))
  expect_equal(sizefactor_config(e_endo)$control, "endogenous")

  e_spk <- estimate_size_factors(dds, list(control = "spike_in"))
  expect_equal(sizefactor_config(e_spk)$control, "spike_in")
  # a different control set generally yields different size factors
  expect_false(isTRUE(all.equal(unname(DESeq2::sizeFactors(e_endo)),
                                unname(DESeq2::sizeFactors(e_spk)))))

  e_cus <- estimate_size_factors(dds, list(control = "custom", custom_ids = endo_ids[1:10]))
  expect_equal(sizefactor_config(e_cus)$control, "custom")
  expect_equal(sizefactor_config(e_cus)$custom_ids, endo_ids[1:10])

  # endogenous wrapper == generalized with control="endogenous"
  expect_equal(unname(DESeq2::sizeFactors(estimate_size_factors_endogenous(dds))),
               unname(DESeq2::sizeFactors(e_endo)))
})

test_that("reestimate_size_factors reuses the stored config; ensure_size_factors respects existing", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 50, n_per_group = 3, n_spike = 6, seed = 4)

  # no size factors yet -> reestimate is a no-op
  expect_null(DESeq2::sizeFactors(reestimate_size_factors(dds)))

  # ensure materializes (provenance auto) when absent
  d_auto <- ensure_size_factors(dds)
  expect_false(is.null(DESeq2::sizeFactors(d_auto)))
  expect_equal(sizefactor_config(d_auto)$provenance, "auto")

  # a spike-in config re-estimates under itself after a (simulated) edit
  d_spk <- estimate_size_factors(dds, list(control = "spike_in"))
  d_re <- reestimate_size_factors(d_spk[, -1])                 # drop a sample -> re-estimate
  expect_equal(sizefactor_config(d_re)$control, "spike_in")
  expect_equal(length(DESeq2::sizeFactors(d_re)), ncol(dds) - 1L)

  # ensure keeps existing factors (provenance loaded)
  d_keep <- ensure_size_factors(d_auto)
  expect_equal(DESeq2::sizeFactors(d_keep), DESeq2::sizeFactors(d_auto))
  expect_equal(sizefactor_config(d_keep)$provenance, "loaded")
})

test_that("all estimator types honor the control set via the row-subset path", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 80, n_per_group = 4, n_spike = 10, seed = 6)
  # iterate ignores DESeq2's controlGenes, but the subset-inherit path makes it
  # respect the control set; factors land on all samples.
  d_it <- estimate_size_factors(dds, list(control = "spike_in", type = "iterate"))
  expect_length(DESeq2::sizeFactors(d_it), ncol(dds))
  expect_equal(sizefactor_config(d_it)$type, "iterate")
  # poscounts on a custom set works too
  d_pc <- estimate_size_factors(dds, list(control = "custom",
                                          custom_ids = rownames(dds)[1:30], type = "poscounts"))
  expect_length(DESeq2::sizeFactors(d_pc), ncol(dds))
})

test_that("empty control sets error; iterate refuses a too-small control set", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 0, seed = 8)  # no spike-ins
  expect_error(estimate_size_factors(dds, list(control = "spike_in")), "No control genes")
  expect_error(estimate_size_factors(dds, list(control = "custom", custom_ids = "nope")),
               "No control genes")
  # iterate with fewer control genes than the design rank is refused
  expect_error(estimate_size_factors(dds, list(control = "custom",
                                               custom_ids = rownames(dds)[1], type = "iterate")),
               "iterate")
})

test_that("reestimate_size_factors keeps externally loaded factors", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_size_factors(make_mock_dds(n_genes = 50, n_per_group = 3, n_spike = 4, seed = 9))
  # simulate a loaded object: mark provenance loaded, keep its factors
  dds <- set_sizefactor_config(dds, utils::modifyList(sizefactor_config(dds),
                                                      list(provenance = "loaded")))
  sf0 <- DESeq2::sizeFactors(dds)
  d_re <- reestimate_size_factors(dds[, -1])            # a structural edit
  expect_equal(unname(DESeq2::sizeFactors(d_re)), unname(sf0[-1]))   # loaded factors kept (subset)
})

test_that("Size-factors UI mounts the three pills + control/estimator inputs", {
  ui <- as.character(mod_sizefactors_ui("sf"))
  # Estimate pill
  expect_match(ui, "sf-sf_control")
  expect_match(ui, "sf-sf_type")
  expect_match(ui, "sf-sf_estimate")
  expect_match(ui, "Custom set")
  expect_match(ui, "All genes", fixed = TRUE)          # the discouraged option
  expect_match(ui, "Estimate using:", fixed = TRUE)    # renamed label
  # Per-sample pill + Compare pill (last)
  expect_match(ui, "sf-pers_render")
  expect_match(ui, "sf-cmp_x_control")
  expect_match(ui, "sf-cmp_y_control")
  expect_match(ui, "sf-cmp_render")
})

test_that("'all_genes' resolves to every row and estimates", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 50, n_per_group = 3, n_spike = 6, seed = 11)
  expect_equal(length(.sf_control_index(dds, list(control = "all_genes"))), nrow(dds))
  e <- estimate_size_factors(dds, list(control = "all_genes"))
  expect_false(is.null(DESeq2::sizeFactors(e)))
  expect_equal(sizefactor_config(e)$control, "all_genes")
})

test_that("Compare pill computes two size-factor vectors read-only (no data_version bump)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::reactiveConsole(TRUE); on.exit(shiny::reactiveConsole(FALSE), add = TRUE)
  shiny::testServer(mod_sizefactors_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 6, seed = 7),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    v0 <- state$data_version
    session$setInputs(cmp_x_control = "endogenous", cmp_x_type = "ratio",
                      cmp_y_control = "spike_in", cmp_y_type = "ratio", cmp_render = 1)
    session$flushReact()
    cv <- compare_shown$value()
    expect_true(isTRUE(cv$ok))
    expect_length(cv$sf_x, ncol(state$working))
    expect_length(cv$sf_y, ncol(state$working))
    # viewing is read-only: the working dds + its size factors are untouched
    expect_identical(state$data_version, v0)
    expect_equal(sizefactor_config(state$working)$control, "endogenous")
    # the built ggplot renders (exercises geom_smooth / lm)
    p <- build_cmp_gg(FALSE)
    expect_s3_class(p, "ggplot")
  })
})

test_that("Compare pill carries a graceful message when a control set is empty", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::reactiveConsole(TRUE); on.exit(shiny::reactiveConsole(FALSE), add = TRUE)
  shiny::testServer(mod_sizefactors_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 0, seed = 9),
               source = "demo", meta = list(feature_type = "gene"))       # no spike-ins
    session$flushReact()
    session$setInputs(cmp_x_control = "endogenous", cmp_x_type = "ratio",
                      cmp_y_control = "spike_in", cmp_y_type = "ratio", cmp_auto = TRUE)
    session$flushReact()
    cv <- compare_shown$value()
    expect_false(isTRUE(cv$ok))
    expect_true(nzchar(cv$msg))
  })
})

test_that("Per-sample pill renders bar (Sample) and grouped (colData) plots", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::reactiveConsole(TRUE); on.exit(shiny::reactiveConsole(FALSE), add = TRUE)
  shiny::testServer(mod_sizefactors_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 4, seed = 5),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    session$setInputs(pers_x = "__sample__", pers_colour = "__none__", pers_auto = TRUE)
    session$flushReact()
    expect_s3_class(build_pers_gg(FALSE), "ggplot")
    session$setInputs(pers_x = "condition")
    session$flushReact()
    expect_s3_class(build_pers_gg(FALSE), "ggplot")
  })
})

test_that("Size-factors server: confirm default commits (auto->user), then value-idempotent", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_sizefactors_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 6, seed = 7),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    # load materialized endogenous defaults, provenance "auto"
    expect_equal(sizefactor_config(state$working)$control, "endogenous")
    expect_equal(sizefactor_config(state$working)$provenance, "auto")
    v0 <- state$data_version

    # confirming the endogenous default commits it (auto -> user): one bump, even
    # though the values are unchanged (the user deliberately set the config).
    session$setInputs(sf_control = "endogenous", sf_type = "ratio", sf_estimate = 1)
    session$flushReact()
    expect_equal(state$data_version, v0 + 1L)
    expect_equal(sizefactor_config(state$working)$provenance, "user")
    v1 <- state$data_version

    # re-running the SAME user config now is a value no-op (no bump)
    session$setInputs(sf_estimate = 2); session$flushReact()
    expect_equal(state$data_version, v1)

    # switch to a custom control set via the shared gene search, then estimate
    gnames <- SummarizedExperiment::rowData(state$working)$gene_name[1:12]
    session$setInputs(sf_control = "custom", sf_searchby = "gene_name",
                      sf_q = paste(gnames, collapse = ", "), sf_type = "ratio")
    session$elapse(400); session$flushReact()
    session$setInputs(sf_estimate = 3); session$flushReact()
    expect_equal(state$data_version, v1 + 1L)                  # a real change bumped it
    expect_equal(sizefactor_config(state$working)$control, "custom")
    expect_equal(sizefactor_config(state$working)$provenance, "user")
    expect_false(is.null(DESeq2::sizeFactors(state$working)))
  })
})
