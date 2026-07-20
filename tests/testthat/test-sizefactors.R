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

test_that("Size-factors UI mounts control + estimator inputs", {
  ui <- as.character(mod_sizefactors_ui("sf"))
  expect_match(ui, "sf-sf_control")
  expect_match(ui, "sf-sf_type")
  expect_match(ui, "sf-sf_estimate")
  expect_match(ui, "Custom set")
})

test_that("Size-factors server estimates on a custom set and is idempotent", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_sizefactors_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 6, seed = 7),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    # load materialized endogenous defaults
    expect_equal(sizefactor_config(state$working)$control, "endogenous")
    v0 <- state$data_version

    # re-estimating the SAME (endogenous) config is a no-op (no version bump)
    session$setInputs(sf_control = "endogenous", sf_type = "ratio", sf_estimate = 1)
    session$flushReact()
    expect_equal(state$data_version, v0)

    # switch to a custom control set via the shared gene search, then estimate
    gnames <- SummarizedExperiment::rowData(state$working)$gene_name[1:12]
    session$setInputs(sf_control = "custom", sf_searchby = "gene_name",
                      sf_q = paste(gnames, collapse = ", "), sf_type = "ratio")
    session$elapse(400); session$flushReact()
    session$setInputs(sf_estimate = 2); session$flushReact()
    expect_equal(state$data_version, v0 + 1L)                  # a real change bumped it
    expect_equal(sizefactor_config(state$working)$control, "custom")
    expect_equal(sizefactor_config(state$working)$provenance, "user")
    expect_false(is.null(DESeq2::sizeFactors(state$working)))
  })
})
