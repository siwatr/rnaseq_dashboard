# The shared design builder sub-module (R/mod_design_builder.R) for P5b.

test_that("mod_design_builder applies a design (design-scoped) + relevels the reference", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_design_builder_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 2, seed = 1)),
               source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    dv0 <- state$data_version

    session$setInputs(primary = "condition", covariates = character(0),
                      ref_condition = "treated")     # dynamic reference-level input
    session$flushReact()
    expect_true(session$returned()$ok)               # ~ condition is full rank

    session$setInputs(apply = 1)
    session$flushReact()
    expect_match(paste(deparse(DESeq2::design(state$working)), collapse = " "), "condition")
    expect_equal(state$design_version, 1L)           # design-scoped bump
    expect_equal(state$data_version, dv0)            # data_version untouched
    expect_equal(levels(SummarizedExperiment::colData(state$working)$condition)[1], "treated")
  })
})

test_that("mod_design_builder refuses to apply a confounded (rank-deficient) design", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_design_builder_server, args = list(state = state), {
    dds <- make_mock_dds(n_genes = 30, n_per_group = 3, n_spike = 0, seed = 2)
    SummarizedExperiment::colData(dds)$batch <-
      SummarizedExperiment::colData(dds)$condition   # perfectly confounded with condition
    state_load(state, dds, source = "demo")
    session$flushReact()

    session$setInputs(primary = "condition", covariates = "batch")   # ~ condition + batch (not full rank)
    session$flushReact()
    session$setInputs(apply = 1)
    session$flushReact()
    # the guard refuses: the design stays ~ condition, no design_version bump
    expect_equal(state$design_version, 0L)
    expect_false(grepl("batch", paste(deparse(DESeq2::design(state$working)), collapse = " ")))
  })
})
