test_that("Feature-info module tags feature_class via state_mutate", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_feature_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 20, n_per_group = 2, n_spike = 1, seed = 1),
               source = "demo")
    session$flushReact()
    target <- rownames(state$working)[1]
    v0 <- state$data_version
    session$setInputs(feat = target, fclass = "exogenous")
    session$setInputs(tag = 1)
    expect_equal(state$data_version, v0 + 1L)
    fc <- SummarizedExperiment::rowData(state$working)$feature_class
    names(fc) <- rownames(state$working)
    expect_equal(as.character(fc[target]), "exogenous")
  })
})
