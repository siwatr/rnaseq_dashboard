# Expression page (R/mod_expression.R) for Phase 7a - Single genes tab.

test_that("Expression UI mounts the tabs + single-gene containers", {
  ui <- as.character(mod_expression_ui("ex"))
  expect_match(ui, "ex-tabs")
  expect_match(ui, "ex-gene_container")
  expect_match(ui, "ex-render")
  expect_match(ui, "Single genes")
  expect_match(ui, "Gene sets")
})

test_that("single-gene plot builds; value matrix caches; gene id drives the sig", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 100, n_per_group = 3, n_spike = 6, seed = 1)),
               source = "demo", meta = list(feature_type = "gene"))
    gname <- SummarizedExperiment::rowData(state$working)$gene_name[20]
    gid   <- rownames(state$working)[20]
    session$setInputs(tabs = "Single genes", val_assay = "logcounts", val_transform = "none",
                      auto = TRUE, x_group = "condition", colour_by = "condition",
                      gene_searchby = "gene_name", gene_q = gname)
    session$elapse(300); session$flushReact()

    mm <- mat_shown$value()
    expect_false(is.null(mm))
    expect_true(is.matrix(mm$mat))
    gv <- gene_values(mm)
    expect_equal(gv$gene_id, gid)
    expect_length(gv$values, ncol(state$working))
    # value matrix cached in derived under expr_value_mat, keyed on the assay
    expect_equal(get("expr_value_mat", envir = state$derived)$params, list("logcounts"))
    expect_s3_class(build_gene_gg(FALSE), "ggplot")
  })
})

test_that("VST value is endogenous-only and rejects a spike-in gene", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 100, n_per_group = 3, n_spike = 8, seed = 4)),
               source = "demo", meta = list(feature_type = "gene"))
    rd <- SummarizedExperiment::rowData(state$working)
    spike_id <- rownames(state$working)[which(rd$feature_class == "spike_in")[1]]
    session$setInputs(tabs = "Single genes", val_assay = "vst", val_transform = "none",
                      auto = TRUE, x_group = "condition", colour_by = "condition",
                      gene_searchby = "__rownames__", gene_q = spike_id)
    session$elapse(300); session$flushReact()
    # VST (endogenous-only) doesn't carry the spike-in row -> a clear validation.
    expect_error(build_gene_gg(FALSE), "endogenous-only")
  })
})

test_that("distribution geoms appear only for large groups; 'Showing' is display-only", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 12, n_spike = 4, seed = 7)),
               source = "demo", meta = list(feature_type = "gene"))
    gname <- SummarizedExperiment::rowData(state$working)$gene_name[10]
    session$setInputs(tabs = "Single genes", val_assay = "logcounts", val_transform = "none",
                      auto = TRUE, x_group = "condition", colour_by = "condition",
                      gene_searchby = "gene_name", gene_q = gname)
    session$elapse(300); session$flushReact()
    expect_equal(group_sizes(), c(12L, 12L))
    expect_true(geom_avail()$dist_shown)                 # groups of 12 >= g1(10)
    expect_false(geom_avail()$dots_default)              # 12 > g1 -> dots off by default
    expect_s3_class(build_gene_gg(FALSE), "ggplot")

    dv <- state$data_version
    session$setInputs(expr_show_by = "condition", expr_show_values = "control")
    session$flushReact()
    expect_equal(state$data_version, dv)                 # display subset, no data bump
    expect_equal(group_sizes(), 12L)                     # only the kept group remains
  })
})
