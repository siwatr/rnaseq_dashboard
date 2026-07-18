# Expression page (R/mod_expression.R) for Phase 7a - Single genes tab.

test_that("Expression UI mounts the tabs + single-gene containers", {
  ui <- as.character(mod_expression_ui("ex"))
  expect_match(ui, "ex-tabs")
  expect_match(ui, "ex-gene_container")
  expect_match(ui, "ex-render")
  expect_match(ui, "Single genes")
  expect_match(ui, "Gene sets")
})

test_that("Expression UI mounts the gene-set aggregate pill", {
  ui <- as.character(mod_expression_ui("ex"))
  expect_match(ui, "ex-geneset_container")
  expect_match(ui, "ex-set_render")
  expect_match(ui, "ex-set_source")
  expect_match(ui, "ex-set_method")
  expect_match(ui, "Aggregate expression")
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

    # per-layer styling controls feed the plot (live re-plot, no deferred re-run)
    session$setInputs(violin_width = 1.2, violin_alpha = 0.5, box_width = 0.3,
                      box_alpha = 0.6, dot_size = 3, dot_alpha = 0.6,
                      dot_method = "quasirandom", dot_width = 0.2, dot_cex = 1.5)
    session$flushReact()
    expect_silent(ggplot2::ggplot_build(build_gene_gg(FALSE)))   # quasirandom: no warning
    session$setInputs(dot_method = "beeswarm"); session$flushReact()
    expect_silent(ggplot2::ggplot_build(build_gene_gg(FALSE)))   # beeswarm + cex: no warning
    session$setInputs(dot_method = "jitter"); session$flushReact()
    expect_s3_class(build_gene_gg(FALSE), "ggplot")              # explicit jitter layout
  })
})

test_that("gene-set aggregate builds from a saved set + reports gene accounting", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 100, n_per_group = 6, n_spike = 6, seed = 2)),
               source = "demo", meta = list(feature_type = "gene"))
    ids <- rownames(state$working)[1:8]
    absent <- c("__ghost_a__", "__ghost_b__")            # authored-but-absent members
    state$gene_sets <- list(SetA = new_gene_set(c(ids, absent)))

    session$setInputs(tabs = "Gene sets", set_source = "saved", set_pick = "SetA",
                      set_val_assay = "logcounts", set_val_transform = "none",
                      set_method = "mean", set_zscore = TRUE, set_only_expr = TRUE,
                      set_auto = TRUE, set_x_group = "condition", set_colour_by = "condition")
    session$elapse(300); session$flushReact()

    mm <- set_mat_shown$value()
    expect_false(is.null(mm))
    red <- set_values(mm)
    expect_equal(red$accounting$n_total, 10L)            # 8 present + 2 absent
    expect_equal(red$accounting$n_present, 8L)
    expect_length(red$values, ncol(state$working))
    expect_match(red$subtitle, "of 10 genes")
    expect_s3_class(build_set_gg(FALSE), "ggplot")

    # median + no z-score is also a valid live re-plot
    session$setInputs(set_method = "median", set_zscore = FALSE)
    session$flushReact()
    expect_s3_class(build_set_gg(FALSE), "ggplot")
  })
})

test_that("gene-set aggregate builds from a quick (uncommitted) search", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 4, n_spike = 4, seed = 5)),
               source = "demo", meta = list(feature_type = "gene"))
    gnames <- SummarizedExperiment::rowData(state$working)$gene_name[1:5]
    session$setInputs(tabs = "Gene sets", set_source = "search",
                      setsearch_searchby = "gene_name",
                      setsearch_q = paste(gnames, collapse = ", "),
                      set_val_assay = "logcounts", set_val_transform = "none",
                      set_method = "mean", set_zscore = TRUE, set_only_expr = TRUE,
                      set_auto = TRUE, set_x_group = "condition", set_colour_by = "condition")
    session$elapse(400); session$flushReact()
    mm <- set_mat_shown$value()
    red <- set_values(mm)
    expect_equal(red$accounting$n_present, 5L)           # searched ids all resolve
    expect_s3_class(build_set_gg(FALSE), "ggplot")
  })
})

test_that("norm_logcounts is selectable and resolves (not silently downgraded)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    dds <- estimate_size_factors_endogenous(
      ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 4, seed = 9)))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    gname <- SummarizedExperiment::rowData(state$working)$gene_name[5]
    # norm_logcounts is a computed key (not in assayNames) - it must still resolve.
    session$setInputs(tabs = "Single genes", val_assay = "norm_logcounts",
                      val_transform = "none", auto = TRUE, x_group = "condition",
                      colour_by = "condition", gene_searchby = "gene_name", gene_q = gname)
    session$elapse(300); session$flushReact()
    mm <- mat_shown$value()
    expect_false(is.null(mm))
    expect_match(mm$label, "normalized")
    expect_equal(get("expr_value_mat", envir = state$derived)$params, list("norm_logcounts"))
    expect_s3_class(build_gene_gg(FALSE), "ggplot")
    # and expr_default_assay picks it when size factors exist
    expect_equal(expr_default_assay(state$working), "norm_logcounts")
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
    expect_true(geom_avail()$dist_shown)                 # groups of 12 >= dist_min(10)
    expect_true(geom_avail()$dots_default)               # 12 < dots_max(100) -> dots on by default
    expect_s3_class(build_gene_gg(FALSE), "ggplot")

    dv <- state$data_version
    session$setInputs(expr_show_by = "condition", expr_show_values = "control")
    session$flushReact()
    expect_equal(state$data_version, dv)                 # display subset, no data bump
    expect_equal(group_sizes(), 12L)                     # only the kept group remains
  })
})
