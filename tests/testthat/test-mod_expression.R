# Expression page (R/mod_expression.R) for Phase 7a - Single genes tab.

test_that("Expression UI mounts the tabs + single-gene containers", {
  ui <- as.character(mod_expression_ui("ex"))
  expect_match(ui, "ex-tabs")
  expect_match(ui, "ex-gene_container")
  expect_match(ui, "ex-render")
  expect_match(ui, "ex-ylim_min")
  expect_match(ui, "ex-ylim_max")
  expect_match(ui, "Single genes")
  expect_match(ui, "Gene sets")
})

test_that("Expression UI mounts the gene-set heatmap pill (P7c)", {
  ui <- as.character(mod_expression_ui("ex"))
  expect_match(ui, "ex-hm_plot_container")       # reactive-size container
  expect_match(ui, "ex-hm_render")
  expect_match(ui, "ex-hm_zscore")
  expect_match(ui, "ex-hm_cluster_rows")
  expect_match(ui, "ex-hm_ramp_cname")          # extracted continuous-palette control
  expect_match(ui, "ex-hm_height")              # plot-size sliders
  expect_match(ui, "ex-hm_width")
  expect_match(ui, "ex-hm_collapse_all")        # collapse / expand all
  expect_match(ui, "ex-hm_expand_all")
  expect_match(ui, "ex-hm_acc")                 # accordion id (collapse/expand target)
  expect_match(ui, "Heatmap")
})

test_that(".expr_single_plot clamps out-of-range points to boundary triangles", {
  df <- data.frame(sample = paste0("S", 1:6),
                   group = factor(rep(c("a", "b"), each = 3)),
                   value = c(1, 2, 50, 1, 2, 3))            # 50 is the outlier
  p <- .expr_single_plot(df, "grp", "val", "t",
                         show_violin = FALSE, show_box = FALSE, show_dots = TRUE,
                         y_range = c(NA, 10))
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  ys <- unlist(lapply(b$data, function(d) d$y))
  expect_lte(max(ys, na.rm = TRUE), 10)                     # clamped to the max
  # the clamped point renders as a triangle (shape 17); in-range as circles (16)
  shapes <- unlist(lapply(b$data, function(d) d$shape))
  expect_true(any(shapes == 17, na.rm = TRUE))
  expect_true(any(shapes == 16, na.rm = TRUE))
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

    got <- gene_out$value()                          # gated: list(mat_cols, red)
    expect_false(is.null(got))
    red <- got$red
    expect_true(red$ok)
    expect_equal(red$gene_id, gid)
    expect_length(red$values, ncol(state$working))
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

    # y-axis limit is a live display input (no re-render needed)
    session$setInputs(ylim_min = 0, ylim_max = 5); session$flushReact()
    expect_s3_class(build_gene_gg(FALSE), "ggplot")
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

    got <- set_out$value()
    expect_false(is.null(got))
    red <- got$red
    expect_true(red$ok)
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
    red <- set_out$value()$red
    expect_true(red$ok)
    expect_equal(red$accounting$n_present, 5L)           # searched ids all resolve
    expect_s3_class(build_set_gg(FALSE), "ggplot")
  })
})

test_that("gene-set aggregate respects the Render gate when auto-render is off", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 4, n_spike = 4, seed = 8)),
               source = "demo", meta = list(feature_type = "gene"))
    rn <- rownames(state$working)
    state$gene_sets <- list(SetA = new_gene_set(rn[1:6]), SetB = new_gene_set(rn[10:20]))
    session$setInputs(tabs = "Gene sets", set_source = "saved", set_pick = "SetA",
                      set_val_assay = "logcounts", set_val_transform = "none",
                      set_method = "mean", set_zscore = TRUE, set_only_expr = TRUE,
                      set_auto = FALSE, set_x_group = "condition", set_colour_by = "condition",
                      set_render = 1)
    session$flushReact()
    expect_equal(set_out$value()$red$accounting$n_present, 6L)   # SetA rendered

    # Switch set with auto-render OFF -> gated: the rendered result does NOT change
    session$setInputs(set_pick = "SetB"); session$flushReact()
    expect_equal(set_out$value()$red$accounting$n_present, 6L)   # still SetA
    expect_true(set_out$stale())                                 # stale banner would show

    # Click Render -> now reflects SetB (rn[10:20] = 11 genes)
    session$setInputs(set_render = 2); session$flushReact()
    expect_equal(set_out$value()$red$accounting$n_present, 11L)
    expect_false(set_out$stale())
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
    red <- gene_out$value()$red
    expect_true(red$ok)
    expect_match(red$y_lab, "normalized")
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

# --- P7c: gene-set heatmap pill --------------------------------------------

test_that("heatmap snapshot builds from a saved set; value matrix caches; draws", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 4, n_spike = 4, seed = 3)),
               source = "demo", meta = list(feature_type = "gene"))
    rn <- rownames(state$working)
    state$gene_sets <- list(SetA = new_gene_set(rn[1:8]))
    session$setInputs(tabs = "Gene sets", hm_source = "saved", hm_pick = "SetA",
                      hm_val_assay = "logcounts", hm_val_transform = "none",
                      hm_zscore = TRUE, hm_only_expr = TRUE, hm_ramp_src = "custom",
                      hm_row_mode = "auto", hm_col_mode = "auto",
                      hm_cluster_rows = TRUE, hm_cluster_cols = TRUE,
                      hm_row_dend = "auto", hm_col_dend = "auto",
                      hm_anno = "condition", hm_render = 1)
    session$flushReact()

    s <- hm_out$value()
    expect_false(is.null(s))
    expect_equal(s$hm$n_present, 8L)
    expect_equal(nrow(s$mat), 8L)
    expect_equal(ncol(s$mat), ncol(state$working))       # all samples (no Showing filter)
    expect_true(s$zscored)
    expect_equal(s$row_mode, "all")                      # 8 rows < row-label threshold
    expect_true(s$cluster_rows && s$cluster_cols)
    expect_false(is.null(s$anno))                        # condition annotation snapshot
    # value matrix cached in derived under expr_hm_value_mat, keyed on the assay
    expect_equal(get("expr_hm_value_mat", envir = state$derived)$params, list("logcounts"))

    expect_match(s$legend, "z-score")                    # z-score keeps the value label
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("circlize")
    cf <- .hm_col_fun(s$mat, s$ramp, s$zscored)          # ramp from the snapshot
    ht <- .hm_build(s, NULL, cf)                          # snapshot -> a real Heatmap
    expect_s4_class(ht, "Heatmap")
  })
})

test_that("heatmap respects the Render gate (no auto-render) + stale banner", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 4, n_spike = 4, seed = 6)),
               source = "demo", meta = list(feature_type = "gene"))
    rn <- rownames(state$working)
    state$gene_sets <- list(SetA = new_gene_set(rn[1:6]), SetB = new_gene_set(rn[10:24]))
    session$setInputs(tabs = "Gene sets", hm_source = "saved", hm_pick = "SetA",
                      hm_val_assay = "logcounts", hm_val_transform = "none",
                      hm_zscore = TRUE, hm_only_expr = TRUE, hm_ramp_src = "custom",
                      hm_row_mode = "auto", hm_col_mode = "auto",
                      hm_cluster_rows = TRUE, hm_cluster_cols = TRUE,
                      hm_row_dend = "auto", hm_col_dend = "auto", hm_render = 1)
    session$flushReact()
    expect_equal(nrow(hm_out$value()$mat), 6L)           # SetA

    # switch set WITHOUT rendering -> snapshot unchanged, stale banner fires
    session$setInputs(hm_pick = "SetB"); session$flushReact()
    expect_equal(nrow(hm_out$value()$mat), 6L)           # still SetA
    expect_true(hm_out$stale())

    session$setInputs(hm_render = 2); session$flushReact()
    expect_equal(nrow(hm_out$value()$mat), 15L)          # SetB (rn[10:24])
    expect_false(hm_out$stale())

    # For this slow static heatmap ALL aesthetics are gated behind Render (unlike
    # the fast dual_plot pages): changing the colour source must STALE the plot,
    # not auto-recolour.
    session$setInputs(hm_ramp_src = "palette"); session$flushReact()
    expect_true(hm_out$stale())
  })
})

test_that("heatmap 'selected' row labels mark searched genes + report coverage", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_expression_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 4, n_spike = 4, seed = 7)),
               source = "demo", meta = list(feature_type = "gene"))
    rn <- rownames(state$working)
    state$gene_sets <- list(SetA = new_gene_set(rn[1:10]))
    gn <- SummarizedExperiment::rowData(state$working)$gene_name
    in_set  <- gn[1:2]     # two genes inside SetA (plotted rows)
    out_set <- gn[20]      # resolves to a real dds id, but NOT in SetA -> unshowable
    session$setInputs(tabs = "Gene sets", hm_source = "saved", hm_pick = "SetA",
                      hm_val_assay = "logcounts", hm_val_transform = "none",
                      hm_zscore = TRUE, hm_only_expr = TRUE, hm_ramp_src = "custom",
                      hm_row_mode = "selected", hm_col_mode = "none",
                      hmrowsel_searchby = "gene_name",
                      hmrowsel_q = paste(c(in_set, out_set), collapse = ", "),
                      hm_cluster_rows = TRUE, hm_cluster_cols = TRUE,
                      hm_row_dend = "auto", hm_col_dend = "auto", hm_render = 1)
    session$flushReact()
    s <- hm_out$value()
    expect_equal(s$row_mode, "selected")
    expect_length(s$row_mark_at, 2L)                     # 2 in-set genes marked
    expect_false(is.null(s$row_cov))
    expect_equal(s$row_cov$n_hidden, 1L)                 # the out-of-set gene can't be shown
  })
})
