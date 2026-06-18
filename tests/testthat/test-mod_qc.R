test_that("QC module caches the metric table keyed on data_version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 2, seed = 1),
               source = "demo")
    session$flushReact()
    session$setInputs(x_axis = "sample", metric = "library_size",
                      group = "condition", sort = "none", auto = FALSE)
    session$setInputs(render = 1)

    # The derived cache holds the metric table at the current data_version.
    expect_true(exists("qc_metrics", envir = state$derived, inherits = FALSE))
    cached <- get("qc_metrics", envir = state$derived)
    expect_equal(cached$version, state$data_version)
    expect_equal(nrow(cached$value), ncol(state$working))
    expect_setequal(colnames(cached$value),
                    c("sample", "library_size", "detected", "pct_mito", "pct_spike"))
  })
})

test_that("QC module re-renders when auto-render is enabled (no button press)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 1, seed = 2),
               source = "demo")
    session$flushReact()
    session$setInputs(x_axis = "pct_spike", metric = "detected",
                      group = "condition", sort = "none", auto = TRUE)
    session$flushReact()
    expect_false(is.null(output$plot))
  })
})

test_that(".qc_metric_plot builds for sample and metric x-axes, both themes", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 2, seed = 1)
  tbl <- qc_per_sample_metrics(dds)
  tbl$group <- factor(rep(c("a", "b"), length.out = nrow(tbl)))

  # Discrete x (per-sample bar), and numeric x (metric-vs-metric scatter).
  p_bar <- ddsdashboard:::.qc_metric_plot(tbl, "sample", "library_size", "condition")
  p_sc  <- ddsdashboard:::.qc_metric_plot(tbl, "pct_spike", "library_size", "condition",
                                          dark_theme = TRUE)
  expect_s3_class(p_bar, "ggplot")
  expect_s3_class(p_sc, "ggplot")
})

test_that(".qc_metric_plot sorts the discrete x-axis by the metric value", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 40, n_per_group = 2, n_spike = 1, seed = 5)
  tbl <- qc_per_sample_metrics(dds)
  tbl$group <- factor(rep("all", nrow(tbl)))
  p <- ddsdashboard:::.qc_metric_plot(tbl, "sample", "detected", sort = "decreasing")
  lvls <- levels(p$data$x)
  # Levels should follow detected, descending.
  ord <- tbl$sample[order(tbl$detected, decreasing = TRUE)]
  expect_equal(lvls, ord)
})

# ---- Dataset diagnostics (P3b) ---------------------------------------------

test_that("dataset-diagnostic ggplot builders return ggplots (light + dark)", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("vsn")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 2, n_spike = 2, seed = 1))
  vst_mat <- SummarizedExperiment::assay(qc_vst(dds))
  # vsn emits a benign upstream aes_string() deprecation; not our concern here.
  expect_s3_class(suppressWarnings(ddsdashboard:::.qc_meansd_plot(vst_mat, FALSE)), "ggplot")
  expect_s3_class(suppressWarnings(ddsdashboard:::.qc_meansd_plot(vst_mat, TRUE)), "ggplot")

  rle <- qc_rle_matrix(dds)
  rle_long <- data.frame(
    sample = factor(rep(colnames(dds), each = nrow(rle)), levels = colnames(dds)),
    group  = rep(factor(SummarizedExperiment::colData(dds)$condition), each = nrow(rle)),
    value  = as.numeric(rle))
  expect_s3_class(ddsdashboard:::.qc_rle_plot(rle_long, FALSE, ncol(dds)), "ggplot")

  expr_long <- qc_expression_long(dds)
  expr_long$group <- rep(factor(SummarizedExperiment::colData(dds)$condition), each = nrow(dds))
  expect_s3_class(ddsdashboard:::.qc_density_plot(expr_long, TRUE, ncol(dds)), "ggplot")
})

test_that(".qc_correlation_heatmap returns a Heatmap", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("ComplexHeatmap")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 1, seed = 2))
  cm <- qc_sample_correlation(dds)
  anno <- stats::setNames(as.character(SummarizedExperiment::colData(dds)$condition), colnames(dds))
  ht <- ddsdashboard:::.qc_correlation_heatmap(cm, anno, "condition", n_samples = ncol(dds))
  expect_s4_class(ht, "Heatmap")
})

test_that("QC module caches VST and sample correlation keyed on data_version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 2,
                                                      n_spike = 2, seed = 1)),
               source = "demo")
    session$flushReact()
    session$setInputs(diag_auto = FALSE, diag_render = 1,
                      cor_method = "spearman", cor_anno = "condition",
                      cor_auto = FALSE, cor_render = 1)
    expect_true(exists("vst", envir = state$derived, inherits = FALSE))
    expect_equal(get("vst", envir = state$derived)$version, state$data_version)
    expect_true(exists("sample_cor", envir = state$derived, inherits = FALSE))
    expect_equal(get("sample_cor", envir = state$derived)$version, state$data_version)
  })
})
