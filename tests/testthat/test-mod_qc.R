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

test_that("dataset/sample ggplot builders return ggplots (light + dark)", {
  skip_if_not_installed("DESeq2")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 2, n_spike = 2, seed = 1))
  # Mean-SD is now a point-density scatter (MASS-based), no vsn/hexbin.
  vst_mat <- SummarizedExperiment::assay(qc_vst(dds))
  expect_s3_class(ddsdashboard:::.qc_meansd_plot(vst_mat, FALSE), "ggplot")
  expect_s3_class(ddsdashboard:::.qc_meansd_plot(vst_mat, TRUE), "ggplot")

  cond <- factor(SummarizedExperiment::colData(dds)$condition)
  rle <- qc_rle_matrix(dds)              # endogenous-only -> nrow(rle) <= nrow(dds)
  rle_long <- data.frame(
    sample = factor(rep(colnames(dds), each = nrow(rle)), levels = colnames(dds)),
    group  = rep(cond, each = nrow(rle)),
    value  = as.numeric(rle))
  expect_s3_class(ddsdashboard:::.qc_rle_plot(rle_long, FALSE, ncol(dds)), "ggplot")

  expr_long <- qc_expression_long(dds)
  expr_long$group <- rep(cond, each = nrow(expr_long) / ncol(dds))
  expect_s3_class(ddsdashboard:::.qc_density_plot(expr_long, TRUE, ncol(dds)), "ggplot")

  wg <- qc_within_group_correlation(dds)
  expect_s3_class(ddsdashboard:::.qc_within_group_plot(wg, TRUE), "ggplot")
})

test_that(".qc_correlation_heatmap handles annotation variants + labels the method", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("ComplexHeatmap")
  dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 1, seed = 2))
  cm <- qc_sample_correlation(dds)
  cd <- as.data.frame(SummarizedExperiment::colData(dds))

  # Multiple annotation columns (e.g. ~ condition + bio_rep); method label set.
  ht_multi <- ddsdashboard:::.qc_correlation_heatmap(
    cm, cd[, c("condition", "bio_rep"), drop = FALSE], n_samples = ncol(dds),
    method = "pearson")
  expect_s4_class(ht_multi, "Heatmap")
  expect_match(as.character(ht_multi@column_title), "Pearson")  # method shown in title

  # No annotation.
  ht_none <- ddsdashboard:::.qc_correlation_heatmap(cm, NULL, n_samples = ncol(dds))
  expect_s4_class(ht_none, "Heatmap")
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
                      cor_method = "spearman",
                      cor_anno = c("condition", "bio_rep"),
                      cor_auto = FALSE, cor_render = 1)
    expect_true(exists("vst", envir = state$derived, inherits = FALSE))
    expect_equal(get("vst", envir = state$derived)$version, state$data_version)
    expect_true(exists("sample_cor", envir = state$derived, inherits = FALSE))
    expect_equal(get("sample_cor", envir = state$derived)$version, state$data_version)
  })
})

# ---- Filtering (P3c) --------------------------------------------------------

test_that("applying a feature removal shrinks the dataset and logs the action", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3,
                                                      n_spike = 2, seed = 1)), source = "demo")
    session$setInputs(feat_use_fbe = FALSE, feat_min_count = 0, feat_use_min_samples = FALSE)
    session$flushReact()
    drop <- rownames(state$working)[1:3]
    feat_pool(drop)
    n0 <- nrow(state$working); v0 <- state$data_version
    session$setInputs(feat_apply_ok = 1)                    # confirm-modal OK
    expect_equal(nrow(state$working), n0 - 3L)
    expect_false(any(drop %in% rownames(state$working)))
    expect_gt(state$data_version, v0)
    last <- state$history[[length(state$history)]]
    expect_equal(last$action, "filter_features")
    expect_equal(last$n_dropped, 3L)
  })
})

test_that("feature pool buttons adopt suggestions (union) and clear", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3,
                                                      n_spike = 2, seed = 2)), source = "demo")
    total <- rowSums(as.matrix(SummarizedExperiment::assay(state$working, "counts")))
    session$setInputs(feat_use_fbe = FALSE, feat_min_count = stats::median(total),
                      feat_use_min_samples = FALSE)
    session$flushReact()
    sugg <- feat_flags()$feature_id[feat_flags()$suggested_drop]
    expect_gt(length(sugg), 0)
    feat_pool(character(0))
    session$setInputs(feat_adopt = 1)
    expect_setequal(feat_pool(), sugg)                      # union with empty = suggestions
    session$setInputs(feat_clear = 1)
    expect_length(feat_pool(), 0)
  })
})

test_that("flagged samples are highlight-only (pool starts empty)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 3,
                                                      n_spike = 1, seed = 3)), source = "demo")
    session$flushReact()
    expect_length(samp_pool(), 0)
    # samp_flags still computes the per-reason schema.
    expect_true(all(c("flagged", "within_group_outlier") %in% colnames(samp_flags())))
  })
})

test_that("samp_flags incorporates spike-in criteria when a rule is set", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 3,
                                                      n_spike = 6, seed = 31)), source = "demo")
    session$flushReact()
    expect_gt(length(spike_ids()), 0)
    expect_true(all(c("low_spike_detected", "bad_slope", "few_spike_points") %in%
                    colnames(samp_flags())))
    # An impossible detected-spike floor flags every sample for that reason; the
    # derive must re-key on the new spike threshold.
    session$setInputs(samp_spike_detected_min = 9999)
    session$flushReact()
    fl <- samp_flags()
    expect_true(all(fl$low_spike_detected))
    expect_true(all(grepl("few detected spikes", fl$reason)))
  })
})

test_that("spike filtering is inert when the dataset has no spike-ins", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 3,
                                                      n_spike = 0, seed = 32)), source = "demo")
    session$flushReact()
    expect_equal(length(spike_ids()), 0)
    session$setInputs(samp_spike_detected_min = 9999)
    session$flushReact()
    fl <- samp_flags()
    expect_false(any(fl$low_spike_detected))
    expect_true(all(is.na(fl$n_spike_detected)))
  })
})

test_that("spike observed assay prefers a length-normalized assay (TPM > FPKM > CPM)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    dds <- add_normalized_assays(
      ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 2, n_spike = 4, seed = 33)),
      c("CPM", "TPM"))
    state_load(state, dds, source = "demo")
    session$flushReact()
    expect_equal(spike_default_assay(), "TPM")
    expect_equal(samp_spike_assay(), "TPM")        # fallback before the spike tab UI renders
  })
})

test_that("spike observed assay falls back to CPM without TPM/FPKM", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 2,
                                                      n_spike = 4, seed = 34)), source = "demo")
    session$flushReact()
    expect_equal(spike_default_assay(), "CPM")
  })
})

test_that("Showing subset filters plotted samples without bumping data_version", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3,
                                                      n_spike = 1, seed = 4)), source = "demo")
    session$flushReact()
    v0 <- state$data_version
    # Edit one tab's control; the canonical state + showing_samples() follow.
    session$setInputs(gen_show_by = "condition")
    session$setInputs(gen_show_values = "control")
    session$flushReact()
    shown <- showing_samples()
    cd <- as.data.frame(SummarizedExperiment::colData(state$working))
    expect_true(all(as.character(cd[shown, "condition"]) == "control"))
    expect_lt(length(shown), ncol(state$working))
    expect_equal(state$data_version, v0)                    # view-only: no bump
  })
})

test_that("Showing controls share one canonical selection across tabs", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 30, n_per_group = 3,
                                                      n_spike = 1, seed = 5)), source = "demo")
    session$flushReact()
    session$setInputs(rle_show_by = "condition")            # edit via the RLE tab
    session$setInputs(rle_show_values = "treated")
    session$flushReact()
    # The General-QC plot honors the same selection (canonical, not per-tab).
    cd <- as.data.frame(SummarizedExperiment::colData(state$working))
    expect_setequal(showing_samples(), rownames(cd)[cd$condition == "treated"])
  })
})

test_that("a deferred plot reports stale when a setting changes before re-render", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 2,
                                                      n_spike = 1, seed = 6)), source = "demo")
    session$setInputs(x_axis = "sample", metric = "library_size", group = "condition",
                      sort = "none", auto = FALSE, render = 1)
    session$flushReact()
    expect_false(isTRUE(gen_shown$stale()))           # just rendered
    session$setInputs(metric = "detected")            # a deferred setting changes
    session$flushReact()
    expect_true(isTRUE(gen_shown$stale()))            # stale until re-render
    session$setInputs(render = 2)
    session$flushReact()
    expect_false(isTRUE(gen_shown$stale()))
  })
})

test_that("Filtering tables use boolean columns; selection moves ids in/out of the pool", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 2,
                                                      n_spike = 1, seed = 7)), source = "demo")
    session$setInputs(feat_use_fbe = FALSE, feat_min_count = 0, feat_use_min_samples = FALSE)
    session$flushReact()
    disp <- feat_display(feat_flags(), character(0))
    expect_type(disp[["Suggested Removal"]], "logical")
    expect_type(disp[["In Removal Pool"]], "logical")
    ids <- feat_flags()$feature_id
    session$setInputs(feat_tbl_rows_selected = 1:3)
    session$setInputs(feat_add = 1)
    expect_setequal(feat_pool(), ids[1:3])
    session$setInputs(feat_remove = 1)
    expect_length(feat_pool(), 0)
  })
})

test_that("QC UI exposes the Filtering tab, pool actions, and per-sidebar Showing", {
  html <- paste(as.character(mod_qc_ui("qc")), collapse = " ")
  expect_match(html, "Filtering")
  expect_match(html, "Removal Pool")                 # pool-button section header
  expect_match(html, "Select all")
  expect_match(html, "Add selected to pool")
  expect_match(html, "Remove Samples")               # renamed apply buttons
  expect_match(html, "Remove Features")
  expect_match(html, "Set auto thresholds")          # scoped general-filters Auto button
  expect_match(html, "Sample QC filters")            # collapsible filter group
  expect_match(html, "Plot Showing")                 # per-sidebar Showing accordion
  expect_match(html, "gen_show_by")                  # one of the synced controls
  expect_match(html, "spike_show_by")                # Showing now on the spike tab too
  expect_match(html, "samp_lib_auto")                # a per-field Auto button
})

test_that("QC per-tab Reset Feature Removal restores removed features and is undoable", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 2,
                                                     n_spike = 2, seed = 8)), source = "demo")
    session$setInputs(feat_use_fbe = FALSE, feat_min_count = 0, feat_use_min_samples = FALSE)
    session$flushReact()
    n_full <- nrow(state$working)
    feat_pool(rownames(state$working)[1:4])
    session$setInputs(feat_apply_ok = 1)
    expect_equal(nrow(state$working), n_full - 4L)
    v <- state$data_version
    session$setInputs(feat_reset = 1)
    expect_equal(nrow(state$working), n_full)          # all features restored
    expect_gt(state$data_version, v)
  })
})

# ---- Spike-in (ERCC) dose-response (P3d) ------------------------------------

test_that("spike-in plot builders return ggplots (light + dark)", {
  skip_if_not_installed("DESeq2")
  dds <- make_mock_dds(n_genes = 60, n_per_group = 2, n_spike = 10, seed = 1)
  dr <- spike_dose_response(dds, assay = "CPM", source = "column")
  long <- dr$long; long$group <- factor(rep("a", nrow(long)))
  ps <- dr$per_sample; ps$group <- factor(rep("a", nrow(ps)))
  expect_s3_class(ddsdashboard:::.qc_spike_dr_plot(long, FALSE), "ggplot")
  expect_s3_class(ddsdashboard:::.qc_spike_dr_plot(long, TRUE), "ggplot")
  expect_s3_class(ddsdashboard:::.qc_spike_summary_plot(ps, "r_squared", TRUE), "ggplot")
  expect_s3_class(ddsdashboard:::.qc_spike_summary_plot(ps, "pct_spike", FALSE), "ggplot")
})

test_that("QC spike-in view caches dose-response keyed on source + assay", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 2,
                                                     n_spike = 8, seed = 2)), source = "demo")
    session$flushReact()
    session$setInputs(spike_source = "mix1", spike_assay = "CPM",
                      spike_group = "condition", spike_auto = FALSE, spike_render = 1)
    expect_true(exists("spike_dr", envir = state$derived, inherits = FALSE))
    cached <- get("spike_dr", envir = state$derived)
    expect_equal(cached$version, state$data_version)
    expect_equal(nrow(cached$value$per_sample), ncol(state$working))
  })
})

test_that("QC 'Remove all spike-in features' drops spikes and round-trips via reset", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_qc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 2,
                                                     n_spike = 6, seed = 3)), source = "demo")
    session$flushReact()
    n0 <- nrow(state$working)
    session$setInputs(feat_drop_spike_ok = 1)               # confirm-modal OK
    expect_equal(nrow(state$working), n0 - 6L)
    expect_equal(state_meta(state)$n_spike_in, 0L)
    session$setInputs(feat_reset = 1)                       # Reset Feature Removal
    expect_equal(nrow(state$working), n0)
    expect_equal(state_meta(state)$n_spike_in, 6L)
  })
})

test_that("QC UI exposes the Spike-in QC Matrix pill (table isolated from summary)", {
  html <- paste(as.character(mod_qc_ui("qc")), collapse = " ")
  expect_match(html, "Spike-in QC Matrix", fixed = TRUE)
  expect_match(html, "Per-sample summary", fixed = TRUE)
})
