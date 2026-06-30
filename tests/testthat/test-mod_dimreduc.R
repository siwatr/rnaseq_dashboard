# PCA page (R/mod_dimreduc.R) for Phase 4a.

test_that("PCA page UI mounts the controls + dual-plot containers", {
  ui <- as.character(mod_dimreduc_ui("dr"))
  expect_match(ui, "dr-pca_container")
  expect_match(ui, "dr-scree_container")
  expect_match(ui, "dr-render")
})

test_that("embedding computes, caches, and is stable across a 'Showing' change", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 150, n_per_group = 3, n_spike = 6, seed = 1)),
               source = "demo", meta = list(feature_type = "gene"))
    session$setInputs(assay = "vst", n_top = 100, auto = TRUE); session$flushReact()
    v1 <- pca_shown$value()
    expect_false(is.null(v1))
    expect_equal(get("pca", envir = state$derived)$params, list("vst", 100L, FALSE))
    # A display-only 'Showing' change must NOT recompute the embedding (stable axes).
    session$setInputs(pca_show_by = "condition", pca_show_values = "control"); session$flushReact()
    expect_identical(pca_shown$value(), v1)
    # Changing n_top *does* recompute (new params).
    session$setInputs(n_top = 60); session$flushReact()
    expect_equal(get("pca", envir = state$derived)$params, list("vst", 60L, FALSE))
    expect_false(identical(pca_shown$value()$scores, v1$scores))
  })
})

test_that("scatter colours by metadata + gene; scree builds; gene-not-found validates", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 4, seed = 2)),
               source = "demo", meta = list(feature_type = "gene"))
    session$setInputs(assay = "vst", n_top = 80, auto = TRUE,
                      pc_x = "PC1", pc_y = "PC2"); session$flushReact()
    session$setInputs(colour_by = "condition"); session$flushReact()
    expect_s3_class(build_pca_gg(FALSE), "ggplot")
    expect_s3_class(build_scree_gg(FALSE), "ggplot")
    # Colour by a real gene resolves (search the gene_name field); a bogus one
    # validates with a clear message. The gene box is debounced -> elapse the timer.
    gname <- SummarizedExperiment::rowData(state$working)$gene_name[1]
    session$setInputs(colour_by = "__gene__", gene_searchby = "gene_name",
                      gene_assay = "logcounts", gene = gname)
    session$elapse(300); session$flushReact()
    expect_s3_class(build_pca_gg(FALSE), "ggplot")
    session$setInputs(gene = "NoSuchGene"); session$elapse(300); session$flushReact()
    expect_error(build_pca_gg(FALSE), "not found")
  })
})

test_that("gene search: field selector, transform, duplicate + suggestion hints", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 4, seed = 11))
    rd <- SummarizedExperiment::rowData(dds)
    rd$gene_name[2] <- rd$gene_name[1]                 # force a duplicate name
    SummarizedExperiment::rowData(dds) <- rd
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$setInputs(assay = "vst", n_top = 80, auto = TRUE, colour_by = "__gene__",
                      gene_searchby = "gene_name", gene_assay = "logcounts"); session$flushReact()
    dupname <- rd$gene_name[1]
    session$setInputs(gene = dupname); session$elapse(300); session$flushReact()
    expect_match(as.character(output$gene_hint$html), "matched")     # "N features matched ...; showing the first"
    expect_s3_class(build_pca_gg(FALSE), "ggplot")                   # first match used

    # A near miss yields a "Did you mean ...?" suggestion (always case-insensitive).
    real <- rd$gene_name[50]
    session$setInputs(gene_ci = FALSE, gene = tolower(substr(real, 1, nchar(real) - 1)))
    session$elapse(300); session$flushReact()
    expect_match(as.character(output$gene_hint$html), "Did you mean")

    # A case-insensitive hit labels the legend + caption with the STORED value
    # (e.g. "Gene45"), not the typed lowercase query.
    session$setInputs(gene_ci = TRUE, gene = tolower(real)); session$elapse(300); session$flushReact()
    expect_match(build_pca_gg(FALSE)$labels$colour, real, fixed = TRUE)
    cap <- as.character(output$gene_caption$html)
    expect_match(cap, "Plotting expression of")
    expect_match(cap, real, fixed = TRUE)                            # matched value
    expect_match(cap, rownames(state$working)[50], fixed = TRUE)     # true unique id (rowname)

    # Searching by Feature ID (rownames) resolves an id directly.
    rid <- rownames(state$working)[3]
    session$setInputs(gene_ci = FALSE, gene_searchby = "__rownames__", gene = rid)
    session$elapse(300); session$flushReact()
    expect_silent(g <- build_pca_gg(FALSE))
    expect_s3_class(g, "ggplot")

    # A log transform applies via expr_transform (colourbar label reflects it).
    session$setInputs(gene_transform = "log2", gene_pseudo = 1); session$elapse(300); session$flushReact()
    g <- build_pca_gg(FALSE)
    expect_match(g$labels$colour, "log2")
  })
})

test_that("colour-by optgroups are ordered General -> This session -> Data metadata", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 50, n_per_group = 3, n_spike = 2, seed = 21)),
               source = "demo"); session$flushReact()
    html <- as.character(output$colour_ui$html)
    pos <- function(lbl) regexpr(sprintf('label="%s"', lbl), html, fixed = TRUE)
    expect_true(pos("General") < pos("This session"))
    expect_true(pos("This session") < pos("Data metadata"))
  })
})

test_that("PCA offers + colours by a spike-in metric when the dataset has spikes", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 90, n_per_group = 3, n_spike = 6, seed = 41)),
               source = "demo"); session$flushReact()
    expect_match(as.character(output$colour_ui$html), "Spike-in")          # optgroup present
    expect_match(as.character(output$colour_ui$html), "__spike__slope", fixed = TRUE)
    session$setInputs(assay = "vst", n_top = 60, auto = TRUE,
                      colour_by = "__spike__n_spike_detected"); session$flushReact()
    g <- build_pca_gg(FALSE)
    expect_s3_class(g, "ggplot")
    expect_equal(g$labels$colour, "Detected spike features")
  })
})

test_that("colour by a per-sample QC metric builds a continuous scale", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 90, n_per_group = 3, n_spike = 4, seed = 12)),
               source = "demo")
    session$setInputs(assay = "vst", n_top = 60, auto = TRUE,
                      colour_by = "__qc__library_size"); session$flushReact()
    g <- build_pca_gg(FALSE)
    expect_s3_class(g, "ggplot")
    expect_equal(g$labels$colour, "Library size")
    # Metrics come from the shared derived cache (keyed on data_version), not a
    # private recompute -> the QC page and PCA page share one frame.
    expect_false(is.null(get0("qc_metrics", envir = state$derived)))
  })
})

test_that("PCA colours/shapes by promoted removal flags + pool (shared state)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 2, seed = 31))
    state_load(state, dds, source = "demo")
    samples <- colnames(state$working)
    # Simulate what the QC Filtering page promotes into shared state.
    state$samp_flags <- data.frame(sample = samples,
                                   flagged = c(TRUE, rep(FALSE, length(samples) - 1L)),
                                   stringsAsFactors = FALSE)
    state$samp_pool <- samples[2]
    session$setInputs(assay = "vst", n_top = 50, auto = TRUE,
                      colour_by = "__removal__", shape_by = "__pool__"); session$flushReact()
    g <- build_pca_gg(FALSE)
    expect_s3_class(g, "ggplot")
    expect_equal(g$labels$colour, "Suggested removal")
    expect_equal(g$labels$shape, "Removal pool")
    # Swap the roles: colour by pool, shape by removal status.
    session$setInputs(colour_by = "__pool__", shape_by = "__removal__"); session$flushReact()
    g2 <- build_pca_gg(FALSE)
    expect_equal(g2$labels$colour, "Removal pool")
    expect_equal(g2$labels$shape, "Suggested removal")
    # The session items are offered in both selectors' "This session" group.
    expect_match(as.character(output$colour_ui$html), "__removal__", fixed = TRUE)
    expect_match(as.character(output$shape_ui$html), "__pool__", fixed = TRUE)
  })
})

test_that("PCA removal colour validates clearly when flags are not yet computed", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 2, seed = 32)),
               source = "demo")                       # state$samp_flags is NULL
    session$setInputs(assay = "vst", n_top = 40, auto = TRUE,
                      colour_by = "__removal__"); session$flushReact()
    expect_error(build_pca_gg(FALSE), "not ready")
  })
})

test_that("shape selector offers only <=6-level discrete columns; legend position applies", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 4, n_spike = 2, seed = 13))
    # A high-cardinality discrete column (one level per sample, 8 > 6) must be excluded.
    SummarizedExperiment::colData(dds)$sample_id <- factor(colnames(dds))
    state_load(state, dds, source = "demo"); session$flushReact()
    sh <- as.character(output$shape_ui$html)
    expect_match(sh, "condition")
    expect_false(grepl('value=\"sample_id\"', sh))
    session$setInputs(assay = "vst", n_top = 40, auto = TRUE, legend_pos = "none"); session$flushReact()
    g <- build_pca_gg(FALSE)
    expect_equal(g$theme$legend.position, "none")
  })
})

test_that("PCA needs >=3 samples (clear validate message)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 1, n_spike = 2, seed = 3)),
               source = "demo")                       # 2 samples
    session$setInputs(assay = "vst", n_top = 50, auto = TRUE); session$flushReact()
    expect_error(pca_spec(), "at least 3 samples")
  })
})

test_that("stale banner fires on an embedding-input change when auto-render is off", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 4, seed = 5)),
               source = "demo")
    session$setInputs(assay = "vst", n_top = 100, auto = FALSE, render = 1); session$flushReact()
    expect_false(is.null(pca_shown$value()))
    expect_false(pca_shown$stale())
    session$setInputs(n_top = 50); session$flushReact()      # changed, not yet re-rendered
    expect_true(pca_shown$stale())
    session$setInputs(render = 2); session$flushReact()
    expect_false(pca_shown$stale())
  })
})

test_that("PC-axis selection is preserved across a re-render (not reset to PC1/PC2)", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 120, n_per_group = 3, n_spike = 4, seed = 6)),
               source = "demo")
    session$setInputs(assay = "vst", n_top = 100, auto = TRUE, pc_x = "PC3"); session$flushReact()
    session$setInputs(n_top = 60); session$flushReact()       # forces pc_ui to re-render
    html <- as.character(output$pc_ui$html)
    expect_match(html, "PC3")                                 # still offered + kept (isolate)
    expect_match(html, 'value=\"PC3\"[^>]*selected', perl = TRUE)
  })
})

test_that("default colour-by follows the DESeq2 design's first variable, with fallbacks", {
  skip_if_not_installed("DESeq2")
  mk <- function() ensure_logcounts(make_mock_dds(n_genes = 60, n_per_group = 3, n_spike = 2, seed = 7))
  # ~ condition (default mock design) -> condition selected.
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, mk(), source = "demo"); session$flushReact()
    expect_match(as.character(output$colour_ui$html), 'value=\"condition\"[^>]*selected', perl = TRUE)
  })
  # ~ batch + condition -> the FIRST design var (batch).
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    dds <- mk()
    SummarizedExperiment::colData(dds)$batch <- factor(rep(c("a", "b"), length.out = ncol(dds)))
    DESeq2::design(dds) <- ~ batch + condition
    state_load(state, dds, source = "demo"); session$flushReact()
    expect_match(as.character(output$colour_ui$html), 'value=\"batch\"[^>]*selected', perl = TRUE)
  })
  # No usable design + no "condition" column -> (none).
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    dds <- mk(); cd <- SummarizedExperiment::colData(dds)
    colnames(cd)[colnames(cd) == "condition"] <- "grp"
    SummarizedExperiment::colData(dds) <- cd
    DESeq2::design(dds) <- ~ 1
    state_load(state, dds, source = "demo"); session$flushReact()
    expect_match(as.character(output$colour_ui$html), 'value=\"__none__\"[^>]*selected', perl = TRUE)
  })
})

test_that("discrete colour uses thematic default (no manual scale) unless a Palette config is set", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 2, seed = 8)),
               source = "demo")
    session$setInputs(assay = "vst", n_top = 60, auto = TRUE, colour_by = "condition"); session$flushReact()
    colour_scales <- function(g) Filter(function(s) "colour" %in% s$aesthetics, g$scales$scales)
    expect_length(colour_scales(build_pca_gg(FALSE)), 0L)         # no config -> thematic default
    # Configure a Palette mapping for condition -> a manual colour scale appears.
    lv <- levels(SummarizedExperiment::colData(state$working)$condition)
    state$palette <- list(colData = list(condition = list(name = "Okabe-Ito",
      colors = palette_discrete(lv, NULL, "Okabe-Ito"))))
    session$flushReact()
    expect_gte(length(colour_scales(build_pca_gg(FALSE))), 1L)
  })
})

test_that("plot-aesthetic controls reach the plot (point size + 1:1 ratio); dark mode builds", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_dimreduc_server, args = list(state = state, dark_mode = reactive(TRUE)), {
    state_load(state, ensure_logcounts(make_mock_dds(n_genes = 80, n_per_group = 3, n_spike = 2, seed = 9)),
               source = "demo")
    session$setInputs(assay = "vst", n_top = 60, auto = TRUE,
                      point_size = 1.5, fixed_ratio = TRUE); session$flushReact()
    g <- build_pca_gg(FALSE)
    expect_equal(g$layers[[1]]$aes_params$size, 1.5)
    expect_equal(g$coordinates$ratio, 1)                          # coord_fixed applied
    expect_s3_class(build_scree_gg(FALSE), "ggplot")              # dark-mode path builds
  })
})
