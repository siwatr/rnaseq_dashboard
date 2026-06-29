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
    # Colour by a real gene resolves; a bogus one validates with a clear message.
    gname <- SummarizedExperiment::rowData(state$working)$gene_name[1]
    session$setInputs(colour_by = "__gene__", gene = gname); session$flushReact()
    expect_s3_class(build_pca_gg(FALSE), "ggplot")
    session$setInputs(gene = "NoSuchGene"); session$flushReact()
    expect_error(build_pca_gg(FALSE), "not found")
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
