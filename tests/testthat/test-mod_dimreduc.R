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
