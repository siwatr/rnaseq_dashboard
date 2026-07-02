# Gene Sets page server (R/mod_geneset.R) -- Manage tab (P6b).

test_that("mod_geneset UI mounts the Manage controls + tables", {
  ui <- as.character(mod_geneset_ui("gs"))
  expect_match(ui, "gs-set_name")
  expect_match(ui, "gs-add_mode")
  expect_match(ui, "gs-sets_table")
  expect_match(ui, "gs-members_table")
  expect_match(ui, "gs-paste_q")            # the shared gene-search box (suffix paste)
})

test_that("paste source adds a set; New auto-suffixes; Append unions; no data_version bump", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 1))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    dv <- state$data_version
    r1 <- rownames(state$working)[1]; r2 <- rownames(state$working)[2]
    r3 <- rownames(state$working)[3]

    session$setInputs(paste_searchby = "__rownames__", paste_mode = "exact", paste_ci = FALSE,
                      set_name = "mySet", add_mode = "new", paste_literal = FALSE,
                      paste_q = paste(r1, r2, sep = ", "))
    session$elapse(300); session$flushReact()
    session$setInputs(add_paste = 1); session$flushReact()
    expect_named(state$gene_sets, "mySet")
    expect_setequal(state$gene_sets$mySet$ids, c(r1, r2))
    expect_equal(state$gene_sets$mySet$source, "paste")
    expect_equal(state$data_version, dv)                 # a gene-set add is NOT a data edit

    # New again with the same name -> auto-suffixed, original untouched.
    session$setInputs(paste_q = r3); session$elapse(300); session$flushReact()
    session$setInputs(add_paste = 2); session$flushReact()
    expect_true("mySet_2" %in% names(state$gene_sets))
    expect_setequal(state$gene_sets$mySet$ids, c(r1, r2))

    # Append into mySet unions.
    session$setInputs(add_mode = "append", set_name = "mySet",
                      paste_q = r3); session$elapse(300); session$flushReact()
    session$setInputs(add_paste = 3); session$flushReact()
    expect_setequal(state$gene_sets$mySet$ids, c(r1, r2, r3))
  })
})

test_that("paste 'add unmatched as literal IDs' authors out-of-dataset members", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 2))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]

    session$setInputs(paste_searchby = "__rownames__", paste_mode = "exact", paste_ci = FALSE,
                      set_name = "lit", add_mode = "new", paste_literal = TRUE,
                      paste_q = paste(r1, "GHOST_ID", sep = ", "))
    session$elapse(300); session$flushReact()
    session$setInputs(add_paste = 1); session$flushReact()
    expect_setequal(state$gene_sets$lit$ids, c(r1, "GHOST_ID"))
    # The literal id is absent from the dataset (non-destructive membership).
    expect_equal(gene_set_absent(state$gene_sets$lit, rownames(state$working)), "GHOST_ID")
  })
})

test_that("DE DEGs source snapshots ids for the chosen direction/thresholds", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 3))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    rn <- rownames(state$working)
    # Fabricate a DE result: rows 1-2 up, 3-4 down, rest no-change.
    df <- data.frame(
      baseMean = 100, padj = c(0.001, 0.001, 0.001, 0.001, rep(0.9, length(rn) - 4)),
      log2FoldChange = c(3, 2, -3, -2, rep(0, length(rn) - 4)), row.names = rn)
    de <- state$de; de$results <- list("C1" = df); de$active <- "C1"; state$de <- de
    session$flushReact()

    session$setInputs(deg_contrast = "C1", deg_dir = "up", deg_padj = 0.05, deg_lfc = 1,
                      deg_shrunk = FALSE, set_name = "up_genes", add_mode = "new")
    session$flushReact()
    session$setInputs(add_deg = 1); session$flushReact()
    expect_setequal(state$gene_sets$up_genes$ids, rn[1:2])
    expect_match(state$gene_sets$up_genes$source, "^DE: C1")

    session$setInputs(deg_dir = "both", set_name = "both_genes"); session$flushReact()
    session$setInputs(add_deg = 2); session$flushReact()
    expect_setequal(state$gene_sets$both_genes$ids, rn[1:4])
  })
})

test_that("top-variable source adds the requested count", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 40, n_per_group = 3, n_spike = 4, seed = 4))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    session$setInputs(topvar_n = 5, set_name = "hv", add_mode = "new"); session$flushReact()
    session$setInputs(add_topvar = 1); session$flushReact()
    expect_length(state$gene_sets$hv$ids, 5L)
    # endogenous only -> no spike-in ids
    expect_false(any(grepl("^ERCC", state$gene_sets$hv$ids)))
  })
})

test_that("membership is non-destructive across a feature drop; load clears sets", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 5))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]; r2 <- rownames(state$working)[2]
    session$setInputs(paste_searchby = "__rownames__", paste_mode = "exact", paste_ci = FALSE,
                      set_name = "s", add_mode = "new", paste_literal = FALSE,
                      paste_q = paste(r1, r2, sep = ", "))
    session$elapse(300); session$flushReact()
    session$setInputs(add_paste = 1); session$flushReact()
    expect_setequal(state$gene_sets$s$ids, c(r1, r2))

    # Drop feature r1 -> the set keeps its full authored membership (non-destructive).
    state_mutate(state, function(d) d[-1, ], action = list(action = "drop_feature"))
    session$flushReact()
    expect_setequal(state$gene_sets$s$ids, c(r1, r2))            # unchanged
    expect_equal(gene_set_absent(state$gene_sets$s, rownames(state$working)), r1)

    # Loading a new object clears gene sets (id-keyed -> meaningless on a new dataset).
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    expect_length(state$gene_sets, 0L)
  })
})

test_that("rename and delete manage the store + log history", {
  skip_if_not_installed("DESeq2")
  state <- new_app_state()
  shiny::testServer(mod_geneset_server, args = list(state = state), {
    dds <- ensure_logcounts(make_mock_dds(n_genes = 20, n_per_group = 3, n_spike = 0, seed = 6))
    state_load(state, dds, source = "demo", meta = list(feature_type = "gene"))
    session$flushReact()
    r1 <- rownames(state$working)[1]
    session$setInputs(paste_searchby = "__rownames__", paste_mode = "exact", paste_ci = FALSE,
                      set_name = "orig", add_mode = "new", paste_literal = FALSE, paste_q = r1)
    session$elapse(300); session$flushReact()
    session$setInputs(add_paste = 1); session$flushReact()
    expect_equal(selected_set(), "orig")               # add selects the set

    # Rename orig -> renamed.
    session$setInputs(rename_to = "renamed"); session$setInputs(rename_ok = 1); session$flushReact()
    expect_named(state$gene_sets, "renamed")
    expect_equal(selected_set(), "renamed")

    # Delete it.
    session$setInputs(delete = 1); session$flushReact()
    expect_length(state$gene_sets, 0L)
    expect_null(selected_set())
  })
})
